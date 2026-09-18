"""Pipeline runner orchestration.

This module contains the core pipeline logic that is shared between
the CLI and web server execution modes.
"""

import json
import logging
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metagomics2 import __version__
from metagomics2.core.aggregation import (
    AggregationResult,
    aggregate_go_taxonomy_combos,
    aggregate_peptide_annotations,
    validate_aggregation_invariants,
)
from metagomics2.core.annotation import (
    PeptideAnnotation,
    SubjectAnnotation,
    annotate_peptide,
    load_subject_annotations_from_dict,
)
from metagomics2.core.diamond import (
    DEFAULT_DIAMOND_EVALUE,
    DEFAULT_MAX_TARGET_SEQS,
    count_queries_at_cap,
    run_diamond,
)
from metagomics2.core.fasta import (
    build_protein_dict,
    parse_fasta,
    write_subset_fasta,
)
from metagomics2.core.filtering import (
    FilterPolicy,
    HomologyHit,
    filter_all_hits,
    filter_all_hits_with_hits,
)
from metagomics2.core.go import GODAG
from metagomics2.core.matching import MatchResult, match_peptides_to_proteins
from metagomics2.core.peptides import Peptide, parse_peptide_list
from metagomics2.core.reference_loader import (
    create_reference_snapshot,
    get_bundled_reference_dir,
    get_reference_metadata,
    load_go_data,
    load_taxonomy_data,
)
from metagomics2.core.reporting import (
    create_manifest,
    write_coverage_csv,
    write_go_taxonomy_combo_csv,
    write_go_terms_csv,
    write_manifest_json,
    write_peptide_mapping_parquet,
    write_taxonomy_nodes_csv,
)
from metagomics2.core.subject_lookup import load_subject_annotations
from metagomics2.core.taxonomy import TaxonomyTree
from metagomics2.logging_setup import format_bytes, peak_rss_bytes

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for a pipeline run."""

    # Input paths
    fasta_path: Path
    peptide_list_paths: list[Path]

    # Output directory
    output_dir: Path

    # Search tool configuration
    search_tool: str = "diamond"  # "diamond" or "blast"
    annotated_db_path: Path | None = None
    threads: int = 1

    # Filter policy
    filter_policy: FilterPolicy = field(default_factory=FilterPolicy)

    # Reference data paths (if None, uses bundled reference data)
    go_data_path: Path | None = None
    taxonomy_data_path: Path | None = None

    # Job directory for snapshots (set automatically for web mode)
    job_dir: Path | None = None

    # Directory for intermediate files (subset FASTA, DIAMOND output, reference
    # snapshot).  None means ``job_dir / "work"`` when job_dir is set (web
    # mode) and ``output_dir / "work"`` otherwise (CLI mode).
    work_dir: Path | None = None

    # GO closure settings
    go_edge_types: set[str] = field(default_factory=lambda: {"is_a", "part_of"})
    go_include_self: bool = True

    # Annotations database path (companion SQLite for taxonomy/GO lookup)
    annotations_db_path: Path | None = None

    # Mock mode for testing (bypasses actual homology search)
    mock_hits_path: Path | None = None
    mock_subject_annotations_path: Path | None = None

    # DIAMOND tuning (None = DIAMOND's defaults). See README "DIAMOND performance
    # and memory" for what these do.
    diamond_block_size: float | None = None
    diamond_index_chunks: int | None = None
    diamond_tmpdir: Path | None = None
    # Per-query hit cap passed to DIAMOND (0 = unlimited).  Raised to top_k if
    # top_k is larger, so the tie-aware top_k filter never sees a list that
    # DIAMOND cut shorter than K.
    diamond_max_target_seqs: int = DEFAULT_MAX_TARGET_SEQS

    def effective_max_target_seqs(self) -> int:
        """The ``--max-target-seqs`` value actually passed to DIAMOND."""
        if self.diamond_max_target_seqs == 0:
            return 0
        return max(self.diamond_max_target_seqs, self.filter_policy.top_k or 0)

    def effective_diamond_evalue(self) -> float:
        """The ``--evalue`` pre-filter actually passed to DIAMOND.

        The filter policy's ``max_evalue`` when set; otherwise DIAMOND still
        needs a threshold, and ``DEFAULT_DIAMOND_EVALUE`` is used.
        """
        if self.filter_policy.max_evalue is not None:
            return self.filter_policy.max_evalue
        return DEFAULT_DIAMOND_EVALUE


# Weighted progress milestones (out of 1000) for each pipeline stage.
# DIAMOND homology search is weighted heaviest as it is typically the
# longest-running step.
_PROGRESS_TOTAL = 1000
_PROGRESS_INITIALIZING = 50
_PROGRESS_PARSING = 80
_PROGRESS_MATCHING = 130
_PROGRESS_SUBSET_FASTA = 150
_PROGRESS_HOMOLOGY = 650
_PROGRESS_FILTERING = 670
_PROGRESS_SUBJECT_ANNOTATIONS = 750
_PROGRESS_PER_LIST_START = 750
_PROGRESS_PER_LIST_END = 1000


@dataclass
class PipelineProgress:
    """Progress tracking for pipeline execution."""

    total_peptide_lists: int = 0
    completed_peptide_lists: int = 0
    current_stage: str = ""
    current_list_id: str = ""
    progress_done: int = 0
    progress_total: int = _PROGRESS_TOTAL


@dataclass
class PeptideListResult:
    """Result for a single peptide list."""

    list_id: str
    filename: str
    n_peptides: int
    n_matched: int
    n_unmatched: int
    aggregation: AggregationResult | None = None


@dataclass
class PipelineResult:
    """Result of a complete pipeline run."""

    success: bool
    error_message: str | None = None
    peptide_list_results: list[PeptideListResult] = field(default_factory=list)
    output_dir: Path | None = None


ProgressCallback = Callable[[PipelineProgress], None]


class PipelineRunner:
    """Orchestrates the metagomics pipeline execution."""

    # Monotonic timestamp of the most recent stage change, for stage timing logs.
    _stage_started: float | None = None
    # The DIAMOND command line actually run, recorded in the manifest.
    diamond_command: str = ""

    def __init__(
        self,
        config: PipelineConfig,
        progress_callback: ProgressCallback | None = None,
    ):
        self.config = config
        self.progress_callback = progress_callback
        self.progress = PipelineProgress()

        # Loaded data
        self.proteins: dict[str, str] = {}
        self.taxonomy_tree: TaxonomyTree | None = None
        self.go_dag: GODAG | None = None
        self.subject_annotations: dict[str, SubjectAnnotation] = {}

        # Shared across peptide lists
        self.protein_to_subjects: dict[str, set[str]] = {}
        self.protein_to_subject_hits: dict[str, dict[str, HomologyHit]] = {}

        # Per-list parsed peptides and match results
        self.parsed_peptide_lists: dict[str, list[Peptide]] = {}
        self.per_list_match_results: dict[str, Any] = {}

        # Union of all hit protein IDs across all peptide lists
        self.all_hit_proteins: set[str] = set()

        # Path to the subset FASTA written for DIAMOND input
        self.subset_fasta_path: Path | None = None

        # Reference snapshot directory
        self.ref_snapshot_dir: Path | None = None
        self.ref_metadata: dict[str, str] = {}

    def _work_dir(self) -> Path:
        """Directory for intermediate files (see ``PipelineConfig.work_dir``)."""
        if self.config.work_dir is not None:
            return self.config.work_dir
        if self.config.job_dir is not None:
            return self.config.job_dir / "work"
        return self.config.output_dir / "work"

    def _update_progress(
        self, stage: str, list_id: str = "", progress_done: int | None = None
    ) -> None:
        """Update and report progress.

        Also logs how long the previous stage took and the process's peak
        resident memory, so a stalled job can be located from the log alone.
        """
        now = time.monotonic()
        previous_stage = self.progress.current_stage
        if previous_stage and self._stage_started is not None:
            logger.info(
                f"Finished stage '{previous_stage}' in {now - self._stage_started:.1f}s "
                f"(peak rss {format_bytes(peak_rss_bytes())})"
            )
        self._stage_started = now
        self.progress.current_stage = stage
        self.progress.current_list_id = list_id
        if progress_done is not None:
            self.progress.progress_done = progress_done
        if self.progress_callback:
            self.progress_callback(self.progress)
        logger.info(f"Stage: {stage}" + (f" (list: {list_id})" if list_id else ""))

    def run(self) -> PipelineResult:
        """Execute the complete pipeline.

        Returns:
            PipelineResult with success status and results
        """
        peptide_list_results: list[PeptideListResult] = []
        pipeline_started = time.monotonic()

        try:
            logger.info(
                f"Pipeline config: {len(self.config.peptide_list_paths)} peptide list(s), "
                f"fasta={self.config.fasta_path}, "
                f"annotated_db={self.config.annotated_db_path}, "
                f"annotations_db={self.config.annotations_db_path}, "
                f"threads={self.config.threads}, "
                f"diamond_block_size={self.config.diamond_block_size}, "
                f"diamond_index_chunks={self.config.diamond_index_chunks}, "
                f"diamond_tmpdir={self.config.diamond_tmpdir}, "
                f"diamond_max_target_seqs={self.config.effective_max_target_seqs()}, "
                f"filters={self.config.filter_policy.to_dict()}, "
                f"go_edge_types={sorted(self.config.go_edge_types)}, "
                f"go_include_self={self.config.go_include_self}, "
                f"output_dir={self.config.output_dir}"
            )

            # Stage 0: Initialize (load FASTA, reference data)
            self._update_progress("Initializing", progress_done=0)
            self._initialize()

            self.progress.total_peptide_lists = len(self.config.peptide_list_paths)

            # Stage 1: Parse all peptide lists
            self._update_progress("Parsing peptide lists", progress_done=_PROGRESS_INITIALIZING)
            for i, peptide_list_path in enumerate(self.config.peptide_list_paths):
                list_id = f"list_{i:03d}"
                peptides = parse_peptide_list(peptide_list_path)
                self.parsed_peptide_lists[list_id] = peptides
                logger.info(f"Parsed {len(peptides)} peptides from {peptide_list_path.name}")

            # Stage 2: Match all peptides against background proteome
            self._update_progress(
                "Matching peptides to background proteome", progress_done=_PROGRESS_PARSING
            )
            self._match_all_peptides()

            # Stage 3: Write subset FASTA of hit proteins (for future DIAMOND search)
            self._update_progress("Writing subset FASTA", progress_done=_PROGRESS_MATCHING)
            self._write_subset_fasta()

            # Stage 4: Homology search (placeholder for DIAMOND)
            if not self.config.mock_hits_path:
                self._run_homology_search()
            else:
                self.progress.progress_done = _PROGRESS_FILTERING

            # Stage 4b: Load subject annotations from companion DB
            if not self.config.mock_subject_annotations_path:
                self._load_subject_annotations()
            else:
                self.progress.progress_done = _PROGRESS_SUBJECT_ANNOTATIONS

            # Stage 5+: Per-list annotation, aggregation, reporting
            n_lists = self.progress.total_peptide_lists
            per_list_budget = _PROGRESS_PER_LIST_END - _PROGRESS_PER_LIST_START
            for i, peptide_list_path in enumerate(self.config.peptide_list_paths):
                list_id = f"list_{i:03d}"
                list_start = _PROGRESS_PER_LIST_START + (per_list_budget * i) // max(n_lists, 1)
                self._update_progress(
                    f"Processing peptide list {i + 1}/{n_lists}",
                    list_id,
                    progress_done=list_start,
                )
                result = self._process_peptide_list(peptide_list_path, list_id)
                peptide_list_results.append(result)
                self.progress.completed_peptide_lists = i + 1
                list_done = _PROGRESS_PER_LIST_START + (per_list_budget * (i + 1)) // max(
                    n_lists, 1
                )
                self._update_progress(
                    f"Completed {i + 1}/{n_lists} peptide lists",
                    progress_done=list_done,
                )

            self._update_progress("Pipeline completed", progress_done=_PROGRESS_TOTAL)
            logger.info(
                f"Pipeline finished in {time.monotonic() - pipeline_started:.0f}s "
                f"(peak rss {format_bytes(peak_rss_bytes())})"
            )

            return PipelineResult(
                success=True,
                peptide_list_results=peptide_list_results,
                output_dir=self.config.output_dir,
            )

        except Exception as e:
            logger.exception(
                f"Pipeline failed after {time.monotonic() - pipeline_started:.0f}s "
                f"in stage '{self.progress.current_stage}'"
            )
            return PipelineResult(
                success=False,
                error_message=str(e),
            )

    def _initialize(self) -> None:
        """Stage 0: Initialize job - prepare directories, load reference data."""
        # Create output directory structure
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Create reference snapshot if job_dir is specified
        if self.config.job_dir:
            self._create_reference_snapshot()

        # Load background FASTA
        fasta_size = self.config.fasta_path.stat().st_size if self.config.fasta_path.exists() else 0
        logger.info(
            f"Loading background FASTA: {self.config.fasta_path} ({format_bytes(fasta_size)})"
        )
        records = parse_fasta(self.config.fasta_path)
        self.proteins = build_protein_dict(records)
        logger.info(f"Loaded {len(self.proteins)} background proteins")

        # Load reference data
        self._load_reference_data()

        # Load mock hits if in mock mode
        if self.config.mock_hits_path:
            self._load_mock_hits()

    def _create_reference_snapshot(self) -> None:
        """Create per-job snapshot of reference data."""
        if not self.config.job_dir:
            return

        self.ref_snapshot_dir = self._work_dir() / "ref_snapshot"

        logger.info(f"Creating reference snapshot in {self.ref_snapshot_dir}")

        # Get bundled reference directory
        bundled_ref = get_bundled_reference_dir()

        if bundled_ref.exists():
            # Create snapshot from bundled reference
            create_reference_snapshot(
                bundled_ref,
                self.ref_snapshot_dir,
                use_hardlinks=True,
            )

            # Load metadata
            self.ref_metadata = get_reference_metadata(bundled_ref)
            logger.info(f"Created reference snapshot with metadata: {self.ref_metadata}")
        else:
            logger.warning(f"Bundled reference directory not found: {bundled_ref}")

    def _load_reference_data(self) -> None:
        """Load GO and taxonomy reference data."""
        # Determine source paths
        taxonomy_source: Path | None
        go_source: Path | None
        if self.ref_snapshot_dir and self.ref_snapshot_dir.exists():
            # Load from snapshot
            taxonomy_source = self.ref_snapshot_dir / "taxonomy"
            go_source = self.ref_snapshot_dir / "go" / "go.obo"
            logger.info("Loading reference data from job snapshot")
        elif self.config.taxonomy_data_path or self.config.go_data_path:
            # Load from specified paths
            taxonomy_source = self.config.taxonomy_data_path
            go_source = self.config.go_data_path
            logger.info("Loading reference data from specified paths")
        else:
            # Load from bundled reference
            bundled_ref = get_bundled_reference_dir()
            taxonomy_source = bundled_ref / "taxonomy"
            go_source = bundled_ref / "go" / "go.obo"
            logger.info("Loading reference data from bundled sources")

        # Load taxonomy
        if taxonomy_source and taxonomy_source.exists():
            logger.info(f"Loading taxonomy: {taxonomy_source}")
            self.taxonomy_tree = load_taxonomy_data(taxonomy_source)
            logger.info(f"Loaded {len(self.taxonomy_tree.nodes)} taxonomy nodes")
        else:
            logger.warning("No taxonomy data available")

        # Load GO
        if go_source and go_source.exists():
            logger.info(f"Loading GO: {go_source}")
            self.go_dag = load_go_data(go_source)
            logger.info(f"Loaded {len(self.go_dag.terms)} GO terms")
        else:
            logger.warning("No GO data available")

        # Load subject annotations
        if self.config.mock_subject_annotations_path:
            logger.info(f"Loading subject annotations: {self.config.mock_subject_annotations_path}")
            with open(self.config.mock_subject_annotations_path) as f:
                data = json.load(f)
            self.subject_annotations = load_subject_annotations_from_dict(data)
            logger.info(f"Loaded {len(self.subject_annotations)} subject annotations")

    def _load_subject_annotations(self) -> None:
        """Load subject annotations from the companion SQLite database.

        Collects all unique subject IDs from the homology search results
        and looks up their taxonomy IDs and GO terms.
        """
        if not self.config.annotations_db_path:
            raise ValueError(
                "No annotations database path configured. Each annotated database "
                "must have a companion .annotations.db file. Build one with: "
                "metagomics2-build-annotations"
            )

        if not self.config.annotations_db_path.exists():
            raise FileNotFoundError(
                f"Annotations database not found: {self.config.annotations_db_path}. "
                f"Build it with: metagomics2-build-annotations"
            )

        if not self.protein_to_subjects:
            logger.warning("No homology search results, skipping subject annotation lookup")
            return

        self._update_progress("Loading subject annotations", progress_done=_PROGRESS_FILTERING)

        # Collect all unique subject IDs across all background proteins
        all_subject_ids: set[str] = set()
        for subjects in self.protein_to_subjects.values():
            all_subject_ids |= subjects

        if not all_subject_ids:
            logger.warning("No subject IDs found in homology results")
            return

        self.subject_annotations = load_subject_annotations(
            self.config.annotations_db_path,
            all_subject_ids,
        )
        logger.info(f"Loaded annotations for {len(self.subject_annotations)} subjects")

    def _load_mock_hits(self) -> None:
        """Load mock hits for testing."""
        logger.info(f"Loading mock hits: {self.config.mock_hits_path}")
        with open(self.config.mock_hits_path) as f:  # type: ignore
            data = json.load(f)

        # Convert to expected format
        self.protein_to_subjects = {
            k: set(v) for k, v in data.get("background_to_subjects", {}).items()
        }
        logger.info(f"Loaded mock hits for {len(self.protein_to_subjects)} proteins")

    def _match_all_peptides(self) -> None:
        """Match peptides from all lists against the background proteome.

        Builds per-list match results and the union of all hit proteins.
        """
        all_peptide_sequences: set[str] = set()
        per_list_sequences: dict[str, set[str]] = {}

        for list_id, peptides in self.parsed_peptide_lists.items():
            seqs = {p.sequence for p in peptides}
            per_list_sequences[list_id] = seqs
            all_peptide_sequences |= seqs

        # Single Aho-Corasick pass over the background proteome
        combined_result = match_peptides_to_proteins(all_peptide_sequences, self.proteins)

        # Partition into per-list results
        for list_id, seqs in per_list_sequences.items():
            list_match = MatchResult(
                peptide_to_proteins={
                    pep: combined_result.peptide_to_proteins.get(pep, set())
                    for pep in seqs
                },
                matched_peptides=combined_result.matched_peptides & seqs,
                unmatched_peptides=combined_result.unmatched_peptides & seqs,
                hit_proteins={
                    pid
                    for pep in seqs
                    for pid in combined_result.peptide_to_proteins.get(pep, set())
                },
            )
            self.per_list_match_results[list_id] = list_match
            logger.info(
                f"{list_id}: {list_match.n_matched} matched, "
                f"{list_match.n_unmatched} unmatched"
            )

        self.all_hit_proteins = combined_result.hit_proteins
        logger.info(
            f"Total unique hit proteins across all lists: {len(self.all_hit_proteins)}"
        )

    def _write_subset_fasta(self) -> None:
        """Write a FASTA file containing only proteins hit by at least one peptide.

        This subset FASTA is the input for the DIAMOND homology search against
        an annotated database (e.g. UniProt).
        """
        if not self.all_hit_proteins:
            logger.warning("No proteins matched any peptides, skipping subset FASTA")
            return

        self.subset_fasta_path = self._work_dir() / "hit_proteins.fasta"
        n_written = write_subset_fasta(
            self.proteins, self.all_hit_proteins, self.subset_fasta_path
        )
        logger.info(
            f"Wrote {n_written} hit proteins to subset FASTA: {self.subset_fasta_path}"
        )

    def _run_homology_search(self) -> None:
        """Run DIAMOND blastp and filter the results.

        Searches the subset FASTA (hit proteins from the background proteome)
        against an annotated database. Populates self.protein_to_subjects:
            background_protein_id -> set of annotated DB accessions
        """
        self._update_progress("Homology search", progress_done=_PROGRESS_SUBSET_FASTA)

        if not self.config.annotated_db_path:
            raise ValueError(
                "No annotated database specified. An annotated database is required "
                "for homology search."
            )

        if not self.subset_fasta_path or not self.subset_fasta_path.exists():
            logger.warning("No subset FASTA available, skipping homology search")
            return

        work_dir = self.subset_fasta_path.parent
        diamond_output = work_dir / "diamond_results.tsv"
        # Keep DIAMOND's console log with the job's logs (which survive cleanup)
        # in web mode; next to the work files in CLI mode.
        log_dir = self.config.job_dir / "logs" if self.config.job_dir else work_dir

        # DIAMOND applies its own per-query hit cap (--max-target-seqs) before
        # the pipeline's tie-aware top_k filter runs, and that cap is not
        # tie-aware.  Pass it explicitly, never below top_k, and report
        # queries that hit it so truncation is visible in the log.
        max_target_seqs = self.config.effective_max_target_seqs()
        diamond_result = run_diamond(
            query_fasta=self.subset_fasta_path,
            db_path=self.config.annotated_db_path,
            output_path=diamond_output,
            evalue=self.config.effective_diamond_evalue(),
            max_target_seqs=max_target_seqs,
            threads=self.config.threads,
            log_path=log_dir / "diamond.log",
            block_size=self.config.diamond_block_size,
            index_chunks=self.config.diamond_index_chunks,
            tmpdir=self.config.diamond_tmpdir,
        )
        self.diamond_command = " ".join(diamond_result.command)

        logger.info(
            f"DIAMOND returned {diamond_result.n_hits} hits "
            f"for {diamond_result.n_queries} query proteins"
        )
        n_capped = count_queries_at_cap(diamond_result.hits_by_query, max_target_seqs)
        if n_capped:
            logger.warning(
                f"{n_capped} of {diamond_result.n_queries} query proteins returned the "
                f"DIAMOND per-query hit cap of {max_target_seqs}; further hits, including "
                "any tied with the last one kept, were discarded by DIAMOND. Raise "
                "METAGOMICS_DIAMOND_MAX_TARGET_SEQS (or --diamond-max-target-seqs; 0 = "
                "unlimited) if those hits matter."
            )

        # Apply filter policy (pident, evalue thresholds, top_k ranking, etc.)
        self._update_progress("Filtering homology hits", progress_done=_PROGRESS_HOMOLOGY)
        self.protein_to_subjects = filter_all_hits(
            diamond_result.hits_by_query, self.config.filter_policy
        )
        self.protein_to_subject_hits = filter_all_hits_with_hits(
            diamond_result.hits_by_query, self.config.filter_policy
        )

        n_with_hits = sum(1 for s in self.protein_to_subjects.values() if s)
        logger.info(
            f"After filtering: {n_with_hits} background proteins "
            f"have at least one accepted subject hit"
        )

    def _process_peptide_list(
        self,
        peptide_list_path: Path,
        list_id: str,
    ) -> PeptideListResult:
        """Process a single peptide list through the pipeline.

        Peptide parsing and matching have already been done in the shared
        stages.  This method handles annotation, aggregation, and reporting.
        """
        peptides = self.parsed_peptide_lists[list_id]
        match_result = self.per_list_match_results[list_id]

        # Annotate peptides
        self._update_progress("Annotating peptides", list_id)
        annotations = self._annotate_peptides(peptides, match_result.peptide_to_proteins)

        # Stage 7: Aggregate
        self._update_progress("Aggregating results", list_id)
        aggregation = aggregate_peptide_annotations(annotations)

        # Validate invariants
        violations = validate_aggregation_invariants(aggregation)
        if violations:
            logger.warning(f"Aggregation invariant violations: {violations}")

        # Stage 8: Write reports
        self._update_progress("Writing reports", list_id)
        self._write_reports(list_id, peptide_list_path, aggregation, annotations)

        return PeptideListResult(
            list_id=list_id,
            filename=peptide_list_path.name,
            n_peptides=len(peptides),
            n_matched=match_result.n_matched,
            n_unmatched=match_result.n_unmatched,
            aggregation=aggregation,
        )

    def _annotate_peptides(
        self,
        peptides: list[Peptide],
        peptide_to_proteins: dict[str, set[str]],
    ) -> list[PeptideAnnotation]:
        """Annotate all peptides."""
        annotations: list[PeptideAnnotation] = []

        for peptide in peptides:
            ann = annotate_peptide(
                peptide=peptide.sequence,
                quantity=peptide.quantity,
                peptide_to_proteins=peptide_to_proteins,
                protein_to_subjects=self.protein_to_subjects,
                subject_annotations=self.subject_annotations,
                taxonomy_tree=self.taxonomy_tree or TaxonomyTree(),
                go_dag=self.go_dag or GODAG(),
                go_edge_types=self.config.go_edge_types,
                go_include_self=self.config.go_include_self,
            )
            annotations.append(ann)

        # Peptides with homology hits but nothing usable behind them are
        # counted as unannotated; say so, because it usually means the
        # annotations database or taxonomy dump is stale or mismatched.
        n_silent = sum(1 for a in annotations if a.implied_subjects and not a.is_annotated)
        if n_silent:
            logger.warning(
                f"{n_silent} of {len(annotations)} peptides matched homology hits but "
                "received no taxonomy or GO annotation (hit accessions missing from the "
                "annotations database, or taxonomy IDs unknown to the taxonomy tree); "
                "they are counted as unannotated"
            )

        return annotations

    def _write_reports(
        self,
        list_id: str,
        peptide_list_path: Path,
        aggregation: AggregationResult,
        annotations: list[PeptideAnnotation] | None = None,
    ) -> None:
        """Write output reports for a peptide list."""
        # Create output directory for this list
        list_output_dir = self.config.output_dir / list_id
        list_output_dir.mkdir(parents=True, exist_ok=True)

        # Write taxonomy nodes CSV
        if self.taxonomy_tree:
            write_taxonomy_nodes_csv(
                aggregation,
                self.taxonomy_tree,
                list_output_dir / "taxonomy_nodes.csv",
            )

        # Write GO terms CSV
        if self.go_dag:
            write_go_terms_csv(
                aggregation,
                self.go_dag,
                list_output_dir / "go_terms.csv",
                edge_types=self.config.go_edge_types,
            )

        # Write GO-taxonomy combo CSV
        if self.taxonomy_tree and self.go_dag and annotations:
            combos = aggregate_go_taxonomy_combos(annotations, aggregation)
            write_go_taxonomy_combo_csv(
                combos,
                self.taxonomy_tree,
                self.go_dag,
                list_output_dir / "go_taxonomy_combo.csv",
                edge_types=self.config.go_edge_types,
            )

        # Write coverage CSV
        write_coverage_csv(
            aggregation.coverage,
            list_output_dir / "coverage.csv",
        )

        # Write peptide mapping Parquet
        if annotations:
            write_peptide_mapping_parquet(
                annotations,
                self.per_list_match_results[list_id].peptide_to_proteins,
                self.protein_to_subjects,
                list_output_dir / "peptide_mapping.parquet",
                self.protein_to_subject_hits,
            )

        # Write manifest
        manifest = create_manifest(
            metagomics2_version=__version__,
            search_tool=self.config.search_tool,
            search_tool_command=self.diamond_command,
            annotated_db_choice=str(self.config.annotated_db_path or "mock"),
            input_fasta_path=self.config.fasta_path,
            peptide_list_path=peptide_list_path,
            parameters={
                **self.config.filter_policy.to_dict(),
                "go_edge_types": sorted(self.config.go_edge_types),
                "go_include_self": self.config.go_include_self,
                "threads": self.config.threads,
                "diamond_block_size": self.config.diamond_block_size,
                "diamond_index_chunks": self.config.diamond_index_chunks,
                "diamond_tmpdir": (
                    str(self.config.diamond_tmpdir) if self.config.diamond_tmpdir else None
                ),
                "diamond_max_target_seqs": self.config.effective_max_target_seqs(),
                # What DIAMOND was actually run with; null when the homology
                # stage was mocked and DIAMOND did not run.
                "diamond_evalue": (
                    None if self.config.mock_hits_path else self.config.effective_diamond_evalue()
                ),
            },
            go_snapshot_dir=self.ref_snapshot_dir / "go" if self.ref_snapshot_dir else None,
            taxonomy_snapshot_dir=(
                self.ref_snapshot_dir / "taxonomy" if self.ref_snapshot_dir else None
            ),
            annotated_db_path=self.config.annotated_db_path,
        )

        # Add reference metadata to manifest
        if self.ref_metadata:
            manifest.parameters["reference_metadata"] = self.ref_metadata

        write_manifest_json(manifest, list_output_dir / "run_manifest.json")

        logger.info(f"Wrote reports to {list_output_dir}")


def run_pipeline(
    config: PipelineConfig,
    progress_callback: ProgressCallback | None = None,
) -> PipelineResult:
    """Convenience function to run the pipeline.

    Args:
        config: Pipeline configuration
        progress_callback: Optional callback for progress updates

    Returns:
        PipelineResult
    """
    runner = PipelineRunner(config, progress_callback)
    return runner.run()
