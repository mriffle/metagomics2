"""Pipeline runner orchestration.

This module contains the core pipeline logic that is shared between
the CLI and web server execution modes.
"""

import json
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from metagomics2 import __version__
from metagomics2.core.aggregation import (
    AggregationResult,
    aggregate_peptide_annotations,
    validate_aggregation_invariants,
)
from metagomics2.core.annotation import (
    PeptideAnnotation,
    SubjectAnnotation,
    annotate_peptide,
    load_subject_annotations_from_dict,
)
from metagomics2.core.fasta import build_protein_dict, compute_file_sha256, parse_fasta
from metagomics2.core.filtering import FilterPolicy, filter_all_hits, parse_blast_tabular
from metagomics2.core.go import GODAG, load_go_from_dict, load_go_from_json
from metagomics2.core.matching import match_peptides_to_proteins
from metagomics2.core.peptides import Peptide, parse_peptide_list
from metagomics2.core.reference_loader import (
    create_reference_snapshot,
    get_bundled_reference_dir,
    get_reference_metadata,
    load_go_data,
    load_taxonomy_data,
)
from metagomics2.core.reporting import (
    ManifestInfo,
    create_manifest,
    write_coverage_csv,
    write_go_terms_csv,
    write_manifest_json,
    write_taxonomy_nodes_csv,
)
from metagomics2.core.taxonomy import TaxonomyTree, load_taxonomy_from_dict, load_taxonomy_from_json

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

    # GO closure settings
    go_edge_types: set[str] = field(default_factory=lambda: {"is_a"})
    go_include_self: bool = True

    # Mock mode for testing (bypasses actual homology search)
    mock_hits_path: Path | None = None
    mock_subject_annotations_path: Path | None = None


@dataclass
class PipelineProgress:
    """Progress tracking for pipeline execution."""

    total_peptide_lists: int = 0
    completed_peptide_lists: int = 0
    current_stage: str = ""
    current_list_id: str = ""


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
        
        # Reference snapshot directory
        self.ref_snapshot_dir: Path | None = None
        self.ref_metadata: dict[str, str] = {}

    def _update_progress(self, stage: str, list_id: str = "") -> None:
        """Update and report progress."""
        self.progress.current_stage = stage
        self.progress.current_list_id = list_id
        if self.progress_callback:
            self.progress_callback(self.progress)
        logger.info(f"Stage: {stage}" + (f" (list: {list_id})" if list_id else ""))

    def run(self) -> PipelineResult:
        """Execute the complete pipeline.

        Returns:
            PipelineResult with success status and results
        """
        peptide_list_results: list[PeptideListResult] = []
        
        try:
            # Stage 0: Initialize
            self._update_progress("Initializing")
            self._initialize()

            # Stage 1-2: Process each peptide list
            self.progress.total_peptide_lists = len(self.config.peptide_list_paths)

            for i, peptide_list_path in enumerate(self.config.peptide_list_paths):
                list_id = f"list_{i:03d}"
                result = self._process_peptide_list(peptide_list_path, list_id)
                peptide_list_results.append(result)
                self.progress.completed_peptide_lists = i + 1
                self._update_progress(
                    f"Completed {i + 1}/{self.progress.total_peptide_lists} peptide lists"
                )

            self._update_progress("Pipeline completed")

            return PipelineResult(
                success=True,
                peptide_list_results=peptide_list_results,
                output_dir=self.config.output_dir,
            )

        except Exception as e:
            logger.exception("Pipeline failed")
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
        logger.info(f"Loading background FASTA: {self.config.fasta_path}")
        records = parse_fasta(self.config.fasta_path)
        self.proteins = build_protein_dict(records)
        logger.info(f"Loaded {len(self.proteins)} background proteins")

        # Load reference data
        self._load_reference_data()

        # Run homology search if not in mock mode
        if self.config.mock_hits_path:
            self._load_mock_hits()
        else:
            self._run_homology_search()

    def _create_reference_snapshot(self) -> None:
        """Create per-job snapshot of reference data."""
        if not self.config.job_dir:
            return
            
        work_dir = self.config.job_dir / "work"
        self.ref_snapshot_dir = work_dir / "ref_snapshot"
        
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

    def _run_homology_search(self) -> None:
        """Run DIAMOND or BLAST homology search."""
        if not self.config.annotated_db_path:
            logger.warning("No annotated database specified, skipping homology search")
            return

        # Get union of all hit proteins across all peptide lists
        # For now, we'll do this per-list in _process_peptide_list
        # This is a placeholder for the full implementation
        logger.info("Homology search would run here (not implemented in mock mode)")

    def _process_peptide_list(
        self,
        peptide_list_path: Path,
        list_id: str,
    ) -> PeptideListResult:
        """Process a single peptide list through the pipeline."""
        self._update_progress("Parsing peptide list", list_id)

        # Stage 1: Parse peptide list
        peptides = parse_peptide_list(peptide_list_path)
        logger.info(f"Parsed {len(peptides)} peptides from {peptide_list_path.name}")

        # Stage 2: Exact matching
        self._update_progress("Matching peptides to proteins", list_id)
        peptide_sequences = {p.sequence for p in peptides}
        match_result = match_peptides_to_proteins(peptide_sequences, self.proteins)
        logger.info(
            f"Matched {match_result.n_matched} peptides, "
            f"{match_result.n_unmatched} unmatched"
        )

        # Stage 3-4: Homology search and filtering
        # (Already done in _initialize for mock mode)

        # Stage 5-6: Annotate peptides
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
        self._write_reports(list_id, peptide_list_path, aggregation)

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

        return annotations

    def _write_reports(
        self,
        list_id: str,
        peptide_list_path: Path,
        aggregation: AggregationResult,
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
            )

        # Write coverage CSV
        write_coverage_csv(
            aggregation.coverage,
            list_output_dir / "coverage.csv",
        )

        # Write manifest
        manifest = create_manifest(
            metagomics2_version=__version__,
            search_tool=self.config.search_tool,
            search_tool_command="",  # TODO: capture actual command
            annotated_db_choice=str(self.config.annotated_db_path or "mock"),
            input_fasta_path=self.config.fasta_path,
            peptide_list_path=peptide_list_path,
            parameters=self.config.filter_policy.to_dict(),
            go_snapshot_dir=self.ref_snapshot_dir / "go" if self.ref_snapshot_dir else None,
            taxonomy_snapshot_dir=self.ref_snapshot_dir / "taxonomy" if self.ref_snapshot_dir else None,
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
