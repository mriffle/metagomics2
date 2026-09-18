"""Report generation: CSV outputs and manifest."""

import csv
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import polars as pl

from metagomics2.core.aggregation import AggregationResult, ComboAggregate, CoverageStats
from metagomics2.core.annotation import PeptideAnnotation
from metagomics2.core.go import GODAG
from metagomics2.core.taxonomy import TaxonomyTree


@dataclass
class ManifestInfo:
    """Information for the run manifest."""

    metagomics2_version: str = ""
    git_sha: str | None = None
    python_version: str = ""
    search_tool: str = ""
    search_tool_version: str = ""
    search_tool_command: str = ""
    annotated_db_choice: str = ""
    annotated_db_version: str = ""
    annotated_db_hash: str = ""
    input_fasta_hash: str = ""
    input_fasta_path: str = ""
    peptide_list_hash: str = ""
    peptide_list_path: str = ""
    go_snapshot_files: dict[str, str] = field(default_factory=dict)  # filename -> hash
    taxonomy_snapshot_files: dict[str, str] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)
    timestamp_utc: str = ""


def format_number(value: float) -> str:
    """Format a quantity or ratio for a CSV cell without losing precision.

    Uses Python's shortest round-trip representation (``repr``), so every
    value reads back as exactly the double that was computed: ``18.0``,
    ``0.5``, ``3e-11``.  Fixed-point formatting with a set number of decimals
    would silently round small normalized abundances to zero.
    """
    return repr(float(value))


def write_taxonomy_nodes_csv(
    result: AggregationResult,
    taxonomy_tree: TaxonomyTree,
    output_path: Path,
) -> None:
    """Write taxonomy_nodes.csv.

    Columns: tax_id, name, rank, parent_tax_id, quantity, ratio_total, ratio_annotated, n_peptides

    Args:
        result: Aggregation result
        taxonomy_tree: Taxonomy tree for node metadata
        output_path: Path to write CSV
    """
    # Sort by quantity descending, then by tax_id for determinism
    sorted_nodes = sorted(
        result.taxonomy_nodes.items(),
        key=lambda x: (-x[1].quantity, x[0]),
    )

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "tax_id",
            "name",
            "rank",
            "parent_tax_id",
            "quantity",
            "ratio_total",
            "ratio_annotated",
            "n_peptides",
        ])

        for tax_id, node in sorted_nodes:
            tax_node = taxonomy_tree.nodes.get(tax_id)
            name = tax_node.name if tax_node else ""
            rank = tax_node.rank if tax_node else ""
            parent_tax_id = tax_node.parent_tax_id if tax_node else None

            ratio_annotated_str = (
                format_number(node.ratio_annotated)
                if node.ratio_annotated is not None
                else ""
            )

            writer.writerow([
                tax_id,
                name,
                rank,
                parent_tax_id if parent_tax_id is not None else "",
                format_number(node.quantity),
                format_number(node.ratio_total),
                ratio_annotated_str,
                node.n_peptides,
            ])


def write_go_terms_csv(
    result: AggregationResult,
    go_dag: GODAG,
    output_path: Path,
    parent_delimiter: str = ";",
    edge_types: set[str] | None = None,
) -> None:
    """Write go_terms.csv.

    Columns: go_id, name, namespace, parent_go_ids, quantity, ratio_total,
    ratio_annotated, n_peptides

    Args:
        result: Aggregation result
        go_dag: GO DAG for term metadata
        output_path: Path to write CSV
        parent_delimiter: Delimiter for parent IDs (default: ";")
        edge_types: Edge types to include in parent_go_ids (default: all)
    """
    # Sort by quantity descending, then by go_id for determinism
    sorted_terms = sorted(
        result.go_terms.items(),
        key=lambda x: (-x[1].quantity, x[0]),
    )

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "go_id",
            "name",
            "namespace",
            "parent_go_ids",
            "quantity",
            "ratio_total",
            "ratio_annotated",
            "n_peptides",
        ])

        for go_id, node in sorted_terms:
            go_term = go_dag.terms.get(go_id)
            if go_term:
                name = go_term.name
                namespace = go_term.namespace
            else:
                obsolete = go_dag.obsolete_terms.get(go_id)
                if obsolete:
                    name = f"OBSOLETE: {obsolete.name}"
                    namespace = obsolete.namespace
                else:
                    name = ""
                    namespace = ""

            # Collect parent IDs for the specified edge types
            parent_ids: set[str] = set()
            if go_term:
                if edge_types is not None:
                    for et in edge_types:
                        parent_ids |= go_term.parents.get(et, set())
                else:
                    for parents in go_term.parents.values():
                        parent_ids |= parents
            parent_ids_str = parent_delimiter.join(sorted(parent_ids))

            ratio_annotated_str = (
                format_number(node.ratio_annotated)
                if node.ratio_annotated is not None
                else ""
            )

            writer.writerow([
                go_id,
                name,
                namespace,
                parent_ids_str,
                format_number(node.quantity),
                format_number(node.ratio_total),
                ratio_annotated_str,
                node.n_peptides,
            ])


def write_coverage_csv(
    coverage: CoverageStats,
    output_path: Path,
) -> None:
    """Write coverage.csv.

    Single row with coverage statistics.

    Args:
        coverage: Coverage statistics
        output_path: Path to write CSV
    """
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "total_peptide_quantity",
            "annotated_peptide_quantity",
            "unannotated_peptide_quantity",
            "annotation_coverage_ratio",
            "n_peptides_total",
            "n_peptides_annotated",
            "n_peptides_unannotated",
        ])
        writer.writerow([
            format_number(coverage.total_peptide_quantity),
            format_number(coverage.annotated_peptide_quantity),
            format_number(coverage.unannotated_peptide_quantity),
            format_number(coverage.annotation_coverage_ratio),
            coverage.n_peptides_total,
            coverage.n_peptides_annotated,
            coverage.n_peptides_unannotated,
        ])


def write_go_taxonomy_combo_csv(
    combos: dict[tuple[int, str], ComboAggregate],
    taxonomy_tree: TaxonomyTree,
    go_dag: GODAG,
    output_path: Path,
    parent_delimiter: str = ";",
    edge_types: set[str] | None = None,
) -> None:
    """Write go_taxonomy_combo.csv.

    Cross-tabulation of taxonomy nodes and GO terms with enough
    parent information to rebuild both trees.

    Columns: tax_id, tax_name, tax_rank, parent_tax_id, go_id, go_name,
             go_namespace, parent_go_ids, quantity, fraction_of_taxon,
             fraction_of_go, ratio_total_taxon, ratio_total_go, n_peptides

    Args:
        combos: Dictionary mapping (tax_id, go_id) to ComboAggregate
        taxonomy_tree: Taxonomy tree for node metadata
        go_dag: GO DAG for term metadata
        output_path: Path to write CSV
        parent_delimiter: Delimiter for parent GO IDs (default: ";")
        edge_types: Edge types to include in parent_go_ids (default: all)
    """
    # Sort by quantity descending, then tax_id, then go_id for determinism
    sorted_combos = sorted(
        combos.values(),
        key=lambda c: (-c.quantity, c.tax_id, c.go_id),
    )

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "tax_id",
            "tax_name",
            "tax_rank",
            "parent_tax_id",
            "go_id",
            "go_name",
            "go_namespace",
            "parent_go_ids",
            "quantity",
            "fraction_of_taxon",
            "fraction_of_go",
            "ratio_total_taxon",
            "ratio_total_go",
            "n_peptides",
        ])

        for combo in sorted_combos:
            # Taxonomy metadata
            tax_node = taxonomy_tree.nodes.get(combo.tax_id)
            tax_name = tax_node.name if tax_node else ""
            tax_rank = tax_node.rank if tax_node else ""
            parent_tax_id = tax_node.parent_tax_id if tax_node else None

            # GO metadata
            go_term = go_dag.terms.get(combo.go_id)
            if go_term:
                go_name = go_term.name
                go_namespace = go_term.namespace
            else:
                obsolete = go_dag.obsolete_terms.get(combo.go_id)
                if obsolete:
                    go_name = f"OBSOLETE: {obsolete.name}"
                    go_namespace = obsolete.namespace
                else:
                    go_name = ""
                    go_namespace = ""

            # Collect parent GO IDs
            parent_go_ids: set[str] = set()
            if go_term:
                if edge_types is not None:
                    for et in edge_types:
                        parent_go_ids |= go_term.parents.get(et, set())
                else:
                    for parents in go_term.parents.values():
                        parent_go_ids |= parents
            parent_go_ids_str = parent_delimiter.join(sorted(parent_go_ids))

            writer.writerow([
                combo.tax_id,
                tax_name,
                tax_rank,
                parent_tax_id if parent_tax_id is not None else "",
                combo.go_id,
                go_name,
                go_namespace,
                parent_go_ids_str,
                format_number(combo.quantity),
                format_number(combo.fraction_of_taxon),
                format_number(combo.fraction_of_go),
                format_number(combo.ratio_total_taxon),
                format_number(combo.ratio_total_go),
                combo.n_peptides,
            ])


def write_peptide_mapping_parquet(
    annotations: list[PeptideAnnotation],
    peptide_to_proteins: dict[str, set[str]],
    protein_to_subjects: dict[str, set[str]],
    output_path: Path,
    protein_to_subject_hits: "dict[str, dict[str, Any]] | None" = None,
) -> None:
    """Write peptide_mapping.parquet.

    Each row represents one (peptide, background_protein, annotated_protein) triple
    for annotated peptides only.

    Schema:
        peptide: Utf8
        peptide_lca_tax_ids: List[Int64]
        peptide_go_terms: List[Utf8]
        background_protein: Utf8
        annotated_protein: Utf8
        evalue: Float64 (nullable)
        pident: Float64 (nullable)

    Args:
        annotations: List of PeptideAnnotation objects
        peptide_to_proteins: Mapping from peptide sequence to background protein IDs
        protein_to_subjects: Mapping from background protein ID to annotated subject IDs
        output_path: Path to write the Parquet file
        protein_to_subject_hits: Optional mapping from background protein ID to
            {subject_id: HomologyHit} for evalue/pident lookup
    """
    rows_peptide: list[str] = []
    rows_lca_tax_ids: list[list[int]] = []
    rows_go_terms: list[list[str]] = []
    rows_bg_protein: list[str] = []
    rows_ann_protein: list[str] = []
    rows_evalue: list[float | None] = []
    rows_pident: list[float | None] = []

    for ann in annotations:
        if not ann.is_annotated:
            continue

        go_terms_sorted = sorted(ann.go_terms)
        lca_tax_ids = sorted(ann.taxonomy_nodes)

        bg_proteins = peptide_to_proteins.get(ann.peptide, set())
        for bg_protein in sorted(bg_proteins):
            subjects = protein_to_subjects.get(bg_protein, set())
            subject_hits = (protein_to_subject_hits or {}).get(bg_protein, {})
            for subject in sorted(subjects):
                hit = subject_hits.get(subject)
                rows_peptide.append(ann.peptide)
                rows_lca_tax_ids.append(lca_tax_ids)
                rows_go_terms.append(go_terms_sorted)
                rows_bg_protein.append(bg_protein)
                rows_ann_protein.append(subject)
                rows_evalue.append(hit.evalue if hit is not None else None)
                rows_pident.append(hit.pident if hit is not None else None)

    schema: dict[str, pl.DataType] = {
        "peptide": pl.Utf8(),
        "peptide_lca_tax_ids": pl.List(pl.Int64),
        "peptide_go_terms": pl.List(pl.Utf8),
        "background_protein": pl.Utf8(),
        "annotated_protein": pl.Utf8(),
        "evalue": pl.Float64(),
        "pident": pl.Float64(),
    }

    if rows_peptide:
        df = pl.DataFrame(
            {
                "peptide": rows_peptide,
                "peptide_lca_tax_ids": rows_lca_tax_ids,
                "peptide_go_terms": rows_go_terms,
                "background_protein": rows_bg_protein,
                "annotated_protein": rows_ann_protein,
                "evalue": rows_evalue,
                "pident": rows_pident,
            },
            schema=schema,
        )
    else:
        df = pl.DataFrame(schema=schema)

    df.write_parquet(output_path)


def get_tool_version(tool: str) -> str:
    """Get version string for a tool.

    Args:
        tool: Tool name (diamond, blastp, etc.)

    Returns:
        Version string or empty if not found
    """
    try:
        if tool == "diamond":
            result = subprocess.run(
                ["diamond", "version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            # Output like "diamond version 2.1.8"
            return result.stdout.strip()
        elif tool in ("blastp", "blast"):
            result = subprocess.run(
                ["blastp", "-version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            # First line contains version
            return result.stdout.strip().split("\n")[0]
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        pass
    return ""


def get_git_sha() -> str | None:
    """Get current git SHA if in a git repository.

    Returns:
        Git SHA or None
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        pass
    return None


def compute_file_hash(file_path: Path) -> str:
    """Compute SHA256 hash of a file.

    Args:
        file_path: Path to file

    Returns:
        Hex-encoded SHA256 hash
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def write_manifest_json(
    manifest: ManifestInfo,
    output_path: Path,
) -> None:
    """Write run_manifest.json.

    Args:
        manifest: Manifest information
        output_path: Path to write JSON
    """
    data = {
        "metagomics2_version": manifest.metagomics2_version,
        "git_sha": manifest.git_sha,
        "python_version": manifest.python_version,
        "search_tool": manifest.search_tool,
        "search_tool_version": manifest.search_tool_version,
        "search_tool_command": manifest.search_tool_command,
        "annotated_db": {
            "choice": manifest.annotated_db_choice,
            "version": manifest.annotated_db_version,
            "hash": manifest.annotated_db_hash,
        },
        "inputs": {
            "fasta": {
                "path": manifest.input_fasta_path,
                "sha256": manifest.input_fasta_hash,
            },
            "peptide_list": {
                "path": manifest.peptide_list_path,
                "sha256": manifest.peptide_list_hash,
            },
        },
        "reference_snapshots": {
            "go_files": manifest.go_snapshot_files,
            "taxonomy_files": manifest.taxonomy_snapshot_files,
        },
        "parameters": manifest.parameters,
        "timestamp_utc": manifest.timestamp_utc,
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def create_manifest(
    metagomics2_version: str,
    search_tool: str,
    search_tool_command: str,
    annotated_db_choice: str,
    input_fasta_path: Path,
    peptide_list_path: Path,
    parameters: dict[str, Any],
    go_snapshot_dir: Path | None = None,
    taxonomy_snapshot_dir: Path | None = None,
    annotated_db_path: Path | None = None,
    annotated_db_hash: str | None = None,
) -> ManifestInfo:
    """Create a manifest with all provenance information.

    Args:
        metagomics2_version: Version of metagomics2
        search_tool: Search tool used (diamond/blast)
        search_tool_command: Full command line used
        annotated_db_choice: Name/identifier of annotated database
        input_fasta_path: Path to input FASTA
        peptide_list_path: Path to peptide list
        parameters: Dictionary of parameters used
        go_snapshot_dir: Directory containing GO snapshot files
        taxonomy_snapshot_dir: Directory containing taxonomy snapshot files
        annotated_db_path: Path to annotated database file, hashed here
            unless ``annotated_db_hash`` is given
        annotated_db_hash: Precomputed SHA256 of the annotated database.
            The pipeline writes one manifest per peptide list and the
            database can be tens of gigabytes, so the runner hashes it once
            and passes the digest in rather than re-reading the file per list.

    Returns:
        ManifestInfo object
    """
    manifest = ManifestInfo(
        metagomics2_version=metagomics2_version,
        git_sha=get_git_sha(),
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        search_tool=search_tool,
        search_tool_version=get_tool_version(search_tool),
        search_tool_command=search_tool_command,
        annotated_db_choice=annotated_db_choice,
        input_fasta_path=str(input_fasta_path),
        input_fasta_hash=compute_file_hash(input_fasta_path) if input_fasta_path.exists() else "",
        peptide_list_path=str(peptide_list_path),
        peptide_list_hash=(
            compute_file_hash(peptide_list_path) if peptide_list_path.exists() else ""
        ),
        parameters=parameters,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
    )

    # Hash GO snapshot files
    if go_snapshot_dir and go_snapshot_dir.exists():
        for file_path in sorted(go_snapshot_dir.iterdir()):
            if file_path.is_file():
                manifest.go_snapshot_files[file_path.name] = compute_file_hash(file_path)

    # Hash taxonomy snapshot files
    if taxonomy_snapshot_dir and taxonomy_snapshot_dir.exists():
        for file_path in sorted(taxonomy_snapshot_dir.iterdir()):
            if file_path.is_file():
                manifest.taxonomy_snapshot_files[file_path.name] = compute_file_hash(file_path)

    # Hash annotated database (or use the digest the caller already computed)
    if annotated_db_hash is not None:
        manifest.annotated_db_hash = annotated_db_hash
    elif annotated_db_path and annotated_db_path.exists():
        manifest.annotated_db_hash = compute_file_hash(annotated_db_path)

    return manifest
