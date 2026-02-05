"""Command-line interface for Metagomics 2."""

import argparse
import json
import logging
import sys
from pathlib import Path

from metagomics2 import __version__
from metagomics2.core.filtering import FilterPolicy
from metagomics2.pipeline.runner import PipelineConfig, PipelineProgress, run_pipeline

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def progress_callback(progress: PipelineProgress) -> None:
    """Print progress updates to stderr."""
    msg = f"[{progress.completed_peptide_lists}/{progress.total_peptide_lists}] {progress.current_stage}"
    if progress.current_list_id:
        msg += f" ({progress.current_list_id})"
    print(msg, file=sys.stderr)


def parse_filter_params(args: argparse.Namespace) -> FilterPolicy:
    """Parse filter parameters from CLI args or params file."""
    if args.params:
        with open(args.params) as f:
            params = json.load(f)
        return FilterPolicy(
            max_evalue=params.get("max_evalue"),
            min_pident=params.get("min_pident"),
            min_qcov=params.get("min_qcov"),
            min_alnlen=params.get("min_alnlen"),
            top_k=params.get("top_k"),
            delta_bitscore=params.get("delta_bitscore"),
            best_hit_only=params.get("best_hit_only", False),
        )

    return FilterPolicy(
        max_evalue=args.max_evalue,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        min_alnlen=args.min_alnlen,
        top_k=args.top_k,
        delta_bitscore=args.delta_bitscore,
        best_hit_only=args.best_hit_only,
    )


def cmd_run(args: argparse.Namespace) -> int:
    """Execute the 'run' command."""
    # Validate inputs
    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        print(f"Error: FASTA file not found: {fasta_path}", file=sys.stderr)
        return 1

    peptide_paths = [Path(p) for p in args.peptides]
    for p in peptide_paths:
        if not p.exists():
            print(f"Error: Peptide file not found: {p}", file=sys.stderr)
            return 1

    output_dir = Path(args.outdir)

    # Parse filter policy
    filter_policy = parse_filter_params(args)

    # Build config
    config = PipelineConfig(
        fasta_path=fasta_path,
        peptide_list_paths=peptide_paths,
        output_dir=output_dir,
        search_tool=args.search_tool,
        annotated_db_path=Path(args.db) if args.db else None,
        threads=args.threads,
        filter_policy=filter_policy,
        go_data_path=Path(args.go) if args.go else None,
        taxonomy_data_path=Path(args.taxonomy) if args.taxonomy else None,
        go_edge_types=set(args.go_edge_types.split(",")) if args.go_edge_types else {"is_a"},
        go_include_self=not args.go_exclude_self,
        mock_hits_path=Path(args.mock_hits) if args.mock_hits else None,
        mock_subject_annotations_path=Path(args.mock_annotations) if args.mock_annotations else None,
    )

    # Run pipeline
    result = run_pipeline(config, progress_callback if not args.quiet else None)

    if result.success:
        print(f"Pipeline completed successfully. Output: {output_dir}", file=sys.stderr)
        return 0
    else:
        print(f"Pipeline failed: {result.error_message}", file=sys.stderr)
        return 1


def cmd_version(args: argparse.Namespace) -> int:
    """Execute the 'version' command."""
    print(f"metagomics2 {__version__}")
    return 0


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        prog="metagomics2",
        description="Metaproteomics annotation and aggregation tool",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # 'run' command
    run_parser = subparsers.add_parser(
        "run",
        help="Run the annotation pipeline",
        description="Run the metagomics annotation pipeline on peptide data",
    )

    # Required arguments
    run_parser.add_argument(
        "--fasta",
        required=True,
        help="Path to background proteome FASTA file",
    )
    run_parser.add_argument(
        "--peptides",
        required=True,
        action="append",
        help="Path to peptide list file (CSV/TSV). Can be specified multiple times.",
    )
    run_parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for results",
    )

    # Search tool options
    run_parser.add_argument(
        "--search-tool",
        choices=["diamond", "blast"],
        default="diamond",
        help="Homology search tool (default: diamond)",
    )
    run_parser.add_argument(
        "--db",
        help="Path to indexed annotated database",
    )
    run_parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for homology search (default: 1)",
    )

    # Filter parameters
    filter_group = run_parser.add_argument_group("Filter parameters")
    filter_group.add_argument(
        "--max-evalue",
        type=float,
        help="Maximum e-value threshold",
    )
    filter_group.add_argument(
        "--min-pident",
        type=float,
        help="Minimum percent identity threshold",
    )
    filter_group.add_argument(
        "--min-qcov",
        type=float,
        help="Minimum query coverage threshold",
    )
    filter_group.add_argument(
        "--min-alnlen",
        type=int,
        help="Minimum alignment length threshold",
    )
    filter_group.add_argument(
        "--top-k",
        type=int,
        help="Keep only top K hits by bitscore",
    )
    filter_group.add_argument(
        "--delta-bitscore",
        type=float,
        help="Keep hits within delta of best bitscore",
    )
    filter_group.add_argument(
        "--best-hit-only",
        action="store_true",
        help="Keep only the single best hit per query",
    )
    filter_group.add_argument(
        "--params",
        help="Path to JSON file with filter parameters",
    )

    # Reference data options
    ref_group = run_parser.add_argument_group("Reference data")
    ref_group.add_argument(
        "--go",
        help="Path to GO data file (JSON format)",
    )
    ref_group.add_argument(
        "--taxonomy",
        help="Path to taxonomy data file (JSON format)",
    )
    ref_group.add_argument(
        "--go-edge-types",
        default="is_a",
        help="Comma-separated GO edge types for closure (default: is_a)",
    )
    ref_group.add_argument(
        "--go-exclude-self",
        action="store_true",
        help="Exclude GO terms themselves from closure",
    )

    # Testing/mock options
    test_group = run_parser.add_argument_group("Testing options")
    test_group.add_argument(
        "--mock-hits",
        help="Path to mock hits JSON file (for testing)",
    )
    test_group.add_argument(
        "--mock-annotations",
        help="Path to mock subject annotations JSON file (for testing)",
    )

    # Output options
    run_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )
    run_parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress output",
    )

    run_parser.set_defaults(func=cmd_run)

    # 'version' command (alternative to --version)
    version_parser = subparsers.add_parser(
        "version",
        help="Show version information",
    )
    version_parser.set_defaults(func=cmd_version)

    return parser


def main() -> int:
    """Main entry point for CLI."""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    # Setup logging
    verbose = getattr(args, "verbose", False)
    setup_logging(verbose)

    # Execute command
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
