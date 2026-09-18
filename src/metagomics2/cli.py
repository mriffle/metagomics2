"""Command-line interface for Metagomics 2."""

import argparse
import json
import logging
import sys
from pathlib import Path

from metagomics2 import __version__
from metagomics2.core.diamond import DEFAULT_MAX_TARGET_SEQS
from metagomics2.core.filtering import FilterPolicy
from metagomics2.core.go import GO_EDGE_TYPES, parse_go_edge_types
from metagomics2.logging_setup import configure_logging
from metagomics2.pipeline.runner import PipelineConfig, PipelineProgress, run_pipeline

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI."""
    configure_logging("cli", None, logging.DEBUG if verbose else logging.INFO)


def progress_callback(progress: PipelineProgress) -> None:
    """Print progress updates to stderr."""
    total = progress.progress_total
    pct = (progress.progress_done * 100) // total if total else 0
    msg = f"[{pct:3d}%] {progress.current_stage}"
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
        )

    return FilterPolicy(
        max_evalue=args.max_evalue,
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        min_alnlen=args.min_alnlen,
        top_k=args.top_k,
    )


def cmd_run(args: argparse.Namespace) -> int:
    """Execute the 'run' command."""
    # Display version
    print(f"Metagomics 2 v{__version__}", file=sys.stderr)

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

    # Validate database paths (required unless in mock mode)
    db_path = None
    annotations_db_path = None

    if args.db:
        db_path = Path(args.db)
        if not db_path.exists():
            print(f"Error: Annotated database not found: {db_path}", file=sys.stderr)
            return 1
    elif not args.mock_hits:
        print("Error: --db is required (or use --mock-hits for testing)", file=sys.stderr)
        return 1

    if args.annotations_db:
        annotations_db_path = Path(args.annotations_db)
        if not annotations_db_path.exists():
            print(f"Error: Annotations database not found: {annotations_db_path}", file=sys.stderr)
            print("Build it with: python scripts/build_annotations_db.py", file=sys.stderr)
            return 1
    elif not args.mock_annotations:
        print(
            "Error: --annotations-db is required (or use --mock-annotations for testing)",
            file=sys.stderr,
        )
        return 1

    if args.diamond_block_size is not None and args.diamond_block_size <= 0:
        print("Error: --diamond-block-size must be greater than 0", file=sys.stderr)
        return 1
    if args.diamond_index_chunks is not None and args.diamond_index_chunks < 1:
        print("Error: --diamond-index-chunks must be at least 1", file=sys.stderr)
        return 1
    if args.diamond_max_target_seqs is not None and args.diamond_max_target_seqs < 0:
        print(
            "Error: --diamond-max-target-seqs must be 0 (unlimited) or a positive integer",
            file=sys.stderr,
        )
        return 1

    try:
        go_edge_types = parse_go_edge_types(args.go_edge_types)
    except ValueError as e:
        print(f"Error: --go-edge-types: {e}", file=sys.stderr)
        return 1

    # Build config
    config = PipelineConfig(
        fasta_path=fasta_path,
        peptide_list_paths=peptide_paths,
        output_dir=output_dir,
        search_tool=args.search_tool,
        annotated_db_path=db_path,
        annotations_db_path=annotations_db_path,
        threads=args.threads,
        filter_policy=filter_policy,
        go_data_path=Path(args.go) if args.go else None,
        taxonomy_data_path=Path(args.taxonomy) if args.taxonomy else None,
        go_edge_types=go_edge_types,
        go_include_self=not args.go_exclude_self,
        mock_hits_path=Path(args.mock_hits) if args.mock_hits else None,
        mock_subject_annotations_path=(
            Path(args.mock_annotations) if args.mock_annotations else None
        ),
        diamond_block_size=args.diamond_block_size,
        diamond_index_chunks=args.diamond_index_chunks,
        diamond_tmpdir=Path(args.diamond_tmpdir) if args.diamond_tmpdir else None,
        diamond_max_target_seqs=(
            args.diamond_max_target_seqs
            if args.diamond_max_target_seqs is not None
            else DEFAULT_MAX_TARGET_SEQS
        ),
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
        help="Path to DIAMOND-formatted annotated database (.dmnd). "
             "Required unless --mock-hits is used.",
    )
    run_parser.add_argument(
        "--annotations-db",
        help="Path to companion annotations SQLite database (.annotations.db). "
             "Required unless --mock-annotations is used.",
    )
    run_parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for homology search (default: 1)",
    )

    # DIAMOND tuning
    diamond_group = run_parser.add_argument_group(
        "DIAMOND tuning",
        "Control DIAMOND's memory use and speed. Defaults are DIAMOND's own. "
        "Expect up to about 6 GB of memory per unit of block size.",
    )
    diamond_group.add_argument(
        "--diamond-block-size",
        type=float,
        default=None,
        metavar="GIGALETTERS",
        help="DIAMOND --block-size: billions of database letters processed per block. "
             "Larger is faster but uses more memory (default: DIAMOND's 2.0).",
    )
    diamond_group.add_argument(
        "--diamond-index-chunks",
        type=int,
        default=None,
        metavar="N",
        help="DIAMOND --index-chunks: 1 is fastest and uses the most memory "
             "(default: DIAMOND's 4).",
    )
    diamond_group.add_argument(
        "--diamond-tmpdir",
        default=None,
        metavar="DIR",
        help="DIAMOND --tmpdir for intermediate files, e.g. /dev/shm to keep them in RAM "
             "(default: the work directory).",
    )
    diamond_group.add_argument(
        "--diamond-max-target-seqs",
        type=int,
        default=None,
        metavar="N",
        help="DIAMOND --max-target-seqs: hits kept per query protein before filtering; "
             "raised to --top-k if that is larger, 0 means unlimited "
             f"(default: {DEFAULT_MAX_TARGET_SEQS}). DIAMOND's own default of 25 would "
             "silently drop hits, including ties, before the tie-aware top-k filter.",
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
        default="is_a,part_of",
        help="Comma-separated GO edge types for closure (default: is_a,part_of). "
             f"Allowed: {', '.join(sorted(GO_EDGE_TYPES))}",
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
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
