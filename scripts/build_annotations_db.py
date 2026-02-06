#!/usr/bin/env python3
"""Build a companion annotations SQLite database for a UniProt DIAMOND database.

Usage:
    python -m scripts.build_annotations_db \
        --fasta uniprot_sprot.fasta.gz \
        --gaf goa_uniprot_all.gaf.gz \
        --output uniprot_sprot.annotations.db

Or via the CLI entry point:
    metagomics2-build-annotations \
        --fasta uniprot_sprot.fasta.gz \
        --gaf goa_uniprot_all.gaf.gz \
        --output uniprot_sprot.annotations.db
"""

import argparse
import logging
import sqlite3
import sys
import time
from pathlib import Path

# Add src to path so we can import metagomics2 modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from metagomics2.core.gaf_parser import parse_gaf_file
from metagomics2.core.uniprot_fasta import parse_uniprot_fasta_annotations

logger = logging.getLogger(__name__)

BATCH_SIZE = 50_000


def create_schema(conn: sqlite3.Connection) -> None:
    """Create the annotations database schema."""
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS taxonomy (
            accession TEXT PRIMARY KEY,
            tax_id INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS go_annotations (
            accession TEXT NOT NULL,
            go_id TEXT NOT NULL,
            aspect TEXT NOT NULL,
            UNIQUE(accession, go_id)
        );

        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
    """)


def build_taxonomy(conn: sqlite3.Connection, fasta_path: Path) -> tuple[int, set[str]]:
    """Parse UniProt FASTA and insert taxonomy mappings.

    Returns:
        Tuple of (number of rows inserted, set of accessions found)
    """
    logger.info(f"Parsing taxonomy from FASTA: {fasta_path}")
    count = 0
    accessions: set[str] = set()
    batch: list[tuple[str, int]] = []

    for accession, tax_id in parse_uniprot_fasta_annotations(fasta_path):
        accessions.add(accession)
        batch.append((accession, tax_id))
        if len(batch) >= BATCH_SIZE:
            conn.executemany(
                "INSERT OR REPLACE INTO taxonomy (accession, tax_id) VALUES (?, ?)",
                batch,
            )
            count += len(batch)
            batch.clear()
            if count % 500_000 == 0:
                logger.info(f"  taxonomy: {count:,} rows inserted...")

    if batch:
        conn.executemany(
            "INSERT OR REPLACE INTO taxonomy (accession, tax_id) VALUES (?, ?)",
            batch,
        )
        count += len(batch)

    conn.commit()
    logger.info(f"Taxonomy: {count:,} total rows inserted")
    return count, accessions


def build_go_annotations(
    conn: sqlite3.Connection,
    gaf_path: Path,
    fasta_accessions: set[str],
) -> int:
    """Parse GAF file and insert GO annotations.

    Only includes rows whose accession is present in ``fasta_accessions``
    (i.e., accessions that appear in the FASTA used to build the DIAMOND DB).
    Also excludes rows where the Qualifier contains "NOT" (negative
    annotations) and rows where the evidence code is "ND" (No biological
    Data).

    Returns:
        Number of rows inserted
    """
    logger.info(f"Parsing GO annotations from GAF: {gaf_path}")
    logger.info(f"Filtering to {len(fasta_accessions):,} accessions from FASTA")
    count = 0
    skipped = 0
    batch: list[tuple[str, str, str]] = []

    for record in parse_gaf_file(gaf_path):
        if record.accession not in fasta_accessions:
            skipped += 1
            if skipped % 10_000_000 == 0:
                logger.info(f"  skipped {skipped:,} GAF rows (accession not in FASTA)...")
            continue
        batch.append((record.accession, record.go_id, record.aspect))
        if len(batch) >= BATCH_SIZE:
            conn.executemany(
                "INSERT OR IGNORE INTO go_annotations (accession, go_id, aspect) "
                "VALUES (?, ?, ?)",
                batch,
            )
            count += len(batch)
            batch.clear()
            if count % 500_000 == 0:
                logger.info(f"  go_annotations: {count:,} rows inserted...")

    if batch:
        conn.executemany(
            "INSERT OR IGNORE INTO go_annotations (accession, go_id, aspect) "
            "VALUES (?, ?, ?)",
            batch,
        )
        count += len(batch)

    conn.commit()

    # Get actual row count (INSERT OR IGNORE silently skips duplicates)
    actual_count = conn.execute("SELECT COUNT(*) FROM go_annotations").fetchone()[0]
    duplicates = count - actual_count
    logger.info(
        f"GO annotations: {actual_count:,} unique rows inserted, "
        f"{duplicates:,} duplicates ignored, "
        f"{skipped:,} rows skipped (accession not in FASTA)"
    )
    return actual_count


def create_indexes(conn: sqlite3.Connection) -> None:
    """Create indexes after bulk inserts for performance."""
    logger.info("Creating indexes...")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_go_accession ON go_annotations(accession)")
    conn.commit()
    logger.info("Indexes created")


def build_annotations_db(
    fasta_path: Path,
    gaf_path: Path,
    output_path: Path,
) -> None:
    """Build the complete annotations SQLite database.

    Args:
        fasta_path: Path to UniProt FASTA file (plain or .gz)
        gaf_path: Path to GOA GAF file (plain or .gz)
        output_path: Path for the output SQLite database
    """
    start = time.time()

    # Remove existing DB if present
    if output_path.exists():
        output_path.unlink()

    output_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(str(output_path))
    conn.execute("PRAGMA journal_mode=DELETE")
    conn.execute("PRAGMA synchronous=NORMAL")

    try:
        create_schema(conn)

        tax_count, fasta_accessions = build_taxonomy(conn, fasta_path)
        go_count = build_go_annotations(conn, gaf_path, fasta_accessions)

        create_indexes(conn)

        # Store metadata
        conn.execute(
            "INSERT INTO metadata (key, value) VALUES (?, ?)",
            ("fasta_source", str(fasta_path.name)),
        )
        conn.execute(
            "INSERT INTO metadata (key, value) VALUES (?, ?)",
            ("gaf_source", str(gaf_path.name)),
        )
        conn.execute(
            "INSERT INTO metadata (key, value) VALUES (?, ?)",
            ("taxonomy_count", str(tax_count)),
        )
        conn.execute(
            "INSERT INTO metadata (key, value) VALUES (?, ?)",
            ("go_annotation_count", str(go_count)),
        )
        conn.execute(
            "INSERT INTO metadata (key, value) VALUES (?, ?)",
            ("gaf_filters", "exclude NOT qualifier, exclude ND evidence"),
        )
        conn.commit()

        elapsed = time.time() - start
        logger.info(
            f"Done. {tax_count:,} taxonomy + {go_count:,} GO annotations "
            f"in {elapsed:.1f}s -> {output_path}"
        )

    finally:
        conn.close()


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build a companion annotations SQLite database for a UniProt DIAMOND database."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Path to UniProt FASTA file (plain or .gz)",
    )
    parser.add_argument(
        "--gaf",
        required=True,
        type=Path,
        help="Path to GOA GAF file (plain or .gz)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output path for the annotations SQLite database",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    build_annotations_db(
        fasta_path=args.fasta,
        gaf_path=args.gaf,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
