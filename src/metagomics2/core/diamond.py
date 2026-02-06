"""DIAMOND homology search execution and result parsing."""

import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from metagomics2.core.filtering import HomologyHit, parse_blast_tabular

logger = logging.getLogger(__name__)

# Matches UniProt-style subject IDs: db|ACCESSION|ENTRY_NAME
_UNIPROT_ID_RE = re.compile(r"^[a-z]{2}\|([A-Za-z0-9_-]+)\|")


def parse_uniprot_accession(subject_id: str) -> str:
    """Extract the bare UniProt accession from a DIAMOND subject ID.

    Handles formats like:
        sp|Q21HH2|RS2_SACD2  -> Q21HH2
        tr|A0A0A0MQG0|...    -> A0A0A0MQG0
        P12345                -> P12345  (bare accession, returned as-is)

    Args:
        subject_id: Full subject ID string from DIAMOND output

    Returns:
        Bare UniProt accession string
    """
    m = _UNIPROT_ID_RE.match(subject_id)
    if m:
        return m.group(1)
    # Assume it's already a bare accession
    return subject_id


class DiamondError(Exception):
    """Raised when DIAMOND execution fails."""

    pass


@dataclass
class DiamondResult:
    """Result of a DIAMOND search."""

    hits_by_query: dict[str, list[HomologyHit]]
    output_path: Path
    n_queries: int
    n_hits: int


def run_diamond(
    query_fasta: Path,
    db_path: Path,
    output_path: Path,
    evalue: float = 1e-10,
    max_target_seqs: int = 1,
    threads: int = 4,
) -> DiamondResult:
    """Run DIAMOND blastp and parse the results.

    Args:
        query_fasta: Path to the query FASTA file (subset of background proteome)
        db_path: Path to the DIAMOND-formatted database (.dmnd)
        output_path: Path to write the tabular output
        evalue: Maximum e-value threshold for DIAMOND search
        max_target_seqs: Maximum number of target sequences per query
        threads: Number of CPU threads to use

    Returns:
        DiamondResult with parsed hits

    Raises:
        DiamondError: If DIAMOND execution fails
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(db_path),
        "--outfmt", "6",
        "--evalue", str(evalue),
        "--max-target-seqs", str(max_target_seqs),
        "--threads", str(threads),
        "--out", str(output_path),
    ]

    logger.info(f"Running DIAMOND: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        raise DiamondError(
            "DIAMOND executable not found. Ensure 'diamond' is installed and on PATH."
        )

    if result.returncode != 0:
        stderr = result.stderr.strip()
        raise DiamondError(
            f"DIAMOND exited with code {result.returncode}: {stderr}"
        )

    logger.info(f"DIAMOND completed. Output: {output_path}")

    # Parse results
    return parse_diamond_output(output_path)


def parse_diamond_output(output_path: Path) -> DiamondResult:
    """Parse DIAMOND outfmt 6 tabular output.

    Args:
        output_path: Path to the DIAMOND output file

    Returns:
        DiamondResult with parsed hits
    """
    if not output_path.exists():
        return DiamondResult(
            hits_by_query={},
            output_path=output_path,
            n_queries=0,
            n_hits=0,
        )

    with open(output_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    hits_by_query = parse_blast_tabular(lines)

    n_hits = sum(len(hits) for hits in hits_by_query.values())

    logger.info(
        f"Parsed {n_hits} DIAMOND hits for {len(hits_by_query)} query proteins"
    )

    return DiamondResult(
        hits_by_query=hits_by_query,
        output_path=output_path,
        n_queries=len(hits_by_query),
        n_hits=n_hits,
    )
