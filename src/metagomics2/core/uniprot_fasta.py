"""Parser for UniProt FASTA headers to extract accession and taxonomy ID."""

import gzip
import logging
import re
from pathlib import Path
from typing import IO, Iterator

logger = logging.getLogger(__name__)

# Matches UniProt FASTA header: >db|ACCESSION|ENTRY_NAME ... OX=TAXID ...
_HEADER_RE = re.compile(r"^>([a-z]{2})\|([A-Za-z0-9_-]+)\|")
_OX_RE = re.compile(r"\bOX=(\d+)\b")


def parse_uniprot_fasta_annotations_stream(
    fileobj: IO[str],
) -> Iterator[tuple[str, int]]:
    """Yield (accession, tax_id) from UniProt FASTA header lines.

    Only processes header lines (starting with '>').  Sequence lines are
    skipped.  Headers without a parseable OX= field are skipped with a
    warning.

    Args:
        fileobj: Text file handle for the FASTA file

    Yields:
        Tuple of (bare_accession, taxonomy_id)
    """
    for line in fileobj:
        if not line.startswith(">"):
            continue

        header_match = _HEADER_RE.match(line)
        if not header_match:
            continue

        accession = header_match.group(2)

        ox_match = _OX_RE.search(line)
        if not ox_match:
            logger.debug(f"No OX= field in header for {accession}, skipping")
            continue

        tax_id = int(ox_match.group(1))
        yield accession, tax_id


def parse_uniprot_fasta_annotations(
    fasta_path: Path | str,
) -> Iterator[tuple[str, int]]:
    """Parse a UniProt FASTA file and yield (accession, tax_id) tuples.

    Supports plain text and gzipped (.gz) FASTA files.

    Args:
        fasta_path: Path to the UniProt FASTA file

    Yields:
        Tuple of (bare_accession, taxonomy_id)
    """
    fasta_path = Path(fasta_path)

    if fasta_path.name.endswith(".gz"):
        with gzip.open(fasta_path, "rt", encoding="utf-8") as f:
            yield from parse_uniprot_fasta_annotations_stream(f)
    else:
        with open(fasta_path, "r", encoding="utf-8") as f:
            yield from parse_uniprot_fasta_annotations_stream(f)
