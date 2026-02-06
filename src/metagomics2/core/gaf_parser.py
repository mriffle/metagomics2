"""Parser for Gene Ontology Annotation (GAF) 2.2 files."""

import gzip
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterator

logger = logging.getLogger(__name__)

# GAF 2.2 column indices (0-based)
_COL_DB_OBJECT_ID = 1       # UniProt accession
_COL_QUALIFIER = 3          # Qualifier (e.g., "enables", "NOT|enables")
_COL_GO_ID = 4              # GO term ID (e.g., GO:0005634)
_COL_EVIDENCE_CODE = 6      # Evidence code (e.g., IEA, IDA)
_COL_ASPECT = 8             # Aspect: F (function), P (process), C (component)
_MIN_COLUMNS = 9            # Minimum columns required


@dataclass(frozen=True)
class GOARecord:
    """A single GO annotation record from a GAF file."""

    accession: str
    go_id: str
    aspect: str
    evidence_code: str


def parse_gaf_stream(
    fileobj: IO[str],
    exclude_not_qualifier: bool = True,
    exclude_nd_evidence: bool = True,
) -> Iterator[GOARecord]:
    """Yield GO annotation records from a GAF file stream.

    Args:
        fileobj: Text file handle (or gzip-decoded stream) for the GAF file
        exclude_not_qualifier: If True, skip rows where the Qualifier
            column contains "NOT" (e.g., "NOT|enables").
        exclude_nd_evidence: If True, skip rows where the evidence
            code is "ND" (No biological Data available).

    Yields:
        GOARecord for each qualifying annotation line
    """
    for line in fileobj:
        # Skip comment/header lines
        if line.startswith("!"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < _MIN_COLUMNS:
            continue

        # Skip negative annotations (qualifier contains "NOT")
        if exclude_not_qualifier and "NOT" in parts[_COL_QUALIFIER].upper():
            continue

        # Skip ND (No biological Data) evidence code
        if exclude_nd_evidence and parts[_COL_EVIDENCE_CODE] == "ND":
            continue

        yield GOARecord(
            accession=parts[_COL_DB_OBJECT_ID],
            go_id=parts[_COL_GO_ID],
            aspect=parts[_COL_ASPECT],
            evidence_code=parts[_COL_EVIDENCE_CODE],
        )


def parse_gaf_file(
    gaf_path: Path | str,
    exclude_not_qualifier: bool = True,
    exclude_nd_evidence: bool = True,
) -> Iterator[GOARecord]:
    """Parse a GAF file (plain text or gzipped).

    Args:
        gaf_path: Path to the GAF file (.gaf or .gaf.gz)
        exclude_not_qualifier: If True, skip rows where Qualifier contains "NOT"
        exclude_nd_evidence: If True, skip rows where evidence code is "ND"

    Yields:
        GOARecord for each qualifying annotation line
    """
    gaf_path = Path(gaf_path)

    if gaf_path.name.endswith(".gz"):
        with gzip.open(gaf_path, "rt", encoding="utf-8") as f:
            yield from parse_gaf_stream(f, exclude_not_qualifier, exclude_nd_evidence)
    else:
        with open(gaf_path, "r", encoding="utf-8") as f:
            yield from parse_gaf_stream(f, exclude_not_qualifier, exclude_nd_evidence)
