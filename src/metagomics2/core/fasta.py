"""FASTA file parsing and hashing."""

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""

    id: str
    description: str
    sequence: str


class FastaParsingError(Exception):
    """Raised when FASTA parsing fails."""

    pass


def parse_fasta_header(line: str) -> tuple[str, str]:
    """Parse a FASTA header line.

    Args:
        line: Header line starting with '>'

    Returns:
        Tuple of (id, description)

    Raises:
        FastaParsingError: If header is invalid
    """
    if not line.startswith(">"):
        raise FastaParsingError(f"Invalid header line (must start with '>'): {line}")

    header = line[1:].strip()
    if not header:
        raise FastaParsingError("Empty header")

    # Split on first whitespace
    parts = header.split(None, 1)
    record_id = parts[0]
    description = parts[1] if len(parts) > 1 else ""

    return record_id, description


def parse_fasta(file_path: Path | str) -> list[FastaRecord]:
    """Parse a FASTA file.

    Args:
        file_path: Path to the FASTA file

    Returns:
        List of FastaRecord objects

    Raises:
        FastaParsingError: If parsing fails
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FastaParsingError(f"File not found: {file_path}")

    with open(file_path, "r", encoding="utf-8") as f:
        return list(parse_fasta_from_handle(f))


def parse_fasta_from_handle(handle: TextIO) -> Iterator[FastaRecord]:
    """Parse FASTA records from a file handle.

    Args:
        handle: File handle to read from

    Yields:
        FastaRecord objects

    Raises:
        FastaParsingError: If parsing fails
    """
    current_id: str | None = None
    current_description: str = ""
    current_sequence_parts: list[str] = []
    line_num = 0

    for line in handle:
        line_num += 1
        line = line.rstrip("\n\r")

        if not line:
            continue  # Skip empty lines

        if line.startswith(">"):
            # Yield previous record if exists
            if current_id is not None:
                yield FastaRecord(
                    id=current_id,
                    description=current_description,
                    sequence="".join(current_sequence_parts),
                )

            try:
                current_id, current_description = parse_fasta_header(line)
            except FastaParsingError as e:
                raise FastaParsingError(f"Line {line_num}: {e}")

            current_sequence_parts = []
        else:
            if current_id is None:
                raise FastaParsingError(
                    f"Line {line_num}: Sequence data before first header"
                )
            # Accumulate sequence (handles line wrapping)
            current_sequence_parts.append(line.strip())

    # Yield last record
    if current_id is not None:
        yield FastaRecord(
            id=current_id,
            description=current_description,
            sequence="".join(current_sequence_parts),
        )


def compute_file_sha256(file_path: Path | str) -> str:
    """Compute SHA256 hash of a file.

    Args:
        file_path: Path to the file

    Returns:
        Hex-encoded SHA256 hash
    """
    file_path = Path(file_path)
    sha256 = hashlib.sha256()

    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)

    return sha256.hexdigest()


def compute_string_sha256(content: str) -> str:
    """Compute SHA256 hash of a string.

    Args:
        content: String content

    Returns:
        Hex-encoded SHA256 hash
    """
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def build_protein_dict(records: list[FastaRecord]) -> dict[str, str]:
    """Build a dictionary mapping protein IDs to sequences.

    Args:
        records: List of FastaRecord objects

    Returns:
        Dictionary mapping protein ID to sequence
    """
    return {record.id: record.sequence for record in records}
