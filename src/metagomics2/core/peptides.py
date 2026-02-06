"""Peptide list parsing and normalization."""

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

# Standard amino acid alphabet
STANDARD_AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
# Extended alphabet including ambiguous codes
EXTENDED_AA_ALPHABET = STANDARD_AA_ALPHABET | set("BJOUXZ")


@dataclass(frozen=True)
class Peptide:
    """A peptide with its associated quantity."""

    sequence: str
    quantity: float


class PeptideParsingError(Exception):
    """Raised when peptide parsing fails."""

    pass


_NON_UPPER_RE = re.compile(r"[^A-Z]")


def normalize_sequence(
    sequence: str,
    allowed_alphabet: set[str] | None = None,
) -> str:
    """Normalize a peptide sequence.

    - Strip whitespace
    - Convert to uppercase
    - Remove all non-uppercase-letter characters (e.g. modification
      annotations like ``[+80]``, dots, dashes, etc.)
    - Validate against allowed alphabet (if provided)

    Args:
        sequence: Raw peptide sequence
        allowed_alphabet: Set of allowed characters. If None, uses EXTENDED_AA_ALPHABET.

    Returns:
        Normalized sequence containing only uppercase amino-acid letters

    Raises:
        PeptideParsingError: If the resulting sequence is empty
    """
    if allowed_alphabet is None:
        allowed_alphabet = EXTENDED_AA_ALPHABET

    # Uppercase first, then strip everything that isn't A-Z
    normalized = _NON_UPPER_RE.sub("", sequence.strip().upper())

    if not normalized:
        raise PeptideParsingError("Empty peptide sequence after removing non-letter characters")

    invalid_chars = set(normalized) - allowed_alphabet
    if invalid_chars:
        raise PeptideParsingError(
            f"Invalid characters in sequence '{normalized}': {sorted(invalid_chars)}"
        )

    return normalized


def parse_quantity(value: str) -> float:
    """Parse a quantity value.

    Args:
        value: String representation of quantity

    Returns:
        Parsed quantity as float

    Raises:
        PeptideParsingError: If value is not a valid non-negative number
    """
    value = value.strip()

    if not value:
        raise PeptideParsingError("Empty quantity value")

    try:
        quantity = float(value)
    except ValueError:
        raise PeptideParsingError(f"Invalid quantity value: '{value}'")

    if quantity < 0:
        raise PeptideParsingError(f"Negative quantity not allowed: {quantity}")

    if not (quantity == quantity):  # NaN check
        raise PeptideParsingError("NaN quantity not allowed")

    return quantity


def detect_delimiter(line: str) -> str:
    """Detect the delimiter used in a CSV/TSV line.

    Args:
        line: First line of the file

    Returns:
        Detected delimiter (tab or comma)
    """
    if "\t" in line:
        return "\t"
    return ","


def _is_numeric(value: str) -> bool:
    """Check if a string represents a numeric value."""
    try:
        float(value.strip())
        return True
    except ValueError:
        return False


def parse_peptide_list(
    file_path: Path | str,
    allowed_alphabet: set[str] | None = None,
) -> list[Peptide]:
    """Parse a peptide list from a CSV/TSV file.

    The file should have two columns: peptide sequence and count/abundance.
    A header row is auto-detected: if the second column of the first row is
    numeric it is treated as data; otherwise it is skipped as a header.

    Args:
        file_path: Path to the peptide list file
        allowed_alphabet: Set of allowed amino acid characters

    Returns:
        List of Peptide objects

    Raises:
        PeptideParsingError: If parsing fails
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise PeptideParsingError(f"File not found: {file_path}")

    with open(file_path, "r", newline="", encoding="utf-8") as f:
        return parse_peptide_list_from_handle(
            f,
            allowed_alphabet=allowed_alphabet,
        )


def parse_peptide_list_from_handle(
    handle: TextIO,
    allowed_alphabet: set[str] | None = None,
) -> list[Peptide]:
    """Parse a peptide list from a file handle.

    The file should have two columns: peptide sequence and count/abundance.
    A header row is auto-detected: if the second column of the first row is
    numeric it is treated as data; otherwise it is skipped as a header.

    Args:
        handle: File handle to read from
        allowed_alphabet: Set of allowed amino acid characters

    Returns:
        List of Peptide objects

    Raises:
        PeptideParsingError: If parsing fails
    """
    first_line = handle.readline()
    if not first_line:
        raise PeptideParsingError("Empty file")

    delimiter = detect_delimiter(first_line)

    # Reset to beginning
    handle.seek(0)

    reader = csv.reader(handle, delimiter=delimiter)

    # Read first row and auto-detect header
    try:
        first_row = next(reader)
    except StopIteration:
        raise PeptideParsingError("Empty file")

    # Always use positional columns: first=sequence, second=quantity
    seq_idx = 0
    qty_idx = 1

    # Auto-detect header: if second column is numeric, first row is data
    has_header = len(first_row) > 1 and not _is_numeric(first_row[qty_idx])

    peptides = []

    if not has_header:
        # First row is data, process it
        if len(first_row) <= max(seq_idx, qty_idx):
            raise PeptideParsingError(
                f"Line 1: Not enough columns (expected at least 2)"
            )
        try:
            sequence = normalize_sequence(first_row[seq_idx], allowed_alphabet)
        except PeptideParsingError as e:
            raise PeptideParsingError(f"Line 1: {e}")
        try:
            quantity = parse_quantity(first_row[qty_idx])
        except PeptideParsingError as e:
            raise PeptideParsingError(f"Line 1: {e}")
        peptides.append(Peptide(sequence=sequence, quantity=quantity))

    for line_num, row in enumerate(reader, start=2):
        if not row or all(cell.strip() == "" for cell in row):
            continue  # Skip empty rows

        if len(row) <= max(seq_idx, qty_idx):
            raise PeptideParsingError(
                f"Line {line_num}: Not enough columns (expected at least 2)"
            )

        try:
            sequence = normalize_sequence(row[seq_idx], allowed_alphabet)
        except PeptideParsingError as e:
            raise PeptideParsingError(f"Line {line_num}: {e}")

        try:
            quantity = parse_quantity(row[qty_idx])
        except PeptideParsingError as e:
            raise PeptideParsingError(f"Line {line_num}: {e}")

        peptides.append(Peptide(sequence=sequence, quantity=quantity))

    return _aggregate_peptides(peptides)


def _aggregate_peptides(peptides: list[Peptide]) -> list[Peptide]:
    """Merge peptides that share the same sequence by summing quantities."""
    totals: dict[str, float] = {}
    for p in peptides:
        totals[p.sequence] = totals.get(p.sequence, 0.0) + p.quantity
    return [Peptide(sequence=seq, quantity=qty) for seq, qty in totals.items()]
