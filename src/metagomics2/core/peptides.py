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


def normalize_sequence(
    sequence: str,
    allowed_alphabet: set[str] | None = None,
) -> str:
    """Normalize a peptide sequence.

    - Strip whitespace
    - Convert to uppercase
    - Validate against allowed alphabet

    Args:
        sequence: Raw peptide sequence
        allowed_alphabet: Set of allowed characters. If None, uses EXTENDED_AA_ALPHABET.

    Returns:
        Normalized sequence

    Raises:
        PeptideParsingError: If sequence contains invalid characters
    """
    if allowed_alphabet is None:
        allowed_alphabet = EXTENDED_AA_ALPHABET

    normalized = sequence.strip().upper()

    if not normalized:
        raise PeptideParsingError("Empty peptide sequence")

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


def parse_peptide_list(
    file_path: Path | str,
    sequence_column: str = "peptide_sequence",
    quantity_column: str = "quantity",
    allowed_alphabet: set[str] | None = None,
    has_header: bool = True,
) -> list[Peptide]:
    """Parse a peptide list from a CSV/TSV file.

    Args:
        file_path: Path to the peptide list file
        sequence_column: Name of the sequence column (if has_header=True)
        quantity_column: Name of the quantity column (if has_header=True)
        allowed_alphabet: Set of allowed amino acid characters
        has_header: Whether the file has a header row

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
            sequence_column=sequence_column,
            quantity_column=quantity_column,
            allowed_alphabet=allowed_alphabet,
            has_header=has_header,
        )


def parse_peptide_list_from_handle(
    handle: TextIO,
    sequence_column: str = "peptide_sequence",
    quantity_column: str = "quantity",
    allowed_alphabet: set[str] | None = None,
    has_header: bool = True,
) -> list[Peptide]:
    """Parse a peptide list from a file handle.

    Args:
        handle: File handle to read from
        sequence_column: Name of the sequence column (if has_header=True)
        quantity_column: Name of the quantity column (if has_header=True)
        allowed_alphabet: Set of allowed amino acid characters
        has_header: Whether the file has a header row

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

    if has_header:
        try:
            header = next(reader)
        except StopIteration:
            raise PeptideParsingError("Empty file")

        # Normalize header names
        header = [col.strip().lower() for col in header]

        try:
            seq_idx = header.index(sequence_column.lower())
        except ValueError:
            raise PeptideParsingError(
                f"Sequence column '{sequence_column}' not found in header: {header}"
            )

        try:
            qty_idx = header.index(quantity_column.lower())
        except ValueError:
            raise PeptideParsingError(
                f"Quantity column '{quantity_column}' not found in header: {header}"
            )
    else:
        # Without header, assume first column is sequence, second is quantity
        seq_idx = 0
        qty_idx = 1

    peptides = []
    for line_num, row in enumerate(reader, start=2 if has_header else 1):
        if not row or all(cell.strip() == "" for cell in row):
            continue  # Skip empty rows

        if len(row) <= max(seq_idx, qty_idx):
            raise PeptideParsingError(
                f"Line {line_num}: Not enough columns (expected at least {max(seq_idx, qty_idx) + 1})"
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

    return peptides
