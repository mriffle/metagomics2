"""Peptide list parsing and normalization."""

import csv
import math
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

# One innermost bracketed group: ``[...]``, ``(...)`` or ``{...}`` with no
# bracket of the same kind inside.  Applied repeatedly so nested annotations
# such as ``(Oxidation (M))`` or ``[Phospho (STY)]`` are removed from the
# inside out.  Search engines put modification names and mass deltas in these
# groups, and the letters in a name must never leak into the sequence.
_BRACKET_GROUP_RE = re.compile(r"\[[^\[\]]*\]|\([^()]*\)|\{[^{}]*\}")

# Comet/MSFragger terminal-modification markers: a lowercase ``n`` at the
# start, or a lowercase ``c`` at the end, immediately attached to a bracketed
# group (``n[42.0106]PEPTIDE``, ``PEPTIDEKc[-0.98]``), optionally inside
# flanking-residue notation (``K.n[42.0106]PEPTIDE.R``).  Anchored to the ends
# so a lowercase ``c`` anywhere else is still cysteine.
_NTERM_MARKER_RE = re.compile(r"^((?:[A-Za-z]|-)\.)?n(?=[\[({])")
_CTERM_MARKER_RE = re.compile(r"c(?=[\[({][^\[\](){}]*[\])}](?:\.(?:[A-Za-z]|-))?$)")

# Flanking-residue notation ``K.PEPTIDE.R`` (``-`` for a protein terminus):
# a single letter or dash, a dot, the peptide, a dot, a single letter or dash.
# The peptide part must contain no dot, so a string of dot-separated letters
# is not mistaken for a flanked peptide.
_FLANKED_RE = re.compile(r"^(?:[A-Za-z]|-)\.([^.]+)\.(?:[A-Za-z]|-)$")


def normalize_sequence(
    sequence: str,
    allowed_alphabet: set[str] | None = None,
) -> str:
    """Normalize a peptide sequence.

    Rules, applied in this order:

    1. Strip surrounding whitespace.
    2. Remove modification annotations: every bracketed group ``[...]``,
       ``(...)`` or ``{...}``, including nested ones, whether it holds a mass
       delta (``PEPT[+79.966]IDE``) or a name (``C[Carbamidomethyl]``,
       ``M(Oxidation (M))``, ``(UniMod:4)``).  A lowercase ``n`` at the start or
       ``c`` at the end that is attached to such a group is a terminal-
       modification marker and is removed with it.
    3. Strip flanking residues: ``K.PEPTIDE.R`` or ``-.PEPTIDE.K`` becomes
       ``PEPTIDE``.
    4. Convert to uppercase and remove every remaining character that is not a
       letter (``_PEPTIDE_``, ``PEP*TIDE``, ``PEPT-IDE``, charge suffixes).
    5. Validate against the allowed alphabet.

    Args:
        sequence: Raw peptide sequence
        allowed_alphabet: Set of allowed characters. If None, uses EXTENDED_AA_ALPHABET.

    Returns:
        Normalized sequence containing only uppercase amino-acid letters

    Raises:
        PeptideParsingError: If the resulting sequence is empty or contains
            characters outside the allowed alphabet
    """
    if allowed_alphabet is None:
        allowed_alphabet = EXTENDED_AA_ALPHABET

    stripped = sequence.strip()

    # Terminal markers first: they are only meaningful next to a bracket group
    stripped = _NTERM_MARKER_RE.sub(r"\1", stripped)
    stripped = _CTERM_MARKER_RE.sub("", stripped)

    # Remove bracketed groups from the inside out until none are left
    while True:
        without_groups = _BRACKET_GROUP_RE.sub("", stripped)
        if without_groups == stripped:
            break
        stripped = without_groups

    flanked = _FLANKED_RE.match(stripped)
    if flanked:
        stripped = flanked.group(1)

    # Uppercase, then strip everything that isn't A-Z
    normalized = _NON_UPPER_RE.sub("", stripped.upper())

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
        PeptideParsingError: If value is not a finite, non-negative number
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

    if math.isnan(quantity):
        raise PeptideParsingError("NaN quantity not allowed")

    if math.isinf(quantity):
        raise PeptideParsingError(f"Infinite quantity not allowed: '{value}'")

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


def _check_no_extra_values(row: list[str], line_num: int) -> None:
    """Reject a data row that carries a value beyond the two expected columns.

    Trailing empty cells (``PEPTIDE,10,,`` as spreadsheets export them) are
    fine.  A non-empty third value is not: it usually means the quantity was
    written with a thousands separator (``PEPTIDE,1,000``) and has just been
    split into two columns, so reading only the second one would silently
    record the wrong number.
    """
    extra = [cell.strip() for cell in row[2:] if cell.strip()]
    if extra:
        raise PeptideParsingError(
            f"Line {line_num}: expected 2 columns (sequence, quantity) but found "
            f"an extra value '{extra[0]}'. Check for thousands separators in the "
            "quantity column or remove the additional columns."
        )


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

    with open(file_path, newline="", encoding="utf-8") as f:
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
    seen_raw: dict[str, int] = {}  # raw sequence -> first line number

    if not has_header:
        # First row is data, process it
        if len(first_row) <= max(seq_idx, qty_idx):
            raise PeptideParsingError(
                "Line 1: Not enough columns (expected at least 2)"
            )
        _check_no_extra_values(first_row, 1)
        raw_seq = first_row[seq_idx].strip()
        if raw_seq in seen_raw:
            raise PeptideParsingError(
                f"Line 1: Duplicate peptide '{raw_seq}' "
                f"(first seen on line {seen_raw[raw_seq]})"
            )
        seen_raw[raw_seq] = 1
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
        _check_no_extra_values(row, line_num)

        raw_seq = row[seq_idx].strip()
        if raw_seq in seen_raw:
            raise PeptideParsingError(
                f"Line {line_num}: Duplicate peptide '{raw_seq}' "
                f"(first seen on line {seen_raw[raw_seq]})"
            )
        seen_raw[raw_seq] = line_num

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
