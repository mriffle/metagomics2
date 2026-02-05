"""Unit tests for peptide list parsing and normalization."""

import io
from pathlib import Path

import pytest

from metagomics2.core.peptides import (
    EXTENDED_AA_ALPHABET,
    STANDARD_AA_ALPHABET,
    Peptide,
    PeptideParsingError,
    normalize_sequence,
    parse_peptide_list,
    parse_peptide_list_from_handle,
    parse_quantity,
)


class TestNormalizeSequence:
    """Tests for sequence normalization."""

    def test_strips_whitespace(self):
        assert normalize_sequence("  PEPTIDE  ") == "PEPTIDE"
        assert normalize_sequence("\tPEPTIDE\n") == "PEPTIDE"

    def test_converts_to_uppercase(self):
        assert normalize_sequence("peptide") == "PEPTIDE"
        assert normalize_sequence("PePtIdE") == "PEPTIDE"

    def test_combined_normalization(self):
        assert normalize_sequence("  pepTide  ") == "PEPTIDE"

    def test_rejects_invalid_characters_default_alphabet(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            normalize_sequence("PEP*TIDE")
        assert "*" in str(exc_info.value)

    def test_rejects_invalid_characters_shows_all(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            normalize_sequence("PEP*TI#DE")
        error_msg = str(exc_info.value)
        assert "*" in error_msg
        assert "#" in error_msg

    def test_accepts_custom_alphabet(self):
        # Extended alphabet allows more characters
        result = normalize_sequence("PEPTIDEX", allowed_alphabet=EXTENDED_AA_ALPHABET)
        assert result == "PEPTIDEX"

    def test_rejects_empty_sequence(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            normalize_sequence("")
        assert "Empty" in str(exc_info.value)

    def test_rejects_whitespace_only_sequence(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            normalize_sequence("   ")
        assert "Empty" in str(exc_info.value)

    def test_valid_standard_amino_acids(self):
        # All standard amino acids should be accepted
        all_aa = "ACDEFGHIKLMNPQRSTVWY"
        assert normalize_sequence(all_aa) == all_aa


class TestParseQuantity:
    """Tests for quantity parsing."""

    def test_parses_integer(self):
        assert parse_quantity("10") == 10.0

    def test_parses_float(self):
        assert parse_quantity("10.5") == 10.5

    def test_parses_scientific_notation(self):
        assert parse_quantity("1e5") == 100000.0

    def test_strips_whitespace(self):
        assert parse_quantity("  10  ") == 10.0

    def test_rejects_empty(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            parse_quantity("")
        assert "Empty" in str(exc_info.value)

    def test_rejects_non_numeric(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            parse_quantity("abc")
        assert "Invalid" in str(exc_info.value)

    def test_rejects_negative(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            parse_quantity("-5")
        assert "Negative" in str(exc_info.value)

    def test_accepts_zero(self):
        assert parse_quantity("0") == 0.0

    def test_rejects_nan(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            parse_quantity("nan")
        assert "NaN" in str(exc_info.value)


class TestParsePeptideListTSV:
    """Tests for parsing TSV peptide lists."""

    def test_parse_tsv_with_header(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t10\nABC\t5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)
        assert peptides[1] == Peptide(sequence="ABC", quantity=5.0)

    def test_parse_tsv_normalizes_sequences(self):
        content = "peptide_sequence\tquantity\n  peptide  \t10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert peptides[0].sequence == "PEPTIDE"

    def test_parse_tsv_without_header(self):
        content = "PEPTIDE\t10\nABC\t5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle, has_header=False)

        assert len(peptides) == 2
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)

    def test_parse_tsv_skips_empty_rows(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t10\n\nABC\t5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2

    def test_parse_tsv_case_insensitive_header(self):
        content = "Peptide_Sequence\tQuantity\nPEPTIDE\t10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1

    def test_parse_tsv_missing_sequence_column(self):
        content = "wrong_column\tquantity\nPEPTIDE\t10\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "peptide_sequence" in str(exc_info.value).lower()

    def test_parse_tsv_missing_quantity_column(self):
        content = "peptide_sequence\twrong\nPEPTIDE\t10\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "quantity" in str(exc_info.value).lower()

    def test_parse_tsv_invalid_sequence_reports_line(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t10\nINVALID*\t5\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Line 3" in str(exc_info.value)

    def test_parse_tsv_invalid_quantity_reports_line(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\tabc\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Line 2" in str(exc_info.value)


class TestParsePeptideListCSV:
    """Tests for parsing CSV peptide lists."""

    def test_parse_csv_with_header(self):
        content = "peptide_sequence,quantity\nPEPTIDE,10\nABC,5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)

    def test_auto_detects_csv_delimiter(self):
        content = "peptide_sequence,quantity\nPEPTIDE,10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1


class TestParsePeptideListFile:
    """Tests for parsing peptide lists from files."""

    def test_parse_fixture_file(self, small_peptides_tsv: Path):
        peptides = parse_peptide_list(small_peptides_tsv)

        assert len(peptides) == 3
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)
        assert peptides[1] == Peptide(sequence="ABC", quantity=5.0)
        assert peptides[2] == Peptide(sequence="NOMATCH", quantity=3.0)

    def test_file_not_found(self, tmp_path: Path):
        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list(tmp_path / "nonexistent.tsv")
        assert "not found" in str(exc_info.value).lower()

    def test_empty_file(self, tmp_path: Path):
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("")

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list(empty_file)
        assert "Empty" in str(exc_info.value)


class TestPeptideDataclass:
    """Tests for the Peptide dataclass."""

    def test_peptide_is_frozen(self):
        p = Peptide(sequence="PEPTIDE", quantity=10.0)
        with pytest.raises(AttributeError):
            p.sequence = "OTHER"  # type: ignore

    def test_peptide_equality(self):
        p1 = Peptide(sequence="PEPTIDE", quantity=10.0)
        p2 = Peptide(sequence="PEPTIDE", quantity=10.0)
        assert p1 == p2

    def test_peptide_hashable(self):
        p = Peptide(sequence="PEPTIDE", quantity=10.0)
        # Should be usable in sets/dicts
        s = {p}
        assert p in s
