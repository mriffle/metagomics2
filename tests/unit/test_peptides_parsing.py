"""Unit tests for peptide list parsing and normalization."""

import io
from pathlib import Path

import pytest

from metagomics2.core.peptides import (
    EXTENDED_AA_ALPHABET,
    STANDARD_AA_ALPHABET,
    Peptide,
    PeptideParsingError,
    _aggregate_peptides,
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

    def test_strips_non_letter_characters(self):
        assert normalize_sequence("PEP*TIDE") == "PEPTIDE"
        assert normalize_sequence("PEP[+80]TIDE") == "PEPTIDE"
        assert normalize_sequence("P.E.P.T.I.D.E") == "PEPTIDE"
        assert normalize_sequence("PEPT-IDE") == "PEPTIDE"

    def test_strips_modification_annotations(self):
        assert normalize_sequence("PEPT[+79.966]IDE") == "PEPTIDE"
        assert normalize_sequence("[+42]PEPTIDE") == "PEPTIDE"
        assert normalize_sequence("PEPTIDE[+80]") == "PEPTIDE"

    def test_rejects_empty_after_stripping(self):
        with pytest.raises(PeptideParsingError) as exc_info:
            normalize_sequence("[+80]")
        assert "Empty" in str(exc_info.value)

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
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)

    def test_parse_tsv_skips_empty_rows(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t10\n\nABC\t5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2

    def test_parse_tsv_any_header_names(self):
        content = "seq\tcount\nPEPTIDE\t10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)

    def test_parse_tsv_strips_modifications(self):
        content = "peptide_sequence\tquantity\nPEPT[+80]IDE\t10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1
        assert peptides[0] == Peptide(sequence="PEPTIDE", quantity=10.0)

    def test_parse_tsv_empty_after_strip_reports_line_with_header(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t10\n[+80]\t5\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Line 3" in str(exc_info.value)

    def test_parse_tsv_empty_after_strip_reports_line_without_header(self):
        content = "[+80]\t5\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Line 1" in str(exc_info.value)

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


class TestAggregatePeptides:
    """Tests for merging peptides with identical sequences."""

    def test_no_duplicates(self):
        peptides = [
            Peptide(sequence="PEPTIDE", quantity=10.0),
            Peptide(sequence="ABC", quantity=5.0),
        ]
        result = _aggregate_peptides(peptides)
        assert len(result) == 2

    def test_sums_duplicates(self):
        peptides = [
            Peptide(sequence="PEPTIDE", quantity=10.0),
            Peptide(sequence="PEPTIDE", quantity=5.0),
        ]
        result = _aggregate_peptides(peptides)
        assert len(result) == 1
        assert result[0].sequence == "PEPTIDE"
        assert result[0].quantity == 15.0

    def test_sums_multiple_duplicates(self):
        peptides = [
            Peptide(sequence="PEPTIDE", quantity=10.0),
            Peptide(sequence="ABC", quantity=3.0),
            Peptide(sequence="PEPTIDE", quantity=5.0),
            Peptide(sequence="ABC", quantity=7.0),
        ]
        result = _aggregate_peptides(peptides)
        assert len(result) == 2
        by_seq = {p.sequence: p.quantity for p in result}
        assert by_seq["PEPTIDE"] == 15.0
        assert by_seq["ABC"] == 10.0

    def test_empty_list(self):
        assert _aggregate_peptides([]) == []


class TestDuplicateRawPeptideDetection:
    """Tests that exact duplicate raw peptide strings are rejected."""

    def test_exact_duplicate_raises_error_with_header(self):
        content = "seq\tqty\nPEPTIDE\t10\nPEPTIDE\t5\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Duplicate peptide 'PEPTIDE'" in str(exc_info.value)
        assert "line 2" in str(exc_info.value)

    def test_exact_duplicate_raises_error_without_header(self):
        content = "PEPTIDE\t10\nPEPTIDE\t5\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Duplicate peptide 'PEPTIDE'" in str(exc_info.value)
        assert "Line 2" in str(exc_info.value)
        assert "line 1" in str(exc_info.value)

    def test_exact_duplicate_reports_first_occurrence_line(self):
        content = "seq\tqty\nAAA\t1\nBBB\t2\nCCC\t3\nBBB\t4\n"
        handle = io.StringIO(content)

        with pytest.raises(PeptideParsingError) as exc_info:
            parse_peptide_list_from_handle(handle)
        assert "Line 5" in str(exc_info.value)
        assert "line 3" in str(exc_info.value)

    def test_different_raw_forms_same_normalized_allowed(self):
        content = "seq\tqty\nPEPTIDE\t10\nPEPT[+80]IDE\t5\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1
        assert peptides[0].sequence == "PEPTIDE"
        assert peptides[0].quantity == 15.0


class TestModificationStrippingAndAggregation:
    """End-to-end tests: modifications are stripped and duplicates are summed."""

    def test_modified_and_unmodified_are_merged(self):
        content = "peptide_sequence\tquantity\nPEPTIDE\t132\nPEPT[+80]IDE\t10\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1
        assert peptides[0].sequence == "PEPTIDE"
        assert peptides[0].quantity == 142.0

    def test_multiple_modifications_merged(self):
        content = (
            "seq\tqty\n"
            "PEPTIDE\t100\n"
            "PEPT[+80]IDE\t10\n"
            "[+42]PEPTIDE\t5\n"
            "PEPT[+79.966]IDE[+16]\t3\n"
        )
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 1
        assert peptides[0].sequence == "PEPTIDE"
        assert peptides[0].quantity == 118.0

    def test_different_sequences_stay_separate(self):
        content = "seq\tqty\nPEPTIDE\t10\nANOTHER\t5\nPEPT[+80]IDE\t3\n"
        handle = io.StringIO(content)
        peptides = parse_peptide_list_from_handle(handle)

        assert len(peptides) == 2
        by_seq = {p.sequence: p.quantity for p in peptides}
        assert by_seq["PEPTIDE"] == 13.0
        assert by_seq["ANOTHER"] == 5.0
