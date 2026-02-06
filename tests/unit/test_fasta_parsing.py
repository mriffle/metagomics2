"""Unit tests for FASTA file parsing and hashing."""

import io
from pathlib import Path

import pytest

from metagomics2.core.fasta import (
    FastaParsingError,
    FastaRecord,
    build_protein_dict,
    compute_file_sha256,
    compute_string_sha256,
    parse_fasta,
    parse_fasta_from_handle,
    parse_fasta_header,
    write_subset_fasta,
)


class TestParseFastaHeader:
    """Tests for FASTA header parsing."""

    def test_simple_header(self):
        record_id, description = parse_fasta_header(">protein1")
        assert record_id == "protein1"
        assert description == ""

    def test_header_with_description(self):
        record_id, description = parse_fasta_header(">protein1 some description here")
        assert record_id == "protein1"
        assert description == "some description here"

    def test_header_with_tabs(self):
        record_id, description = parse_fasta_header(">protein1\tdescription")
        assert record_id == "protein1"
        assert description == "description"

    def test_invalid_header_no_gt(self):
        with pytest.raises(FastaParsingError) as exc_info:
            parse_fasta_header("protein1")
        assert "must start with '>'" in str(exc_info.value)

    def test_empty_header(self):
        with pytest.raises(FastaParsingError) as exc_info:
            parse_fasta_header(">")
        assert "Empty header" in str(exc_info.value)

    def test_whitespace_only_header(self):
        with pytest.raises(FastaParsingError) as exc_info:
            parse_fasta_header(">   ")
        assert "Empty header" in str(exc_info.value)


class TestParseFastaFromHandle:
    """Tests for parsing FASTA from file handles."""

    def test_single_record(self):
        content = ">protein1\nMPEPTIDEKAAA\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 1
        assert records[0].id == "protein1"
        assert records[0].sequence == "MPEPTIDEKAAA"

    def test_multiple_records(self):
        content = ">protein1\nMPEPTIDEKAAA\n>protein2\nTTTPEPTIDETTT\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 2
        assert records[0].id == "protein1"
        assert records[1].id == "protein2"

    def test_multiline_sequence(self):
        content = ">protein1\nMPEPTIDE\nKAAA\nMORE\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 1
        assert records[0].sequence == "MPEPTIDEKAAAMORE"

    def test_empty_lines_ignored(self):
        content = ">protein1\n\nMPEPTIDEKAAA\n\n>protein2\nTTT\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 2
        assert records[0].sequence == "MPEPTIDEKAAA"

    def test_sequence_before_header_raises(self):
        content = "MPEPTIDEKAAA\n>protein1\nSEQ\n"
        handle = io.StringIO(content)

        with pytest.raises(FastaParsingError) as exc_info:
            list(parse_fasta_from_handle(handle))
        assert "before first header" in str(exc_info.value)

    def test_preserves_description(self):
        content = ">protein1 this is a description\nMPEPTIDEKAAA\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert records[0].description == "this is a description"

    def test_empty_file(self):
        content = ""
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 0

    def test_only_empty_lines(self):
        content = "\n\n\n"
        handle = io.StringIO(content)
        records = list(parse_fasta_from_handle(handle))

        assert len(records) == 0


class TestParseFastaFile:
    """Tests for parsing FASTA files."""

    def test_parse_fixture_file(self, small_background_fasta: Path):
        records = parse_fasta(small_background_fasta)

        assert len(records) == 3
        assert records[0].id == "B1"
        assert records[0].sequence == "MPEPTIDEKAAA"
        assert records[1].id == "B2"
        assert records[1].sequence == "TTTPEPTIDETTT"
        assert records[2].id == "B3"
        assert records[2].sequence == "XXXXABCYYYY"

    def test_file_not_found(self, tmp_path: Path):
        with pytest.raises(FastaParsingError) as exc_info:
            parse_fasta(tmp_path / "nonexistent.fasta")
        assert "not found" in str(exc_info.value).lower()


class TestComputeHash:
    """Tests for hash computation."""

    def test_string_sha256_deterministic(self):
        content = "test content"
        hash1 = compute_string_sha256(content)
        hash2 = compute_string_sha256(content)
        assert hash1 == hash2

    def test_string_sha256_known_value(self):
        # Known SHA256 for "hello"
        expected = "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        assert compute_string_sha256("hello") == expected

    def test_file_sha256_deterministic(self, small_background_fasta: Path):
        hash1 = compute_file_sha256(small_background_fasta)
        hash2 = compute_file_sha256(small_background_fasta)
        assert hash1 == hash2

    def test_file_sha256_is_hex(self, small_background_fasta: Path):
        file_hash = compute_file_sha256(small_background_fasta)
        # SHA256 produces 64 hex characters
        assert len(file_hash) == 64
        assert all(c in "0123456789abcdef" for c in file_hash)


class TestBuildProteinDict:
    """Tests for building protein dictionary."""

    def test_builds_dict(self):
        records = [
            FastaRecord(id="P1", description="", sequence="AAA"),
            FastaRecord(id="P2", description="", sequence="BBB"),
        ]
        protein_dict = build_protein_dict(records)

        assert protein_dict == {"P1": "AAA", "P2": "BBB"}

    def test_empty_list(self):
        protein_dict = build_protein_dict([])
        assert protein_dict == {}


class TestFastaRecord:
    """Tests for the FastaRecord dataclass."""

    def test_record_is_frozen(self):
        r = FastaRecord(id="P1", description="desc", sequence="AAA")
        with pytest.raises(AttributeError):
            r.id = "P2"  # type: ignore

    def test_record_equality(self):
        r1 = FastaRecord(id="P1", description="desc", sequence="AAA")
        r2 = FastaRecord(id="P1", description="desc", sequence="AAA")
        assert r1 == r2


class TestWriteSubsetFasta:
    """Tests for writing a subset FASTA from selected protein IDs."""

    def test_writes_selected_proteins(self, tmp_path):
        proteins = {"P1": "ACDE", "P2": "FGHI", "P3": "KLMN"}
        out = tmp_path / "subset.fasta"
        n = write_subset_fasta(proteins, {"P1", "P3"}, out)
        assert n == 2
        content = out.read_text()
        assert ">P1" in content
        assert "ACDE" in content
        assert ">P3" in content
        assert "KLMN" in content
        assert ">P2" not in content

    def test_skips_missing_ids(self, tmp_path):
        proteins = {"P1": "ACDE"}
        out = tmp_path / "subset.fasta"
        n = write_subset_fasta(proteins, {"P1", "MISSING"}, out)
        assert n == 1

    def test_empty_protein_ids(self, tmp_path):
        proteins = {"P1": "ACDE"}
        out = tmp_path / "subset.fasta"
        n = write_subset_fasta(proteins, set(), out)
        assert n == 0
        assert out.read_text() == ""

    def test_line_wrapping(self, tmp_path):
        proteins = {"P1": "A" * 100}
        out = tmp_path / "subset.fasta"
        write_subset_fasta(proteins, {"P1"}, out, line_width=40)
        lines = out.read_text().strip().split("\n")
        assert lines[0] == ">P1"
        assert lines[1] == "A" * 40
        assert lines[2] == "A" * 40
        assert lines[3] == "A" * 20

    def test_creates_parent_dirs(self, tmp_path):
        proteins = {"P1": "ACDE"}
        out = tmp_path / "deep" / "nested" / "subset.fasta"
        n = write_subset_fasta(proteins, {"P1"}, out)
        assert n == 1
        assert out.exists()

    def test_roundtrip(self, tmp_path):
        """Subset FASTA can be parsed back correctly."""
        proteins = {"P1": "ACDEFGHIKLMN", "P2": "PQRSTVWY"}
        out = tmp_path / "subset.fasta"
        write_subset_fasta(proteins, {"P1", "P2"}, out)
        records = parse_fasta(out)
        result = build_protein_dict(records)
        assert result == proteins
