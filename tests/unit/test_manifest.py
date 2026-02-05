"""Unit tests for manifest generation and provenance."""

import json
from pathlib import Path

import pytest

from metagomics2.core.reporting import (
    ManifestInfo,
    compute_file_hash,
    create_manifest,
    write_manifest_json,
)


class TestManifestRequiredKeys:
    """Tests for manifest required keys."""

    def test_manifest_contains_metagomics2_version(self, tmp_path: Path):
        manifest = ManifestInfo(metagomics2_version="0.1.0")
        output_path = tmp_path / "manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["metagomics2_version"] == "0.1.0"

    def test_manifest_contains_python_version(self, tmp_path: Path):
        manifest = ManifestInfo(python_version="3.11.0")
        output_path = tmp_path / "manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["python_version"] == "3.11.0"

    def test_manifest_contains_tool_versions(self, tmp_path: Path):
        manifest = ManifestInfo(
            search_tool="diamond",
            search_tool_version="diamond version 2.1.8",
            search_tool_command="diamond blastp -d db -q query.fasta",
        )
        output_path = tmp_path / "manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["search_tool"] == "diamond"
        assert data["search_tool_version"] == "diamond version 2.1.8"
        assert "diamond blastp" in data["search_tool_command"]

    def test_manifest_contains_file_hashes(self, tmp_path: Path):
        # Create test files
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">P1\nMPEPTIDE\n")
        peptide_file = tmp_path / "peptides.tsv"
        peptide_file.write_text("peptide_sequence\tquantity\nPEPTIDE\t10\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=fasta_file,
            peptide_list_path=peptide_file,
            parameters={},
        )

        assert manifest.input_fasta_hash != ""
        assert manifest.peptide_list_hash != ""
        assert len(manifest.input_fasta_hash) == 64
        assert len(manifest.peptide_list_hash) == 64

    def test_manifest_contains_parameters(self, tmp_path: Path):
        manifest = ManifestInfo(
            parameters={
                "max_evalue": 1e-5,
                "min_pident": 80.0,
                "top_k": 10,
            }
        )
        output_path = tmp_path / "manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["parameters"]["max_evalue"] == 1e-5
        assert data["parameters"]["min_pident"] == 80.0
        assert data["parameters"]["top_k"] == 10


class TestManifestFileHashCorrectness:
    """Tests for file hash correctness in manifest."""

    def test_fasta_hash_matches_computed(self, small_background_fasta: Path, tmp_path: Path):
        expected_hash = compute_file_hash(small_background_fasta)

        # Create a dummy peptide file
        peptide_file = tmp_path / "peptides.tsv"
        peptide_file.write_text("peptide_sequence\tquantity\nPEPTIDE\t10\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=small_background_fasta,
            peptide_list_path=peptide_file,
            parameters={},
        )

        assert manifest.input_fasta_hash == expected_hash

    def test_peptide_hash_matches_computed(self, small_peptides_tsv: Path, tmp_path: Path):
        expected_hash = compute_file_hash(small_peptides_tsv)

        # Create a dummy fasta file
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">P1\nMPEPTIDE\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=fasta_file,
            peptide_list_path=small_peptides_tsv,
            parameters={},
        )

        assert manifest.peptide_list_hash == expected_hash


class TestManifestTimestamp:
    """Tests for timestamp format."""

    def test_timestamp_is_iso8601_utc(self, tmp_path: Path):
        # Create test files
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">P1\nMPEPTIDE\n")
        peptide_file = tmp_path / "peptides.tsv"
        peptide_file.write_text("peptide_sequence\tquantity\nPEPTIDE\t10\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=fasta_file,
            peptide_list_path=peptide_file,
            parameters={},
        )

        # Should be ISO 8601 format with timezone
        timestamp = manifest.timestamp_utc
        assert "T" in timestamp
        assert "+" in timestamp or "Z" in timestamp


class TestManifestSnapshotHashes:
    """Tests for reference snapshot file hashes."""

    def test_go_snapshot_files_hashed(self, tmp_path: Path):
        # Create snapshot directory with files
        go_dir = tmp_path / "go_snapshot"
        go_dir.mkdir()
        (go_dir / "go.obo").write_text("test go content")
        (go_dir / "go-basic.obo").write_text("test go basic content")

        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">P1\nMPEPTIDE\n")
        peptide_file = tmp_path / "peptides.tsv"
        peptide_file.write_text("peptide_sequence\tquantity\nPEPTIDE\t10\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=fasta_file,
            peptide_list_path=peptide_file,
            parameters={},
            go_snapshot_dir=go_dir,
        )

        assert "go.obo" in manifest.go_snapshot_files
        assert "go-basic.obo" in manifest.go_snapshot_files
        assert len(manifest.go_snapshot_files["go.obo"]) == 64

    def test_taxonomy_snapshot_files_hashed(self, tmp_path: Path):
        # Create snapshot directory with files
        tax_dir = tmp_path / "taxonomy_snapshot"
        tax_dir.mkdir()
        (tax_dir / "nodes.dmp").write_text("test nodes content")
        (tax_dir / "names.dmp").write_text("test names content")

        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">P1\nMPEPTIDE\n")
        peptide_file = tmp_path / "peptides.tsv"
        peptide_file.write_text("peptide_sequence\tquantity\nPEPTIDE\t10\n")

        manifest = create_manifest(
            metagomics2_version="0.1.0",
            search_tool="diamond",
            search_tool_command="",
            annotated_db_choice="test_db",
            input_fasta_path=fasta_file,
            peptide_list_path=peptide_file,
            parameters={},
            taxonomy_snapshot_dir=tax_dir,
        )

        assert "nodes.dmp" in manifest.taxonomy_snapshot_files
        assert "names.dmp" in manifest.taxonomy_snapshot_files
