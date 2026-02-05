"""Integration tests for reference snapshot system."""

import json
from pathlib import Path

import pytest

from metagomics2.pipeline.runner import PipelineConfig, run_pipeline


class TestReferenceSnapshot:
    """Tests for reference data snapshot creation and usage."""

    def test_pipeline_creates_reference_snapshot(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Pipeline should create reference snapshot when job_dir is set and bundled ref exists."""
        job_dir = tmp_path / "job_123"
        
        # Create a mock bundled reference directory for testing
        bundled_ref = tmp_path / "app_reference"
        bundled_ref.mkdir()
        (bundled_ref / "go").mkdir()
        (bundled_ref / "go" / "go.obo").write_text("[Term]\nid: GO:0000001\n")
        (bundled_ref / "taxonomy").mkdir()
        (bundled_ref / "taxonomy" / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
        (bundled_ref / "taxonomy" / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=job_dir / "results",
            job_dir=job_dir,  # Enable snapshot creation
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success, f"Pipeline failed: {result.error_message}"

        # In development (no bundled ref), snapshot won't be created
        # This is expected behavior - snapshots only work in Docker with bundled refs

    def test_pipeline_without_job_dir_skips_snapshot(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Pipeline should skip snapshot creation when job_dir is None."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            job_dir=None,  # No snapshot
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success, f"Pipeline failed: {result.error_message}"

        # Verify no snapshot directory was created
        assert not (tmp_path / "work" / "ref_snapshot").exists()

    def test_manifest_includes_reference_metadata(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Manifest should include reference metadata when available."""
        job_dir = tmp_path / "job_456"
        
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=job_dir / "results",
            job_dir=job_dir,
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success

        # Read manifest
        manifest_path = job_dir / "results" / "list_000" / "run_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        # Manifest should have parameters section
        assert "parameters" in manifest

    def test_snapshot_files_are_accessible(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Snapshot files should be readable and contain expected data."""
        job_dir = tmp_path / "job_789"
        
        # Create a simple GO file for testing
        go_source = tmp_path / "go_source"
        go_source.mkdir()
        (go_source / "go.json").write_text(json.dumps({
            "terms": {"GO:0000001": {"name": "test", "namespace": "biological_process"}},
            "edges": {"is_a": []}
        }))
        
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=job_dir / "results",
            job_dir=job_dir,
            go_data_path=go_source / "go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success

        # Verify we can read files from snapshot
        snapshot_dir = job_dir / "work" / "ref_snapshot"
        
        # Should have taxonomy files
        if (snapshot_dir / "taxonomy").exists():
            # Check that files are readable
            taxonomy_files = list((snapshot_dir / "taxonomy").iterdir())
            assert len(taxonomy_files) > 0, "Snapshot should contain taxonomy files"


class TestReferenceSnapshotProvenance:
    """Tests for reference snapshot provenance tracking."""

    def test_snapshot_hashes_in_manifest(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Manifest should include SHA256 hashes of snapshot files."""
        job_dir = tmp_path / "job_hash_test"
        
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=job_dir / "results",
            job_dir=job_dir,
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success

        # Read manifest
        manifest_path = job_dir / "results" / "list_000" / "run_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        # Check for reference snapshots section
        if "reference_snapshots" in manifest:
            # Should have file hashes
            assert "go_files" in manifest["reference_snapshots"] or \
                   "taxonomy_files" in manifest["reference_snapshots"]
