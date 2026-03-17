"""Integration tests for pipeline with mocked homology search."""

import csv
import json
from pathlib import Path

import pytest

from metagomics2.pipeline.runner import (
    PipelineConfig,
    PipelineProgress,
    _PROGRESS_PER_LIST_END,
    _PROGRESS_PER_LIST_START,
    _PROGRESS_TOTAL,
    run_pipeline,
)


class TestPipelineMockedHomology:
    """End-to-end pipeline tests with mocked homology output."""

    def test_full_pipeline_produces_expected_outputs(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Run full pipeline and verify output files are created."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)

        assert result.success, f"Pipeline failed: {result.error_message}"

        # Check output files exist
        list_dir = tmp_path / "results" / "list_000"
        assert list_dir.exists()
        assert (list_dir / "taxonomy_nodes.csv").exists()
        assert (list_dir / "go_terms.csv").exists()
        assert (list_dir / "coverage.csv").exists()
        assert (list_dir / "run_manifest.json").exists()

    def test_coverage_csv_values(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify coverage.csv has correct values."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)
        assert result.success

        coverage_path = tmp_path / "results" / "list_000" / "coverage.csv"
        with open(coverage_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        # PEPTIDE=10 (annotated), ABC=5 (unannotated), NOMATCH=3 (unannotated)
        assert float(row["total_peptide_quantity"]) == pytest.approx(18.0)
        assert float(row["annotated_peptide_quantity"]) == pytest.approx(10.0)
        assert float(row["unannotated_peptide_quantity"]) == pytest.approx(8.0)
        assert int(row["n_peptides_total"]) == 3
        assert int(row["n_peptides_annotated"]) == 1
        assert int(row["n_peptides_unannotated"]) == 2

    def test_taxonomy_nodes_csv_values(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify taxonomy_nodes.csv has correct values."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)
        assert result.success

        taxonomy_path = tmp_path / "results" / "list_000" / "taxonomy_nodes.csv"
        with open(taxonomy_path) as f:
            reader = csv.DictReader(f)
            rows = {int(row["tax_id"]): row for row in reader}

        # PEPTIDE maps to U1(70), U2(71), U3(72)
        # LCA of 70, 71, 72 is 30 (ClassA)
        # Lineage: 30, 20, 10, 1
        assert 30 in rows
        assert 20 in rows
        assert 10 in rows
        assert 1 in rows

        # Each node should have quantity 10 (only PEPTIDE is annotated)
        for tax_id in [30, 20, 10, 1]:
            assert float(rows[tax_id]["quantity"]) == pytest.approx(10.0)
            assert int(rows[tax_id]["n_peptides"]) == 1

        # Verify names
        assert rows[30]["name"] == "ClassA"
        assert rows[1]["name"] == "root"

    def test_go_terms_csv_values(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify go_terms.csv has correct values."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)
        assert result.success

        go_path = tmp_path / "results" / "list_000" / "go_terms.csv"
        with open(go_path) as f:
            reader = csv.DictReader(f)
            rows = {row["go_id"]: row for row in reader}

        # PEPTIDE maps to U1(C), U2(D), U3(E)
        # Closures (is_a only):
        # C -> A, B -> root_BP
        # D -> A -> root_BP
        # E -> B -> root_BP
        # Union: C, D, E, A, B, root_BP
        expected_terms = {
            "GO:0000001",  # root_BP
            "GO:0000002",  # A
            "GO:0000003",  # B
            "GO:0000004",  # C
            "GO:0000005",  # D
            "GO:0000006",  # E
        }
        assert set(rows.keys()) == expected_terms

        # Each term should have quantity 10
        for go_id in expected_terms:
            assert float(rows[go_id]["quantity"]) == pytest.approx(10.0)

    def test_manifest_contains_required_fields(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify run_manifest.json has required fields."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)
        assert result.success

        manifest_path = tmp_path / "results" / "list_000" / "run_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        assert "metagomics2_version" in manifest
        assert "python_version" in manifest
        assert "inputs" in manifest
        assert "fasta" in manifest["inputs"]
        assert "sha256" in manifest["inputs"]["fasta"]
        assert "parameters" in manifest
        assert "timestamp_utc" in manifest

    def test_progress_callback_called(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify progress callback is invoked with weighted progress."""
        progress_updates: list[PipelineProgress] = []

        def callback(progress: PipelineProgress) -> None:
            progress_updates.append(
                PipelineProgress(
                    total_peptide_lists=progress.total_peptide_lists,
                    completed_peptide_lists=progress.completed_peptide_lists,
                    current_stage=progress.current_stage,
                    current_list_id=progress.current_list_id,
                    progress_done=progress.progress_done,
                    progress_total=progress.progress_total,
                )
            )

        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config, progress_callback=callback)
        assert result.success

        # Should have received multiple progress updates
        assert len(progress_updates) > 0

        # Should include initialization and completion stages
        stages = [p.current_stage for p in progress_updates]
        assert any("Initializing" in s for s in stages)
        assert any("completed" in s.lower() for s in stages)

        # progress_total should always be _PROGRESS_TOTAL (1000)
        for p in progress_updates:
            assert p.progress_total == _PROGRESS_TOTAL

        # First update should start at 0, last should reach _PROGRESS_TOTAL
        assert progress_updates[0].progress_done == 0
        assert progress_updates[-1].progress_done == _PROGRESS_TOTAL

        # progress_done should be monotonically non-decreasing
        for i in range(1, len(progress_updates)):
            assert progress_updates[i].progress_done >= progress_updates[i - 1].progress_done, (
                f"progress_done decreased from {progress_updates[i - 1].progress_done} "
                f"to {progress_updates[i].progress_done} at stage "
                f"'{progress_updates[i].current_stage}'"
            )

    def test_progress_covers_all_stages(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify progress updates cover pre-list and per-list stages."""
        progress_updates: list[PipelineProgress] = []

        def callback(progress: PipelineProgress) -> None:
            progress_updates.append(
                PipelineProgress(
                    total_peptide_lists=progress.total_peptide_lists,
                    completed_peptide_lists=progress.completed_peptide_lists,
                    current_stage=progress.current_stage,
                    current_list_id=progress.current_list_id,
                    progress_done=progress.progress_done,
                    progress_total=progress.progress_total,
                )
            )

        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[fixtures_dir / "peptides" / "small_peptides.tsv"],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config, progress_callback=callback)
        assert result.success

        stages = [p.current_stage for p in progress_updates]

        # Pre-list stages should all appear
        assert any("Initializing" in s for s in stages)
        assert any("Parsing" in s for s in stages)
        assert any("Matching" in s for s in stages)
        assert any("subset FASTA" in s for s in stages)

        # Per-list stages should appear
        assert any("Processing peptide list" in s for s in stages)
        assert any("Completed" in s for s in stages)

        # Pre-list stages should have progress_done < _PROGRESS_PER_LIST_START
        pre_list = [p for p in progress_updates if "Initializing" in p.current_stage]
        assert all(p.progress_done < _PROGRESS_PER_LIST_START for p in pre_list)

    def test_progress_multiple_lists_divides_evenly(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Verify per-list progress is divided evenly among multiple lists."""
        progress_updates: list[PipelineProgress] = []

        def callback(progress: PipelineProgress) -> None:
            progress_updates.append(
                PipelineProgress(
                    total_peptide_lists=progress.total_peptide_lists,
                    completed_peptide_lists=progress.completed_peptide_lists,
                    current_stage=progress.current_stage,
                    current_list_id=progress.current_list_id,
                    progress_done=progress.progress_done,
                    progress_total=progress.progress_total,
                )
            )

        peptide_path = fixtures_dir / "peptides" / "small_peptides.tsv"
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[peptide_path, peptide_path, peptide_path],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config, progress_callback=callback)
        assert result.success

        # progress_done must be monotonically non-decreasing
        for i in range(1, len(progress_updates)):
            assert progress_updates[i].progress_done >= progress_updates[i - 1].progress_done

        # Final progress should reach _PROGRESS_TOTAL
        assert progress_updates[-1].progress_done == _PROGRESS_TOTAL

        # Each "Completed" update should be strictly higher than the previous one
        completed = [p for p in progress_updates if "Completed" in p.current_stage]
        assert len(completed) == 3
        for i in range(1, len(completed)):
            assert completed[i].progress_done > completed[i - 1].progress_done

        # The per-list progress should span from _PROGRESS_PER_LIST_START to _PROGRESS_PER_LIST_END
        per_list_starts = [p for p in progress_updates if "Processing peptide list" in p.current_stage]
        assert per_list_starts[0].progress_done == _PROGRESS_PER_LIST_START
        assert completed[-1].progress_done == _PROGRESS_PER_LIST_END

    def test_multiple_peptide_lists(
        self,
        fixtures_dir: Path,
        tmp_path: Path,
    ):
        """Test processing multiple peptide lists."""
        # Use the same peptide list twice for simplicity
        peptide_path = fixtures_dir / "peptides" / "small_peptides.tsv"

        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[peptide_path, peptide_path],
            output_dir=tmp_path / "results",
            go_data_path=fixtures_dir / "go" / "small_go.json",
            taxonomy_data_path=fixtures_dir / "taxonomy" / "small_taxonomy.json",
            mock_hits_path=fixtures_dir / "hits" / "accepted_hits.json",
            mock_subject_annotations_path=fixtures_dir / "annotations" / "subjects.json",
        )

        result = run_pipeline(config)
        assert result.success

        # Should have output for both lists
        assert (tmp_path / "results" / "list_000").exists()
        assert (tmp_path / "results" / "list_001").exists()


class TestPipelineErrorHandling:
    """Tests for pipeline error handling."""

    def test_missing_fasta_file(self, tmp_path: Path):
        """Pipeline should fail gracefully with missing FASTA."""
        config = PipelineConfig(
            fasta_path=tmp_path / "nonexistent.fasta",
            peptide_list_paths=[tmp_path / "peptides.tsv"],
            output_dir=tmp_path / "results",
        )

        result = run_pipeline(config)

        assert not result.success
        assert result.error_message is not None

    def test_missing_peptide_file(self, fixtures_dir: Path, tmp_path: Path):
        """Pipeline should fail gracefully with missing peptide file."""
        config = PipelineConfig(
            fasta_path=fixtures_dir / "fasta" / "small_background.fasta",
            peptide_list_paths=[tmp_path / "nonexistent.tsv"],
            output_dir=tmp_path / "results",
        )

        result = run_pipeline(config)

        assert not result.success
        assert result.error_message is not None
