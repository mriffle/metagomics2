"""Unit tests for DIAMOND execution and result parsing."""

import io
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from metagomics2.core.diamond import (
    DiamondError,
    DiamondResult,
    parse_diamond_output,
    run_diamond,
)


class TestParseDiamondOutput:
    """Tests for parsing DIAMOND outfmt 6 output."""

    def test_parse_standard_output(self, tmp_path):
        output = tmp_path / "results.tsv"
        output.write_text(
            "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
            "protA\tsp|P67890|UNIPROT\t85.0\t90\t10\t1\t1\t90\t1\t90\t1e-30\t150.0\n"
            "protB\tsp|P11111|UNIPROT\t99.0\t200\t2\t0\t1\t200\t1\t200\t1e-80\t400.0\n"
        )
        result = parse_diamond_output(output)
        assert result.n_queries == 2
        assert result.n_hits == 3
        assert len(result.hits_by_query["protA"]) == 2
        assert len(result.hits_by_query["protB"]) == 1
        assert result.hits_by_query["protA"][0].subject_id == "sp|P12345|UNIPROT"
        assert result.hits_by_query["protA"][0].pident == 95.0
        assert result.hits_by_query["protB"][0].evalue == 1e-80

    def test_parse_empty_output(self, tmp_path):
        output = tmp_path / "results.tsv"
        output.write_text("")
        result = parse_diamond_output(output)
        assert result.n_queries == 0
        assert result.n_hits == 0

    def test_parse_missing_file(self, tmp_path):
        output = tmp_path / "nonexistent.tsv"
        result = parse_diamond_output(output)
        assert result.n_queries == 0
        assert result.n_hits == 0

    def test_parse_skips_comments(self, tmp_path):
        output = tmp_path / "results.tsv"
        output.write_text(
            "# comment line\n"
            "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
        )
        result = parse_diamond_output(output)
        assert result.n_queries == 1
        assert result.n_hits == 1


class TestRunDiamond:
    """Tests for DIAMOND execution."""

    @patch("metagomics2.core.diamond.subprocess.run")
    def test_successful_run(self, mock_run, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        db.write_text("")  # dummy
        output = tmp_path / "results.tsv"

        # DIAMOND writes output file; simulate it
        def side_effect(*args, **kwargs):
            output.write_text(
                "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
            )
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stderr = ""
            return mock_result

        mock_run.side_effect = side_effect

        result = run_diamond(query, db, output, threads=2)
        assert result.n_queries == 1
        assert result.n_hits == 1

        # Verify command was called correctly
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert cmd[0] == "diamond"
        assert cmd[1] == "blastp"
        assert "--threads" in cmd
        assert "2" in cmd

    @patch("metagomics2.core.diamond.subprocess.run")
    def test_diamond_failure(self, mock_run, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "Error: database not found"
        mock_run.return_value = mock_result

        with pytest.raises(DiamondError, match="database not found"):
            run_diamond(query, db, output)

    @patch("metagomics2.core.diamond.subprocess.run")
    def test_diamond_not_found(self, mock_run, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_run.side_effect = FileNotFoundError()

        with pytest.raises(DiamondError, match="not found"):
            run_diamond(query, db, output)

    @patch("metagomics2.core.diamond.subprocess.run")
    def test_passes_evalue_and_max_target_seqs(self, mock_run, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        mock_run.return_value = mock_result
        output.write_text("")

        run_diamond(query, db, output, evalue=1e-5, max_target_seqs=3, threads=8)

        cmd = mock_run.call_args[0][0]
        assert "1e-05" in cmd
        assert "3" in cmd
        assert "8" in cmd
