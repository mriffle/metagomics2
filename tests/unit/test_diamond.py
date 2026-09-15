"""Unit tests for DIAMOND execution and result parsing."""

import logging
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from metagomics2.core.diamond import (
    DiamondError,
    _iter_lines_with_progress,
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


def _fake_popen(returncode: int = 0, console: str = "", output_text: str | None = None,
                output_path=None, timeouts: int = 0):
    """Build a Popen side effect that writes DIAMOND-like console/output files."""

    def side_effect(cmd, stdout=None, stderr=None, **kwargs):
        if console:
            stdout.write(console.encode())
            stdout.flush()
        if output_text is not None and output_path is not None:
            output_path.write_text(output_text)
        proc = MagicMock()
        proc.pid = 4242
        remaining = [timeouts]

        def wait(timeout=None):
            if remaining[0] > 0:
                remaining[0] -= 1
                raise subprocess.TimeoutExpired(cmd, timeout)
            return returncode

        proc.wait.side_effect = wait
        return proc

    return side_effect


class TestRunDiamond:
    """Tests for DIAMOND execution."""

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_successful_run(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        db.write_text("")  # dummy
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = _fake_popen(
            console="Processing query block 1, reference block 1/1\nTotal time = 1.0s\n",
            output_text="protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\n",
            output_path=output,
        )

        result = run_diamond(query, db, output, threads=2)
        assert result.n_queries == 1
        assert result.n_hits == 1

        # Verify command was called correctly
        cmd = mock_popen.call_args[0][0]
        assert cmd[0] == "diamond"
        assert cmd[1] == "blastp"
        assert "--threads" in cmd
        assert "2" in cmd

        # DIAMOND's console output is preserved in a log file next to the output
        log_text = (tmp_path / "diamond.log").read_text()
        assert log_text.startswith("# ")
        assert "command: diamond blastp" in log_text
        assert "Processing query block 1" in log_text

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_custom_log_path(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "work" / "results.tsv"
        log_path = tmp_path / "logs" / "diamond.log"

        mock_popen.side_effect = _fake_popen(console="hello\n", output_text="", output_path=output)
        run_diamond(query, db, output, log_path=log_path)

        assert log_path.exists()
        assert "hello" in log_path.read_text()
        # stdout and stderr both go to the log file
        kwargs = mock_popen.call_args[1]
        assert kwargs["stderr"] is subprocess.STDOUT

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_diamond_failure_includes_log_tail(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = _fake_popen(returncode=1, console="Error: database not found\n")

        with pytest.raises(DiamondError, match="exited with code 1") as exc_info:
            run_diamond(query, db, output)
        assert "database not found" in str(exc_info.value)

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_diamond_killed_by_signal(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = _fake_popen(returncode=-9, console="Processing query block 3\n")

        with pytest.raises(DiamondError, match="killed by signal 9 \\(SIGKILL\\)") as exc_info:
            run_diamond(query, db, output)
        message = str(exc_info.value)
        assert "out-of-memory" in message
        assert "Processing query block 3" in message

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_diamond_not_found(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = FileNotFoundError()

        with pytest.raises(DiamondError, match="not found"):
            run_diamond(query, db, output)

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_passes_evalue_and_max_target_seqs(self, mock_popen, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = _fake_popen(output_text="", output_path=output)

        run_diamond(query, db, output, evalue=1e-5, max_target_seqs=3, threads=8)

        cmd = mock_popen.call_args[0][0]
        assert "1e-05" in cmd
        assert "3" in cmd
        assert "8" in cmd

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_heartbeat_logged_while_running(self, mock_popen, tmp_path, caplog):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "results.tsv"

        mock_popen.side_effect = _fake_popen(
            console="Processing query block 1, reference block 2/5\n",
            output_text="",
            output_path=output,
            timeouts=2,
        )

        with caplog.at_level(logging.INFO, logger="metagomics2.core.diamond"):
            run_diamond(query, db, output, heartbeat_seconds=0.01)

        heartbeats = [r for r in caplog.records if "DIAMOND still running" in r.message]
        assert len(heartbeats) == 2
        assert "pid 4242" in heartbeats[0].message
        assert "reference block 2/5" in heartbeats[0].message


class TestStreamingParse:
    """The DIAMOND output is streamed rather than read into memory."""

    def test_progress_logged_every_n_lines(self, tmp_path, caplog):
        output = tmp_path / "results.tsv"
        row = "q{i}\ts\t90\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n"
        output.write_text("".join(row.format(i=i) for i in range(5)))

        with caplog.at_level(logging.INFO, logger="metagomics2.core.diamond"):
            lines = list(_iter_lines_with_progress(output, every=2))

        assert len(lines) == 5
        progress = [r.message for r in caplog.records if "lines read so far" in r.message]
        assert progress == [
            "Parsing DIAMOND output: 2 lines read so far",
            "Parsing DIAMOND output: 4 lines read so far",
        ]
