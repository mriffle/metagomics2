"""Unit tests for DIAMOND execution and result parsing."""

import logging
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from metagomics2.core.diamond import (
    DIAMOND_OUTFMT_COLUMNS,
    DiamondError,
    _iter_lines_with_progress,
    count_queries_at_cap,
    estimate_diamond_memory_bytes,
    parse_diamond_output,
    run_diamond,
)


class TestParseDiamondOutput:
    """Tests for parsing DIAMOND outfmt 6 output."""

    def test_parse_standard_output(self, tmp_path):
        output = tmp_path / "results.tsv"
        output.write_text(
            "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\t100.0\n"
            "protA\tsp|P67890|UNIPROT\t85.0\t90\t10\t1\t1\t90\t1\t90\t1e-30\t150.0\t90.0\n"
            "protB\tsp|P11111|UNIPROT\t99.0\t200\t2\t0\t1\t200\t1\t200\t1e-80\t400.0\t66.7\n"
        )
        result = parse_diamond_output(output)
        assert result.n_queries == 2
        assert result.n_hits == 3
        assert len(result.hits_by_query["protA"]) == 2
        assert len(result.hits_by_query["protB"]) == 1
        assert result.hits_by_query["protA"][0].subject_id == "sp|P12345|UNIPROT"
        assert result.hits_by_query["protA"][0].pident == 95.0
        assert result.hits_by_query["protB"][0].evalue == 1e-80
        # Query coverage comes from the qcovhsp column DIAMOND is asked for
        assert result.hits_by_query["protA"][1].qcov == 90.0
        assert result.hits_by_query["protB"][0].qcov == 66.7

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
            "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\t100.0\n"
        )
        result = parse_diamond_output(output)
        assert result.n_queries == 1
        assert result.n_hits == 1

    def test_parse_rejects_rows_without_coverage_column(self, tmp_path):
        """A file written with plain --outfmt 6 (12 columns) is a column mismatch."""
        output = tmp_path / "results.tsv"
        output.write_text(
            "protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
        )
        with pytest.raises(ValueError, match="expected 13 tab-separated columns"):
            parse_diamond_output(output)


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
            output_text="protA\tsp|P12345|UNIPROT\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.0\t100.0\n",
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


class TestDiamondTuningOptions:
    """Block size, index chunks and tmpdir are passed through to DIAMOND."""

    def _run(self, mock_popen, tmp_path, **kwargs):
        query = tmp_path / "query.fasta"
        query.write_text(">protA\nACDE\n")
        db = tmp_path / "db.dmnd"
        output = tmp_path / "work" / "results.tsv"
        mock_popen.side_effect = _fake_popen(output_text="", output_path=output)
        result = run_diamond(query, db, output, **kwargs)
        return mock_popen.call_args[0][0], result

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_defaults_add_no_tuning_flags(self, mock_popen, tmp_path):
        cmd, result = self._run(mock_popen, tmp_path)
        assert "--block-size" not in cmd
        assert "--index-chunks" not in cmd
        assert "--tmpdir" not in cmd
        assert result.command == cmd

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_flags_passed_through(self, mock_popen, tmp_path):
        tmpdir = tmp_path / "shm"
        cmd, _ = self._run(
            mock_popen, tmp_path, block_size=8.0, index_chunks=1, tmpdir=tmpdir
        )
        assert cmd[cmd.index("--block-size") + 1] == "8"
        assert cmd[cmd.index("--index-chunks") + 1] == "1"
        assert cmd[cmd.index("--tmpdir") + 1] == str(tmpdir)
        # tmpdir is created so DIAMOND does not fail on a missing directory
        assert tmpdir.is_dir()

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_outfmt_requests_query_coverage_column(self, mock_popen, tmp_path):
        """--outfmt names every column explicitly, ending with qcovhsp."""
        cmd, _ = self._run(mock_popen, tmp_path)
        start = cmd.index("--outfmt")
        assert cmd[start + 1] == "6"
        assert cmd[start + 2 : start + 2 + len(DIAMOND_OUTFMT_COLUMNS)] == DIAMOND_OUTFMT_COLUMNS
        assert DIAMOND_OUTFMT_COLUMNS[-1] == "qcovhsp"
        # The next argument after the column list is another option, not a column
        assert cmd[start + 2 + len(DIAMOND_OUTFMT_COLUMNS)].startswith("--")

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_max_target_seqs_flag_pairing(self, mock_popen, tmp_path):
        cmd, _ = self._run(mock_popen, tmp_path, max_target_seqs=500)
        assert cmd[cmd.index("--max-target-seqs") + 1] == "500"
        cmd, _ = self._run(mock_popen, tmp_path, max_target_seqs=0)
        assert cmd[cmd.index("--max-target-seqs") + 1] == "0"
        cmd, _ = self._run(mock_popen, tmp_path)
        assert "--max-target-seqs" not in cmd

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_fractional_block_size_formatting(self, mock_popen, tmp_path):
        cmd, _ = self._run(mock_popen, tmp_path, block_size=0.5)
        assert cmd[cmd.index("--block-size") + 1] == "0.5"

    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_memory_budget_logged(self, mock_popen, tmp_path, caplog):
        with caplog.at_level(logging.INFO, logger="metagomics2.core.diamond"):
            self._run(mock_popen, tmp_path, block_size=4.0)
        budget = [r.message for r in caplog.records if "expect up to about" in r.message]
        assert budget == [
            "DIAMOND block size 4 billion letters: expect up to about 24.0 GB of memory"
        ]

    @patch("metagomics2.core.diamond.get_total_memory", return_value=8 * 1024**3)
    @patch("metagomics2.core.diamond.get_cgroup_memory_limit", return_value=None)
    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_warns_when_budget_exceeds_machine(self, mock_popen, _limit, _total, tmp_path, caplog):
        with caplog.at_level(logging.WARNING, logger="metagomics2.core.diamond"):
            self._run(mock_popen, tmp_path, block_size=4.0)
        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "machine has 8.0 GB" in warnings[0]
        assert "lower METAGOMICS_DIAMOND_BLOCK_SIZE" in warnings[0]

    @patch("metagomics2.core.diamond.get_total_memory", return_value=700 * 1024**3)
    @patch("metagomics2.core.diamond.get_cgroup_memory_limit", return_value=16 * 1024**3)
    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_warns_when_budget_exceeds_container_limit(
        self, mock_popen, _limit, _total, tmp_path, caplog
    ):
        with caplog.at_level(logging.WARNING, logger="metagomics2.core.diamond"):
            self._run(mock_popen, tmp_path, block_size=4.0)
        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert len(warnings) == 1
        assert "container memory limit is 16.0 GB" in warnings[0]

    @patch("metagomics2.core.diamond.get_total_memory", return_value=700 * 1024**3)
    @patch("metagomics2.core.diamond.get_cgroup_memory_limit", return_value=None)
    @patch("metagomics2.core.diamond.subprocess.Popen")
    def test_no_warning_when_budget_fits(self, mock_popen, _limit, _total, tmp_path, caplog):
        with caplog.at_level(logging.WARNING, logger="metagomics2.core.diamond"):
            self._run(mock_popen, tmp_path, block_size=20.0)
        assert not [r for r in caplog.records if r.levelno == logging.WARNING]

    def test_estimate_uses_diamond_default_when_unset(self):
        assert estimate_diamond_memory_bytes(None) == estimate_diamond_memory_bytes(2.0)
        assert estimate_diamond_memory_bytes(1.0) == 6 * 1024**3


class TestStreamingParse:
    """The DIAMOND output is streamed rather than read into memory."""

    def test_progress_logged_every_n_lines(self, tmp_path, caplog):
        output = tmp_path / "results.tsv"
        row = "q{i}\ts\t90\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\t100\n"
        output.write_text("".join(row.format(i=i) for i in range(5)))

        with caplog.at_level(logging.INFO, logger="metagomics2.core.diamond"):
            lines = list(_iter_lines_with_progress(output, every=2))

        assert len(lines) == 5
        progress = [r.message for r in caplog.records if "lines read so far" in r.message]
        assert progress == [
            "Parsing DIAMOND output: 2 lines read so far",
            "Parsing DIAMOND output: 4 lines read so far",
        ]


class TestCountQueriesAtCap:
    def _hits(self, n):
        from metagomics2.core.filtering import HomologyHit
        return [
            HomologyHit("q", f"s{i}", evalue=1e-9, bitscore=50.0, pident=90.0, qcov=90.0, alnlen=10)
            for i in range(n)
        ]

    def test_counts_queries_at_or_above_cap(self):
        hits = {"a": self._hits(3), "b": self._hits(2), "c": self._hits(4)}
        assert count_queries_at_cap(hits, 3) == 2

    def test_zero_cap_means_unlimited(self):
        assert count_queries_at_cap({"a": self._hits(3)}, 0) == 0

    def test_empty(self):
        assert count_queries_at_cap({}, 25) == 0
