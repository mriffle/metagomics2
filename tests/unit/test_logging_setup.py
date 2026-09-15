"""Unit tests for the shared logging configuration."""

import logging
from pathlib import Path

from metagomics2 import logging_setup
from metagomics2.logging_setup import (
    attach_file_handler,
    configure_logging,
    detach_handler,
    format_bytes,
    get_cgroup_memory_limit,
    log_system_info,
    peak_rss_bytes,
    process_rss_bytes,
)


def _installed_handlers() -> list[logging.Handler]:
    return list(logging_setup._installed)


class TestConfigureLogging:
    def test_writes_to_stderr_and_file(self, tmp_path: Path, capsys):
        configure_logging("worker", tmp_path / "logs", "INFO")
        logging.getLogger("test.x").info("hello file and stderr")

        log_file = tmp_path / "logs" / "worker.log"
        assert log_file.exists()
        text = log_file.read_text()
        assert "hello file and stderr" in text
        assert "test.x" in text
        assert "INFO" in text
        assert "hello file and stderr" in capsys.readouterr().err

    def test_no_file_when_log_dir_is_none(self, tmp_path: Path):
        configure_logging("cli", None, "DEBUG")
        assert len(_installed_handlers()) == 1
        assert logging.getLogger().level == logging.DEBUG

    def test_reconfigure_replaces_handlers(self, tmp_path: Path):
        configure_logging("worker", tmp_path / "a", "INFO")
        first = _installed_handlers()
        configure_logging("worker", tmp_path / "b", "WARNING")
        second = _installed_handlers()

        root = logging.getLogger()
        for handler in first:
            assert handler not in root.handlers
        assert len(second) == 2
        assert root.level == logging.WARNING

    def test_level_accepts_int_or_name(self, tmp_path: Path):
        configure_logging("cli", None, logging.ERROR)
        assert logging.getLogger().level == logging.ERROR
        configure_logging("cli", None, "debug")
        assert logging.getLogger().level == logging.DEBUG

    def test_unwritable_log_dir_falls_back_to_stderr(self, tmp_path: Path):
        blocker = tmp_path / "file"
        blocker.write_text("")
        # A directory cannot be created under a regular file
        configure_logging("worker", blocker / "logs", "INFO")
        assert len(_installed_handlers()) == 1


class TestJobLogHandler:
    def test_attach_and_detach(self, tmp_path: Path):
        configure_logging("worker", None, "INFO")
        job_log = tmp_path / "jobs" / "abc" / "logs" / "pipeline.log"
        handler = attach_file_handler(job_log)
        assert handler is not None
        logging.getLogger("metagomics2.pipeline").info("stage one")
        detach_handler(handler)
        logging.getLogger("metagomics2.pipeline").info("after detach")

        text = job_log.read_text()
        assert "stage one" in text
        assert "after detach" not in text
        assert handler not in logging.getLogger().handlers

    def test_attach_unwritable_returns_none(self, tmp_path: Path):
        blocker = tmp_path / "file"
        blocker.write_text("")
        assert attach_file_handler(blocker / "x" / "pipeline.log") is None
        detach_handler(None)  # must be a no-op


class TestDiagnostics:
    def test_format_bytes(self):
        assert format_bytes(None) == "unknown"
        assert format_bytes(512) == "512 B"
        assert format_bytes(4 * 1024**3) == "4.0 GB"
        assert format_bytes(1536) == "1.5 KB"

    def test_cgroup_v2_limit(self, tmp_path: Path):
        v2 = tmp_path / "memory.max"
        v2.write_text("4294967296\n")
        assert get_cgroup_memory_limit(v2, tmp_path / "missing") == 4 * 1024**3

    def test_cgroup_v2_unlimited(self, tmp_path: Path):
        v2 = tmp_path / "memory.max"
        v2.write_text("max\n")
        assert get_cgroup_memory_limit(v2, tmp_path / "missing") is None

    def test_cgroup_v1_limit_and_unlimited(self, tmp_path: Path):
        v1 = tmp_path / "memory.limit_in_bytes"
        v1.write_text("2147483648\n")
        assert get_cgroup_memory_limit(tmp_path / "missing", v1) == 2 * 1024**3
        v1.write_text("9223372036854771712\n")
        assert get_cgroup_memory_limit(tmp_path / "missing", v1) is None

    def test_cgroup_files_missing(self, tmp_path: Path):
        assert get_cgroup_memory_limit(tmp_path / "a", tmp_path / "b") is None

    def test_peak_rss_positive(self):
        assert peak_rss_bytes() > 0

    def test_process_rss_missing_pid(self):
        assert process_rss_bytes(2**22 + 12345) is None

    def test_log_system_info(self, caplog):
        logger = logging.getLogger("test.sysinfo")
        with caplog.at_level(logging.INFO, logger="test.sysinfo"):
            log_system_info(logger)
        messages = [r.message for r in caplog.records]
        assert any(m.startswith("System:") and "CPUs" in m for m in messages)
        assert any(m.startswith("DIAMOND executable:") for m in messages)
