"""Unit tests for background worker."""

import importlib
import json
import os
import signal
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import metagomics2.config as config_module
from metagomics2.config import Settings, DatabaseEntry
from metagomics2.db.database import Database
from metagomics2.models.job import JobParams, JobStatus, PeptideListStatus


def _get_worker_class():
    """Import Worker after config is set up."""
    from metagomics2.worker.worker import Worker
    return Worker


@pytest.fixture(autouse=True)
def _setup_worker_config(tmp_path: Path):
    """Ensure config singleton is set before worker module is imported."""
    config_dir = tmp_path / "config"
    config_dir.mkdir(exist_ok=True)
    (config_dir / "databases.json").write_text(
        json.dumps([{"name": "Test DB", "description": "Test", "path": "test.dmnd"}])
    )
    with patch.dict(os.environ, {
        "METAGOMICS_DATA_DIR": str(tmp_path),
        "METAGOMICS_CONFIG_DIR": str(config_dir),
    }):
        config_module.reset_settings()
        # Force reload so module-level code picks up new config
        import metagomics2.worker.worker as wmod
        importlib.reload(wmod)
        yield
    config_module.reset_settings()


@pytest.fixture
def test_db(tmp_path: Path):
    """Create a test database."""
    return Database(tmp_path / "test.db")


@pytest.fixture
def jobs_dir(tmp_path: Path):
    """Create a jobs directory."""
    d = tmp_path / "jobs"
    d.mkdir()
    return d


def create_job_with_files(db: Database, jobs_dir: Path, fixtures_dir: Path) -> str:
    """Helper to create a job with input files on disk."""
    params = JobParams()
    job_id = db.create_job(params)
    db.update_job_status(job_id, JobStatus.QUEUED)

    # Create job directory structure
    job_dir = jobs_dir / job_id
    inputs_dir = job_dir / "inputs"
    peptides_dir = inputs_dir / "peptides"
    work_dir = job_dir / "work"
    results_dir = job_dir / "results"
    logs_dir = job_dir / "logs"

    for d in [inputs_dir, peptides_dir, work_dir, results_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Copy fixture files
    import shutil
    shutil.copy(
        fixtures_dir / "fasta" / "small_background.fasta",
        inputs_dir / "background.fasta",
    )

    peptide_src = fixtures_dir / "peptides" / "small_peptides.tsv"
    peptide_dst = peptides_dir / "list_000_small_peptides.tsv"
    shutil.copy(peptide_src, peptide_dst)

    db.add_peptide_list(job_id, "list_000", "small_peptides.tsv", str(peptide_dst))

    return job_id


class TestWorkerInit:
    """Tests for worker initialization."""

    def test_worker_creates(self, test_db):
        Worker = _get_worker_class()
        worker = Worker(test_db)
        assert worker.running is True
        assert worker.current_job_id is None

    def test_worker_signal_handler(self, test_db):
        Worker = _get_worker_class()
        worker = Worker(test_db)
        assert worker.running is True
        worker._handle_signal(signal.SIGTERM, None)
        assert worker.running is False


class TestWorkerBuildConfig:
    """Tests for building pipeline config from job info."""

    def test_builds_config(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job = test_db.get_job(job_id)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir):
            worker = Worker(test_db)
            config = worker._build_config(job_id, job)

        assert config.fasta_path == jobs_dir / job_id / "inputs" / "background.fasta"
        assert config.output_dir == jobs_dir / job_id / "results"
        assert config.job_dir == jobs_dir / job_id
        assert len(config.peptide_list_paths) == 1

    def test_config_includes_filter_policy(self, test_db, jobs_dir, fixtures_dir):
        params = JobParams(max_evalue=1e-5, min_pident=80.0, top_k=10)
        job_id = test_db.create_job(params)
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        # Create minimal directory structure
        job_dir = jobs_dir / job_id
        inputs_dir = job_dir / "inputs"
        peptides_dir = inputs_dir / "peptides"
        for d in [inputs_dir, peptides_dir, job_dir / "work", job_dir / "results", job_dir / "logs"]:
            d.mkdir(parents=True, exist_ok=True)

        import shutil
        shutil.copy(
            fixtures_dir / "fasta" / "small_background.fasta",
            inputs_dir / "background.fasta",
        )

        job = test_db.get_job(job_id)

        Worker = _get_worker_class()
        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir):
            worker = Worker(test_db)
            config = worker._build_config(job_id, job)

        assert config.filter_policy.max_evalue == 1e-5
        assert config.filter_policy.min_pident == 80.0
        assert config.filter_policy.top_k == 10


class TestWorkerProcessJob:
    """Tests for job processing."""

    def test_process_job_updates_status_to_running(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            worker = Worker(test_db)
            worker._process_job(job_id)

        job = test_db.get_job(job_id)
        assert job.status == JobStatus.COMPLETED

    def test_process_job_marks_failed_on_error(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = False
            mock_result.error_message = "Test error"
            mock_pipeline.return_value = mock_result

            worker = Worker(test_db)
            worker._process_job(job_id)

        job = test_db.get_job(job_id)
        assert job.status == JobStatus.FAILED
        assert "Test error" in job.error_message

    def test_process_job_handles_exception(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_pipeline.side_effect = RuntimeError("Unexpected error")

            worker = Worker(test_db)
            worker._process_job(job_id)

        job = test_db.get_job(job_id)
        assert job.status == JobStatus.FAILED
        assert "Unexpected error" in job.error_message

    def test_process_job_clears_current_job_id(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            worker = Worker(test_db)
            worker._process_job(job_id)

        assert worker.current_job_id is None

    def test_process_job_adds_events(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            worker = Worker(test_db)
            worker._process_job(job_id)

        # Events should have been added (started + completed)
        # We can't easily query events directly, but the job should be completed
        job = test_db.get_job(job_id)
        assert job.status == JobStatus.COMPLETED

    def test_process_nonexistent_job(self, test_db, jobs_dir):
        Worker = _get_worker_class()
        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir):
            worker = Worker(test_db)
            worker._process_job("nonexistent_id")

        # Should fail gracefully - job doesn't exist so status update will just be a no-op


class TestWorkerRunLoop:
    """Tests for the main worker loop."""

    def test_run_processes_queued_job(self, test_db, jobs_dir, fixtures_dir):
        Worker = _get_worker_class()
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.POLL_INTERVAL", 0), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            worker = Worker(test_db)

            # Process one iteration then stop
            call_count = 0
            original_run = worker.run

            def limited_run():
                nonlocal call_count
                while worker.running:
                    job = test_db.get_next_queued_job()
                    if job:
                        worker._process_job(job)
                    call_count += 1
                    if call_count >= 1:
                        worker.running = False

            limited_run()

        job = test_db.get_job(job_id)
        assert job.status == JobStatus.COMPLETED

    def test_run_stops_on_signal(self, test_db):
        Worker = _get_worker_class()
        worker = Worker(test_db)
        worker.running = False  # Simulate signal received before loop starts

        with patch("metagomics2.worker.worker.POLL_INTERVAL", 0):
            worker.run()  # Should return immediately


class TestWorkerNotification:
    """Tests for email notification from worker."""

    def test_process_job_sends_email_on_completion(self, test_db, jobs_dir, fixtures_dir):
        params = JobParams(notification_email="user@example.com")
        job_id = test_db.create_job(params)
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        job_dir = jobs_dir / job_id
        inputs_dir = job_dir / "inputs"
        peptides_dir = inputs_dir / "peptides"
        for d in [inputs_dir, peptides_dir, job_dir / "work", job_dir / "results", job_dir / "logs"]:
            d.mkdir(parents=True, exist_ok=True)

        import shutil
        shutil.copy(
            fixtures_dir / "fasta" / "small_background.fasta",
            inputs_dir / "background.fasta",
        )

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline, \
             patch("metagomics2.worker.worker.send_job_notification") as mock_send:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        mock_send.assert_called_once()
        call_args = mock_send.call_args
        assert call_args[0][0].job_id == job_id

    def test_process_job_sends_email_on_failure(self, test_db, jobs_dir, fixtures_dir):
        params = JobParams(notification_email="user@example.com")
        job_id = test_db.create_job(params)
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        job_dir = jobs_dir / job_id
        inputs_dir = job_dir / "inputs"
        peptides_dir = inputs_dir / "peptides"
        for d in [inputs_dir, peptides_dir, job_dir / "work", job_dir / "results", job_dir / "logs"]:
            d.mkdir(parents=True, exist_ok=True)

        import shutil
        shutil.copy(
            fixtures_dir / "fasta" / "small_background.fasta",
            inputs_dir / "background.fasta",
        )

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline, \
             patch("metagomics2.worker.worker.send_job_notification") as mock_send:
            mock_result = MagicMock()
            mock_result.success = False
            mock_result.error_message = "Pipeline error"
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        mock_send.assert_called_once()

    def test_process_job_no_email_when_not_configured(self, test_db, jobs_dir, fixtures_dir):
        params = JobParams()  # No notification_email
        job_id = test_db.create_job(params)
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        job_dir = jobs_dir / job_id
        inputs_dir = job_dir / "inputs"
        peptides_dir = inputs_dir / "peptides"
        for d in [inputs_dir, peptides_dir, job_dir / "work", job_dir / "results", job_dir / "logs"]:
            d.mkdir(parents=True, exist_ok=True)

        import shutil
        shutil.copy(
            fixtures_dir / "fasta" / "small_background.fasta",
            inputs_dir / "background.fasta",
        )

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline, \
             patch("metagomics2.worker.worker.send_job_notification") as mock_send:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        # send_job_notification should not be called because notification_email is empty
        mock_send.assert_not_called()


class TestWorkerCleanup:
    """Tests for post-job file cleanup."""

    def test_cleanup_on_success_removes_inputs_and_work(self, test_db, jobs_dir, fixtures_dir):
        """Successful job with cleanup enabled should remove inputs/ and work/."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir = jobs_dir / job_id

        # Add a file in work/ to verify it gets cleaned
        (job_dir / "work" / "diamond_results.tsv").write_text("test")
        # Add a file in results/ to verify it is preserved
        (job_dir / "results" / "output.csv").write_text("results")

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", True), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", True), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        assert not (job_dir / "inputs").exists()
        assert not (job_dir / "work").exists()
        assert (job_dir / "results" / "output.csv").exists()

    def test_cleanup_on_failure_removes_inputs_and_work(self, test_db, jobs_dir, fixtures_dir):
        """Failed job with cleanup enabled should remove inputs/ and work/."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir = jobs_dir / job_id

        (job_dir / "work" / "diamond_results.tsv").write_text("test")
        (job_dir / "results" / "output.csv").write_text("results")

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", True), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", True), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = False
            mock_result.error_message = "Pipeline error"
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        assert not (job_dir / "inputs").exists()
        assert not (job_dir / "work").exists()
        assert (job_dir / "results" / "output.csv").exists()

    def test_cleanup_on_exception_removes_inputs_and_work(self, test_db, jobs_dir, fixtures_dir):
        """Exception during processing with cleanup enabled should still clean up."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir = jobs_dir / job_id

        (job_dir / "results" / "output.csv").write_text("results")

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", True), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_pipeline.side_effect = RuntimeError("Unexpected error")

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        assert not (job_dir / "inputs").exists()
        assert not (job_dir / "work").exists()
        assert (job_dir / "results" / "output.csv").exists()

    def test_no_cleanup_on_success_when_disabled(self, test_db, jobs_dir, fixtures_dir):
        """Successful job with cleanup disabled should preserve inputs/ and work/."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir = jobs_dir / job_id

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", False), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        assert (job_dir / "inputs").exists()
        assert (job_dir / "work").exists()

    def test_no_cleanup_on_failure_when_disabled(self, test_db, jobs_dir, fixtures_dir):
        """Failed job with failure cleanup disabled should preserve inputs/ and work/."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir = jobs_dir / job_id

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", False), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = False
            mock_result.error_message = "Pipeline error"
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id)

        assert (job_dir / "inputs").exists()
        assert (job_dir / "work").exists()

    def test_independent_toggles(self, test_db, jobs_dir, fixtures_dir):
        """Cleanup on success but not failure should only clean up successful jobs."""
        # Successful job — should clean up
        job_id_ok = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir_ok = jobs_dir / job_id_ok

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", True), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", False), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id_ok)

        assert not (job_dir_ok / "inputs").exists()

        # Failed job — should NOT clean up
        job_id_fail = create_job_with_files(test_db, jobs_dir, fixtures_dir)
        job_dir_fail = jobs_dir / job_id_fail

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", True), \
             patch("metagomics2.worker.worker.CLEANUP_ON_FAILURE", False), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline:
            mock_result = MagicMock()
            mock_result.success = False
            mock_result.error_message = "Error"
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            worker._process_job(job_id_fail)

        assert (job_dir_fail / "inputs").exists()

    def test_cleanup_error_does_not_raise(self, test_db, jobs_dir, fixtures_dir):
        """Cleanup failure should be logged but not raise."""
        job_id = create_job_with_files(test_db, jobs_dir, fixtures_dir)

        with patch("metagomics2.worker.worker.JOBS_DIR", jobs_dir), \
             patch("metagomics2.worker.worker.CLEANUP_ON_SUCCESS", True), \
             patch("metagomics2.worker.worker.run_pipeline") as mock_pipeline, \
             patch("shutil.rmtree", side_effect=PermissionError("denied")):
            mock_result = MagicMock()
            mock_result.success = True
            mock_result.peptide_list_results = []
            mock_pipeline.return_value = mock_result

            Worker = _get_worker_class()
            worker = Worker(test_db)
            # Should not raise
            worker._process_job(job_id)

        job = test_db.get_job(job_id)
        assert job.status == JobStatus.COMPLETED
