"""Unit tests for SQLite database operations."""

import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from metagomics2.db.database import Database, generate_job_id
from metagomics2.models.job import JobParams, JobStatus, PeptideListStatus


class TestGenerateJobId:
    """Tests for job ID generation."""

    def test_generates_string(self):
        job_id = generate_job_id()
        assert isinstance(job_id, str)

    def test_generates_unique_ids(self):
        ids = {generate_job_id() for _ in range(100)}
        assert len(ids) == 100

    def test_url_safe(self):
        job_id = generate_job_id()
        # URL-safe base64 uses only these characters
        allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_=")
        assert all(c in allowed for c in job_id)

    def test_sufficient_length(self):
        job_id = generate_job_id()
        # 16 bytes -> ~22 base64 chars
        assert len(job_id) >= 20


class TestDatabaseInit:
    """Tests for database initialization."""

    def test_creates_db_file(self, tmp_path: Path):
        db_path = tmp_path / "test.db"
        db = Database(db_path)
        assert db_path.exists()

    def test_creates_parent_directories(self, tmp_path: Path):
        db_path = tmp_path / "subdir" / "nested" / "test.db"
        db = Database(db_path)
        assert db_path.exists()

    def test_idempotent_init(self, tmp_path: Path):
        db_path = tmp_path / "test.db"
        db1 = Database(db_path)
        db2 = Database(db_path)  # Should not raise
        assert db_path.exists()


class TestCreateJob:
    """Tests for job creation."""

    def test_creates_job_returns_id(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams()
        job_id = db.create_job(params)
        assert isinstance(job_id, str)
        assert len(job_id) > 0

    def test_job_initial_status_is_uploaded(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams()
        job_id = db.create_job(params)
        job = db.get_job(job_id)
        assert job is not None
        assert job.status == JobStatus.UPLOADED

    def test_job_stores_params(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams(
            search_tool="blast",
            max_evalue=1e-5,
            min_pident=80.0,
            top_k=10,
        )
        job_id = db.create_job(params)
        job = db.get_job(job_id)
        assert job.params.search_tool == "blast"
        assert job.params.max_evalue == 1e-5
        assert job.params.min_pident == 80.0
        assert job.params.top_k == 10

    def test_job_has_created_at(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams()
        job_id = db.create_job(params)
        job = db.get_job(job_id)
        assert job.created_at is not None

    def test_multiple_jobs_unique_ids(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams()
        ids = [db.create_job(params) for _ in range(10)]
        assert len(set(ids)) == 10


class TestGetJob:
    """Tests for job retrieval."""

    def test_returns_none_for_nonexistent(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        assert db.get_job("nonexistent_id") is None

    def test_returns_job_info(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams(search_tool="diamond")
        job_id = db.create_job(params)
        job = db.get_job(job_id)
        assert job is not None
        assert job.job_id == job_id
        assert job.params.search_tool == "diamond"

    def test_includes_peptide_lists(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        params = JobParams()
        job_id = db.create_job(params)
        db.add_peptide_list(job_id, "list_000", "peptides.tsv", "/path/to/file")
        job = db.get_job(job_id)
        assert len(job.peptide_lists) == 1
        assert job.peptide_lists[0].list_id == "list_000"
        assert job.peptide_lists[0].filename == "peptides.tsv"


class TestUpdateJobStatus:
    """Tests for job status updates."""

    def test_update_to_queued(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_status(job_id, JobStatus.QUEUED)
        job = db.get_job(job_id)
        assert job.status == JobStatus.QUEUED

    def test_update_to_running(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_status(job_id, JobStatus.RUNNING)
        job = db.get_job(job_id)
        assert job.status == JobStatus.RUNNING

    def test_update_to_completed(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_status(job_id, JobStatus.COMPLETED)
        job = db.get_job(job_id)
        assert job.status == JobStatus.COMPLETED

    def test_update_to_failed_with_error(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_status(job_id, JobStatus.FAILED, "Something went wrong")
        job = db.get_job(job_id)
        assert job.status == JobStatus.FAILED
        assert job.error_message == "Something went wrong"


class TestUpdateJobProgress:
    """Tests for job progress updates."""

    def test_update_progress(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_progress(job_id, 2, 5, "Processing list 2")
        job = db.get_job(job_id)
        assert job.progress_done == 2
        assert job.progress_total == 5
        assert job.current_step == "Processing list 2"

    def test_update_progress_without_total(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_progress(job_id, 0, 3)
        db.update_job_progress(job_id, 1, current_step="Step 1")
        job = db.get_job(job_id)
        assert job.progress_done == 1
        assert job.progress_total == 3
        assert job.current_step == "Step 1"


class TestPeptideLists:
    """Tests for peptide list management."""

    def test_add_peptide_list(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.add_peptide_list(job_id, "list_000", "peptides_1.tsv", "/data/jobs/x/peptides_1.tsv")
        job = db.get_job(job_id)
        assert len(job.peptide_lists) == 1
        assert job.peptide_lists[0].status == PeptideListStatus.PENDING

    def test_add_multiple_peptide_lists(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.add_peptide_list(job_id, "list_000", "pep1.tsv", "/path/pep1.tsv")
        db.add_peptide_list(job_id, "list_001", "pep2.tsv", "/path/pep2.tsv")
        job = db.get_job(job_id)
        assert len(job.peptide_lists) == 2
        assert job.peptide_lists[0].list_id == "list_000"
        assert job.peptide_lists[1].list_id == "list_001"

    def test_update_peptide_list_status(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.add_peptide_list(job_id, "list_000", "pep.tsv", "/path/pep.tsv")
        db.update_peptide_list_status(
            job_id, "list_000", PeptideListStatus.DONE,
            n_peptides=100, n_matched=80, n_unmatched=20,
        )
        job = db.get_job(job_id)
        pl = job.peptide_lists[0]
        assert pl.status == PeptideListStatus.DONE
        assert pl.n_peptides == 100
        assert pl.n_matched == 80
        assert pl.n_unmatched == 20

    def test_unique_list_id_per_job(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.add_peptide_list(job_id, "list_000", "pep.tsv", "/path/pep.tsv")
        with pytest.raises(Exception):
            db.add_peptide_list(job_id, "list_000", "pep2.tsv", "/path/pep2.tsv")


class TestGetNextQueuedJob:
    """Tests for job queue retrieval."""

    def test_returns_none_when_empty(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        assert db.get_next_queued_job() is None

    def test_returns_none_when_no_queued(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        # Status is 'uploaded', not 'queued'
        assert db.get_next_queued_job() is None

    def test_returns_queued_job(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        db.update_job_status(job_id, JobStatus.QUEUED)
        assert db.get_next_queued_job() == job_id

    def test_returns_oldest_queued_job(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        id1 = db.create_job(JobParams())
        id2 = db.create_job(JobParams())
        db.update_job_status(id1, JobStatus.QUEUED)
        db.update_job_status(id2, JobStatus.QUEUED)
        assert db.get_next_queued_job() == id1

    def test_skips_running_jobs(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        id1 = db.create_job(JobParams())
        id2 = db.create_job(JobParams())
        db.update_job_status(id1, JobStatus.RUNNING)
        db.update_job_status(id2, JobStatus.QUEUED)
        assert db.get_next_queued_job() == id2


class TestAddEvent:
    """Tests for event logging."""

    def test_add_event(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        job_id = db.create_job(JobParams())
        # Should not raise
        db.add_event(job_id, "started", "Job processing started")
        db.add_event(job_id, "completed", "Job completed successfully")


class TestListJobs:
    """Tests for listing jobs."""

    def test_list_empty(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        jobs = db.list_jobs()
        assert jobs == []

    def test_list_returns_jobs(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        id1 = db.create_job(JobParams())
        id2 = db.create_job(JobParams())
        jobs = db.list_jobs()
        assert len(jobs) == 2

    def test_list_ordered_by_created_at_desc(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        id1 = db.create_job(JobParams())
        id2 = db.create_job(JobParams())
        jobs = db.list_jobs()
        # Most recent first
        assert jobs[0].job_id == id2
        assert jobs[1].job_id == id1

    def test_list_respects_limit(self, tmp_path: Path):
        db = Database(tmp_path / "test.db")
        for _ in range(5):
            db.create_job(JobParams())
        jobs = db.list_jobs(limit=3)
        assert len(jobs) == 3
