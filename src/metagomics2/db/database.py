"""SQLite database operations for job tracking."""

import json
import secrets
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Generator

from metagomics2.models.job import (
    JobInfo,
    JobParams,
    JobStatus,
    PeptideListInfo,
    PeptideListStatus,
)

# Token length for job IDs (128-bit entropy = 16 bytes = 32 hex chars)
JOB_ID_BYTES = 16


def generate_job_id() -> str:
    """Generate a cryptographically secure job ID.

    Returns:
        URL-safe token with at least 128-bit entropy
    """
    return secrets.token_urlsafe(JOB_ID_BYTES)


class Database:
    """SQLite database for job tracking."""

    def __init__(self, db_path: Path | str):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _init_db(self) -> None:
        """Initialize database schema."""
        with self._get_connection() as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    created_at TEXT NOT NULL,
                    status TEXT NOT NULL DEFAULT 'uploaded',
                    params_json TEXT,
                    db_choice TEXT,
                    search_tool TEXT,
                    progress_total INTEGER DEFAULT 0,
                    progress_done INTEGER DEFAULT 0,
                    current_step TEXT,
                    error_message TEXT
                );

                CREATE TABLE IF NOT EXISTS peptide_lists (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    job_id TEXT NOT NULL,
                    list_id TEXT NOT NULL,
                    filename TEXT NOT NULL,
                    path TEXT NOT NULL,
                    status TEXT NOT NULL DEFAULT 'pending',
                    n_peptides INTEGER,
                    n_matched INTEGER,
                    n_unmatched INTEGER,
                    FOREIGN KEY (job_id) REFERENCES jobs(job_id),
                    UNIQUE (job_id, list_id)
                );

                CREATE TABLE IF NOT EXISTS job_events (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    job_id TEXT NOT NULL,
                    timestamp TEXT NOT NULL,
                    event_type TEXT NOT NULL,
                    message TEXT,
                    FOREIGN KEY (job_id) REFERENCES jobs(job_id)
                );

                CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
                CREATE INDEX IF NOT EXISTS idx_peptide_lists_job ON peptide_lists(job_id);
                CREATE INDEX IF NOT EXISTS idx_events_job ON job_events(job_id);
                """
            )

    @contextmanager
    def _get_connection(self) -> Generator[sqlite3.Connection, None, None]:
        """Get a database connection."""
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    def create_job(self, params: JobParams) -> str:
        """Create a new job.

        Args:
            params: Job parameters

        Returns:
            Generated job ID
        """
        job_id = generate_job_id()
        created_at = datetime.now(timezone.utc).isoformat()

        with self._get_connection() as conn:
            conn.execute(
                """
                INSERT INTO jobs (job_id, created_at, status, params_json, db_choice, search_tool)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    job_id,
                    created_at,
                    JobStatus.UPLOADED.value,
                    params.model_dump_json(),
                    params.db_choice,
                    params.search_tool,
                ),
            )

        return job_id

    def get_job(self, job_id: str) -> JobInfo | None:
        """Get job information.

        Args:
            job_id: Job ID

        Returns:
            JobInfo or None if not found
        """
        with self._get_connection() as conn:
            row = conn.execute(
                "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()

            if not row:
                return None

            # Get peptide lists
            lists_rows = conn.execute(
                "SELECT * FROM peptide_lists WHERE job_id = ? ORDER BY list_id",
                (job_id,),
            ).fetchall()

            peptide_lists = [
                PeptideListInfo(
                    list_id=r["list_id"],
                    filename=r["filename"],
                    status=PeptideListStatus(r["status"]),
                    n_peptides=r["n_peptides"],
                    n_matched=r["n_matched"],
                    n_unmatched=r["n_unmatched"],
                )
                for r in lists_rows
            ]

            params = JobParams.model_validate_json(row["params_json"]) if row["params_json"] else JobParams()

            return JobInfo(
                job_id=row["job_id"],
                created_at=datetime.fromisoformat(row["created_at"]),
                status=JobStatus(row["status"]),
                params=params,
                progress_total=row["progress_total"] or 0,
                progress_done=row["progress_done"] or 0,
                current_step=row["current_step"],
                error_message=row["error_message"],
                peptide_lists=peptide_lists,
            )

    def update_job_status(
        self,
        job_id: str,
        status: JobStatus,
        error_message: str | None = None,
    ) -> None:
        """Update job status.

        Args:
            job_id: Job ID
            status: New status
            error_message: Optional error message
        """
        with self._get_connection() as conn:
            conn.execute(
                "UPDATE jobs SET status = ?, error_message = ? WHERE job_id = ?",
                (status.value, error_message, job_id),
            )

    def update_job_progress(
        self,
        job_id: str,
        progress_done: int,
        progress_total: int | None = None,
        current_step: str | None = None,
    ) -> None:
        """Update job progress.

        Args:
            job_id: Job ID
            progress_done: Number of completed items
            progress_total: Total number of items (optional)
            current_step: Current step description (optional)
        """
        with self._get_connection() as conn:
            if progress_total is not None:
                conn.execute(
                    """
                    UPDATE jobs
                    SET progress_done = ?, progress_total = ?, current_step = ?
                    WHERE job_id = ?
                    """,
                    (progress_done, progress_total, current_step, job_id),
                )
            else:
                conn.execute(
                    "UPDATE jobs SET progress_done = ?, current_step = ? WHERE job_id = ?",
                    (progress_done, current_step, job_id),
                )

    def add_peptide_list(
        self,
        job_id: str,
        list_id: str,
        filename: str,
        path: str,
    ) -> None:
        """Add a peptide list to a job.

        Args:
            job_id: Job ID
            list_id: List identifier
            filename: Original filename
            path: Path to stored file
        """
        with self._get_connection() as conn:
            conn.execute(
                """
                INSERT INTO peptide_lists (job_id, list_id, filename, path, status)
                VALUES (?, ?, ?, ?, ?)
                """,
                (job_id, list_id, filename, path, PeptideListStatus.PENDING.value),
            )

    def update_peptide_list_status(
        self,
        job_id: str,
        list_id: str,
        status: PeptideListStatus,
        n_peptides: int | None = None,
        n_matched: int | None = None,
        n_unmatched: int | None = None,
    ) -> None:
        """Update peptide list status.

        Args:
            job_id: Job ID
            list_id: List identifier
            status: New status
            n_peptides: Number of peptides (optional)
            n_matched: Number of matched peptides (optional)
            n_unmatched: Number of unmatched peptides (optional)
        """
        with self._get_connection() as conn:
            conn.execute(
                """
                UPDATE peptide_lists
                SET status = ?, n_peptides = ?, n_matched = ?, n_unmatched = ?
                WHERE job_id = ? AND list_id = ?
                """,
                (status.value, n_peptides, n_matched, n_unmatched, job_id, list_id),
            )

    def get_next_queued_job(self) -> str | None:
        """Get the next queued job ID.

        Returns:
            Job ID or None if no queued jobs
        """
        with self._get_connection() as conn:
            row = conn.execute(
                """
                SELECT job_id FROM jobs
                WHERE status = ?
                ORDER BY created_at ASC
                LIMIT 1
                """,
                (JobStatus.QUEUED.value,),
            ).fetchone()

            return row["job_id"] if row else None

    def add_event(self, job_id: str, event_type: str, message: str) -> None:
        """Add an event log entry.

        Args:
            job_id: Job ID
            event_type: Event type
            message: Event message
        """
        timestamp = datetime.now(timezone.utc).isoformat()
        with self._get_connection() as conn:
            conn.execute(
                """
                INSERT INTO job_events (job_id, timestamp, event_type, message)
                VALUES (?, ?, ?, ?)
                """,
                (job_id, timestamp, event_type, message),
            )

    def regenerate_job_id(self, old_job_id: str, jobs_dir: Path) -> str:
        """Generate a new job ID and migrate all references.

        Atomically updates the job_id in all database tables and renames
        the job directory on disk.

        Args:
            old_job_id: Current job ID
            jobs_dir: Base directory containing job directories

        Returns:
            The new job ID

        Raises:
            ValueError: If the old job ID does not exist
            OSError: If the directory rename fails
        """
        new_job_id = generate_job_id()

        with self._get_connection() as conn:
            # Verify old job exists
            row = conn.execute(
                "SELECT job_id FROM jobs WHERE job_id = ?", (old_job_id,)
            ).fetchone()
            if not row:
                raise ValueError(f"Job not found: {old_job_id}")

            # Update all tables
            conn.execute(
                "UPDATE jobs SET job_id = ? WHERE job_id = ?",
                (new_job_id, old_job_id),
            )
            conn.execute(
                "UPDATE peptide_lists SET job_id = ? WHERE job_id = ?",
                (new_job_id, old_job_id),
            )
            conn.execute(
                "UPDATE job_events SET job_id = ? WHERE job_id = ?",
                (new_job_id, old_job_id),
            )

        # Rename directory on disk (if it exists)
        old_dir = jobs_dir / old_job_id
        if old_dir.exists():
            new_dir = jobs_dir / new_job_id
            old_dir.rename(new_dir)

        return new_job_id

    def list_jobs(self, limit: int = 100) -> list[JobInfo]:
        """List recent jobs.

        Args:
            limit: Maximum number of jobs to return

        Returns:
            List of JobInfo objects
        """
        jobs = []
        with self._get_connection() as conn:
            rows = conn.execute(
                "SELECT job_id FROM jobs ORDER BY created_at DESC LIMIT ?",
                (limit,),
            ).fetchall()

        for row in rows:
            job = self.get_job(row["job_id"])
            if job:
                jobs.append(job)

        return jobs
