"""Background worker for processing jobs."""

import logging
import shutil
import signal
import time
from pathlib import Path

from metagomics2.config import get_settings
from metagomics2.core.filtering import FilterPolicy
from metagomics2.db.database import Database
from metagomics2.models.job import JobStatus, PeptideListStatus
from metagomics2.notifications.email import SmtpConfig, send_job_notification
from metagomics2.pipeline.runner import PipelineConfig, PipelineProgress, run_pipeline

logger = logging.getLogger(__name__)

# Load validated settings from centralized config
_cfg = get_settings()

JOBS_DIR = _cfg.jobs_dir
DB_PATH = _cfg.db_path
POLL_INTERVAL = _cfg.poll_interval
THREADS = _cfg.threads
DATABASES_DIR = _cfg.databases_dir
DATABASES: list[dict] = _cfg.databases_as_dicts

# Email notification settings (from centralized config)
SMTP_CONFIG = SmtpConfig(
    host=_cfg.smtp.host,
    port=_cfg.smtp.port,
    username=_cfg.smtp.username,
    password=_cfg.smtp.password,
    from_address=_cfg.smtp.from_address,
)
SITE_URL = _cfg.site_url

# Cleanup settings
CLEANUP_ON_SUCCESS = _cfg.cleanup_on_success
CLEANUP_ON_FAILURE = _cfg.cleanup_on_failure


class Worker:
    """Background worker that processes queued jobs."""

    def __init__(self, db: Database):
        self.db = db
        self.running = True
        self.current_job_id: str | None = None

        # Setup signal handlers
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)

    def _handle_signal(self, signum, frame):
        """Handle shutdown signals."""
        logger.info(f"Received signal {signum}, shutting down...")
        self.running = False

    def run(self) -> None:
        """Main worker loop."""
        logger.info("Worker started")

        while self.running:
            try:
                job_id = self.db.get_next_queued_job()

                if job_id:
                    self._process_job(job_id)
                else:
                    time.sleep(POLL_INTERVAL)

            except Exception as e:
                logger.exception(f"Error in worker loop: {e}")
                time.sleep(POLL_INTERVAL)

        logger.info("Worker stopped")

    def _process_job(self, job_id: str) -> None:
        """Process a single job."""
        self.current_job_id = job_id
        logger.info(f"Processing job {job_id}")

        try:
            # Mark as running
            self.db.update_job_status(job_id, JobStatus.RUNNING)
            self.db.add_event(job_id, "started", "Job processing started")

            # Get job info
            job = self.db.get_job(job_id)
            if not job:
                raise ValueError(f"Job {job_id} not found")

            # Build pipeline config
            job_dir = JOBS_DIR / job_id
            config = self._build_config(job_id, job)

            # Create progress callback
            def progress_callback(progress: PipelineProgress) -> None:
                self.db.update_job_progress(
                    job_id,
                    progress.progress_done,
                    progress.progress_total,
                    progress.current_stage,
                )

            # Run pipeline
            result = run_pipeline(config, progress_callback)

            if result.success:
                # Update per-list status with results
                for pl_result in result.peptide_list_results:
                    self.db.update_peptide_list_status(
                        job_id,
                        pl_result.list_id,
                        PeptideListStatus.DONE,
                        n_peptides=pl_result.n_peptides,
                        n_matched=pl_result.n_matched,
                        n_unmatched=pl_result.n_unmatched,
                    )

                self.db.update_job_status(job_id, JobStatus.COMPLETED)
                self.db.add_event(job_id, "completed", "Job completed successfully")
                logger.info(f"Job {job_id} completed successfully")
            else:
                self.db.update_job_status(
                    job_id, JobStatus.FAILED, result.error_message
                )
                self.db.add_event(job_id, "failed", f"Job failed: {result.error_message}")
                logger.error(f"Job {job_id} failed: {result.error_message}")

            # Send email notification (re-fetch job to get final status)
            self._send_notification(job_id)

            # Clean up intermediate files
            if result.success and CLEANUP_ON_SUCCESS:
                self._cleanup_job_files(job_id)
            elif not result.success and CLEANUP_ON_FAILURE:
                self._cleanup_job_files(job_id)

        except Exception as e:
            logger.exception(f"Error processing job {job_id}")
            self.db.update_job_status(job_id, JobStatus.FAILED, str(e))
            self.db.add_event(job_id, "error", str(e))
            self._send_notification(job_id)
            if CLEANUP_ON_FAILURE:
                self._cleanup_job_files(job_id)

        finally:
            self.current_job_id = None

    def _cleanup_job_files(self, job_id: str) -> None:
        """Remove inputs/ and work/ directories to free disk space."""
        try:
            job_dir = JOBS_DIR / job_id
            for subdir in ("inputs", "work"):
                path = job_dir / subdir
                if path.exists():
                    shutil.rmtree(path)
                    logger.info(f"Cleaned up {path}")
        except Exception:
            logger.exception(f"Error cleaning up files for job {job_id}")

    def _send_notification(self, job_id: str) -> None:
        """Send email notification for a finished job."""
        try:
            job = self.db.get_job(job_id)
            if job and job.params.notification_email:
                send_job_notification(job, SITE_URL, SMTP_CONFIG)
        except Exception:
            logger.exception(f"Error sending notification for job {job_id}")

    def _build_config(self, job_id: str, job) -> PipelineConfig:
        """Build pipeline configuration from job info."""
        job_dir = JOBS_DIR / job_id
        params = job.params

        # Get peptide list paths from database
        peptide_paths = []
        for pl in job.peptide_lists:
            # The path is stored in the database, but we need to reconstruct it
            peptide_path = job_dir / "inputs" / "peptides" / f"{pl.list_id}_{pl.filename}"
            if peptide_path.exists():
                peptide_paths.append(peptide_path)
            else:
                # Try alternative path format
                for f in (job_dir / "inputs" / "peptides").iterdir():
                    if f.name.startswith(pl.list_id):
                        peptide_paths.append(f)
                        break

        # Build filter policy
        filter_policy = FilterPolicy(
            max_evalue=params.max_evalue,
            min_pident=params.min_pident,
            min_qcov=params.min_qcov,
            min_alnlen=params.min_alnlen,
            top_k=params.top_k,
        )

        # Resolve database path: db_choice is relative to DATABASES_DIR
        annotated_db_path = None
        annotations_db_path = None
        if params.db_choice:
            annotated_db_path = DATABASES_DIR / params.db_choice
            # Look up companion annotations DB from database config
            for db_entry in DATABASES:
                if db_entry.get("path") == params.db_choice:
                    ann_path = db_entry.get("annotations")
                    if ann_path:
                        annotations_db_path = DATABASES_DIR / ann_path
                    break

        return PipelineConfig(
            fasta_path=job_dir / "inputs" / "background.fasta",
            peptide_list_paths=peptide_paths,
            output_dir=job_dir / "results",
            search_tool=params.search_tool,
            annotated_db_path=annotated_db_path,
            annotations_db_path=annotations_db_path,
            threads=THREADS,
            filter_policy=filter_policy,
            job_dir=job_dir,  # Enable reference snapshot creation
            go_edge_types=set(params.go_edge_types.split(",")),
            go_include_self=params.go_include_self,
            compute_enrichment_pvalues=params.compute_enrichment_pvalues,
            enrichment_iterations=params.enrichment_iterations,
        )


def main() -> None:
    """Main entry point for worker."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    db = Database(DB_PATH)
    worker = Worker(db)
    worker.run()


if __name__ == "__main__":
    main()
