"""Background worker for processing jobs."""
import logging
import shutil
import signal
import time
from types import FrameType
from typing import Any

from metagomics2 import __version__
from metagomics2.config import get_settings
from metagomics2.core.filtering import FilterPolicy
from metagomics2.db.database import Database
from metagomics2.logging_setup import (
    attach_file_handler,
    configure_logging,
    detach_handler,
    format_bytes,
    log_system_info,
)
from metagomics2.models.job import JobInfo, JobStatus, PeptideListStatus
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
DATABASES: list[dict[str, Any]] = _cfg.databases_as_dicts

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

    def _handle_signal(self, signum: int, frame: FrameType | None) -> None:
        """Handle shutdown signals."""
        logger.info(f"Received signal {signum}, shutting down...")
        self.running = False

    def run(self) -> None:
        """Main worker loop."""
        logger.info("Worker started")
        self._recover_orphaned_jobs()

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

    def _recover_orphaned_jobs(self) -> None:
        """Fail any job left in 'running' state by a worker that died mid-job.

        Without this, a job whose worker was killed (for example by the kernel
        out-of-memory killer) would stay 'running' forever with no error.
        """
        for job_id in self.db.list_job_ids_by_status(JobStatus.RUNNING):
            message = (
                "Worker restarted while this job was running. The previous worker "
                "process probably died before finishing (killed by the out-of-memory "
                "killer, a container memory limit, or a container restart). Check "
                "the worker log and the job's logs/ directory for the last stage reached."
            )
            logger.warning(f"Job {job_id} was left in 'running' state; marking it failed")
            self.db.update_job_status(job_id, JobStatus.FAILED, message)
            self.db.add_event(job_id, "error", message)
            self._send_notification(job_id)

    def _process_job(self, job_id: str) -> None:
        """Process a single job."""
        self.current_job_id = job_id
        job_log_path = JOBS_DIR / job_id / "logs" / "pipeline.log"
        job_log_handler = attach_file_handler(job_log_path)
        logger.info(f"Processing job {job_id} (job log: {job_log_path})")
        started = time.monotonic()

        try:
            # Mark as running
            self.db.update_job_status(job_id, JobStatus.RUNNING)
            self.db.add_event(job_id, "started", "Job processing started")

            # Get job info
            job = self.db.get_job(job_id)
            if not job:
                raise ValueError(f"Job {job_id} not found")

            # Build pipeline config
            config = self._build_config(job_id, job)
            fasta_size = config.fasta_path.stat().st_size if config.fasta_path.exists() else 0
            logger.info(
                f"Job {job_id}: {len(config.peptide_list_paths)} peptide list(s), "
                f"FASTA {format_bytes(fasta_size)}, "
                f"database {job.params.db_name or job.params.db_choice!r}, "
                f"filters {config.filter_policy.to_dict()}, "
                f"notify {job.params.notification_email or '(none)'}"
            )

            # Create progress callback. Every stage change is also recorded as a
            # job event so the job's history shows a timestamped timeline.
            last_stage = ""

            def progress_callback(progress: PipelineProgress) -> None:
                nonlocal last_stage
                self.db.update_job_progress(
                    job_id,
                    progress.progress_done,
                    progress.progress_total,
                    progress.current_stage,
                )
                stage_label = progress.current_stage
                if progress.current_list_id:
                    stage_label += f" ({progress.current_list_id})"
                if stage_label != last_stage:
                    last_stage = stage_label
                    self.db.add_event(job_id, "stage", stage_label)

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
                logger.info(
                    f"Job {job_id} completed successfully in {time.monotonic() - started:.0f}s"
                )
            else:
                self.db.update_job_status(
                    job_id, JobStatus.FAILED, result.error_message
                )
                self.db.add_event(job_id, "failed", f"Job failed: {result.error_message}")
                logger.error(
                    f"Job {job_id} failed after {time.monotonic() - started:.0f}s: "
                    f"{result.error_message}"
                )

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
            detach_handler(job_log_handler)

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

    def _build_config(self, job_id: str, job: JobInfo) -> PipelineConfig:
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
            work_dir=job_dir / "work",
            go_edge_types=set(params.go_edge_types.split(",")),
            go_include_self=params.go_include_self,
            diamond_block_size=_cfg.diamond_block_size,
            diamond_index_chunks=_cfg.diamond_index_chunks,
            diamond_tmpdir=_cfg.diamond_tmpdir,
            diamond_max_target_seqs=_cfg.diamond_max_target_seqs,
        )


def main() -> None:
    """Main entry point for worker."""
    configure_logging("worker", _cfg.logs_dir, _cfg.log_level)
    logger.info(f"Metagomics 2 worker v{__version__} starting")
    log_system_info(logger)
    logger.info(
        f"Config: data_dir={_cfg.data_dir}, jobs_dir={JOBS_DIR}, databases_dir={DATABASES_DIR}, "
        f"threads={THREADS}, poll_interval={POLL_INTERVAL}s, log_level={_cfg.log_level}, "
        f"cleanup_on_success={CLEANUP_ON_SUCCESS}, cleanup_on_failure={CLEANUP_ON_FAILURE}"
    )
    logger.info(f"Configured databases: {[d.get('name') for d in DATABASES]}")
    block_size = _cfg.diamond_block_size
    index_chunks = _cfg.diamond_index_chunks
    logger.info(
        "DIAMOND tuning: "
        f"block_size={block_size if block_size is not None else 'default (2.0)'}, "
        f"index_chunks={index_chunks if index_chunks is not None else 'default (4)'}, "
        f"tmpdir={_cfg.diamond_tmpdir or 'default (job work dir)'}, "
        f"max_target_seqs={_cfg.diamond_max_target_seqs or 'unlimited'}"
    )
    logger.info(
        f"Worker log: {_cfg.logs_dir / 'worker.log'}; per-job logs: {JOBS_DIR}/<job_id>/logs/"
    )

    db = Database(DB_PATH)
    worker = Worker(db)
    try:
        worker.run()
    except BaseException:
        logger.exception("Worker exiting because of an unhandled error")
        raise


if __name__ == "__main__":
    main()
