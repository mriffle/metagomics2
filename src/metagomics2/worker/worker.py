"""Background worker for processing jobs."""

import logging
import os
import signal
import time
from pathlib import Path

from metagomics2.core.filtering import FilterPolicy
from metagomics2.db.database import Database
from metagomics2.models.job import JobStatus, PeptideListStatus
from metagomics2.pipeline.runner import PipelineConfig, PipelineProgress, run_pipeline

logger = logging.getLogger(__name__)

# Configuration from environment
DATA_DIR = Path(os.environ.get("METAGOMICS_DATA_DIR", "/data"))
JOBS_DIR = DATA_DIR / "jobs"
DB_PATH = DATA_DIR / "metagomics2.db"
REFERENCE_DIR = DATA_DIR / "reference"

# Worker settings
POLL_INTERVAL = int(os.environ.get("METAGOMICS_POLL_INTERVAL", "5"))


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
                    progress.completed_peptide_lists,
                    progress.total_peptide_lists,
                    progress.current_stage,
                )

            # Run pipeline
            result = run_pipeline(config, progress_callback)

            if result.success:
                self.db.update_job_status(job_id, JobStatus.COMPLETED)
                self.db.add_event(job_id, "completed", "Job completed successfully")
                logger.info(f"Job {job_id} completed successfully")
            else:
                self.db.update_job_status(
                    job_id, JobStatus.FAILED, result.error_message
                )
                self.db.add_event(job_id, "failed", f"Job failed: {result.error_message}")
                logger.error(f"Job {job_id} failed: {result.error_message}")

        except Exception as e:
            logger.exception(f"Error processing job {job_id}")
            self.db.update_job_status(job_id, JobStatus.FAILED, str(e))
            self.db.add_event(job_id, "error", str(e))

        finally:
            self.current_job_id = None

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
            delta_bitscore=params.delta_bitscore,
            best_hit_only=params.best_hit_only,
        )

        return PipelineConfig(
            fasta_path=job_dir / "inputs" / "background.fasta",
            peptide_list_paths=peptide_paths,
            output_dir=job_dir / "results",
            search_tool=params.search_tool,
            annotated_db_path=Path(params.db_choice) if params.db_choice else None,
            filter_policy=filter_policy,
            job_dir=job_dir,  # Enable reference snapshot creation
            go_edge_types=set(params.go_edge_types.split(",")),
            go_include_self=params.go_include_self,
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
