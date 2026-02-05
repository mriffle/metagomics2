"""FastAPI server application."""

import os
import shutil
from pathlib import Path
from typing import Annotated

import aiofiles
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from metagomics2 import __version__
from metagomics2.db.database import Database
from metagomics2.models.job import (
    JobCreate,
    JobCreateResponse,
    JobInfo,
    JobListResponse,
    JobParams,
    JobStatus,
)

# Configuration from environment
DATA_DIR = Path(os.environ.get("METAGOMICS_DATA_DIR", "/data"))
JOBS_DIR = DATA_DIR / "jobs"
DB_PATH = DATA_DIR / "metagomics2.db"

# Initialize database
db = Database(DB_PATH)

# Create FastAPI app
app = FastAPI(
    title="Metagomics 2",
    description="Metaproteomics annotation and aggregation tool",
    version=__version__,
)

# CORS middleware for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, restrict this
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "version": __version__}


@app.post("/api/jobs", response_model=JobCreateResponse)
async def create_job(
    fasta: Annotated[UploadFile, File(description="Background proteome FASTA file")],
    peptides: Annotated[list[UploadFile], File(description="Peptide list files")],
    params: Annotated[str, Form()] = "{}",
):
    """Create a new job.

    Upload a background FASTA and one or more peptide list files.
    """
    # Parse parameters
    try:
        import json
        params_dict = json.loads(params)
        job_params = JobParams(**params_dict)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid parameters: {e}")

    # Create job in database
    job_id = db.create_job(job_params)

    # Create job directory structure
    job_dir = JOBS_DIR / job_id
    inputs_dir = job_dir / "inputs"
    peptides_dir = inputs_dir / "peptides"
    work_dir = job_dir / "work"
    results_dir = job_dir / "results"
    logs_dir = job_dir / "logs"

    for d in [inputs_dir, peptides_dir, work_dir, results_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Save FASTA file
    fasta_path = inputs_dir / "background.fasta"
    async with aiofiles.open(fasta_path, "wb") as f:
        content = await fasta.read()
        await f.write(content)

    # Save peptide files
    for i, peptide_file in enumerate(peptides):
        list_id = f"list_{i:03d}"
        filename = peptide_file.filename or f"peptides_{i}.tsv"
        safe_filename = f"{list_id}_{filename}"
        peptide_path = peptides_dir / safe_filename

        async with aiofiles.open(peptide_path, "wb") as f:
            content = await peptide_file.read()
            await f.write(content)

        # Register in database
        db.add_peptide_list(job_id, list_id, filename, str(peptide_path))

    # Queue the job
    db.update_job_status(job_id, JobStatus.QUEUED)
    db.update_job_progress(job_id, 0, len(peptides))

    return JobCreateResponse(job_id=job_id, status=JobStatus.QUEUED)


@app.get("/api/jobs/{job_id}", response_model=JobInfo)
async def get_job(job_id: str):
    """Get job status and information."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@app.get("/api/jobs", response_model=JobListResponse)
async def list_jobs(limit: int = 100):
    """List recent jobs."""
    jobs = db.list_jobs(limit)
    return JobListResponse(jobs=jobs)


@app.get("/api/jobs/{job_id}/peptide-lists")
async def get_peptide_lists(job_id: str):
    """Get peptide lists for a job."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"peptide_lists": job.peptide_lists}


@app.get("/api/jobs/{job_id}/results/{list_id}/{filename}")
async def download_result(job_id: str, list_id: str, filename: str):
    """Download a result file."""
    # Validate job exists
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Validate filename (prevent path traversal)
    allowed_files = [
        "taxonomy_nodes.csv",
        "go_terms.csv",
        "coverage.csv",
        "run_manifest.json",
        "peptides_annotated.csv",
    ]
    if filename not in allowed_files:
        raise HTTPException(status_code=400, detail="Invalid filename")

    file_path = JOBS_DIR / job_id / "results" / list_id / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(
        file_path,
        filename=filename,
        media_type="application/octet-stream",
    )


@app.get("/api/jobs/{job_id}/results/all_results.zip")
async def download_all_results(job_id: str):
    """Download all results as a ZIP file."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed")

    results_dir = JOBS_DIR / job_id / "results"
    zip_path = results_dir / "all_results.zip"

    # Create ZIP if it doesn't exist
    if not zip_path.exists():
        shutil.make_archive(
            str(zip_path.with_suffix("")),
            "zip",
            results_dir,
        )

    return FileResponse(
        zip_path,
        filename=f"metagomics2_results_{job_id[:8]}.zip",
        media_type="application/zip",
    )


# Mount static files for frontend (if exists)
FRONTEND_DIR = Path(__file__).parent.parent.parent.parent / "frontend" / "dist"
if FRONTEND_DIR.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIR), html=True), name="frontend")
