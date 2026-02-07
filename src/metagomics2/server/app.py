"""FastAPI server application."""

import hashlib
import json
import os
import secrets
import shutil
from pathlib import Path
from typing import Annotated

import aiofiles
from fastapi import Depends, FastAPI, File, Form, Header, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

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
ADMIN_PASSWORD = os.environ.get("METAGOMICS_ADMIN_PASSWORD", "")
DIAMOND_VERSION = os.environ.get("DIAMOND_VERSION", "")
THREADS = int(os.environ.get("METAGOMICS_THREADS", "4"))
DATABASES_DIR = Path(os.environ.get("METAGOMICS_DATABASES_DIR", "/databases"))

# Parse annotated databases configuration from JSON env var
_databases_raw = os.environ.get("METAGOMICS_DATABASES", "[]")
try:
    DATABASES: list[dict] = json.loads(_databases_raw)
except (json.JSONDecodeError, TypeError):
    DATABASES = []

if not DATABASES:
    raise RuntimeError(
        "No annotated databases configured. Set METAGOMICS_DATABASES in your "
        "environment to a JSON array with at least one database entry. "
        'Example: METAGOMICS_DATABASES=\'[{"name": "UniProt SwissProt", '
        '"description": "Reviewed entries", "path": "uniprot_sprot.dmnd", '
        '"annotations": "uniprot_sprot.annotations.db"}]\''
    )

# Admin session tokens (in-memory, cleared on restart)
_admin_tokens: set[str] = set()

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


# --- Admin auth ---

class AdminAuthRequest(BaseModel):
    password: str


class AdminAuthResponse(BaseModel):
    token: str


def require_admin(authorization: str = Header(default="")):
    """Dependency that validates admin token from Authorization header."""
    token = authorization.replace("Bearer ", "") if authorization.startswith("Bearer ") else authorization
    if not token or token not in _admin_tokens:
        raise HTTPException(status_code=401, detail="Unauthorized")
    return token


@app.post("/api/admin/auth", response_model=AdminAuthResponse)
async def admin_login(body: AdminAuthRequest):
    """Authenticate with admin password and receive a session token."""
    if not ADMIN_PASSWORD:
        raise HTTPException(status_code=403, detail="Admin access is not configured")
    if not secrets.compare_digest(body.password, ADMIN_PASSWORD):
        raise HTTPException(status_code=401, detail="Invalid password")
    token = secrets.token_urlsafe(32)
    _admin_tokens.add(token)
    return AdminAuthResponse(token=token)


@app.get("/api/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "version": __version__}


@app.get("/api/version")
async def get_version():
    """Get application version."""
    return {"version": __version__}


@app.get("/api/config")
async def get_config():
    """Get public application configuration."""
    return {
        "diamond_version": DIAMOND_VERSION,
        "databases": DATABASES,
    }


def _validate_fasta_content(text: str) -> None:
    """Validate that text looks like a FASTA file.

    Checks that the first non-empty line starts with '>' and that at least
    one sequence line follows.  Raises HTTPException(400) on failure.
    """
    lines = [l for l in text.splitlines() if l.strip()]

    if not lines:
        raise HTTPException(
            status_code=400,
            detail="The uploaded FASTA file is empty.",
        )

    if not lines[0].startswith(">"):
        raise HTTPException(
            status_code=400,
            detail=(
                "The uploaded file does not appear to be a valid FASTA file. "
                "FASTA files must begin with a header line starting with '>'."
            ),
        )

    if len(lines) < 2:
        raise HTTPException(
            status_code=400,
            detail=(
                "The uploaded FASTA file contains a header but no sequence data. "
                "Each header line (starting with '>') must be followed by one or "
                "more lines of amino acid sequence."
            ),
        )

    # Check that the second non-empty line is not another header (i.e. it's sequence)
    if lines[1].startswith(">"):
        raise HTTPException(
            status_code=400,
            detail=(
                "The uploaded FASTA file has consecutive header lines with no "
                "sequence data between them. Each '>' header must be followed "
                "by at least one sequence line."
            ),
        )


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
        params_dict = json.loads(params)
        job_params = JobParams(**params_dict)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid parameters: {e}")

    # Store original FASTA filename
    job_params.fasta_filename = fasta.filename or "background.fasta"

    # Validate db_choice against configured databases
    if job_params.db_choice:
        valid_paths = {db_entry.get("path") for db_entry in DATABASES}
        if job_params.db_choice not in valid_paths:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown database: {job_params.db_choice}",
            )

    # Resolve database name from config
    if job_params.db_choice and not job_params.db_name:
        for db_entry in DATABASES:
            if db_entry.get("path") == job_params.db_choice:
                job_params.db_name = db_entry.get("name", "")
                break

    # Validate FASTA file
    fasta_content = await fasta.read()
    try:
        fasta_text = fasta_content.decode("utf-8", errors="replace")
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="The uploaded FASTA file could not be read as text.",
        )

    _validate_fasta_content(fasta_text)
    await fasta.seek(0)

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


@app.get("/api/admin/jobs", response_model=JobListResponse)
async def list_jobs(limit: int = 100, _token: str = Depends(require_admin)):
    """List recent jobs (admin only)."""
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


# Frontend SPA support
FRONTEND_DIR = Path(__file__).parent.parent.parent.parent / "frontend" / "dist"
if FRONTEND_DIR.exists():
    # Serve static assets (JS, CSS, images, etc.)
    app.mount("/assets", StaticFiles(directory=str(FRONTEND_DIR / "assets")), name="assets")

    @app.get("/{full_path:path}")
    async def serve_spa(full_path: str):
        """Serve the SPA index.html for all non-API routes."""
        # Try to serve the exact file first
        file_path = FRONTEND_DIR / full_path
        if full_path and file_path.is_file():
            return FileResponse(file_path)
        # Fall back to index.html for SPA routing
        return FileResponse(FRONTEND_DIR / "index.html")
