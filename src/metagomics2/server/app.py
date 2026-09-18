"""FastAPI server application."""

import json
import logging
import os
import secrets
import shutil
import tempfile
import time
import zipfile
from pathlib import Path
from typing import Annotated, Any

import aiofiles
from fastapi import Depends, FastAPI, File, Form, Header, HTTPException, UploadFile
from fastapi.concurrency import run_in_threadpool
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from metagomics2 import __version__
from metagomics2.config import get_settings
from metagomics2.db.database import Database
from metagomics2.logging_setup import configure_logging
from metagomics2.models.job import (
    JobCreateResponse,
    JobInfo,
    JobListResponse,
    JobParams,
    JobStatus,
)

# Load validated settings from centralized config
_cfg = get_settings()
configure_logging("server", _cfg.logs_dir, _cfg.log_level)
logger = logging.getLogger(__name__)
logger.info(f"Metagomics 2 server v{__version__} starting (log: {_cfg.logs_dir / 'server.log'})")

DATA_DIR = _cfg.data_dir
JOBS_DIR = _cfg.jobs_dir
DB_PATH = _cfg.db_path
ADMIN_PASSWORD = _cfg.admin_password
DIAMOND_VERSION = _cfg.diamond_version
THREADS = _cfg.threads
DATABASES_DIR = _cfg.databases_dir
MAX_UPLOAD_MB = _cfg.max_upload_mb
MAX_UPLOAD_BYTES = _cfg.max_upload_bytes
DATABASES: list[dict[str, Any]] = _cfg.databases_as_dicts

# Chunk size for streaming file writes (1 MB)
_WRITE_CHUNK_SIZE = 1024 * 1024

# Allowed CORS origins from config
_ALLOWED_ORIGINS = _cfg.allowed_origins

# Admin session tokens (in-memory, cleared on restart): token -> expiry as a
# time.monotonic() timestamp.  Tokens expire so a leaked one is not valid
# forever, and the table is capped so repeated logins cannot grow it without
# bound.
_admin_tokens: dict[str, float] = {}
_ADMIN_TOKEN_TTL_SECONDS = 12 * 60 * 60
_ADMIN_TOKEN_LIMIT = 100


def _prune_admin_tokens(now: float) -> None:
    """Drop expired admin tokens."""
    for token, expires_at in list(_admin_tokens.items()):
        if expires_at <= now:
            del _admin_tokens[token]

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
    allow_origins=_ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# --- Admin auth ---

class AdminAuthRequest(BaseModel):
    password: str


class AdminAuthResponse(BaseModel):
    token: str


def require_admin(authorization: str = Header(default="")) -> str:
    """Dependency that validates admin token from Authorization header."""
    token = (
        authorization.replace("Bearer ", "")
        if authorization.startswith("Bearer ")
        else authorization
    )
    _prune_admin_tokens(time.monotonic())
    if not token or token not in _admin_tokens:
        raise HTTPException(status_code=401, detail="Unauthorized")
    return token


@app.post("/api/admin/auth", response_model=AdminAuthResponse)
async def admin_login(body: AdminAuthRequest) -> AdminAuthResponse:
    """Authenticate with admin password and receive a session token."""
    if not ADMIN_PASSWORD:
        raise HTTPException(status_code=403, detail="Admin access is not configured")
    if not secrets.compare_digest(body.password, ADMIN_PASSWORD):
        raise HTTPException(status_code=401, detail="Invalid password")
    now = time.monotonic()
    _prune_admin_tokens(now)
    # Evict the tokens closest to expiry if the table is full
    while len(_admin_tokens) >= _ADMIN_TOKEN_LIMIT:
        oldest = min(_admin_tokens, key=lambda t: _admin_tokens[t])
        del _admin_tokens[oldest]
    token = secrets.token_urlsafe(32)
    _admin_tokens[token] = now + _ADMIN_TOKEN_TTL_SECONDS
    return AdminAuthResponse(token=token)


@app.get("/api/health")
async def health_check() -> dict[str, str]:
    """Health check endpoint."""
    return {"status": "healthy", "version": __version__}


@app.get("/api/version")
async def get_version() -> dict[str, str]:
    """Get application version."""
    return {"version": __version__}


@app.get("/api/config")
async def get_config() -> dict[str, Any]:
    """Get public application configuration."""
    return {
        "diamond_version": DIAMOND_VERSION,
        "databases": DATABASES,
    }


class UploadTooLargeError(Exception):
    """Raised by :func:`_save_upload_streamed` when an upload passes its byte limit."""


async def _save_upload_streamed(
    upload: UploadFile, dest: Path, max_bytes: int | None = None
) -> int:
    """Stream an uploaded file to disk in chunks, returning total bytes written.

    When ``max_bytes`` is given the stream is abandoned as soon as the running
    total passes it, so a client cannot fill the disk before the size check.
    At most one extra chunk (1 MB) beyond the limit reaches the file, and the
    caller is expected to delete it.
    """
    total = 0
    async with aiofiles.open(dest, "wb") as f:
        while True:
            chunk = await upload.read(_WRITE_CHUNK_SIZE)
            if not chunk:
                break
            total += len(chunk)
            if max_bytes is not None and total > max_bytes:
                raise UploadTooLargeError(total)
            await f.write(chunk)
    return total


def _validate_fasta_content(text: str) -> None:
    """Validate that text looks like a FASTA file.

    Checks that the first non-empty line starts with '>' and that at least
    one sequence line follows.  Raises HTTPException(400) on failure.
    """
    lines = [line for line in text.splitlines() if line.strip()]

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
) -> JobCreateResponse:
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

    # Validate db_choice against configured databases.  An empty choice would
    # be accepted here and only fail in the worker after parsing and matching.
    if not job_params.db_choice:
        raise HTTPException(
            status_code=400,
            detail="No annotated database selected (db_choice is required).",
        )
    valid_paths = {db_entry.get("path") for db_entry in DATABASES}
    if job_params.db_choice not in valid_paths:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown database: {job_params.db_choice}",
        )

    # Resolve database name from config
    if not job_params.db_name:
        for db_entry in DATABASES:
            if db_entry.get("path") == job_params.db_choice:
                job_params.db_name = db_entry.get("name", "")
                break

    # Validate FASTA file (only read first 8 KB for header check)
    fasta_header = await fasta.read(8192)
    if not fasta_header:
        raise HTTPException(status_code=400, detail="The uploaded FASTA file is empty.")
    try:
        fasta_header_text = fasta_header.decode("utf-8", errors="replace")
    except Exception:
        raise HTTPException(
            status_code=400,
            detail="The uploaded FASTA file could not be read as text.",
        )

    _validate_fasta_content(fasta_header_text)
    await fasta.seek(0)

    # Create job in database
    job_id = db.create_job(job_params)
    job_dir = JOBS_DIR / job_id

    # Anything that goes wrong from here until the job is queued (an oversized
    # upload, a client that disconnects mid-stream, a disk error) must not
    # leave a half-written job directory or an orphan "uploaded" row behind.
    try:
        await _store_job_inputs(job_id, job_dir, fasta, peptides)
    except UploadTooLargeError as e:
        shutil.rmtree(job_dir, ignore_errors=True)
        db.delete_job(job_id)
        raise HTTPException(status_code=413, detail=str(e))
    except BaseException:
        shutil.rmtree(job_dir, ignore_errors=True)
        db.delete_job(job_id)
        raise

    # Queue the job
    db.update_job_status(job_id, JobStatus.QUEUED)
    db.update_job_progress(job_id, 0, 1000)

    return JobCreateResponse(job_id=job_id, status=JobStatus.QUEUED)


async def _store_job_inputs(
    job_id: str, job_dir: Path, fasta: UploadFile, peptides: list[UploadFile]
) -> None:
    """Create the job directory, stream the uploads into it and register the lists.

    Raises:
        UploadTooLargeError: with a user-facing message, when the FASTA or the
            peptide files together exceed the configured upload limit.
    """
    inputs_dir = job_dir / "inputs"
    peptides_dir = inputs_dir / "peptides"
    work_dir = job_dir / "work"
    results_dir = job_dir / "results"
    logs_dir = job_dir / "logs"

    for d in [inputs_dir, peptides_dir, work_dir, results_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Save FASTA file (streamed in chunks, abandoned once it passes the limit)
    fasta_path = inputs_dir / "background.fasta"
    try:
        await _save_upload_streamed(fasta, fasta_path, MAX_UPLOAD_BYTES)
    except UploadTooLargeError:
        raise UploadTooLargeError(
            f"FASTA file exceeds the maximum upload size of {MAX_UPLOAD_MB} MB."
        )

    # Save peptide files; the limit applies to their combined size
    total_peptide_size = 0
    for i, peptide_file in enumerate(peptides):
        list_id = f"list_{i:03d}"
        filename = peptide_file.filename or f"peptides_{i}.tsv"
        safe_filename = f"{list_id}_{filename}"
        peptide_path = peptides_dir / safe_filename

        try:
            file_size = await _save_upload_streamed(
                peptide_file, peptide_path, MAX_UPLOAD_BYTES - total_peptide_size
            )
        except UploadTooLargeError:
            raise UploadTooLargeError(
                "Total peptide file size exceeds the maximum upload size of "
                f"{MAX_UPLOAD_MB} MB."
            )
        total_peptide_size += file_size

        # Register in database
        db.add_peptide_list(job_id, list_id, filename, str(peptide_path))


@app.get("/api/jobs/{job_id}", response_model=JobInfo)
async def get_job(job_id: str) -> JobInfo:
    """Get job status and information."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


class RegenerateIdResponse(BaseModel):
    new_job_id: str


@app.post("/api/jobs/{job_id}/regenerate-id", response_model=RegenerateIdResponse)
async def regenerate_job_id(job_id: str) -> RegenerateIdResponse:
    """Regenerate the job ID (URL hash) for a job.

    This changes the URL used to access the job, invalidating the old one.
    Useful when a user has shared a link and wants to revoke access.
    """
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # The worker addresses a running job by its ID and directory; renaming
    # either underneath it would fail the pipeline and strand the job.
    if job.status not in (JobStatus.COMPLETED, JobStatus.FAILED):
        raise HTTPException(
            status_code=409,
            detail=(
                "The job ID can only be regenerated once the job has completed or failed "
                f"(current status: {job.status.value})."
            ),
        )

    try:
        new_job_id = db.regenerate_job_id(job_id, JOBS_DIR)
    except ValueError:
        raise HTTPException(status_code=404, detail="Job not found")
    except OSError as e:
        raise HTTPException(status_code=500, detail=f"Failed to rename job directory: {e}")

    return RegenerateIdResponse(new_job_id=new_job_id)


@app.get("/api/admin/jobs", response_model=JobListResponse)
async def list_jobs(limit: int = 100, _token: str = Depends(require_admin)) -> JobListResponse:
    """List recent jobs (admin only)."""
    jobs = db.list_jobs(limit)
    return JobListResponse(jobs=jobs)


@app.get("/api/jobs/{job_id}/peptide-lists")
async def get_peptide_lists(job_id: str) -> dict[str, Any]:
    """Get peptide lists for a job."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"peptide_lists": job.peptide_lists}


@app.get("/api/jobs/{job_id}/results/{list_id}/{filename}")
async def download_result(job_id: str, list_id: str, filename: str) -> FileResponse:
    """Download a result file."""
    # Validate job exists
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Validate filename (prevent path traversal)
    allowed_files = [
        "taxonomy_nodes.csv",
        "go_terms.csv",
        "go_taxonomy_combo.csv",
        "coverage.csv",
        "run_manifest.json",
        "peptides_annotated.csv",
        "peptide_mapping.parquet",
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
async def download_all_results(job_id: str) -> FileResponse:
    """Download all results as a ZIP file."""
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not completed")

    results_dir = JOBS_DIR / job_id / "results"
    zip_path = results_dir / "all_results.zip"

    # Create ZIP if it doesn't exist. Building it is blocking I/O, so it runs
    # in a worker thread rather than stalling every other request.
    if not zip_path.exists():
        await run_in_threadpool(_build_results_zip, results_dir, zip_path)

    return FileResponse(
        zip_path,
        filename=f"metagomics2_results_{job_id[:8]}.zip",
        media_type="application/zip",
    )


def _build_results_zip(results_dir: Path, zip_path: Path) -> None:
    """Archive ``results_dir`` into ``zip_path`` without ever exposing a partial file.

    The archive is written to a temporary file in the job directory (outside
    ``results_dir``, so it cannot include itself) and moved into place with an
    atomic rename once complete.  Any earlier archive of the same name inside
    ``results_dir`` is skipped rather than nested.  Two concurrent builds each
    produce a complete file and the last rename wins.
    """
    fd, tmp_name = tempfile.mkstemp(
        prefix=".all_results.", suffix=".zip.part", dir=zip_path.parent.parent
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        with zipfile.ZipFile(tmp_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in sorted(results_dir.rglob("*")):
                if not file_path.is_file() or file_path.name == zip_path.name:
                    continue
                zf.write(file_path, file_path.relative_to(results_dir))
        os.replace(tmp_path, zip_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


# Frontend SPA support
_DEFAULT_FRONTEND_DIR = Path(__file__).parent.parent.parent.parent / "frontend" / "dist"
FRONTEND_DIR = _cfg.frontend_dir or _DEFAULT_FRONTEND_DIR


def _resolve_frontend_file(frontend_dir: Path, full_path: str) -> Path | None:
    """Return the file inside ``frontend_dir`` that ``full_path`` names, or None.

    The URL path is untrusted. Joining it onto the directory and resolving the
    result lets ``..`` segments and absolute paths escape the directory, so the
    resolved candidate must still lie inside the resolved directory and must be
    a regular file. Anything else falls back to ``index.html``.
    """
    base = frontend_dir.resolve()
    try:
        candidate = (base / full_path).resolve()
    except (OSError, RuntimeError):
        return None
    if not candidate.is_relative_to(base) or not candidate.is_file():
        return None
    return candidate


if FRONTEND_DIR.exists():
    # Serve static assets (JS, CSS, images, etc.)
    app.mount("/assets", StaticFiles(directory=str(FRONTEND_DIR / "assets")), name="assets")

    @app.get("/{full_path:path}")
    async def serve_spa(full_path: str) -> FileResponse:
        """Serve the SPA index.html for all non-API routes."""
        # Serve a real file from the dist directory if the path names one
        file_path = _resolve_frontend_file(FRONTEND_DIR, full_path) if full_path else None
        if file_path is not None:
            return FileResponse(file_path)
        # Fall back to index.html for SPA routing
        return FileResponse(FRONTEND_DIR / "index.html")
