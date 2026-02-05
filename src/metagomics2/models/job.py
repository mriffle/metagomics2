"""Job data models and schemas."""

from datetime import datetime
from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    """Job status enumeration."""

    UPLOADED = "uploaded"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class PeptideListStatus(str, Enum):
    """Peptide list processing status."""

    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"


class JobParams(BaseModel):
    """Parameters for a job."""

    search_tool: str = "diamond"
    db_choice: str = ""

    # Filter parameters
    max_evalue: float | None = None
    min_pident: float | None = None
    min_qcov: float | None = None
    min_alnlen: int | None = None
    top_k: int | None = None
    delta_bitscore: float | None = None
    best_hit_only: bool = False

    # GO settings
    go_edge_types: str = "is_a"
    go_include_self: bool = True


class JobCreate(BaseModel):
    """Request model for creating a job."""

    params: JobParams = Field(default_factory=JobParams)


class PeptideListInfo(BaseModel):
    """Information about a peptide list."""

    list_id: str
    filename: str
    status: PeptideListStatus
    n_peptides: int | None = None
    n_matched: int | None = None
    n_unmatched: int | None = None


class JobInfo(BaseModel):
    """Response model for job information."""

    job_id: str
    created_at: datetime
    status: JobStatus
    params: JobParams
    progress_total: int = 0
    progress_done: int = 0
    current_step: str | None = None
    error_message: str | None = None
    peptide_lists: list[PeptideListInfo] = []


class JobCreateResponse(BaseModel):
    """Response model for job creation."""

    job_id: str
    status: JobStatus


class JobListResponse(BaseModel):
    """Response model for listing jobs."""

    jobs: list[JobInfo]
