"""Job data models and schemas."""

from datetime import datetime
from enum import Enum
from typing import Any, Literal

import math
import re

from pydantic import BaseModel, Field, field_validator

ALLOWED_SEARCH_TOOLS = {"diamond"}


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

    search_tool: Literal["diamond"] = "diamond"
    db_choice: str = ""
    db_name: str = ""

    # Filter parameters
    max_evalue: float | None = None
    min_pident: float | None = None
    min_qcov: float | None = None
    min_alnlen: int | None = None
    top_k: int | None = None

    # Notification
    notification_email: str = ""
    fasta_filename: str = ""

    # GO settings
    go_edge_types: str = "is_a,part_of"
    go_include_self: bool = True
    compute_enrichment_pvalues: bool = False

    @field_validator("max_evalue")
    @classmethod
    def validate_max_evalue(cls, v: float | None) -> float | None:
        if v is None:
            return v
        if not math.isfinite(v):
            raise ValueError("max_evalue must be a finite number")
        if v <= 0:
            raise ValueError("max_evalue must be greater than 0")
        if v > 1000:
            raise ValueError("max_evalue must be at most 1000")
        return v

    @field_validator("min_pident")
    @classmethod
    def validate_min_pident(cls, v: float | None) -> float | None:
        if v is None:
            return v
        if not math.isfinite(v):
            raise ValueError("min_pident must be a finite number")
        if v < 0 or v > 100:
            raise ValueError("min_pident must be between 0 and 100")
        return v

    @field_validator("min_qcov")
    @classmethod
    def validate_min_qcov(cls, v: float | None) -> float | None:
        if v is None:
            return v
        if not math.isfinite(v):
            raise ValueError("min_qcov must be a finite number")
        if v < 0 or v > 100:
            raise ValueError("min_qcov must be between 0 and 100")
        return v

    @field_validator("min_alnlen")
    @classmethod
    def validate_min_alnlen(cls, v: int | None) -> int | None:
        if v is None:
            return v
        if v < 1:
            raise ValueError("min_alnlen must be at least 1")
        return v

    @field_validator("top_k")
    @classmethod
    def validate_top_k(cls, v: int | None) -> int | None:
        if v is None:
            return v
        if v < 1:
            raise ValueError("top_k must be at least 1")
        return v

    @field_validator("db_choice")
    @classmethod
    def validate_db_choice(cls, v: str) -> str:
        if not v:
            return v
        if ".." in v or "/" in v or "\\" in v:
            raise ValueError("db_choice must be a plain filename, not a path")
        return v

    @field_validator("notification_email")
    @classmethod
    def validate_email(cls, v: str) -> str:
        v = v.strip()
        if not v:
            return v
        if not re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", v):
            raise ValueError("Invalid email address")
        return v


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
