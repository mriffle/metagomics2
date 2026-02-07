"""Email notification for job completion."""

import logging
import smtplib
from dataclasses import dataclass
from email.message import EmailMessage

from metagomics2.models.job import JobInfo, JobStatus

logger = logging.getLogger(__name__)


@dataclass
class SmtpConfig:
    """SMTP server configuration."""

    host: str
    port: int = 587
    username: str = ""
    password: str = ""
    from_address: str = ""


def send_job_notification(
    job: JobInfo,
    site_url: str,
    smtp_config: SmtpConfig,
) -> None:
    """Send an email notification about job completion or failure.

    This function never raises; SMTP errors are logged and swallowed so
    that email failures cannot break the pipeline.

    Args:
        job: The completed/failed job info.
        site_url: Base URL of the site (e.g. ``https://metagomics.example.com``).
        smtp_config: SMTP connection details.
    """
    recipient = job.params.notification_email
    if not recipient:
        return
    if not smtp_config.host:
        logger.debug("SMTP not configured, skipping email notification")
        return

    try:
        subject, body = _build_message(job, site_url)
        msg = EmailMessage()
        msg["Subject"] = subject
        msg["From"] = smtp_config.from_address or smtp_config.username
        msg["To"] = recipient
        msg.set_content(body)

        with smtplib.SMTP(smtp_config.host, smtp_config.port, timeout=30) as server:
            server.ehlo()
            server.starttls()
            server.ehlo()
            if smtp_config.username:
                server.login(smtp_config.username, smtp_config.password)
            server.send_message(msg)

        logger.info(f"Notification email sent to {recipient} for job {job.job_id}")

    except Exception:
        logger.exception(
            f"Failed to send notification email to {recipient} for job {job.job_id}"
        )


def _build_message(job: JobInfo, site_url: str) -> tuple[str, str]:
    """Build the email subject and plain-text body.

    Returns:
        A (subject, body) tuple.
    """
    status_label = "completed" if job.status == JobStatus.COMPLETED else "failed"
    subject = f"Metagomics 2 — Job {status_label}"

    job_url = f"{site_url.rstrip('/')}/job/{job.job_id}" if site_url else ""

    lines: list[str] = []
    lines.append(f"Your Metagomics 2 annotation job has {status_label}.")
    lines.append("")

    if job.status == JobStatus.FAILED and job.error_message:
        lines.append(f"Error: {job.error_message}")
        lines.append("")

    # Uploaded files
    lines.append("--- Uploaded Files ---")
    if job.params.fasta_filename:
        lines.append(f"  FASTA file: {job.params.fasta_filename}")
    if job.peptide_lists:
        lines.append(f"  Peptide lists:")
        for pl in job.peptide_lists:
            lines.append(f"    - {pl.filename}")
    lines.append("")

    # Parameters
    lines.append("--- Parameters ---")
    if job.params.db_name:
        lines.append(f"  Database: {job.params.db_name} ({job.params.db_choice})")
    elif job.params.db_choice:
        lines.append(f"  Database: {job.params.db_choice}")
    if job.params.max_evalue is not None:
        lines.append(f"  Max E-value: {job.params.max_evalue}")
    if job.params.min_pident is not None:
        lines.append(f"  Min % Identity: {job.params.min_pident}")
    if job.params.top_k is not None:
        lines.append(f"  Top K Hits: {job.params.top_k}")
    if job.params.min_qcov is not None:
        lines.append(f"  Min Query Coverage: {job.params.min_qcov}")
    if job.params.min_alnlen is not None:
        lines.append(f"  Min Alignment Length: {job.params.min_alnlen}")
    lines.append("")

    # Job link
    if job_url:
        lines.append(f"View your results: {job_url}")
        lines.append("")

    lines.append(f"Job ID: {job.job_id}")

    return subject, "\n".join(lines)
