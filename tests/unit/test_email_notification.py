"""Tests for email notification sending."""

from datetime import datetime, timezone
from unittest.mock import MagicMock, patch

import pytest

from metagomics2.models.job import (
    JobInfo,
    JobParams,
    JobStatus,
    PeptideListInfo,
    PeptideListStatus,
)
from metagomics2.notifications.email import SmtpConfig, send_job_notification


def _make_job(
    status: JobStatus = JobStatus.COMPLETED,
    notification_email: str = "user@example.com",
    fasta_filename: str = "proteome.fasta",
    error_message: str | None = None,
) -> JobInfo:
    """Create a minimal JobInfo for testing."""
    return JobInfo(
        job_id="test-job-123",
        created_at=datetime(2025, 1, 15, 12, 0, 0, tzinfo=timezone.utc),
        status=status,
        params=JobParams(
            db_choice="uniprot_sprot.dmnd",
            db_name="UniProt SwissProt",
            max_evalue=1e-10,
            min_pident=80.0,
            top_k=1,
            notification_email=notification_email,
            fasta_filename=fasta_filename,
        ),
        error_message=error_message,
        peptide_lists=[
            PeptideListInfo(
                list_id="list_000",
                filename="sample_peptides.tsv",
                status=PeptideListStatus.DONE,
            ),
            PeptideListInfo(
                list_id="list_001",
                filename="control_peptides.tsv",
                status=PeptideListStatus.DONE,
            ),
        ],
    )


SMTP_CONFIG = SmtpConfig(
    host="smtp.example.com",
    port=587,
    username="sender@example.com",
    password="secret",
    from_address="noreply@example.com",
)

SITE_URL = "https://metagomics.example.com"


class TestSendJobNotification:
    """Tests for send_job_notification."""

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_send_email_success(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        mock_smtp_cls.assert_called_once_with("smtp.example.com", 587, timeout=30)
        mock_server.send_message.assert_called_once()

        msg = mock_server.send_message.call_args[0][0]
        assert msg["To"] == "user@example.com"
        assert msg["From"] == "noreply@example.com"
        assert "completed" in msg["Subject"]

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_send_email_skipped_no_smtp_host(self, mock_smtp_cls):
        empty_config = SmtpConfig(host="")
        job = _make_job()
        send_job_notification(job, SITE_URL, empty_config)

        mock_smtp_cls.assert_not_called()

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_send_email_skipped_no_recipient(self, mock_smtp_cls):
        job = _make_job(notification_email="")
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        mock_smtp_cls.assert_not_called()

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_send_email_smtp_error_logged_not_raised(self, mock_smtp_cls):
        mock_smtp_cls.side_effect = ConnectionRefusedError("Connection refused")

        job = _make_job()
        # Should not raise
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_email_body_contains_fasta_filename(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        body = msg.get_content()
        assert "proteome.fasta" in body

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_email_body_contains_peptide_filenames(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        body = msg.get_content()
        assert "sample_peptides.tsv" in body
        assert "control_peptides.tsv" in body

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_email_body_contains_parameters(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        body = msg.get_content()
        assert "UniProt SwissProt" in body
        assert "1e-10" in body
        assert "80" in body

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_email_body_contains_job_link(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        body = msg.get_content()
        assert "https://metagomics.example.com/job/test-job-123" in body

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_failed_job_email_contains_error(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job(
            status=JobStatus.FAILED,
            error_message="DIAMOND crashed",
        )
        send_job_notification(job, SITE_URL, SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        assert "failed" in msg["Subject"]
        body = msg.get_content()
        assert "DIAMOND crashed" in body

    @patch("metagomics2.notifications.email.smtplib.SMTP")
    def test_email_no_job_link_when_no_site_url(self, mock_smtp_cls):
        mock_server = MagicMock()
        mock_smtp_cls.return_value.__enter__ = MagicMock(return_value=mock_server)
        mock_smtp_cls.return_value.__exit__ = MagicMock(return_value=False)

        job = _make_job()
        send_job_notification(job, "", SMTP_CONFIG)

        msg = mock_server.send_message.call_args[0][0]
        body = msg.get_content()
        assert "/job/" not in body
