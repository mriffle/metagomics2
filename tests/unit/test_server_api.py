"""Unit tests for FastAPI server endpoints."""

import importlib
import io
import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

import metagomics2.config as config_module
from metagomics2.db.database import Database
from metagomics2.models.job import JobParams, JobStatus


@pytest.fixture
def test_db(tmp_path: Path):
    """Create a test database."""
    return Database(tmp_path / "test.db")


def _setup_config_dir(tmp_path: Path, databases=None):
    """Write a databases.json config file and return the config dir."""
    config_dir = tmp_path / "config"
    config_dir.mkdir(exist_ok=True)
    if databases is None:
        databases = [
            {"name": "Test DB", "description": "Test", "path": "test.dmnd",
             "annotations": "test.annotations.db"}
        ]
    (config_dir / "databases.json").write_text(json.dumps(databases))
    return config_dir


@pytest.fixture
def client(tmp_path: Path, test_db):
    """Create a test client with patched paths.

    We must patch the module-level globals AFTER import, since the module
    creates a Database at import time. We set METAGOMICS_DATA_DIR env var
    and reload the module to get a fresh app with tmp_path-based paths.
    """
    (tmp_path / "jobs").mkdir(exist_ok=True)
    config_dir = _setup_config_dir(tmp_path)

    with patch.dict(os.environ, {
        "METAGOMICS_DATA_DIR": str(tmp_path),
        "METAGOMICS_ADMIN_PASSWORD": "testpass",
        "METAGOMICS_CONFIG_DIR": str(config_dir),
    }):
        config_module.reset_settings()
        import metagomics2.server.app as app_module
        importlib.reload(app_module)

        # Override the db with our test db
        app_module.db = test_db
        app_module.JOBS_DIR = tmp_path / "jobs"

        from fastapi.testclient import TestClient
        yield TestClient(app_module.app)

    config_module.reset_settings()


class TestHealthCheck:
    """Tests for health check endpoint."""

    def test_health_returns_200(self, client):
        response = client.get("/api/health")
        assert response.status_code == 200

    def test_health_returns_status(self, client):
        data = client.get("/api/health").json()
        assert data["status"] == "healthy"

    def test_health_includes_version(self, client):
        data = client.get("/api/health").json()
        assert "version" in data


class TestVersionEndpoint:
    """Tests for version endpoint."""

    def test_version_returns_200(self, client):
        response = client.get("/api/version")
        assert response.status_code == 200

    def test_version_returns_version_string(self, client):
        data = client.get("/api/version").json()
        assert "version" in data
        assert isinstance(data["version"], str)
        assert len(data["version"]) > 0


class TestCreateJob:
    """Tests for job creation endpoint."""

    def test_create_job_success(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "queued"

    def test_create_job_with_params(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"
        params = json.dumps({"search_tool": "diamond", "max_evalue": 1e-5})

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": params},
        )

        assert response.status_code == 200

    def test_create_job_multiple_peptide_files(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        pep1 = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"
        pep2 = b"peptide_sequence\tquantity\nABC\t5\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("pep1.tsv", io.BytesIO(pep1), "text/tab-separated-values"),
                ),
                (
                    "peptides",
                    ("pep2.tsv", io.BytesIO(pep2), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 200

    def test_create_job_invalid_params(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": "not valid json"},
        )

        assert response.status_code == 400

    def test_create_job_saves_files(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        job_id = response.json()["job_id"]
        job_dir = tmp_path / "jobs" / job_id

        assert (job_dir / "inputs" / "background.fasta").exists()
        assert (job_dir / "inputs" / "peptides").exists()
        assert (job_dir / "work").exists()
        assert (job_dir / "results").exists()
        assert (job_dir / "logs").exists()

    def test_create_job_rejects_non_fasta(self, client, tmp_path: Path):
        """A plain text file that isn't FASTA should be rejected."""
        bad_content = b"This is not a FASTA file\nJust some text\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("bad.txt", io.BytesIO(bad_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 400
        assert "does not appear to be a valid FASTA" in response.json()["detail"]

    def test_create_job_rejects_empty_fasta(self, client, tmp_path: Path):
        """An empty file should be rejected."""
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                ("fasta", ("empty.fasta", io.BytesIO(b""), "application/octet-stream")),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 400
        assert "empty" in response.json()["detail"].lower()

    def test_create_job_rejects_header_only_fasta(self, client, tmp_path: Path):
        """A FASTA with a header but no sequence should be rejected."""
        bad_content = b">protein1\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("header_only.fasta", io.BytesIO(bad_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 400
        assert "no sequence data" in response.json()["detail"].lower()

    def test_create_job_rejects_consecutive_headers(self, client, tmp_path: Path):
        """A FASTA with two headers in a row should be rejected."""
        bad_content = b">protein1\n>protein2\nACDE\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                (
                    "fasta",
                    ("bad.fasta", io.BytesIO(bad_content), "application/octet-stream"),
                ),
                (
                    "peptides",
                    ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                ),
            ],
            data={"params": json.dumps({"db_choice": "test.dmnd"})},
        )

        assert response.status_code == 400
        assert "consecutive header" in response.json()["detail"].lower()


class TestUploadSizeLimit:
    """Tests for upload size limit enforcement."""

    def test_fasta_exceeding_limit_rejected(self, tmp_path: Path, test_db):
        """A FASTA file exceeding MAX_UPLOAD_MB should be rejected with 413."""
        (tmp_path / "jobs").mkdir(exist_ok=True)
        config_dir = _setup_config_dir(tmp_path)

        with patch.dict(os.environ, {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_ADMIN_PASSWORD": "testpass",
            "METAGOMICS_CONFIG_DIR": str(config_dir),
            "METAGOMICS_MAX_UPLOAD_MB": "1",  # 1 MB limit
        }):
            config_module.reset_settings()
            import metagomics2.server.app as app_module
            importlib.reload(app_module)
            app_module.db = test_db
            app_module.JOBS_DIR = tmp_path / "jobs"

            from fastapi.testclient import TestClient
            client = TestClient(app_module.app)

            # Create a FASTA file just over 1 MB
            fasta_content = b">P1\n" + b"M" * (1024 * 1024 + 100) + b"\n"
            peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

            response = client.post(
                "/api/jobs",
                files=[
                    (
                        "fasta",
                        ("big.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                    ),
                    (
                        "peptides",
                        ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                    ),
                ],
                data={"params": json.dumps({"db_choice": "test.dmnd"})},
            )

            assert response.status_code == 413
            assert "maximum upload size" in response.json()["detail"].lower()
            # No orphan job row or directory is left behind
            assert test_db.list_jobs() == []
            assert list((tmp_path / "jobs").iterdir()) == []

    def test_peptides_exceeding_limit_rejected(self, tmp_path: Path, test_db):
        """Peptide files are limited by their combined size, and nothing is left behind."""
        (tmp_path / "jobs").mkdir(exist_ok=True)
        config_dir = _setup_config_dir(tmp_path)

        with patch.dict(os.environ, {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_ADMIN_PASSWORD": "testpass",
            "METAGOMICS_CONFIG_DIR": str(config_dir),
            "METAGOMICS_MAX_UPLOAD_MB": "1",
        }):
            config_module.reset_settings()
            import metagomics2.server.app as app_module
            importlib.reload(app_module)
            app_module.db = test_db
            app_module.JOBS_DIR = tmp_path / "jobs"

            from fastapi.testclient import TestClient
            client = TestClient(app_module.app)

            fasta_content = b">P1\nMPEPTIDEK\n"
            # Two peptide files of 600 KB each: fine alone, over 1 MB together
            big = b"peptide_sequence\tquantity\n" + b"PEPTIDE\t1\n" * 60_000
            response = client.post(
                "/api/jobs",
                files=[
                    ("fasta", ("s.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                    ("peptides", ("a.tsv", io.BytesIO(big), "text/tab-separated-values")),
                    ("peptides", ("b.tsv", io.BytesIO(big), "text/tab-separated-values")),
                ],
                data={"params": json.dumps({"db_choice": "test.dmnd"})},
            )

            assert response.status_code == 413
            assert "peptide" in response.json()["detail"].lower()
            assert test_db.list_jobs() == []
            assert list((tmp_path / "jobs").iterdir()) == []

    def test_streaming_aborts_past_limit(self, tmp_path: Path):
        """The stream is abandoned once it passes the limit; the disk is not filled first."""
        import asyncio

        from fastapi import UploadFile

        from metagomics2.server.app import (
            _WRITE_CHUNK_SIZE,
            UploadTooLargeError,
            _save_upload_streamed,
        )

        payload = b"x" * (10 * _WRITE_CHUNK_SIZE)  # 10 MB
        limit = 2 * _WRITE_CHUNK_SIZE
        dest = tmp_path / "upload.bin"

        async def run():
            upload = UploadFile(file=io.BytesIO(payload), filename="big.bin")
            await _save_upload_streamed(upload, dest, limit)

        with pytest.raises(UploadTooLargeError):
            asyncio.run(run())
        # At most the chunks up to the limit reached the disk, not the whole payload
        assert dest.stat().st_size <= limit

    def test_failed_upload_leaves_no_job(self, client, test_db, tmp_path: Path):
        """An error while storing inputs rolls the job back entirely."""
        import metagomics2.server.app as app_module

        with patch.object(app_module, "_save_upload_streamed", side_effect=OSError("disk full")):
            with pytest.raises(OSError):
                client.post(
                    "/api/jobs",
                    files=[
                        ("fasta", ("s.fasta", io.BytesIO(b">P1\nMPEPTIDEK\n"), "text/plain")),
                        ("peptides", ("p.tsv", io.BytesIO(b"PEPTIDE\t1\n"), "text/plain")),
                    ],
                    data={"params": json.dumps({"db_choice": "test.dmnd"})},
                )
        assert test_db.list_jobs() == []
        assert list((tmp_path / "jobs").iterdir()) == []

    def test_fasta_within_limit_accepted(self, tmp_path: Path, test_db):
        """A FASTA file within MAX_UPLOAD_MB should be accepted."""
        (tmp_path / "jobs").mkdir(exist_ok=True)
        config_dir = _setup_config_dir(tmp_path)

        with patch.dict(os.environ, {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_ADMIN_PASSWORD": "testpass",
            "METAGOMICS_CONFIG_DIR": str(config_dir),
            "METAGOMICS_MAX_UPLOAD_MB": "1",  # 1 MB limit
        }):
            config_module.reset_settings()
            import metagomics2.server.app as app_module
            importlib.reload(app_module)
            app_module.db = test_db
            app_module.JOBS_DIR = tmp_path / "jobs"

            from fastapi.testclient import TestClient
            client = TestClient(app_module.app)

            # Small valid FASTA
            fasta_content = b">P1\nMPEPTIDEK\n"
            peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

            response = client.post(
                "/api/jobs",
                files=[
                    (
                        "fasta",
                        ("small.fasta", io.BytesIO(fasta_content), "application/octet-stream"),
                    ),
                    (
                        "peptides",
                        ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values"),
                    ),
                ],
                data={"params": json.dumps({"db_choice": "test.dmnd"})},
            )

            assert response.status_code == 200


class TestGetJob:
    """Tests for job retrieval endpoint."""

    def test_get_existing_job(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        response = client.get(f"/api/jobs/{job_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == job_id
        assert data["status"] == "queued"

    def test_get_nonexistent_job(self, client):
        response = client.get("/api/jobs/nonexistent_id")
        assert response.status_code == 404

    def test_get_job_includes_peptide_lists(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.add_peptide_list(job_id, "list_000", "pep.tsv", "/path/pep.tsv")

        response = client.get(f"/api/jobs/{job_id}")
        data = response.json()
        assert len(data["peptide_lists"]) == 1
        assert data["peptide_lists"][0]["list_id"] == "list_000"


class TestAdminAuth:
    """Tests for admin authentication."""

    def test_login_success(self, client):
        response = client.post("/api/admin/auth", json={"password": "testpass"})
        assert response.status_code == 200
        data = response.json()
        assert "token" in data
        assert len(data["token"]) > 0

    def test_login_wrong_password(self, client):
        response = client.post("/api/admin/auth", json={"password": "wrong"})
        assert response.status_code == 401

    def test_admin_endpoint_without_token(self, client):
        response = client.get("/api/admin/jobs")
        assert response.status_code == 401

    def test_admin_endpoint_with_bad_token(self, client):
        response = client.get("/api/admin/jobs", headers={"Authorization": "Bearer badtoken"})
        assert response.status_code == 401

    def test_token_expires(self, client):
        import metagomics2.server.app as app_module

        token = client.post("/api/admin/auth", json={"password": "testpass"}).json()["token"]
        headers = {"Authorization": f"Bearer {token}"}
        assert client.get("/api/admin/jobs", headers=headers).status_code == 200

        later = app_module.time.monotonic() + app_module._ADMIN_TOKEN_TTL_SECONDS + 1
        with patch("metagomics2.server.app.time.monotonic", return_value=later):
            assert client.get("/api/admin/jobs", headers=headers).status_code == 401
        # The expired token was dropped from the table
        assert token not in app_module._admin_tokens

    def test_token_table_is_capped(self, client):
        import metagomics2.server.app as app_module

        with patch("metagomics2.server.app._ADMIN_TOKEN_LIMIT", 3):
            tokens = [
                client.post("/api/admin/auth", json={"password": "testpass"}).json()["token"]
                for _ in range(5)
            ]
        assert len(app_module._admin_tokens) == 3
        for token in tokens[:2]:
            headers = {"Authorization": f"Bearer {token}"}
            assert client.get("/api/admin/jobs", headers=headers).status_code == 401
        for token in tokens[2:]:
            headers = {"Authorization": f"Bearer {token}"}
            assert client.get("/api/admin/jobs", headers=headers).status_code == 200


def _get_admin_token(client) -> str:
    """Helper to get a valid admin token."""
    response = client.post("/api/admin/auth", json={"password": "testpass"})
    return response.json()["token"]


class TestListJobs:
    """Tests for admin job listing endpoint."""

    def test_list_empty(self, client):
        token = _get_admin_token(client)
        response = client.get("/api/admin/jobs", headers={"Authorization": f"Bearer {token}"})
        assert response.status_code == 200
        data = response.json()
        assert data["jobs"] == []

    def test_list_returns_jobs(self, client, test_db):
        test_db.create_job(JobParams())
        test_db.create_job(JobParams())
        token = _get_admin_token(client)

        response = client.get("/api/admin/jobs", headers={"Authorization": f"Bearer {token}"})
        data = response.json()
        assert len(data["jobs"]) == 2

    def test_list_respects_limit(self, client, test_db):
        for _ in range(5):
            test_db.create_job(JobParams())
        token = _get_admin_token(client)

        response = client.get(
            "/api/admin/jobs?limit=3", headers={"Authorization": f"Bearer {token}"}
        )
        data = response.json()
        assert len(data["jobs"]) == 3


class TestGetPeptideLists:
    """Tests for peptide lists endpoint."""

    def test_get_peptide_lists(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.add_peptide_list(job_id, "list_000", "pep1.tsv", "/path/pep1.tsv")
        test_db.add_peptide_list(job_id, "list_001", "pep2.tsv", "/path/pep2.tsv")

        response = client.get(f"/api/jobs/{job_id}/peptide-lists")
        assert response.status_code == 200
        data = response.json()
        assert len(data["peptide_lists"]) == 2

    def test_get_peptide_lists_nonexistent_job(self, client):
        response = client.get("/api/jobs/nonexistent/peptide-lists")
        assert response.status_code == 404


class TestDownloadResult:
    """Tests for result file download endpoint."""

    def test_download_existing_file(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())

        # Create result file
        result_dir = tmp_path / "jobs" / job_id / "results" / "list_000"
        result_dir.mkdir(parents=True)
        (result_dir / "coverage.csv").write_text("col1,col2\n1,2\n")

        response = client.get(f"/api/jobs/{job_id}/results/list_000/coverage.csv")
        assert response.status_code == 200

    def test_download_nonexistent_job(self, client):
        response = client.get("/api/jobs/nonexistent/results/list_000/coverage.csv")
        assert response.status_code == 404

    def test_download_nonexistent_file(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        response = client.get(f"/api/jobs/{job_id}/results/list_000/coverage.csv")
        assert response.status_code == 404

    def test_download_disallowed_filename(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        response = client.get(f"/api/jobs/{job_id}/results/list_000/secret_data.txt")
        assert response.status_code == 400

    def test_allowed_filenames(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        result_dir = tmp_path / "jobs" / job_id / "results" / "list_000"
        result_dir.mkdir(parents=True)

        allowed = ["taxonomy_nodes.csv", "go_terms.csv", "coverage.csv",
                    "run_manifest.json", "peptides_annotated.csv"]
        for fname in allowed:
            (result_dir / fname).write_text("test content")
            response = client.get(f"/api/jobs/{job_id}/results/list_000/{fname}")
            assert response.status_code == 200, f"Failed for {fname}"


class TestRegenerateJobId:
    """Tests for revoking a shared job link."""

    def test_regenerates_completed_job(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.COMPLETED)
        job_dir = tmp_path / "jobs" / job_id
        (job_dir / "results").mkdir(parents=True)

        response = client.post(f"/api/jobs/{job_id}/regenerate-id")

        assert response.status_code == 200
        new_id = response.json()["new_job_id"]
        assert new_id != job_id
        assert client.get(f"/api/jobs/{new_id}").status_code == 200
        assert client.get(f"/api/jobs/{job_id}").status_code == 404
        assert (tmp_path / "jobs" / new_id / "results").is_dir()
        assert not job_dir.exists()

    @pytest.mark.parametrize("status", [JobStatus.UPLOADED, JobStatus.QUEUED, JobStatus.RUNNING])
    def test_refuses_unfinished_job(self, client, test_db, status: JobStatus):
        """Renaming a job the worker may be processing would strand it."""
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, status)

        response = client.post(f"/api/jobs/{job_id}/regenerate-id")

        assert response.status_code == 409
        assert status.value in response.json()["detail"]
        assert client.get(f"/api/jobs/{job_id}").status_code == 200

    def test_failed_job_can_be_regenerated(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.FAILED, "boom")

        assert client.post(f"/api/jobs/{job_id}/regenerate-id").status_code == 200

    def test_nonexistent_job(self, client):
        assert client.post("/api/jobs/nope/regenerate-id").status_code == 404


class TestDownloadAllResults:
    """Tests for ZIP download endpoint."""

    def test_download_all_completed_job(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.COMPLETED)

        result_dir = tmp_path / "jobs" / job_id / "results"
        result_dir.mkdir(parents=True)
        (result_dir / "list_000").mkdir()
        (result_dir / "list_000" / "coverage.csv").write_text("test")

        response = client.get(f"/api/jobs/{job_id}/results/all_results.zip")
        assert response.status_code == 200
        assert "zip" in response.headers.get("content-type", "")

        import zipfile
        with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
            assert zf.namelist() == ["list_000/coverage.csv"]
            assert zf.read("list_000/coverage.csv") == b"test"
        # No stray temporary file is left in the job directory
        job_dir = tmp_path / "jobs" / job_id
        assert [p.name for p in job_dir.iterdir()] == ["results"]

    def test_zip_never_nests_a_previous_archive(self, tmp_path: Path):
        """A stale all_results.zip inside results/ is skipped, not zipped into itself."""
        import zipfile

        from metagomics2.server.app import _build_results_zip

        results_dir = tmp_path / "job" / "results"
        (results_dir / "list_000").mkdir(parents=True)
        (results_dir / "list_000" / "go_terms.csv").write_text("go")
        (results_dir / "list_001").mkdir()
        (results_dir / "list_001" / "coverage.csv").write_text("cov")
        zip_path = results_dir / "all_results.zip"
        zip_path.write_bytes(b"stale partial archive")

        _build_results_zip(results_dir, zip_path)

        with zipfile.ZipFile(zip_path) as zf:
            assert zf.namelist() == ["list_000/go_terms.csv", "list_001/coverage.csv"]
        assert [p.name for p in (tmp_path / "job").iterdir()] == ["results"]

    def test_zip_build_failure_leaves_no_partial(self, tmp_path: Path):
        from metagomics2.server.app import _build_results_zip

        results_dir = tmp_path / "job" / "results"
        results_dir.mkdir(parents=True)
        (results_dir / "coverage.csv").write_text("x")
        zip_path = results_dir / "all_results.zip"

        with patch("metagomics2.server.app.os.replace", side_effect=OSError("boom")):
            with pytest.raises(OSError):
                _build_results_zip(results_dir, zip_path)

        assert not zip_path.exists()
        assert [p.name for p in (tmp_path / "job").iterdir()] == ["results"]

    def test_download_all_incomplete_job(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.RUNNING)

        response = client.get(f"/api/jobs/{job_id}/results/all_results.zip")
        assert response.status_code == 400

    def test_download_all_nonexistent_job(self, client):
        response = client.get("/api/jobs/nonexistent/results/all_results.zip")
        assert response.status_code == 404


class TestFrontendPathTraversal:
    """The SPA catch-all must never serve files outside the frontend directory."""

    @staticmethod
    def _make_dist(tmp_path: Path) -> Path:
        dist = tmp_path / "dist"
        (dist / "assets").mkdir(parents=True)
        (dist / "index.html").write_text("<html>index</html>")
        (dist / "app.js").write_text("console.log('app')")
        (tmp_path / "secret.txt").write_text("SECRET")
        return dist

    def test_resolve_file_inside_dist(self, tmp_path: Path):
        from metagomics2.server.app import _resolve_frontend_file

        dist = self._make_dist(tmp_path)
        assert _resolve_frontend_file(dist, "app.js") == (dist / "app.js").resolve()

    @pytest.mark.parametrize(
        "full_path",
        ["../secret.txt", "foo/../../secret.txt", "/etc/passwd", "", "missing.js", "assets"],
    )
    def test_resolve_rejects_escapes_and_non_files(self, tmp_path: Path, full_path: str):
        from metagomics2.server.app import _resolve_frontend_file

        dist = self._make_dist(tmp_path)
        assert _resolve_frontend_file(dist, full_path) is None

    def test_http_requests_cannot_escape_dist(self, tmp_path: Path, test_db):
        """Percent-encoded dot segments and absolute paths fall back to index.html.

        The test client normalises a raw ``..`` before sending, but delivers
        ``%2e%2e`` decoded, which is exactly what a hand-crafted request does.
        """
        (tmp_path / "jobs").mkdir(exist_ok=True)
        config_dir = _setup_config_dir(tmp_path)
        dist = self._make_dist(tmp_path)

        with patch.dict(os.environ, {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_ADMIN_PASSWORD": "testpass",
            "METAGOMICS_CONFIG_DIR": str(config_dir),
            "METAGOMICS_FRONTEND_DIR": str(dist),
        }):
            config_module.reset_settings()
            import metagomics2.server.app as app_module
            importlib.reload(app_module)
            app_module.db = test_db
            app_module.JOBS_DIR = tmp_path / "jobs"

            from fastapi.testclient import TestClient
            client = TestClient(app_module.app)

            assert client.get("/app.js").text == "console.log('app')"
            assert client.get("/job/abc").text == "<html>index</html>"

            for path in ("/%2e%2e/secret.txt", "/foo/%2e%2e/%2e%2e/secret.txt", "/%2Fetc/passwd"):
                response = client.get(path)
                assert response.status_code == 200, path
                assert response.text == "<html>index</html>", path

        config_module.reset_settings()
