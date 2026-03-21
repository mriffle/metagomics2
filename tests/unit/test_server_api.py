"""Unit tests for FastAPI server endpoints."""

import asyncio
import importlib
import io
import json
import os
from pathlib import Path
from unittest.mock import patch

import httpx
import pytest
from fastapi import HTTPException

import metagomics2.config as config_module
from metagomics2.db.database import Database
from metagomics2.models.job import JobParams, JobStatus, PeptideListStatus


@pytest.fixture
def test_db(tmp_path: Path):
    """Create a test database."""
    return Database(tmp_path / "test.db")


def _setup_config_dir(tmp_path: Path, databases=None):
    """Write a databases.json config file and return the config dir."""
    config_dir = tmp_path / "config"
    config_dir.mkdir(exist_ok=True)
    if databases is None:
        databases = [{"name": "Test DB", "description": "Test", "path": "test.dmnd"}]
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
        app_module.aiofiles.open = fake_aiofiles_open
        yield SyncASGIClient(app_module.app)

    config_module.reset_settings()


class SyncASGIClient:
    """Synchronous wrapper around httpx.AsyncClient for ASGI apps."""

    def __init__(self, app):
        self.app = app

    def request(self, method: str, url: str, **kwargs):
        async def _request():
            transport = httpx.ASGITransport(app=self.app)
            async with httpx.AsyncClient(
                transport=transport,
                base_url="http://testserver",
            ) as client:
                return await client.request(method, url, **kwargs)

        return asyncio.run(_request())

    def get(self, url: str, **kwargs):
        return self.request("GET", url, **kwargs)

    def post(self, url: str, **kwargs):
        return self.request("POST", url, **kwargs)


class _AsyncFileWrapper:
    """Minimal async file wrapper for tests."""

    def __init__(self, path: Path, mode: str):
        self._file = open(path, mode)

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc, tb):
        self._file.close()

    async def write(self, data):
        self._file.write(data)


def fake_aiofiles_open(path: str | Path, mode: str = "r", *args, **kwargs):
    """Replacement for aiofiles.open used in ASGI transport tests."""
    return _AsyncFileWrapper(Path(path), mode)


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
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
            ],
            data={"params": params},
        )

        assert response.status_code == 200

    def test_create_job_persists_enrichment_toggle(self, client):
        fasta_content = b">P1\nMPEPTIDEK\n"
        peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

        response = client.post(
            "/api/jobs",
            files=[
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
            ],
            data={"params": json.dumps({"compute_enrichment_pvalues": True})},
        )

        assert response.status_code == 200
        job_id = response.json()["job_id"]
        job = client.get(f"/api/jobs/{job_id}").json()
        assert job["params"]["compute_enrichment_pvalues"] is True

    def test_create_job_multiple_peptide_files(self, client, tmp_path: Path):
        fasta_content = b">P1\nMPEPTIDEK\n"
        pep1 = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"
        pep2 = b"peptide_sequence\tquantity\nABC\t5\n"

        response = client.post(
            "/api/jobs",
            files=[
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("pep1.tsv", io.BytesIO(pep1), "text/tab-separated-values")),
                ("peptides", ("pep2.tsv", io.BytesIO(pep2), "text/tab-separated-values")),
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
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("fasta", ("background.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("fasta", ("bad.txt", io.BytesIO(bad_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("fasta", ("header_only.fasta", io.BytesIO(bad_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
                ("fasta", ("bad.fasta", io.BytesIO(bad_content), "application/octet-stream")),
                ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
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
            app_module.aiofiles.open = fake_aiofiles_open
            async def fake_save_upload_streamed(upload, dest):
                dest.write_bytes(b"x")
                return app_module.MAX_UPLOAD_BYTES + 1
            app_module._save_upload_streamed = fake_save_upload_streamed

            client = SyncASGIClient(app_module.app)

            fasta_content = b">P1\nMPEPTIDEK\n"
            peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

            response = client.post(
                "/api/jobs",
                files=[
                    ("fasta", ("big.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                    ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
                ],
                data={"params": json.dumps({"db_choice": "test.dmnd"})},
            )

            assert response.status_code == 413
            assert "maximum upload size" in response.json()["detail"].lower()

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
            app_module.aiofiles.open = fake_aiofiles_open
            async def fake_save_upload_streamed(upload, dest):
                dest.write_bytes(b"x")
                return app_module.MAX_UPLOAD_BYTES - 1
            app_module._save_upload_streamed = fake_save_upload_streamed

            client = SyncASGIClient(app_module.app)

            # Small valid FASTA
            fasta_content = b">P1\nMPEPTIDEK\n"
            peptide_content = b"peptide_sequence\tquantity\nPEPTIDE\t10\n"

            response = client.post(
                "/api/jobs",
                files=[
                    ("fasta", ("small.fasta", io.BytesIO(fasta_content), "application/octet-stream")),
                    ("peptides", ("peptides.tsv", io.BytesIO(peptide_content), "text/tab-separated-values")),
                ],
                data={"params": json.dumps({"db_choice": "test.dmnd"})},
            )

            assert response.status_code == 200


class TestGetJob:
    """Tests for job retrieval endpoint."""

    def test_get_existing_job(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.QUEUED)

        import metagomics2.server.app as app_module

        data = asyncio.run(app_module.get_job(job_id))
        assert data.job_id == job_id
        assert data.status == "queued"

    def test_get_nonexistent_job(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.get_job("nonexistent_id"))

        assert exc_info.value.status_code == 404

    def test_get_job_includes_peptide_lists(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.add_peptide_list(job_id, "list_000", "pep.tsv", "/path/pep.tsv")

        import metagomics2.server.app as app_module

        data = asyncio.run(app_module.get_job(job_id))
        assert len(data.peptide_lists) == 1
        assert data.peptide_lists[0].list_id == "list_000"


class TestAdminAuth:
    """Tests for admin authentication."""

    def test_login_success(self, client):
        response = client.post("/api/admin/auth", json={"password": "testpass"})
        assert response.status_code == 200
        data = response.json()
        assert "token" in data
        assert len(data["token"]) > 0

    def test_login_wrong_password(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.admin_login(app_module.AdminAuthRequest(password="wrong")))

        assert exc_info.value.status_code == 401

    def test_admin_endpoint_without_token(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            app_module.require_admin("")

        assert exc_info.value.status_code == 401

    def test_admin_endpoint_with_bad_token(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            app_module.require_admin("Bearer badtoken")

        assert exc_info.value.status_code == 401


def _get_admin_token(client) -> str:
    """Helper to get a valid admin token."""
    response = client.post("/api/admin/auth", json={"password": "testpass"})
    return response.json()["token"]


class TestListJobs:
    """Tests for admin job listing endpoint."""

    def test_list_empty(self, client):
        import metagomics2.server.app as app_module

        token = asyncio.run(app_module.admin_login(app_module.AdminAuthRequest(password="testpass"))).token
        data = asyncio.run(app_module.list_jobs(100, token))
        assert data.jobs == []

    def test_list_returns_jobs(self, client, test_db):
        test_db.create_job(JobParams())
        test_db.create_job(JobParams())
        import metagomics2.server.app as app_module

        token = asyncio.run(app_module.admin_login(app_module.AdminAuthRequest(password="testpass"))).token
        data = asyncio.run(app_module.list_jobs(100, token))
        assert len(data.jobs) == 2

    def test_list_respects_limit(self, client, test_db):
        for _ in range(5):
            test_db.create_job(JobParams())
        import metagomics2.server.app as app_module

        token = asyncio.run(app_module.admin_login(app_module.AdminAuthRequest(password="testpass"))).token
        data = asyncio.run(app_module.list_jobs(3, token))
        assert len(data.jobs) == 3


class TestGetPeptideLists:
    """Tests for peptide lists endpoint."""

    def test_get_peptide_lists(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.add_peptide_list(job_id, "list_000", "pep1.tsv", "/path/pep1.tsv")
        test_db.add_peptide_list(job_id, "list_001", "pep2.tsv", "/path/pep2.tsv")

        import metagomics2.server.app as app_module

        data = asyncio.run(app_module.get_peptide_lists(job_id))
        assert len(data["peptide_lists"]) == 2

    def test_get_peptide_lists_nonexistent_job(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.get_peptide_lists("nonexistent"))

        assert exc_info.value.status_code == 404


class TestDownloadResult:
    """Tests for result file download endpoint."""

    def test_download_existing_file(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())

        # Create result file
        result_dir = tmp_path / "jobs" / job_id / "results" / "list_000"
        result_dir.mkdir(parents=True)
        (result_dir / "coverage.csv").write_text("col1,col2\n1,2\n")

        import metagomics2.server.app as app_module

        response = asyncio.run(app_module.download_result(job_id, "list_000", "coverage.csv"))
        assert response.status_code == 200

    def test_download_nonexistent_job(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.download_result("nonexistent", "list_000", "coverage.csv"))

        assert exc_info.value.status_code == 404

    def test_download_nonexistent_file(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.download_result(job_id, "list_000", "coverage.csv"))

        assert exc_info.value.status_code == 404

    def test_download_disallowed_filename(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.download_result(job_id, "list_000", "secret_data.txt"))

        assert exc_info.value.status_code == 400

    def test_allowed_filenames(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        result_dir = tmp_path / "jobs" / job_id / "results" / "list_000"
        result_dir.mkdir(parents=True)

        allowed = ["taxonomy_nodes.csv", "go_terms.csv", "coverage.csv",
                    "run_manifest.json", "peptides_annotated.csv"]
        for fname in allowed:
            (result_dir / fname).write_text("test content")
            import metagomics2.server.app as app_module

            response = asyncio.run(app_module.download_result(job_id, "list_000", fname))
            assert response.status_code == 200, f"Failed for {fname}"


class TestDownloadAllResults:
    """Tests for ZIP download endpoint."""

    def test_download_all_completed_job(self, client, test_db, tmp_path: Path):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.COMPLETED)

        result_dir = tmp_path / "jobs" / job_id / "results"
        result_dir.mkdir(parents=True)
        (result_dir / "list_000").mkdir()
        (result_dir / "list_000" / "coverage.csv").write_text("test")

        import metagomics2.server.app as app_module

        response = asyncio.run(app_module.download_all_results(job_id))
        assert response.status_code == 200
        assert response.media_type == "application/zip"

    def test_download_all_incomplete_job(self, client, test_db):
        job_id = test_db.create_job(JobParams())
        test_db.update_job_status(job_id, JobStatus.RUNNING)

        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.download_all_results(job_id))

        assert exc_info.value.status_code == 400

    def test_download_all_nonexistent_job(self, client):
        import metagomics2.server.app as app_module

        with pytest.raises(HTTPException) as exc_info:
            asyncio.run(app_module.download_all_results("nonexistent"))

        assert exc_info.value.status_code == 404
