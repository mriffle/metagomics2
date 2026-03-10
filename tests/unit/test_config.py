"""Unit tests for the centralized configuration module."""

import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

from metagomics2.config import (
    DatabaseEntry,
    Settings,
    SmtpSettings,
    _load_databases_json,
    _load_server_json,
    _parse_bool,
    load_settings,
    reset_settings,
)


# ---------------------------------------------------------------------------
# _parse_bool
# ---------------------------------------------------------------------------

class TestParseBool:
    def test_true_values(self):
        for v in ("true", "True", "TRUE", "1", "yes", "YES"):
            assert _parse_bool(v) is True

    def test_false_values(self):
        for v in ("false", "False", "0", "no", "", "anything"):
            assert _parse_bool(v) is False


# ---------------------------------------------------------------------------
# _load_databases_json
# ---------------------------------------------------------------------------

class TestLoadDatabasesJson:
    def test_valid_file(self, tmp_path: Path):
        db_file = tmp_path / "databases.json"
        db_file.write_text(json.dumps([
            {
                "name": "SwissProt",
                "description": "Reviewed entries",
                "path": "sp.dmnd",
                "annotations": "sp.annotations.db",
            }
        ]))
        entries = _load_databases_json(db_file)
        assert len(entries) == 1
        assert entries[0].name == "SwissProt"
        assert entries[0].path == "sp.dmnd"
        assert entries[0].annotations == "sp.annotations.db"

    def test_missing_file(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            _load_databases_json(tmp_path / "nope.json")

    def test_invalid_json(self, tmp_path: Path):
        f = tmp_path / "bad.json"
        f.write_text("not json {{{")
        with pytest.raises(ValueError, match="Invalid JSON"):
            _load_databases_json(f)

    def test_not_a_list(self, tmp_path: Path):
        f = tmp_path / "obj.json"
        f.write_text('{"name": "bad"}')
        with pytest.raises(ValueError, match="JSON array"):
            _load_databases_json(f)

    def test_missing_required_field(self, tmp_path: Path):
        f = tmp_path / "incomplete.json"
        f.write_text(json.dumps([{"name": "X"}]))
        with pytest.raises(ValueError, match="missing required"):
            _load_databases_json(f)

    def test_annotations_defaults_to_empty(self, tmp_path: Path):
        f = tmp_path / "databases.json"
        f.write_text(json.dumps([
            {"name": "DB", "description": "desc", "path": "db.dmnd"}
        ]))
        entries = _load_databases_json(f)
        assert entries[0].annotations == ""

    def test_multiple_entries(self, tmp_path: Path):
        f = tmp_path / "databases.json"
        f.write_text(json.dumps([
            {"name": "A", "description": "a", "path": "a.dmnd"},
            {"name": "B", "description": "b", "path": "b.dmnd", "annotations": "b.ann.db"},
        ]))
        entries = _load_databases_json(f)
        assert len(entries) == 2


# ---------------------------------------------------------------------------
# _load_server_json
# ---------------------------------------------------------------------------

class TestLoadServerJson:
    def test_valid_file(self, tmp_path: Path):
        f = tmp_path / "server.json"
        f.write_text(json.dumps({"allowed_origins": ["https://example.com"]}))
        data = _load_server_json(f)
        assert data["allowed_origins"] == ["https://example.com"]

    def test_missing_file_returns_empty(self, tmp_path: Path):
        assert _load_server_json(tmp_path / "nope.json") == {}

    def test_invalid_json(self, tmp_path: Path):
        f = tmp_path / "bad.json"
        f.write_text("{bad")
        with pytest.raises(ValueError, match="Invalid JSON"):
            _load_server_json(f)

    def test_not_a_dict(self, tmp_path: Path):
        f = tmp_path / "arr.json"
        f.write_text("[]")
        with pytest.raises(ValueError, match="JSON object"):
            _load_server_json(f)


# ---------------------------------------------------------------------------
# load_settings
# ---------------------------------------------------------------------------

class TestLoadSettings:
    def _make_config_dir(self, tmp_path: Path, databases=None, server=None):
        """Helper: write config files and return the config dir."""
        config_dir = tmp_path / "config"
        config_dir.mkdir()
        if databases is not None:
            (config_dir / "databases.json").write_text(json.dumps(databases))
        if server is not None:
            (config_dir / "server.json").write_text(json.dumps(server))
        return config_dir

    def test_loads_from_json_files(self, tmp_path: Path):
        config_dir = self._make_config_dir(
            tmp_path,
            databases=[{"name": "DB", "description": "d", "path": "x.dmnd"}],
        )
        with patch.dict(os.environ, {"METAGOMICS_DATA_DIR": str(tmp_path)}, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert len(settings.databases) == 1
        assert settings.databases[0].name == "DB"

    def test_fails_when_no_databases_and_required(self, tmp_path: Path):
        config_dir = self._make_config_dir(tmp_path, databases=[])
        with patch.dict(os.environ, {"METAGOMICS_DATA_DIR": str(tmp_path)}, clear=False):
            with pytest.raises(RuntimeError, match="No annotated databases"):
                load_settings(config_dir=config_dir)

    def test_no_databases_ok_when_not_required(self, tmp_path: Path):
        config_dir = self._make_config_dir(tmp_path, databases=[])
        with patch.dict(os.environ, {"METAGOMICS_DATA_DIR": str(tmp_path)}, clear=False):
            settings = load_settings(config_dir=config_dir, require_databases=False)
        assert settings.databases == []

    def test_scalar_env_vars(self, tmp_path: Path):
        config_dir = self._make_config_dir(
            tmp_path,
            databases=[{"name": "DB", "description": "d", "path": "x.dmnd"}],
        )
        env = {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_THREADS": "8",
            "METAGOMICS_ADMIN_PASSWORD": "secret",
            "METAGOMICS_MAX_UPLOAD_MB": "512",
            "METAGOMICS_CLEANUP_ON_SUCCESS": "false",
        }
        with patch.dict(os.environ, env, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert settings.threads == 8
        assert settings.admin_password == "secret"
        assert settings.max_upload_mb == 512
        assert settings.max_upload_bytes == 512 * 1024 * 1024
        assert settings.cleanup_on_success is False

    def test_server_json_allowed_origins(self, tmp_path: Path):
        config_dir = self._make_config_dir(
            tmp_path,
            databases=[{"name": "DB", "description": "d", "path": "x.dmnd"}],
            server={"allowed_origins": ["https://a.com", "https://b.com"]},
        )
        with patch.dict(os.environ, {"METAGOMICS_DATA_DIR": str(tmp_path)}, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert settings.allowed_origins == ["https://a.com", "https://b.com"]

    def test_legacy_env_var_fallback(self, tmp_path: Path):
        """When no databases.json exists, fall back to METAGOMICS_DATABASES env var."""
        config_dir = tmp_path / "config"
        config_dir.mkdir()
        # No databases.json file
        env = {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "METAGOMICS_DATABASES": json.dumps([
                {"name": "Legacy", "description": "d", "path": "old.dmnd"}
            ]),
        }
        with patch.dict(os.environ, env, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert len(settings.databases) == 1
        assert settings.databases[0].name == "Legacy"

    def test_smtp_settings(self, tmp_path: Path):
        config_dir = self._make_config_dir(
            tmp_path,
            databases=[{"name": "DB", "description": "d", "path": "x.dmnd"}],
        )
        env = {
            "METAGOMICS_DATA_DIR": str(tmp_path),
            "SMTP_HOST": "smtp.example.com",
            "SMTP_PORT": "465",
            "SMTP_USERNAME": "user",
            "SMTP_PASSWORD": "pass",
            "SMTP_FROM": "noreply@example.com",
        }
        with patch.dict(os.environ, env, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert settings.smtp.host == "smtp.example.com"
        assert settings.smtp.port == 465
        assert settings.smtp.enabled is True

    def test_derived_paths(self, tmp_path: Path):
        config_dir = self._make_config_dir(
            tmp_path,
            databases=[{"name": "DB", "description": "d", "path": "x.dmnd"}],
        )
        with patch.dict(os.environ, {"METAGOMICS_DATA_DIR": str(tmp_path)}, clear=False):
            settings = load_settings(config_dir=config_dir)
        assert settings.jobs_dir == tmp_path / "jobs"
        assert settings.db_path == tmp_path / "metagomics2.db"


# ---------------------------------------------------------------------------
# Settings properties
# ---------------------------------------------------------------------------

class TestSettingsProperties:
    def test_databases_as_dicts(self):
        s = Settings(
            databases=[
                DatabaseEntry(name="A", description="a", path="a.dmnd", annotations="a.ann"),
            ]
        )
        dicts = s.databases_as_dicts
        assert len(dicts) == 1
        assert dicts[0] == {
            "name": "A",
            "description": "a",
            "path": "a.dmnd",
            "annotations": "a.ann",
        }

    def test_max_upload_bytes(self):
        s = Settings(max_upload_mb=10)
        assert s.max_upload_bytes == 10 * 1024 * 1024


class TestSmtpSettings:
    def test_enabled_when_host_set(self):
        s = SmtpSettings(host="smtp.example.com")
        assert s.enabled is True

    def test_disabled_when_host_empty(self):
        s = SmtpSettings()
        assert s.enabled is False
