"""Centralized configuration loading and validation for Metagomics 2.

This module is the single source of truth for all runtime configuration.
It reads scalar values from environment variables / .env files and structured
data (e.g. database definitions) from JSON config files.  The rest of the
application should import ``get_settings()`` and work exclusively with the
validated :class:`Settings` object it returns.

Precedence (highest → lowest):
    1. Explicit function arguments (used by CLI overrides)
    2. Environment variables / .env
    3. JSON config file values
    4. Built-in defaults
"""

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatabaseEntry:
    """A single searchable annotated database."""

    name: str
    description: str
    path: str
    annotations: str = ""


@dataclass(frozen=True)
class SmtpSettings:
    """SMTP configuration for email notifications."""

    host: str = ""
    port: int = 587
    username: str = ""
    password: str = ""
    from_address: str = ""

    @property
    def enabled(self) -> bool:
        return bool(self.host)


@dataclass(frozen=True)
class Settings:
    """Validated, immutable application settings.

    Every component of the application should receive this object rather than
    reading environment variables or config files directly.
    """

    # --- Filesystem locations ---
    data_dir: Path = Path("/data")
    jobs_dir: Path = field(default=Path("/data/jobs"))
    db_path: Path = field(default=Path("/data/metagomics2.db"))
    databases_dir: Path = Path("/databases")

    # --- Annotated databases (loaded from JSON config) ---
    databases: list[DatabaseEntry] = field(default_factory=list)

    # --- Server / runtime ---
    threads: int = 4
    admin_password: str = ""
    max_upload_mb: int = 1024
    diamond_version: str = ""
    site_url: str = ""
    allowed_origins: list[str] = field(default_factory=lambda: ["*"])

    # --- Worker ---
    poll_interval: int = 5
    cleanup_on_success: bool = True
    cleanup_on_failure: bool = True

    # --- Email ---
    smtp: SmtpSettings = field(default_factory=SmtpSettings)

    @property
    def max_upload_bytes(self) -> int:
        return self.max_upload_mb * 1024 * 1024

    @property
    def databases_as_dicts(self) -> list[dict]:
        """Return database entries as plain dicts (for JSON API responses)."""
        return [
            {
                "name": db.name,
                "description": db.description,
                "path": db.path,
                "annotations": db.annotations,
            }
            for db in self.databases
        ]


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def _resolve_path(raw: str, base: Path) -> Path:
    """Resolve *raw* against *base* if it is relative, otherwise return as-is."""
    p = Path(raw)
    if p.is_absolute():
        return p
    return (base / p).resolve()


def _load_databases_json(path: Path) -> list[DatabaseEntry]:
    """Load and validate a databases JSON config file.

    The file must contain a JSON array of objects.  Each object must have at
    least ``name``, ``description``, and ``path`` keys.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If the JSON is malformed or entries are invalid.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Databases config file not found: {path}")

    text = path.read_text(encoding="utf-8")
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc

    if not isinstance(raw, list):
        raise ValueError(
            f"Databases config must be a JSON array, got {type(raw).__name__} in {path}"
        )

    entries: list[DatabaseEntry] = []
    for i, item in enumerate(raw):
        if not isinstance(item, dict):
            raise ValueError(f"Database entry {i} must be an object, got {type(item).__name__}")

        missing = [k for k in ("name", "description", "path") if k not in item]
        if missing:
            raise ValueError(
                f"Database entry {i} is missing required field(s): {', '.join(missing)}"
            )

        entries.append(
            DatabaseEntry(
                name=item["name"],
                description=item["description"],
                path=item["path"],
                annotations=item.get("annotations", ""),
            )
        )

    return entries


def _load_server_json(path: Path) -> dict:
    """Load an optional server JSON config file.

    Returns an empty dict if the file does not exist.

    Raises:
        ValueError: If the file exists but contains invalid JSON.
    """
    if not path.is_file():
        return {}

    text = path.read_text(encoding="utf-8")
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc

    if not isinstance(data, dict):
        raise ValueError(
            f"Server config must be a JSON object, got {type(data).__name__} in {path}"
        )
    return data


def _parse_bool(value: str) -> bool:
    return value.strip().lower() in ("1", "true", "yes")


# ---------------------------------------------------------------------------
# Main loader
# ---------------------------------------------------------------------------

def load_settings(
    *,
    config_dir: Path | None = None,
    databases_json: Path | None = None,
    server_json: Path | None = None,
    require_databases: bool = True,
) -> Settings:
    """Load, merge, and validate all configuration into a :class:`Settings`.

    Call this **once** at application startup.  The returned object should be
    threaded through the rest of the application.

    Args:
        config_dir: Base directory for config files.  Defaults to the
            ``METAGOMICS_CONFIG_DIR`` env var, or ``./config``.
        databases_json: Explicit path to a databases JSON file.  Overrides
            ``METAGOMICS_DATABASES_JSON`` and the default
            ``<config_dir>/databases.json``.
        server_json: Explicit path to a server JSON file.  Overrides
            ``METAGOMICS_SERVER_JSON`` and the default
            ``<config_dir>/server.json``.
        require_databases: If True (default for server/worker mode), raise
            if no databases are configured.  CLI mode sets this to False
            because databases are provided via CLI args.

    Returns:
        A fully validated :class:`Settings` instance.

    Raises:
        RuntimeError: On configuration errors (missing files, bad JSON, etc.).
    """

    # --- Resolve config directory ---
    if config_dir is None:
        config_dir = Path(os.environ.get("METAGOMICS_CONFIG_DIR", "./config"))
    config_dir = config_dir.resolve()

    # --- Scalar env vars ---
    data_dir = Path(os.environ.get("METAGOMICS_DATA_DIR", "/data"))
    databases_dir = Path(os.environ.get("METAGOMICS_DATABASES_DIR", "/databases"))
    threads = int(os.environ.get("METAGOMICS_THREADS", "4"))
    admin_password = os.environ.get("METAGOMICS_ADMIN_PASSWORD", "")
    max_upload_mb = int(os.environ.get("METAGOMICS_MAX_UPLOAD_MB", "1024"))
    diamond_version = os.environ.get("DIAMOND_VERSION", "")
    site_url = os.environ.get("SITE_URL", "")
    poll_interval = int(os.environ.get("METAGOMICS_POLL_INTERVAL", "5"))
    cleanup_on_success = _parse_bool(os.environ.get("METAGOMICS_CLEANUP_ON_SUCCESS", "true"))
    cleanup_on_failure = _parse_bool(os.environ.get("METAGOMICS_CLEANUP_ON_FAILURE", "true"))

    # --- SMTP ---
    smtp = SmtpSettings(
        host=os.environ.get("SMTP_HOST", ""),
        port=int(os.environ.get("SMTP_PORT", "587")),
        username=os.environ.get("SMTP_USERNAME", ""),
        password=os.environ.get("SMTP_PASSWORD", ""),
        from_address=os.environ.get("SMTP_FROM", ""),
    )

    # --- Databases (from JSON config file) ---
    databases: list[DatabaseEntry] = []
    errors: list[str] = []

    # Determine databases JSON path
    if databases_json is None:
        databases_json_env = os.environ.get("METAGOMICS_DATABASES_JSON", "")
        if databases_json_env:
            databases_json = _resolve_path(databases_json_env, config_dir)
        else:
            databases_json = config_dir / "databases.json"

    if databases_json.is_file():
        try:
            databases = _load_databases_json(databases_json)
            logger.info("Loaded %d database(s) from %s", len(databases), databases_json)
        except (ValueError, FileNotFoundError) as exc:
            errors.append(str(exc))
    else:
        # Fall back: check legacy METAGOMICS_DATABASES env var (JSON string)
        legacy_raw = os.environ.get("METAGOMICS_DATABASES", "")
        if legacy_raw and legacy_raw.strip() not in ("", "[]"):
            try:
                legacy_list = json.loads(legacy_raw)
                if isinstance(legacy_list, list) and legacy_list:
                    for item in legacy_list:
                        databases.append(
                            DatabaseEntry(
                                name=item.get("name", ""),
                                description=item.get("description", ""),
                                path=item.get("path", ""),
                                annotations=item.get("annotations", ""),
                            )
                        )
                    logger.info(
                        "Loaded %d database(s) from legacy METAGOMICS_DATABASES env var",
                        len(databases),
                    )
            except (json.JSONDecodeError, TypeError, AttributeError):
                errors.append(
                    "METAGOMICS_DATABASES env var contains invalid JSON. "
                    "Consider migrating to a databases.json config file."
                )

    if require_databases and not databases:
        errors.append(
            "No annotated databases configured. "
            "Create a databases.json config file or set the METAGOMICS_DATABASES_JSON "
            "environment variable. See config/databases.example.json for the expected format."
        )

    # --- Server config (optional JSON file) ---
    allowed_origins: list[str] = ["*"]

    if server_json is None:
        server_json_env = os.environ.get("METAGOMICS_SERVER_JSON", "")
        if server_json_env:
            server_json = _resolve_path(server_json_env, config_dir)
        else:
            server_json = config_dir / "server.json"

    if server_json.is_file():
        try:
            server_data = _load_server_json(server_json)
            if "allowed_origins" in server_data:
                origins = server_data["allowed_origins"]
                if isinstance(origins, list) and all(isinstance(o, str) for o in origins):
                    allowed_origins = origins
                else:
                    errors.append(
                        "server.json 'allowed_origins' must be a list of strings"
                    )
            logger.info("Loaded server config from %s", server_json)
        except ValueError as exc:
            errors.append(str(exc))

    # --- Fail fast on errors ---
    if errors:
        msg = "Configuration error(s):\n" + "\n".join(f"  - {e}" for e in errors)
        raise RuntimeError(msg)

    # --- Build final settings ---
    jobs_dir = data_dir / "jobs"
    db_path = data_dir / "metagomics2.db"

    return Settings(
        data_dir=data_dir,
        jobs_dir=jobs_dir,
        db_path=db_path,
        databases_dir=databases_dir,
        databases=databases,
        threads=threads,
        admin_password=admin_password,
        max_upload_mb=max_upload_mb,
        diamond_version=diamond_version,
        site_url=site_url,
        allowed_origins=allowed_origins,
        poll_interval=poll_interval,
        cleanup_on_success=cleanup_on_success,
        cleanup_on_failure=cleanup_on_failure,
        smtp=smtp,
    )


# ---------------------------------------------------------------------------
# Singleton accessor
# ---------------------------------------------------------------------------

_settings: Settings | None = None


def get_settings() -> Settings:
    """Return the application-wide settings, loading them on first call.

    The settings are loaded once and cached.  To force a reload (e.g. in
    tests), call :func:`reset_settings` first.
    """
    global _settings
    if _settings is None:
        _settings = load_settings()
    return _settings


def set_settings(settings: Settings) -> None:
    """Replace the cached settings (useful for tests and CLI overrides)."""
    global _settings
    _settings = settings


def reset_settings() -> None:
    """Clear cached settings so the next :func:`get_settings` reloads."""
    global _settings
    _settings = None
