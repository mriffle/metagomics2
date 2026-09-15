"""Shared logging configuration for the worker, web server and CLI.

Every entry point calls :func:`configure_logging` once at startup.  Log lines
always go to stderr (so ``docker logs`` shows them) and, when a log directory
is given, to a rotating file under it as well.  Per-job log files are attached
with :func:`attach_file_handler` for the duration of a job.

The module also exposes a few small helpers for diagnostic logging: the
container's cgroup memory limit, the process's peak resident memory, and a
human-readable byte formatter.
"""

import logging
import logging.handlers
import os
import platform
import resource
import shutil
import sys
from pathlib import Path

LOG_FORMAT = "%(asctime)s %(process)d %(name)s %(levelname)s %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_MAX_BYTES = 20 * 1024 * 1024
_BACKUP_COUNT = 5

# Handlers installed by configure_logging(), so a second call can replace them.
_installed: list[logging.Handler] = []

_CGROUP_V2_LIMIT = Path("/sys/fs/cgroup/memory.max")
_CGROUP_V1_LIMIT = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")

# cgroup v1 reports "no limit" as a very large number close to 2**63.
_CGROUP_V1_UNLIMITED_THRESHOLD = 2**62


def _coerce_level(level: str | int) -> int:
    if isinstance(level, int):
        return level
    resolved = logging.getLevelName(level.upper())
    if isinstance(resolved, int):
        return resolved
    raise ValueError(f"Unknown log level: {level!r}")


def configure_logging(
    component: str,
    log_dir: Path | None = None,
    level: str | int = "INFO",
) -> logging.Logger:
    """Configure the root logger for one process.

    Args:
        component: Short name used for the log file (``worker``, ``server``, ``cli``).
        log_dir: Directory for the rotating ``<component>.log`` file.  ``None``
            disables file logging.
        level: Root log level, as a name or numeric level.

    Returns:
        The root logger.

    Calling this again replaces the handlers installed by the previous call, so
    tests and reloads do not accumulate duplicate output.
    """
    root = logging.getLogger()
    root.setLevel(_coerce_level(level))

    for handler in _installed:
        root.removeHandler(handler)
        handler.close()
    _installed.clear()

    formatter = logging.Formatter(LOG_FORMAT, DATE_FORMAT)

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    root.addHandler(stream_handler)
    _installed.append(stream_handler)

    if log_dir is not None:
        try:
            log_dir.mkdir(parents=True, exist_ok=True)
            file_handler = logging.handlers.RotatingFileHandler(
                log_dir / f"{component}.log",
                maxBytes=_MAX_BYTES,
                backupCount=_BACKUP_COUNT,
                encoding="utf-8",
            )
            file_handler.setFormatter(formatter)
            root.addHandler(file_handler)
            _installed.append(file_handler)
        except OSError as e:
            root.warning(f"File logging disabled: cannot open log file in {log_dir}: {e}")

    return root


def attach_file_handler(path: Path, level: str | int = "INFO") -> logging.Handler | None:
    """Add a plain file handler to the root logger, e.g. for one job's log.

    Returns the handler so the caller can pass it to :func:`detach_handler`,
    or ``None`` if the file could not be opened (a warning is logged).
    """
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(path, encoding="utf-8")
    except OSError as e:
        logging.getLogger(__name__).warning(f"Cannot open log file {path}: {e}")
        return None
    handler.setLevel(_coerce_level(level))
    handler.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))
    logging.getLogger().addHandler(handler)
    return handler


def detach_handler(handler: logging.Handler | None) -> None:
    """Remove a handler added by :func:`attach_file_handler` and close it."""
    if handler is None:
        return
    logging.getLogger().removeHandler(handler)
    handler.close()


def format_bytes(n: int | float | None) -> str:
    """Format a byte count as a short human-readable string."""
    if n is None:
        return "unknown"
    value = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(value) < 1024 or unit == "TB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TB"


def get_cgroup_memory_limit(
    v2_path: Path = _CGROUP_V2_LIMIT,
    v1_path: Path = _CGROUP_V1_LIMIT,
) -> int | None:
    """Return the container's memory limit in bytes, or ``None`` if unlimited/unknown."""
    for path in (v2_path, v1_path):
        try:
            text = path.read_text().strip()
        except OSError:
            continue
        if text == "max":
            return None
        try:
            value = int(text)
        except ValueError:
            continue
        if value >= _CGROUP_V1_UNLIMITED_THRESHOLD:
            return None
        return value
    return None


def get_total_memory() -> int | None:
    """Return total physical memory in bytes, or ``None`` if unavailable."""
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError, AttributeError):
        return None
    if pages <= 0 or page_size <= 0:
        return None
    return int(pages) * int(page_size)


def peak_rss_bytes() -> int:
    """Peak resident set size of this process in bytes."""
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports kilobytes; macOS reports bytes.
    if sys.platform == "darwin":
        return int(maxrss)
    return int(maxrss) * 1024


def process_rss_bytes(pid: int) -> int | None:
    """Current resident set size of another process (Linux only), or ``None``."""
    try:
        with open(f"/proc/{pid}/status", encoding="utf-8") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def log_system_info(logger: logging.Logger) -> None:
    """Log the facts most often needed when diagnosing a stuck or dead job."""
    try:
        cpus: int | None = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count()
    limit = get_cgroup_memory_limit()
    logger.info(
        f"System: python {platform.python_version()} on {platform.platform()}; "
        f"{cpus} CPUs available to this process; "
        f"total memory {format_bytes(get_total_memory())}; "
        f"container memory limit {format_bytes(limit) if limit else 'none'}"
    )
    diamond = shutil.which("diamond")
    logger.info(f"DIAMOND executable: {diamond or 'NOT FOUND on PATH'}")
