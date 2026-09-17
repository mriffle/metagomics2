"""DIAMOND homology search execution and result parsing."""

import logging
import re
import signal
import subprocess
import time
from collections.abc import Iterator
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from metagomics2.core.filtering import HomologyHit, parse_blast_tabular
from metagomics2.logging_setup import (
    format_bytes,
    get_cgroup_memory_limit,
    get_total_memory,
    process_rss_bytes,
)

logger = logging.getLogger(__name__)

# Matches UniProt-style subject IDs: db|ACCESSION|ENTRY_NAME
_UNIPROT_ID_RE = re.compile(r"^[a-z]{2}\|([A-Za-z0-9_-]+)\|")

# Columns requested from DIAMOND with ``--outfmt 6``.  The first twelve are
# the BLAST tabular defaults; ``qcovhsp`` (percent of the query covered by
# the HSP) is added because the ``min_qcov`` filter needs it and plain
# ``--outfmt 6`` does not include any coverage column.
DIAMOND_OUTFMT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp",
]

# Default for ``--max-target-seqs``.  DIAMOND's own default is 25, which
# silently truncates the hit list before the pipeline's tie-aware top_k filter
# ever sees it.  500 is far above any realistic top_k while still bounding the
# output size on TrEMBL-scale databases.  0 means unlimited.
DEFAULT_MAX_TARGET_SEQS = 500


def parse_uniprot_accession(subject_id: str) -> str:
    """Extract the bare UniProt accession from a DIAMOND subject ID.

    Handles formats like:
        sp|Q21HH2|RS2_SACD2  -> Q21HH2
        tr|A0A0A0MQG0|...    -> A0A0A0MQG0
        P12345                -> P12345  (bare accession, returned as-is)

    Args:
        subject_id: Full subject ID string from DIAMOND output

    Returns:
        Bare UniProt accession string
    """
    m = _UNIPROT_ID_RE.match(subject_id)
    if m:
        return m.group(1)
    # Assume it's already a bare accession
    return subject_id


class DiamondError(Exception):
    """Raised when DIAMOND execution fails."""

    pass


@dataclass
class DiamondResult:
    """Result of a DIAMOND search."""

    hits_by_query: dict[str, list[HomologyHit]]
    output_path: Path
    n_queries: int
    n_hits: int
    command: list[str] = field(default_factory=list)


def _tail(path: Path, n: int = 30) -> str:
    """Return the last ``n`` non-empty lines of a text file (best effort)."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = [line.rstrip() for line in f if line.strip()]
    except OSError:
        return ""
    return "\n".join(lines[-n:])


def _last_line(path: Path) -> str:
    tail = _tail(path, 1)
    return tail.splitlines()[-1] if tail else ""


def _file_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _wait_with_heartbeat(
    proc: subprocess.Popen[bytes],
    log_path: Path,
    output_path: Path,
    heartbeat_seconds: float,
    start: float,
) -> int:
    """Wait for DIAMOND to exit, logging a heartbeat while it runs.

    Each heartbeat reports elapsed time, the size of the output file, the
    DIAMOND process's resident memory, and the last line DIAMOND wrote to its
    console log, so ``docker logs`` shows which block it is working on.
    """
    while True:
        try:
            return proc.wait(timeout=heartbeat_seconds)
        except subprocess.TimeoutExpired:
            elapsed = time.monotonic() - start
            logger.info(
                f"DIAMOND still running (pid {proc.pid}): elapsed {elapsed:.0f}s, "
                f"rss {format_bytes(process_rss_bytes(proc.pid))}, "
                f"output {format_bytes(_file_size(output_path))}, "
                f"last console line: {_last_line(log_path) or '(none yet)'}"
            )


def run_diamond(
    query_fasta: Path,
    db_path: Path,
    output_path: Path,
    evalue: float = 1e-10,
    max_target_seqs: int | None = None,
    threads: int = 4,
    log_path: Path | None = None,
    heartbeat_seconds: float = 60.0,
    block_size: float | None = None,
    index_chunks: int | None = None,
    tmpdir: Path | None = None,
) -> DiamondResult:
    """Run DIAMOND blastp and parse the results.

    DIAMOND's console output (progress per query/reference block, timings,
    errors) is appended to ``log_path`` rather than captured in memory, so it
    survives a crash and can be inspected while the search is running.

    Args:
        query_fasta: Path to the query FASTA file (subset of background proteome)
        db_path: Path to the DIAMOND-formatted database (.dmnd)
        output_path: Path to write the tabular output
        evalue: Maximum e-value threshold for DIAMOND search
        max_target_seqs: Maximum number of target sequences per query
            (``--max-target-seqs``; 0 means unlimited).  If None the flag is
            not passed and DIAMOND uses its own default of 25 per query,
            which is too few for the pipeline; callers should pass
            ``DEFAULT_MAX_TARGET_SEQS`` or a larger value.
        threads: Number of CPU threads to use
        log_path: File that receives DIAMOND's console output.  Defaults to
            ``diamond.log`` next to the output file.
        heartbeat_seconds: Interval between progress log lines while waiting.
        block_size: DIAMOND ``--block-size`` in billions of letters.  The
            main control over DIAMOND's memory use and speed; ``None`` keeps
            DIAMOND's default (2.0).
        index_chunks: DIAMOND ``--index-chunks``; ``None`` keeps the default (4).
        tmpdir: DIAMOND ``--tmpdir`` for its intermediate files; ``None`` lets
            DIAMOND use the output file's directory.

    Returns:
        DiamondResult with parsed hits

    Raises:
        DiamondError: If DIAMOND cannot be started, exits non-zero, or is
            killed by a signal.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if log_path is None:
        log_path = output_path.parent / "diamond.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "diamond", "blastp",
        "--query", str(query_fasta),
        "--db", str(db_path),
        "--outfmt", "6", *DIAMOND_OUTFMT_COLUMNS,
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--out", str(output_path),
    ]

    if max_target_seqs is not None:
        cmd.extend(["--max-target-seqs", str(max_target_seqs)])
    if block_size is not None:
        cmd.extend(["--block-size", f"{block_size:g}"])
    if index_chunks is not None:
        cmd.extend(["--index-chunks", str(index_chunks)])
    if tmpdir is not None:
        tmpdir.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--tmpdir", str(tmpdir)])

    logger.info(f"Running DIAMOND: {' '.join(cmd)}")
    _log_memory_expectation(block_size)
    logger.info(
        f"DIAMOND inputs: query {format_bytes(_file_size(query_fasta))}, "
        f"database {format_bytes(_file_size(db_path))}; console output -> {log_path}"
    )

    start = time.monotonic()
    try:
        with open(log_path, "ab") as log_file:
            started_at = datetime.now().isoformat(timespec="seconds")
            log_file.write(f"# {started_at} command: {' '.join(cmd)}\n".encode())
            log_file.flush()
            proc = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT)
            returncode = _wait_with_heartbeat(
                proc, log_path, output_path, heartbeat_seconds, start
            )
    except FileNotFoundError:
        raise DiamondError(
            "DIAMOND executable not found. Ensure 'diamond' is installed and on PATH."
        )
    elapsed = time.monotonic() - start

    if returncode != 0:
        tail = _tail(log_path)
        if returncode < 0:
            signum = -returncode
            try:
                sig_name = signal.Signals(signum).name
            except ValueError:
                sig_name = "unknown"
            message = (
                f"DIAMOND was killed by signal {signum} ({sig_name}) after {elapsed:.0f}s. "
                "This usually means it was killed externally, for example by the kernel "
                "out-of-memory killer or a container memory limit."
            )
        else:
            message = f"DIAMOND exited with code {returncode} after {elapsed:.0f}s"
        logger.error(message)
        raise DiamondError(f"{message}. Last lines of {log_path}:\n{tail}")

    logger.info(
        f"DIAMOND completed in {elapsed:.0f}s (exit code 0). "
        f"Output: {output_path} ({format_bytes(_file_size(output_path))})"
    )

    # Parse results
    result = parse_diamond_output(output_path)
    result.command = cmd
    return result


def count_queries_at_cap(hits_by_query: dict[str, list[HomologyHit]], cap: int) -> int:
    """Number of queries whose hit count reached DIAMOND's per-query cap.

    A query with exactly ``cap`` hits may have had further hits, including
    ties at the cutoff, discarded by DIAMOND.  Returns 0 when ``cap`` is 0
    (unlimited).
    """
    if cap <= 0:
        return 0
    return sum(1 for hits in hits_by_query.values() if len(hits) >= cap)


# DIAMOND's documentation says to expect roughly six times the block size in
# gigabytes of memory.  In practice it is often less, so this is a budget, not
# a prediction.
_DIAMOND_GB_PER_BLOCK_UNIT = 6.0
_DIAMOND_DEFAULT_BLOCK_SIZE = 2.0


def estimate_diamond_memory_bytes(block_size: float | None) -> int:
    """Upper-bound estimate of DIAMOND's memory use for a block size."""
    effective = block_size if block_size is not None else _DIAMOND_DEFAULT_BLOCK_SIZE
    return int(effective * _DIAMOND_GB_PER_BLOCK_UNIT * 1024**3)


def _log_memory_expectation(block_size: float | None) -> None:
    """Log the expected DIAMOND memory budget and warn if it exceeds what is available."""
    expected = estimate_diamond_memory_bytes(block_size)
    effective = block_size if block_size is not None else _DIAMOND_DEFAULT_BLOCK_SIZE
    logger.info(
        f"DIAMOND block size {effective:g} billion letters: expect up to about "
        f"{format_bytes(expected)} of memory"
    )
    limit = get_cgroup_memory_limit()
    total = get_total_memory()
    if limit is not None and expected > limit:
        logger.warning(
            f"DIAMOND may need {format_bytes(expected)} but the container memory limit is "
            f"{format_bytes(limit)}; lower METAGOMICS_DIAMOND_BLOCK_SIZE or raise the limit"
        )
    elif total is not None and expected > total:
        logger.warning(
            f"DIAMOND may need {format_bytes(expected)} but the machine has "
            f"{format_bytes(total)}; lower METAGOMICS_DIAMOND_BLOCK_SIZE"
        )


def _iter_lines_with_progress(path: Path, every: int = 1_000_000) -> Iterator[str]:
    """Yield lines from a file, logging a progress line every ``every`` lines."""
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f, start=1):
            if i % every == 0:
                logger.info(f"Parsing DIAMOND output: {i:,} lines read so far")
            yield line


def parse_diamond_output(output_path: Path) -> DiamondResult:
    """Parse DIAMOND tabular output written with ``DIAMOND_OUTFMT_COLUMNS``.

    Args:
        output_path: Path to the DIAMOND output file

    Returns:
        DiamondResult with parsed hits
    """
    if not output_path.exists():
        return DiamondResult(
            hits_by_query={},
            output_path=output_path,
            n_queries=0,
            n_hits=0,
        )

    logger.info(f"Parsing DIAMOND output: {output_path} ({format_bytes(_file_size(output_path))})")
    hits_by_query = parse_blast_tabular(
        _iter_lines_with_progress(output_path), columns=DIAMOND_OUTFMT_COLUMNS
    )

    n_hits = sum(len(hits) for hits in hits_by_query.values())

    logger.info(
        f"Parsed {n_hits} DIAMOND hits for {len(hits_by_query)} query proteins"
    )

    return DiamondResult(
        hits_by_query=hits_by_query,
        output_path=output_path,
        n_queries=len(hits_by_query),
        n_hits=n_hits,
    )
