# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Metagomics 2 is a metaproteomics annotation tool. Given a background proteome FASTA and one or more peptide lists with quantities, it matches peptides to proteins, runs DIAMOND against an annotated database (UniProt), transfers NCBI taxonomy and GO annotations back to each peptide, and aggregates quantities into taxonomy nodes and GO terms. It runs as a CLI or as a web app (FastAPI server + polling worker) in one Docker container.

The two specification documents are authoritative and contain developer directives you must follow:
`docs/SPECIFICATION_BACKEND.md` and `docs/SPECIFICATION_FRONTEND.md`. Read the relevant one before non-trivial changes and update it when you change behaviour it describes.

## Commands

Backend (always use the project venv, never a global Python):

```bash
./venv/bin/python -m pip install -e ".[dev]"          # one-time setup
./venv/bin/python -m pytest tests/                    # full suite (~30 s)
./venv/bin/python -m pytest tests/unit/test_diamond.py                       # one file
./venv/bin/python -m pytest tests/unit/test_diamond.py -k heartbeat          # one test by name
./venv/bin/python -m pytest -m "not slow"             # skip slow tests
./venv/bin/ruff check src tests                       # lint (must be clean)
./venv/bin/mypy src                                   # strict type check (must be clean)
```

CI runs ruff, mypy, then pytest on every push. All three must pass before committing.

DIAMOND is not part of the Python package. In this checkout the pinned binary lives at `venv/bin/diamond`; the Dockerfile pins the same version via `DIAMOND_VERSION`. The homology stage will fail without it on `PATH`.

Frontend: do **not** run npm on the host. Everything runs through the Docker `frontend-builder` stage:

```bash
docker build --target frontend-builder -t metagomics2-frontend-test .
docker run --rm metagomics2-frontend-test npx tsc --noEmit
docker run --rm metagomics2-frontend-test npx vitest run
```

Running locally without Docker (needs `config/databases.json` copied from the example):

```bash
METAGOMICS_CONFIG_DIR=./config METAGOMICS_DATA_DIR=./_data uvicorn metagomics2.server.app:app --port 8000
METAGOMICS_CONFIG_DIR=./config METAGOMICS_DATA_DIR=./_data python -m metagomics2.worker.worker
```

CLI without a real DIAMOND database: pass `--mock-hits tests/fixtures/hits/accepted_hits.json --mock-annotations tests/fixtures/annotations/subjects.json` and the homology stage is skipped. The integration tests do exactly this.

Releases: bump the version in `pyproject.toml`, `Dockerfile` (`ARG VERSION`), `src/metagomics2/__init__.py`, and `frontend/package.json`; tag `vX.Y.Z`; publish a GitHub release. Publishing triggers `.github/workflows/release-docker-image.yml`, which builds and pushes `ghcr.io/mriffle/metagomics2` with the version taken from the tag.

## Architecture

**One pipeline, two front ends.** `pipeline/runner.py` (`PipelineRunner`, `PipelineConfig`) is the only place the analysis is orchestrated. `cli.py` builds a `PipelineConfig` from argparse; `worker/worker.py` builds one from a job row. Never put analysis logic in either front end.

**Pipeline stages** (all in `core/`): parse peptide lists → one Aho-Corasick pass matching every peptide from every list against the FASTA (`matching.py`) → write a subset FASTA of hit proteins → DIAMOND blastp on that subset (`diamond.py`) → threshold filters plus tie-aware top-k by bitscore (`filtering.py`) → look up tax IDs and GO terms for hit subjects in the companion SQLite `.annotations.db` (`subject_lookup.py`) → per list: LCA of all implied subjects for taxonomy and union of GO closures for function (`annotation.py`) → aggregate quantities into nodes (`aggregation.py`) → write CSVs, a per-peptide Parquet, and a provenance manifest (`reporting.py`). The shared stages run once for all lists; only annotation onward is per list.

**Scientific correctness is the top priority.** LCA, GO closure, tie-aware top-k, and aggregation invariants are covered by unit tests and Hypothesis property tests in `tests/property/`. Changes to `core/` need tests at the same level.

**Configuration** is loaded once by `config.py` into a frozen `Settings` (`get_settings()`). Scalars come from env vars, structured data from `config/databases.json` and optional `config/server.json`. Both `server/app.py` and `worker/worker.py` call `get_settings()` at import time and copy values into module constants; tests therefore set env vars, call `reset_settings()`, and `importlib.reload` the module (see the autouse fixture in `tests/unit/test_worker.py`). Startup fails on invalid config rather than falling back silently. New settings must be added to `.env.example` with a comment explaining their effect, and passed through in `docker-compose.example.yml`.

**Web mode.** The server writes uploads to `jobs/<id>/inputs/`, creates the job in SQLite (`db/database.py`) and sets it `queued`. The worker polls for queued jobs, marks `running`, runs the pipeline with a progress callback that updates the job row and records a `stage` event on every stage change, then marks `completed`/`failed`, emails if configured, and deletes `inputs/` and `work/` (cleanup is configurable). `results/` and `logs/` persist. The frontend reads result CSVs and queries the Parquet file in the browser with DuckDB-WASM.

**Process model in Docker.** `docker-entrypoint.sh` runs the worker and uvicorn as siblings and supervises both; if either exits, the container exits so `restart: unless-stopped` restarts it. On startup the worker fails any job left in `running` by a dead predecessor. The web server has no knowledge of the worker beyond the shared SQLite file.

**Logging** goes through `logging_setup.configure_logging`: stderr for `docker logs` plus a rotating file under `$METAGOMICS_DATA_DIR/logs/`. The worker also attaches a per-job `jobs/<id>/logs/pipeline.log` for the duration of a job. DIAMOND's console output is streamed to `jobs/<id>/logs/diamond.log` (or the work dir in CLI mode) and a heartbeat is logged while it runs. Use `logging.getLogger(__name__)` everywhere; never `print` in library code.

**DIAMOND behaviour worth knowing.** It processes the database in blocks with mostly single-threaded per-block setup and writes no output until all blocks are joined, so long silences with an empty results file are normal on large databases. Block size, index chunks, and tmpdir are configurable (`METAGOMICS_DIAMOND_*`) and the memory budget is roughly 6 GB per block-size unit; the wrapper logs and warns about this. DIAMOND's own `--max-target-seqs` default is 25 and its cut is not tie-aware, so the runner always passes `PipelineConfig.effective_max_target_seqs()` (default 500, raised to `top_k`, `0` = unlimited) and warns when a query fills the cap.

**Reference data.** GO OBO and NCBI taxonomy dumps are baked into the image at `/app/reference` with version files. Web jobs hardlink a snapshot into `work/ref_snapshot/` and record hashes in the manifest. `core/reference_loader.py` accepts either the raw formats or the small JSON forms used by tests (`tests/fixtures/`).

**Frontend conventions** (from the frontend spec): TypeScript strict, Tailwind with mandatory `dark:` variants, all colours consumed by Plotly/Cytoscape defined in `frontend/src/utils/colors.ts`, icons only from `lucide-react`, and every `utils/` parser has a test file.

## Working conventions

- Commit after each substantive change or milestone with a descriptive message; push when the task is done.
- Keep behavioural changes and mechanical cleanups (formatting, typing) in separate commits.
- Line length is 100 (ruff). Test fixture strings holding real UniProt or GAF records may exceed it with a `# noqa: E501` on the closing line.
