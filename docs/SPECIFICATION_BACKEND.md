# Metagomics 2 — Backend Specification

## 1. Purpose and Problem Statement

Metagomics 2 is a metaproteomics annotation and aggregation tool. In metaproteomics experiments, researchers identify peptides from complex microbial communities (e.g., gut microbiomes, environmental samples). These peptides come from mass spectrometry and are matched against a **background proteome** (a FASTA file of protein sequences from the sample). The core challenge is: *given a list of peptides with abundances, what organisms and biological functions are represented, and in what proportions?*

Metagomics 2 solves this by:

1. **Matching peptides** to background proteins (exact substring matching)
2. **Searching** those background proteins against a large annotated database (e.g., UniProt) via homology search (DIAMOND blastp)
3. **Transferring annotations** (NCBI taxonomy IDs and Gene Ontology terms) from the annotated database hits back through the chain to each peptide
4. **Aggregating** peptide quantities into taxonomy nodes and GO terms, producing quantitative summaries of "what organisms" and "what functions" are present

The tool supports both a **command-line interface** for batch processing and a **web interface** (FastAPI server + background worker) for interactive use.

---

## 2. Developer Directives

These directives **must** be followed for all backend development:

- **Python virtual environment**: All Python calls and package installations must use the project's virtual environment (`./venv/`). Run tests with `./venv/bin/python -m pytest tests/`. Never use a global `pip` or `python`.
- **Testing is essential**: Examine the existing testing framework and patterns before writing new code. Write tests for all new Python functions and run them before considering a task complete.
- **Docker is the runtime**: This application always runs inside a Docker container. The multi-stage `Dockerfile` handles all build steps. npm installs do not need to happen on the host — update the Dockerfile as needed.
- **Accuracy of analysis is paramount**: The scientific correctness of annotation transfer, LCA computation, GO closure, and aggregation must be verified by tests. Property-based testing with Hypothesis is used to verify mathematical invariants.

---

## 3. High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        Execution Modes                          │
│                                                                 │
│   CLI Mode                          Web Mode                    │
│   (metagomics2 run ...)             (FastAPI + Worker)          │
│         │                                │                      │
│         ▼                                ▼                      │
│   ┌──────────┐                    ┌──────────────┐              │
│   │  cli.py  │                    │  server/app  │◄── HTTP API  │
│   └────┬─────┘                    └──────┬───────┘              │
│        │                                 │ enqueue              │
│        │                                 ▼                      │
│        │                          ┌──────────────┐              │
│        │                          │ worker/worker│ poll loop    │
│        │                          └──────┬───────┘              │
│        │                                 │                      │
│        ▼                                 ▼                      │
│   ┌──────────────────────────────────────────┐                  │
│   │          pipeline/runner.py              │                  │
│   │   (PipelineRunner — shared orchestrator) │                  │
│   └──────────────────┬───────────────────────┘                  │
│                      │                                          │
│                      ▼                                          │
│   ┌──────────────────────────────────────────┐                  │
│   │              core/ modules                │                  │
│   │  fasta, peptides, matching, diamond,     │                  │
│   │  filtering, annotation, aggregation,     │                  │
│   │  reporting, taxonomy, go, ...            │                  │
│   └──────────────────────────────────────────┘                  │
│                                                                 │
│   Supporting layers:                                            │
│   - db/database.py    (SQLite job tracking)                     │
│   - models/job.py     (Pydantic models)                         │
│   - notifications/    (email on job completion)                 │
│   - scripts/          (build_annotations_db.py — thin wrapper)  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. Project Structure (Backend)

```
src/metagomics2/
├── __init__.py                  # Version (__version__ from env or default)
├── config.py                    # Centralized config loader (env vars + JSON files → Settings)
├── cli.py                       # CLI entry point (argparse)
├── core/                        # Core algorithms (pure logic, no I/O side effects beyond file read/write)
│   ├── aggregation.py           # Quantity roll-up into taxonomy/GO nodes
│   ├── annotation.py            # Peptide→annotation logic (LCA, GO union)
│   ├── diamond.py               # DIAMOND execution and output parsing
│   ├── fasta.py                 # FASTA parsing, hashing, subset writing
│   ├── filtering.py             # Homology hit filtering (thresholds + top_k)
│   ├── gaf_parser.py            # GAF 2.2 file parser (GO annotations)
│   ├── go.py                    # GO DAG model and closure computation
│   ├── matching.py              # Aho-Corasick peptide→protein matching
│   ├── ncbi_parser.py           # NCBI taxonomy dump parser
│   ├── obo_parser.py            # OBO format parser (Gene Ontology)
│   ├── peptides.py              # Peptide list parsing and normalization
│   ├── reference_loader.py      # Reference data loading (GO, taxonomy, snapshots)
│   ├── reporting.py             # CSV/Parquet/JSON output generation
│   ├── subject_lookup.py        # SQLite lookup of subject annotations
│   ├── taxonomy.py              # Taxonomy tree model and LCA computation
│   └── uniprot_fasta.py         # UniProt FASTA header parser (accession + OX=taxid)
├── db/
│   └── database.py              # SQLite job-tracking database (web mode)
├── models/
│   └── job.py                   # Pydantic models (JobParams, JobInfo, JobStatus, etc.)
├── notifications/
│   └── email.py                 # SMTP email notifications
├── pipeline/
│   └── runner.py                # Pipeline orchestration (PipelineRunner, PipelineConfig)
├── scripts/
│   └── build_annotations_db.py  # Build companion .annotations.db (installed package module)
├── server/
│   └── app.py                   # FastAPI application (REST API + SPA serving)
└── worker/
    └── worker.py                # Background worker (poll-based job processing)

scripts/
└── build_annotations_db.py      # Thin wrapper — delegates to metagomics2.scripts.build_annotations_db

tests/
├── conftest.py                  # Shared fixtures (fixtures_dir, small_taxonomy, small_go, etc.)
├── unit/                        # 28 unit test files
├── property/                    # 3 property-based test files (Hypothesis)
├── integration/                 # 5 integration test files
└── fixtures/                    # Test data (FASTA, peptides, taxonomy, GO, hits, annotations)
```

---

## 5. The Annotation Pipeline — Complete Workflow

The pipeline is the heart of Metagomics 2. It is orchestrated by `PipelineRunner` in `pipeline/runner.py` and consists of these stages:

### Stage 0: Initialize
- **Load background FASTA** (`core/fasta.py`): Parse the user-provided FASTA file into a `dict[protein_id, sequence]`.
- **Load reference data** (`core/reference_loader.py`): Load GO DAG (from `.obo` or `.json`) and NCBI taxonomy tree (from dump directory or `.json`). In web mode, creates a per-job snapshot of bundled reference data for provenance.
- **Load mock data** (testing only): If `--mock-hits` and `--mock-annotations` are provided, loads pre-computed mappings instead of running DIAMOND.

### Stage 1: Parse Peptide Lists
- **Parse** each peptide list file (`core/peptides.py`): CSV/TSV with two columns (sequence, quantity). Auto-detects delimiter and header. Normalizes sequences (uppercase, strip modifications). Validates against amino acid alphabet. Aggregates duplicate sequences by summing quantities.

### Stage 2: Match Peptides to Background Proteome
- **Aho-Corasick matching** (`core/matching.py`): All peptide sequences from all lists are combined into a single automaton. A single pass over every background protein sequence finds all exact substring matches. Results in `peptide_to_proteins: dict[peptide_seq, set[protein_id]]`.
- Per-list match results are then partitioned from the combined result.
- Produces `all_hit_proteins`: the union of all background proteins that contain at least one peptide.

### Stage 3: Write Subset FASTA
- Write a FASTA file containing only the hit proteins from Stage 2 (`core/fasta.py: write_subset_fasta`). This becomes the query input for DIAMOND.

### Stage 4: Homology Search (DIAMOND)
- **Run DIAMOND blastp** (`core/diamond.py`): Searches the subset FASTA against an annotated database (e.g., UniProt SwissProt `.dmnd`). Output format: BLAST tabular (outfmt 6). The `max_evalue` from filter policy is passed to DIAMOND as a pre-filter.
- **Parse results** (`core/filtering.py: parse_blast_tabular`): Parse the tabular output into `HomologyHit` objects grouped by query protein.
- **Filter hits** (`core/filtering.py: filter_all_hits`):
  1. **Threshold filters** (AND logic): `max_evalue`, `min_pident`, `min_qcov`, `min_alnlen`
  2. **Tie-aware top_k ranking**: Sort by bitscore descending, keep top K, but include *all* hits tied at the Kth-best bitscore. This ensures annotation is never biased by arbitrary tie-breaking.
- Result: `protein_to_subjects: dict[bg_protein_id, set[subject_id]]`

### Stage 4b: Load Subject Annotations
- **SQLite lookup** (`core/subject_lookup.py`): For all unique subject IDs from DIAMOND results, query the companion `.annotations.db` to get taxonomy IDs and GO terms.
- The `.annotations.db` has two tables: `taxonomy(accession, tax_id)` and `go_annotations(accession, go_id, aspect)`.
- Subject IDs from DIAMOND may be in UniProt format (`sp|Q21HH2|RS2_SACD2`); the bare accession (`Q21HH2`) is extracted for database lookup, but the full ID is used as the key in `subject_annotations`.

### Stage 5: Annotate Peptides (per list)
For each peptide in the list (`core/annotation.py`):

1. **Implied subjects**: `S(p) = ⋃_{b∈B(p)} subjects(b)` — union of all annotated-DB subjects across all background proteins that contain this peptide.
2. **Taxonomy annotation (LCA intersection)**: Collect all `tax_id` values from implied subjects. Compute the Lowest Common Ancestor (LCA) using lineage intersection (`core/taxonomy.py: TaxonomyTree.compute_lca`). The LCA is the deepest node that is an ancestor of *all* subject tax_ids. Then `taxonomy_nodes = lineage(LCA → root)`.
3. **GO annotation (non-redundant union)**: Collect all direct GO terms from implied subjects. For each term, compute the transitive closure (ancestors via `is_a` and/or `part_of` edges) using `core/go.py: GODAG.get_closure`. The peptide's GO terms are the union of all closures: `GO(p) = ⋃_{s∈S(p)} closure(go_terms(s))`.

### Stage 6: Aggregate (per list)
- **Taxonomy aggregation** (`core/aggregation.py`): For each taxonomy node `t`: `quantity(t) = Σ_p q(p) * 1[t ∈ TAX(p)]`. Computes `ratio_total = quantity / total_quantity` and `ratio_annotated = quantity / annotated_quantity`.
- **GO aggregation**: Same formula for GO terms.
- **Coverage stats**: Total, annotated, and unannotated peptide quantities and counts.
- **GO-taxonomy cross-tabulation** (`aggregate_go_taxonomy_combos`): For each `(tax_id, go_id)` pair in `TAX(p) × GO(p)`, accumulates quantity. Computes `fraction_of_taxon` and `fraction_of_go`.
- **Invariant validation** (`validate_aggregation_invariants`): Checks that `0 ≤ quantity(n) ≤ total` and `0 ≤ ratio_total ≤ ratio_annotated ≤ 1`. Violations are logged as warnings.

### Stage 7: Write Reports (per list)
Output directory: `<output_dir>/<list_id>/` (e.g., `results/list_000/`)

| File | Format | Description |
|------|--------|-------------|
| `taxonomy_nodes.csv` | CSV | `tax_id, name, rank, parent_tax_id, quantity, ratio_total, ratio_annotated, n_peptides` |
| `go_terms.csv` | CSV | `go_id, name, namespace, parent_go_ids, quantity, ratio_total, ratio_annotated, n_peptides` |
| `go_taxonomy_combo.csv` | CSV | Cross-tab: `tax_id, tax_name, tax_rank, parent_tax_id, go_id, go_name, go_namespace, parent_go_ids, quantity, fraction_of_taxon, fraction_of_go, ratio_total_taxon, ratio_total_go, n_peptides` |
| `coverage.csv` | CSV | Single row: `total_peptide_quantity, annotated_peptide_quantity, unannotated_peptide_quantity, annotation_coverage_ratio, n_peptides_total, n_peptides_annotated, n_peptides_unannotated` |
| `peptide_mapping.parquet` | Parquet | One row per `(peptide, background_protein, annotated_protein)` triple. Columns: `peptide, peptide_lca_tax_ids (List[Int64]), peptide_go_terms (List[Utf8]), background_protein, annotated_protein, evalue, pident` |
| `run_manifest.json` | JSON | Provenance: version, tool versions, input hashes (SHA256), parameters, reference data hashes, timestamp |

---

## 6. Key Data Structures

### `FastaRecord` (`core/fasta.py`)
```python
@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str
```

### `Peptide` (`core/peptides.py`)
```python
@dataclass(frozen=True)
class Peptide:
    sequence: str   # Normalized uppercase amino acid letters only
    quantity: float  # Non-negative abundance/count
```

### `HomologyHit` (`core/filtering.py`)
```python
@dataclass
class HomologyHit:
    query_id: str      # Background protein ID
    subject_id: str    # Annotated DB accession (possibly in sp|ACC|NAME format)
    evalue: float
    bitscore: float
    pident: float      # Percent identity
    qcov: float        # Query coverage
    alnlen: int        # Alignment length
```

### `FilterPolicy` (`core/filtering.py`)
```python
@dataclass
class FilterPolicy:
    max_evalue: float | None = None
    min_pident: float | None = None
    min_qcov: float | None = None
    min_alnlen: int | None = None
    top_k: int | None = None
```
Threshold filters are AND-combined. The `top_k` ranking is tie-aware: all hits tied at the Kth-best bitscore are kept.

### `SubjectAnnotation` (`core/annotation.py`)
```python
@dataclass
class SubjectAnnotation:
    subject_id: str
    tax_id: int | None = None
    go_terms: set[str] = field(default_factory=set)
```

### `PeptideAnnotation` (`core/annotation.py`)
```python
@dataclass
class PeptideAnnotation:
    peptide: str
    quantity: float
    is_annotated: bool = False
    lca_tax_id: int | None = None
    taxonomy_nodes: set[int] = field(default_factory=set)  # Lineage LCA→root
    go_terms: set[str] = field(default_factory=set)         # Union of closures
    implied_subjects: set[str] = field(default_factory=set)
    background_proteins: set[str] = field(default_factory=set)
```

### `TaxonomyTree` / `TaxonNode` (`core/taxonomy.py`)
```python
@dataclass
class TaxonNode:
    tax_id: int
    name: str
    rank: str
    parent_tax_id: int | None

@dataclass
class TaxonomyTree:
    nodes: dict[int, TaxonNode]
    # Methods: get_lineage(), compute_lca(), get_lineage_set()
```
LCA is computed by intersecting lineage sets and finding the deepest common ancestor.

### `GODAG` / `GOTerm` (`core/go.py`)
```python
@dataclass
class GOTerm:
    id: str
    name: str
    namespace: str
    parents: dict[str, set[str]]  # edge_type -> set of parent GO IDs

@dataclass
class GODAG:
    terms: dict[str, GOTerm]
    obsolete_terms: dict[str, GOTerm]
    # Methods: get_closure(), get_closure_union()
```
Closure is computed via DFS/BFS traversal following specified edge types (default: `is_a`).

### `AggregationResult` / `NodeAggregate` (`core/aggregation.py`)
```python
@dataclass
class NodeAggregate:
    node_id: str | int
    quantity: float = 0.0
    n_peptides: int = 0
    contributing_peptides: set[str]
    # Properties: ratio_total, ratio_annotated

@dataclass
class AggregationResult:
    taxonomy_nodes: dict[int, NodeAggregate]
    go_terms: dict[str, NodeAggregate]
    coverage: CoverageStats
```

### `PipelineConfig` (`pipeline/runner.py`)
```python
@dataclass
class PipelineConfig:
    fasta_path: Path
    peptide_list_paths: list[Path]
    output_dir: Path
    search_tool: str = "diamond"
    annotated_db_path: Path | None = None
    annotations_db_path: Path | None = None
    threads: int = 1
    filter_policy: FilterPolicy
    go_data_path: Path | None = None
    taxonomy_data_path: Path | None = None
    job_dir: Path | None = None          # Set for web mode (enables reference snapshots)
    go_edge_types: set[str] = {"is_a", "part_of"}
    go_include_self: bool = True
    mock_hits_path: Path | None = None   # Testing only
    mock_subject_annotations_path: Path | None = None  # Testing only
```

### `PipelineProgress` (`pipeline/runner.py`)
```python
@dataclass
class PipelineProgress:
    total_peptide_lists: int = 0
    completed_peptide_lists: int = 0
    current_stage: str = ""
    current_list_id: str = ""
    progress_done: int = 0          # Weighted progress (0–1000)
    progress_total: int = 1000      # Always 1000
```

Progress uses a **weighted milestone system** (0–1000 scale) so the progress bar reflects overall pipeline progress, not just peptide list completion. Each stage has a fixed milestone value based on its typical runtime weight:

| Stage | Milestone | % |
|-------|-----------|---|
| Initializing | 0 → 50 | 0–5% |
| Parsing peptide lists | 50 → 80 | 5–8% |
| Matching peptides | 80 → 130 | 8–13% |
| Writing subset FASTA | 130 → 150 | 13–15% |
| Homology search (DIAMOND) | 150 → 650 | 15–65% |
| Filtering hits | 650 → 670 | 65–67% |
| Loading subject annotations | 670 → 750 | 67–75% |
| Per-list processing | 750 → 1000 | 75–100% |

DIAMOND is weighted heaviest (~50%) as it is typically the longest-running step. Per-list processing (annotation, aggregation, reporting) is divided evenly among peptide lists within the 750–1000 range.

---

## 7. Core Algorithms

### 7.1 Aho-Corasick Peptide Matching (`core/matching.py`)

Uses the `pyahocorasick` library to build a finite-state automaton from all peptide sequences. A single scan of each background protein sequence finds all peptides that are exact substrings. This is O(N + M + Z) where N = total protein sequence length, M = total peptide length, Z = number of matches. Much faster than searching each peptide individually.

### 7.2 Taxonomy LCA (`core/taxonomy.py`)

Given a set of tax_ids from the implied subjects of a peptide:
1. Compute the lineage (leaf→root) for each tax_id
2. Intersect all lineage sets to find common ancestors
3. The LCA is the deepest node (first one encountered in any lineage) that appears in the intersection

If subjects have tax_ids 70, 71, 72 and their lineages converge at tax_id 30 (ClassA), then LCA = 30, and the peptide's `taxonomy_nodes` = `{30, 20, 10, 1}` (ClassA → PhylumA → KingdomA → root).

### 7.3 GO Closure Union (`core/go.py`)

Given direct GO terms from all implied subjects:
1. For each term, compute transitive closure by traversing `is_a` and/or `part_of` edges upward
2. Union all closures: `GO(p) = ⋃ closure(term)`

This produces a non-redundant set of all GO terms that describe the peptide's function, including inherited ancestor terms.

### 7.4 Tie-Aware Top-K Filtering (`core/filtering.py`)

The top_k filter ranks hits by bitscore descending. At the boundary (Kth position), if multiple hits share the same bitscore, **all tied hits are retained**. Example: `top_k=1` with 5 hits all scoring 200.0 → all 5 kept. This prevents arbitrary bias from tie-breaking.

### 7.5 Aggregation (`core/aggregation.py`)

For each taxonomy node `t`:
```
quantity(t) = Σ_{p annotated} q(p) × 1[t ∈ taxonomy_nodes(p)]
```

Invariants that must always hold:
- `0 ≤ quantity(n) ≤ total_peptide_quantity`
- `0 ≤ ratio_total ≤ 1`
- `0 ≤ ratio_annotated ≤ 1` (when annotated > 0)
- `ratio_total ≤ ratio_annotated` (since `total ≥ annotated`)

These invariants are checked by `validate_aggregation_invariants()` and verified by property-based tests.

---

## 8. CLI (`cli.py`)

Entry point: `metagomics2` (defined in `pyproject.toml` `[project.scripts]`).

### Commands
- `metagomics2 run` — Run the annotation pipeline
- `metagomics2 version` — Show version

### Key `run` Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Background proteome FASTA |
| `--peptides` | Yes (repeatable) | Peptide list CSV/TSV |
| `--outdir` | Yes | Output directory |
| `--db` | Yes* | DIAMOND database (.dmnd) |
| `--annotations-db` | Yes* | Companion SQLite (.annotations.db) |
| `--threads` | No (default: 1) | DIAMOND threads |
| `--search-tool` | No (default: diamond) | Homology search tool |
| `--max-evalue` | No | Maximum e-value |
| `--min-pident` | No | Minimum percent identity |
| `--min-qcov` | No | Minimum query coverage |
| `--min-alnlen` | No | Minimum alignment length |
| `--top-k` | No | Top K hits by bitscore |
| `--params` | No | JSON file with filter params |
| `--go` | No | Path to GO data (OBO/JSON) |
| `--taxonomy` | No | Path to taxonomy data (dir/JSON) |
| `--go-edge-types` | No (default: is_a,part_of) | Comma-separated edge types for GO closure |
| `--go-exclude-self` | No | Exclude terms themselves from closure |
| `--mock-hits` | No | Mock hits JSON (testing) |
| `--mock-annotations` | No | Mock annotations JSON (testing) |

\* Required unless `--mock-hits` / `--mock-annotations` are used for testing.

---

## 9. Web Backend

### 9.1 FastAPI Server (`server/app.py`)

The server provides a REST API and serves the built frontend SPA.

**Configuration**: The server loads all settings through the centralized config module (`config.py`). At import time it calls `get_settings()`, which reads environment variables for scalar values and JSON config files for structured data (database definitions, allowed origins). See Section 16 for the full environment variable reference.

Key settings consumed by the server:
| Setting | Source | Description |
|---------|--------|-------------|
| `data_dir` | `METAGOMICS_DATA_DIR` env var | Persistent data directory |
| `admin_password` | `METAGOMICS_ADMIN_PASSWORD` env var | Admin panel password |
| `databases` | `databases.json` config file | List of annotated database entries |
| `databases_dir` | `METAGOMICS_DATABASES_DIR` env var | Directory containing .dmnd and .annotations.db files |
| `threads` | `METAGOMICS_THREADS` env var | DIAMOND thread count |
| `max_upload_mb` | `METAGOMICS_MAX_UPLOAD_MB` env var | Max upload size |
| `allowed_origins` | `server.json` config file | CORS allowed origins (default: `["*"]`) |

**API Endpoints**:
| Method | Path | Auth | Description |
|--------|------|------|-------------|
| GET | `/api/health` | No | Health check (returns version) |
| GET | `/api/version` | No | Version info |
| GET | `/api/config` | No | Public config (databases, DIAMOND version) |
| POST | `/api/admin/auth` | No | Admin login (returns session token) |
| POST | `/api/jobs` | No | Create job (multipart: FASTA + peptide files + params JSON) |
| GET | `/api/jobs/{job_id}` | No | Get job status and info |
| POST | `/api/jobs/{job_id}/regenerate-id` | No | Regenerate job URL (revoke old link) |
| GET | `/api/admin/jobs` | Admin | List all jobs |
| GET | `/api/jobs/{job_id}/peptide-lists` | No | Get peptide list info for a job |
| GET | `/api/jobs/{job_id}/results/{list_id}/{filename}` | No | Download result file |
| GET | `/api/jobs/{job_id}/results/all_results.zip` | No | Download all results as ZIP |

**Job creation flow**:
1. Validate FASTA content (first 8KB header check)
2. Validate `db_choice` against configured databases
3. Create job record in SQLite
4. Stream-save uploaded files to `<JOBS_DIR>/<job_id>/inputs/`
5. Register peptide lists in database
6. Set job status to `QUEUED`

**Security**:
- Job IDs are cryptographically random URL-safe tokens (128-bit entropy)
- Admin auth uses `secrets.compare_digest` and session tokens stored in memory
- File downloads are restricted to an allowlist of filenames (prevents path traversal)
- Upload size limits enforced per-file

**Frontend serving**: The built SPA is served from `frontend/dist/`. Static assets are mounted at `/assets/`, and all other non-API routes fall back to `index.html` for client-side routing.

### 9.2 Background Worker (`worker/worker.py`)

The worker runs as a separate process (started by `docker-entrypoint.sh`). It uses a **poll-based** design:

1. Poll SQLite for next `QUEUED` job (FIFO by `created_at`)
2. Mark job as `RUNNING`
3. Build `PipelineConfig` from job parameters
4. Run the pipeline with a progress callback that updates the database
5. On success: update per-list status, mark job `COMPLETED`, send email notification
6. On failure: mark job `FAILED`, log error, send email notification
7. Clean up intermediate files (`inputs/`, `work/` directories) based on config
8. Sleep `POLL_INTERVAL` (default 5s) and repeat

**Worker configuration**: The worker also loads settings via `get_settings()` from the centralized config module. Key settings consumed:
| Setting | Source | Description |
|---------|--------|-------------|
| `poll_interval` | `METAGOMICS_POLL_INTERVAL` env var | Seconds between queue polls |
| `cleanup_on_success` | `METAGOMICS_CLEANUP_ON_SUCCESS` env var | Delete inputs/work after success |
| `cleanup_on_failure` | `METAGOMICS_CLEANUP_ON_FAILURE` env var | Delete inputs/work after failure |
| `smtp.*` | `SMTP_HOST`, `SMTP_PORT`, etc. env vars | Email notification config |
| `site_url` | `SITE_URL` env var | Base URL for job links in emails |

### 9.3 Job Database (`db/database.py`)

SQLite database at `<DATA_DIR>/metagomics2.db` with three tables:

**`jobs`**: `job_id (PK), created_at, status, params_json, db_choice, search_tool, progress_total, progress_done, current_step, error_message`

**`peptide_lists`**: `id (PK), job_id (FK), list_id, filename, path, status, n_peptides, n_matched, n_unmatched`

**`job_events`**: `id (PK), job_id (FK), timestamp, event_type, message`

Job statuses: `uploaded → queued → running → completed | failed`

The `regenerate_job_id()` method atomically updates the job_id across all tables and renames the filesystem directory, enabling users to revoke shared links.

### 9.4 Job Models (`models/job.py`)

Pydantic models with field validators:
- `JobParams`: Validates `max_evalue > 0`, `min_pident ∈ [0,100]`, `min_qcov ∈ [0,100]`, `min_alnlen ≥ 1`, `top_k ≥ 1`, `db_choice` is a plain filename (no path traversal), `notification_email` matches basic email regex.
- `JobInfo`: Response model with nested `PeptideListInfo` objects.
- `JobStatus`: Enum (`uploaded`, `queued`, `running`, `completed`, `failed`).

### 9.5 Email Notifications (`notifications/email.py`)

Sends plain-text email on job completion/failure via SMTP with STARTTLS. The `send_job_notification` function never raises — SMTP errors are logged and swallowed to prevent email failures from breaking the pipeline. Email body includes: status, uploaded filenames, parameters, and a link to view results.

---

## 10. Annotated Database System

Metagomics 2 uses a **two-file system** for each annotated database:

1. **`.dmnd`**: DIAMOND-formatted protein database (built with `diamond makedb`)
2. **`.annotations.db`**: Companion SQLite database mapping accessions to taxonomy IDs and GO terms

### Building the Annotations Database (`metagomics2.scripts.build_annotations_db`)

Entry point: `metagomics2-build-annotations` (defined in `pyproject.toml`). The build logic lives in the installed package at `src/metagomics2/scripts/build_annotations_db.py`. A thin wrapper at `scripts/build_annotations_db.py` delegates to this module for convenience when running directly with `python scripts/build_annotations_db.py`.

The entry point is available in the Docker image, so annotations databases can be built without a local Python installation:

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/mriffle/metagomics2:latest \
  metagomics2-build-annotations \
    --fasta /work/uniprot_sprot.fasta.gz \
    --gaf /work/goa_uniprot_all.gaf.gz \
    --output /work/databases/uniprot_sprot.annotations.db
```

**Inputs**:
- UniProt FASTA file (plain or `.gz`): Accession and taxonomy ID extracted from headers (`OX=` field)
- GOA GAF file (plain or `.gz`): GO annotations filtered to accessions present in the FASTA

**Processing**:
1. Parse FASTA headers → `taxonomy(accession, tax_id)` table
2. Parse GAF → `go_annotations(accession, go_id, aspect)` table, filtered to FASTA accessions only
3. Exclude negative annotations (Qualifier contains "NOT") and ND evidence codes
4. Batch inserts (50,000 rows) for performance
5. Create index on `go_annotations.accession`
6. Store metadata (source files, counts, filters applied)

**Schema**:
```sql
CREATE TABLE taxonomy (accession TEXT PRIMARY KEY, tax_id INTEGER NOT NULL);
CREATE TABLE go_annotations (accession TEXT NOT NULL, go_id TEXT NOT NULL, aspect TEXT NOT NULL, UNIQUE(accession, go_id));
CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
```

### Runtime Lookup (`core/subject_lookup.py`)

Opens the SQLite database in immutable mode (`?immutable=1`). Batch queries (respecting SQLite's 999-variable limit) fetch taxonomy and GO data for all subject IDs from DIAMOND results. Full subject IDs (e.g., `sp|Q21HH2|RS2_SACD2`) are mapped to bare accessions (`Q21HH2`) for lookup, but results are keyed by the full ID.

---

## 11. Reference Data System

### Bundled in Docker Image
- **Gene Ontology OBO**: `/app/reference/go/go.obo` (version specified by `GO_VERSION` build arg)
- **NCBI Taxonomy dump**: `/app/reference/taxonomy/` (`nodes.dmp`, `names.dmp`, etc., with version metadata derived from the upstream `Last-Modified` header and a UTC fetch-date fallback during Docker build)

### Format Support
- **GO**: OBO format (parsed by `core/obo_parser.py`) or JSON
- **Taxonomy**: NCBI dump directory (parsed by `core/ncbi_parser.py`) or JSON

### Per-Job Snapshots (Web Mode)
When `job_dir` is set, the pipeline creates a hardlink-based snapshot at `<job_dir>/work/ref_snapshot/`. All snapshot files are SHA256-hashed and recorded in the manifest for provenance.

### CLI Override
The `--go` and `--taxonomy` flags allow using custom reference data files instead of bundled data.

---

## 12. Docker Infrastructure

### Multi-Stage Dockerfile

**Stage 1** (`frontend-builder`): Node.js 20, installs npm dependencies, builds the React frontend.

**Stage 2** (`python:3.11-slim`): Installs system deps (wget), DIAMOND binary, Python package (editable install), copies built frontend, downloads GO and taxonomy reference data, creates data directories.

### Entrypoint (`docker-entrypoint.sh`)

- If arguments are provided:
  - `metagomics2`, `run`, `version`, or flags (`-*`): passes through to the `metagomics2` CLI
  - Any other command (e.g., `metagomics2-build-annotations`): executed directly via `exec "$@"`
- If no arguments: starts the **worker** as a background process, then starts the **web server** (uvicorn) in the foreground. Traps SIGTERM/SIGINT to cleanly shut down both processes.

### Docker Compose (`docker-compose.example.yml`)

Single service with:
- Port mapping: `8000:8000`
- Volumes: data directory, databases directory (read-only), and config directory (read-only)
- Environment variables for scalar configuration (paths, threads, passwords)
- Structured configuration (database definitions, server settings) loaded from JSON files mounted at `/config`

---

## 13. Dependency Management

### Python (`pyproject.toml`)

Build system: **Hatchling**. Python ≥ 3.10 required.

**Runtime dependencies**:
| Package | Purpose |
|---------|---------|
| `fastapi` | Web framework |
| `uvicorn[standard]` | ASGI server |
| `python-multipart` | File upload handling |
| `pydantic` | Data validation and models |
| `aiofiles` | Async file I/O for uploads |
| `pyahocorasick` | Aho-Corasick automaton for peptide matching |
| `httpx` | HTTP client |
| `polars` | DataFrame library (Parquet output) |
| `pyarrow` | Parquet file format support |

**Dev dependencies** (`[dev]` extra):
| Package | Purpose |
|---------|---------|
| `pytest` | Test runner |
| `pytest-cov` | Coverage reporting |
| `pytest-asyncio` | Async test support |
| `hypothesis` | Property-based testing |
| `freezegun` | Time mocking |
| `ruff` | Linter and formatter |
| `mypy` | Static type checking |

### Version Management
- `pyproject.toml` defines the base version (e.g., `0.1.0`)
- Docker `ARG VERSION` overrides at build time
- Runtime: `__version__ = os.getenv("METAGOMICS_VERSION", "0.1.0")`
- Follows semantic versioning. Release process: update pyproject.toml, build Docker with matching version, tag git.

---

## 14. Testing Framework

### Philosophy
- **Correctness is paramount**: Scientific algorithms (LCA, GO closure, aggregation) must produce mathematically correct results. Tests verify both specific scenarios and general invariants.
- **Every module has tests**: Each `core/` module has a corresponding test file. Server, worker, and database have dedicated test suites.
- **Three test tiers**: Unit → Property-based → Integration
- **Test fixtures provide a known scenario**: A small, hand-verified dataset is used throughout to ensure all stages produce expected results.

### Running Tests

```bash
# Always use the project virtual environment
./venv/bin/python -m pytest tests/

# Run with coverage
./venv/bin/python -m pytest tests/ --cov=metagomics2 --cov-report=html

# Run only unit tests
./venv/bin/python -m pytest tests/unit/

# Run only property-based tests
./venv/bin/python -m pytest tests/property/

# Run only integration tests
./venv/bin/python -m pytest tests/integration/

# Skip slow tests (marked with @pytest.mark.slow)
./venv/bin/python -m pytest tests/ -m "not slow"
```

### CI Pipeline (`.github/workflows/ci.yml`)

Runs on every push:
1. **Python tests**: Creates venv, installs `.[dev]`, runs `pytest tests/`
2. **Frontend tests**: Builds `frontend-builder` Docker target, runs `tsc --noEmit` and `vitest run`
3. **Docker smoke test**: Builds full image, runs `--help` command and `metagomics2-build-annotations --help`

### Test Configuration (`pyproject.toml`)

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_functions = ["test_*"]
addopts = "-v --tb=short"
markers = [
    "slow: marks tests as slow",
    "integration: marks tests as integration tests",
]
```

### Unit Tests (`tests/unit/` — 28 files)

Each test file corresponds to a module. Tests use **class-based grouping** with descriptive method names. Helper factory functions (e.g., `make_hit()`, `make_annotation()`) create test objects concisely.

Key test files and what they verify:

| Test File | Module Under Test | What It Verifies |
|-----------|-------------------|------------------|
| `test_fasta_parsing.py` | `core/fasta.py` | FASTA parsing, header extraction, error handling |
| `test_peptides_parsing.py` | `core/peptides.py` | CSV/TSV parsing, normalization, quantity validation, duplicate handling |
| `test_exact_matching.py` | `core/matching.py` | Aho-Corasick matching correctness, edge cases |
| `test_hit_filtering.py` | `core/filtering.py` | Threshold filters, tie-aware top_k, determinism |
| `test_taxonomy_lca.py` | `core/taxonomy.py` | LCA computation, lineage, edge cases |
| `test_go_closure.py` | `core/go.py` | GO closure, union, edge types |
| `test_peptide_annotation_taxonomy.py` | `core/annotation.py` | Taxonomy annotation via LCA |
| `test_peptide_annotation_go.py` | `core/annotation.py` | GO annotation via union of closures |
| `test_aggregation.py` | `core/aggregation.py` | Quantity rollup, ratios, coverage, combo cross-tab |
| `test_reporting.py` | `core/reporting.py` | CSV/Parquet output format and content |
| `test_diamond.py` | `core/diamond.py` | DIAMOND output parsing, accession extraction |
| `test_obo_parser.py` | `core/obo_parser.py` | OBO format parsing |
| `test_ncbi_parser.py` | `core/ncbi_parser.py` | NCBI dump parsing |
| `test_gaf_parser.py` | `core/gaf_parser.py` | GAF 2.2 parsing, filtering |
| `test_uniprot_fasta.py` | `core/uniprot_fasta.py` | UniProt header parsing (accession + OX=taxid) |
| `test_subject_lookup.py` | `core/subject_lookup.py` | SQLite annotation lookup |
| `test_reference_loader.py` | `core/reference_loader.py` | Reference loading, snapshots |
| `test_manifest.py` | `core/reporting.py` | Manifest generation |
| `test_database.py` | `db/database.py` | SQLite operations, job lifecycle |
| `test_server_api.py` | `server/app.py` | FastAPI endpoints (uses `TestClient`) |
| `test_worker.py` | `worker/worker.py` | Worker lifecycle, config building, cleanup |
| `test_input_sanitization.py` | Multiple | Input validation, path traversal prevention |
| `test_email_notification.py` | `notifications/email.py` | Email message building, SMTP mocking |
| `test_email_validation.py` | `models/job.py` | Email field validation |
| `test_build_annotations_db.py` | `metagomics2/scripts/build_annotations_db.py` | Annotations DB build process |
| `test_accession_parsing.py` | `core/diamond.py` | UniProt accession extraction from DIAMOND IDs |
| `test_write_peptide_mapping_parquet.py` | `core/reporting.py` | Parquet output schema and content |

### Property-Based Tests (`tests/property/` — 3 files)

Use **Hypothesis** to generate random inputs and verify mathematical invariants hold for all cases.

| Test File | What It Verifies |
|-----------|------------------|
| `test_aggregation_invariants.py` | Quantity bounds, ratio bounds, ratio ordering, coverage sums, n_peptides bounds |
| `test_go_union_invariants.py` | GO closure union properties |
| `test_taxonomy_lca_invariants.py` | LCA properties (ancestor of all inputs, depth ordering) |

Example: `test_aggregation_invariants.py` generates random `PeptideAnnotation` lists (up to 50 items, 100 examples) and checks that all aggregation invariants hold.

### Integration Tests (`tests/integration/` — 5 files)

| Test File | What It Verifies |
|-----------|------------------|
| `test_pipeline_mocked_homology.py` | Full pipeline with mock DIAMOND output: output files exist, CSV values correct, manifest fields present, progress callbacks, multiple peptide lists |
| `test_pipeline_cli_end_to_end.py` | CLI subprocess execution: `run`, `version`, `--help`, filter params, params file |
| `test_real_go_obo.py` | Parsing real GO OBO files (marked `slow`) |
| `test_real_ncbi_taxonomy.py` | Parsing real NCBI taxonomy dumps (marked `slow`) |
| `test_reference_snapshot.py` | Reference data snapshot creation and loading |

### Test Fixtures (`tests/fixtures/`)

A small, hand-verified dataset used across all test tiers:

**Background FASTA** (`fasta/small_background.fasta`):
- `B1`: contains "PEPTIDE" substring
- `B2`: contains "PEPTIDE" substring
- `B3`: contains "ABC" substring

**Peptide list** (`peptides/small_peptides.tsv`):
- `PEPTIDE` (quantity 10) — matches B1, B2
- `ABC` (quantity 5) — matches B3
- `NOMATCH` (quantity 3) — matches nothing

**Mock hits** (`hits/accepted_hits.json`):
- B1 → {U1, U2}, B2 → {U3}, B3 → {} (no annotated hits)

**Subject annotations** (`annotations/subjects.json`):
- U1: tax_id=70 (SpeciesA), GO={GO:0000004}
- U2: tax_id=71 (SpeciesB), GO={GO:0000005}
- U3: tax_id=72 (SpeciesC), GO={GO:0000006}

**Expected flow for "PEPTIDE"**:
1. Matches B1, B2 → background_proteins = {B1, B2}
2. Implied subjects: B1→{U1,U2}, B2→{U3} → S(PEPTIDE) = {U1, U2, U3}
3. Tax IDs: {70, 71, 72} → LCA = 30 (ClassA) → taxonomy_nodes = {30, 20, 10, 1}
4. GO terms: U1→{C}, U2→{D}, U3→{E} → closures: C→{C,A,B,root}, D→{D,A,root}, E→{E,B,root} → union = {C,D,E,A,B,root} = {GO:0000001..0000006}

**Expected flow for "ABC"**:
1. Matches B3 → background_proteins = {B3}
2. B3 has no accepted subjects → not annotated

**Taxonomy fixture** (`taxonomy/small_taxonomy.json`): A tree with root → kingdoms → phyla → classes → orders → families → genera → species. Two kingdoms with branching taxonomy.

**GO fixture** (`go/small_go.json`): A small DAG with `is_a` and `part_of` edges: root_BP ← A ← C, D; root_BP ← B ← C, E; D `part_of` B.

### Writing New Tests

When adding new functionality:
1. Create unit tests in `tests/unit/test_<module>.py`
2. Use class-based grouping: `class TestFeatureName:`
3. Use helper factory functions to create test objects
4. Test both happy paths and edge cases (empty input, invalid input, boundary values)
5. For mathematical/algorithmic code, add property-based tests in `tests/property/`
6. For end-to-end behavior, add integration tests
7. Use the existing fixtures when applicable; extend them if needed
8. Run `./venv/bin/python -m pytest tests/` and ensure all tests pass before committing

---

## 15. Provenance and Reproducibility

Every pipeline run produces a `run_manifest.json` capturing:
- Metagomics 2 version and git SHA
- Python version
- Search tool and version (e.g., "diamond version 2.1.21")
- Input file SHA256 hashes (FASTA, peptide lists)
- Annotated database hash
- Reference data file hashes (GO OBO, taxonomy dumps)
- All filter parameters
- UTC timestamp

This enables full reproducibility: given the same inputs, same software version, and same reference data, the output should be identical.

---

## 16. Configuration System

Metagomics 2 uses a hybrid configuration model:

- **Environment variables / `.env`** — simple scalar values (paths, thread counts, passwords, feature flags)
- **JSON config files** — structured data (database definitions, server settings)

All configuration is loaded and validated once at startup by the centralized config module (`config.py`). The module exposes a `get_settings()` function that returns an immutable `Settings` dataclass. Both the server and worker call this at module-import time.

### Config Module (`config.py`)

Key components:
- `Settings` — Frozen dataclass holding all validated runtime settings
- `DatabaseEntry` — Frozen dataclass for a single annotated database (`name`, `description`, `path`, `annotations`)
- `SmtpSettings` — Frozen dataclass for SMTP configuration
- `load_settings()` — Reads env vars and JSON files, validates, returns `Settings`
- `get_settings()` / `set_settings()` / `reset_settings()` — Singleton management

### JSON Config Files

Stored in `METAGOMICS_CONFIG_DIR` (default `./config`, `/config` in Docker):

- **`databases.json`** (required) — JSON array of database entries. Each entry must have `name`, `description`, and `path`; `annotations` is optional. The server refuses to start if no databases are configured.
- **`server.json`** (optional) — JSON object with server settings such as `allowed_origins` (list of strings for CORS).

Example files are provided at `config/databases.example.json` and `config/server.example.json`.

### Legacy Fallback

If no `databases.json` file exists, the config loader falls back to the `METAGOMICS_DATABASES` environment variable (a JSON string). This preserves backward compatibility but is deprecated in favor of the JSON config file.

### Environment Variables Reference

| Variable | Default | Used By | Description |
|----------|---------|---------|-------------|
| `METAGOMICS_DATA_DIR` | `/data` | Server, Worker | Base directory for persistent storage |
| `METAGOMICS_DATABASES_DIR` | `/databases` | Server, Worker | Directory containing .dmnd and .annotations.db files |
| `METAGOMICS_CONFIG_DIR` | `./config` | Server, Worker | Directory containing JSON config files (`databases.json`, `server.json`) |
| `METAGOMICS_ADMIN_PASSWORD` | *(empty)* | Server | Admin panel password |
| `METAGOMICS_THREADS` | `4` | Server, Worker | DIAMOND thread count |
| `METAGOMICS_MAX_UPLOAD_MB` | `1024` | Server | Maximum upload file size in MB |
| `METAGOMICS_POLL_INTERVAL` | `5` | Worker | Seconds between job queue polls |
| `METAGOMICS_CLEANUP_ON_SUCCESS` | `true` | Worker | Delete inputs/work after successful job |
| `METAGOMICS_CLEANUP_ON_FAILURE` | `true` | Worker | Delete inputs/work after failed job |
| `METAGOMICS_VERSION` | `0.1.0` | All | Runtime version (set by Docker build) |
| `DIAMOND_VERSION` | *(set at build)* | Server | DIAMOND version string for config API |
| `SMTP_HOST` | *(empty)* | Worker | SMTP server for notifications |
| `SMTP_PORT` | `587` | Worker | SMTP port |
| `SMTP_USERNAME` | *(empty)* | Worker | SMTP username |
| `SMTP_PASSWORD` | *(empty)* | Worker | SMTP password |
| `SMTP_FROM` | *(empty)* | Worker | Sender email address |
| `SITE_URL` | *(empty)* | Worker | Base URL for job links in emails |

---

## 17. Job Directory Structure (Web Mode)

```
/data/jobs/<job_id>/
├── inputs/
│   ├── background.fasta           # Uploaded FASTA
│   └── peptides/
│       ├── list_000_filename.tsv   # Uploaded peptide files
│       └── list_001_filename.tsv
├── work/
│   ├── hit_proteins.fasta          # Subset FASTA for DIAMOND
│   ├── diamond_results.tsv         # DIAMOND output
│   └── ref_snapshot/               # Reference data snapshot
│       ├── go/
│       │   ├── go.obo
│       │   └── VERSION
│       └── taxonomy/
│           ├── nodes.dmp
│           ├── names.dmp
│           └── VERSION
├── results/
│   ├── list_000/
│   │   ├── taxonomy_nodes.csv
│   │   ├── go_terms.csv
│   │   ├── go_taxonomy_combo.csv
│   │   ├── coverage.csv
│   │   ├── peptide_mapping.parquet
│   │   └── run_manifest.json
│   └── list_001/
│       └── (same structure)
└── logs/
```

After successful completion (with cleanup enabled), `inputs/` and `work/` are deleted to save disk space. Only `results/` persists.

---

## 18. Release and CI/CD

### GitHub Actions CI (`.github/workflows/ci.yml`)
Runs on every push:
1. Python tests: venv → install → `pytest tests/`
2. Frontend tests: Docker `frontend-builder` target → `tsc --noEmit` → `vitest run`
3. Docker smoke test: full image build → `--help` command, `metagomics2-build-annotations --help`

### Docker Image Release (`.github/workflows/release-docker-image.yml`)
Triggered on GitHub release publish:
1. Validate semantic version from tag (e.g., `v0.1.0`)
2. Build and push to `ghcr.io/<repo>` with tags: `0.1.0`, `0.1`, `0`, `latest`
3. Sets `VERSION` build arg from release tag

---

## 19. Glossary

| Term | Definition |
|------|------------|
| **Background proteome** | The user's unannotated FASTA file of protein sequences from their sample |
| **Peptide list** | CSV/TSV of peptide sequences with abundances/counts |
| **Annotated database** | A curated protein database (e.g., UniProt SwissProt) with known taxonomy and GO annotations, formatted as a DIAMOND `.dmnd` file |
| **Companion annotations DB** | SQLite file (`.annotations.db`) mapping accessions to tax_ids and GO terms |
| **Hit protein** | A background protein that contains at least one peptide as a substring |
| **Subject** | An annotated database protein that a hit protein matches via homology search |
| **Implied subjects** | All subjects reachable from a peptide through: peptide → background proteins → homology → subjects |
| **LCA** | Lowest Common Ancestor — the deepest taxonomy node ancestral to all of a peptide's implied subjects |
| **GO closure** | The transitive closure of a GO term following is_a/part_of edges to root |
| **Tie-aware top_k** | Ranking filter that keeps all hits tied at the Kth-best bitscore |
| **Quantity** | The peptide abundance/count value from the input file |
| **ratio_total** | `node_quantity / total_peptide_quantity` |
| **ratio_annotated** | `node_quantity / annotated_peptide_quantity` |
