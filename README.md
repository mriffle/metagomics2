# Metagomics 2

A Dockerized metaproteomics annotation and aggregation tool.

## Overview

Metagomics 2 maps peptides to background proteins, performs homology searches against annotated databases, and transfers taxonomy and Gene Ontology (GO) annotations from reference proteins to peptides.

**Key Features:**
- Exact peptide→protein matching using Aho-Corasick algorithm
- Homology search via DIAMOND or BLAST
- Taxonomy annotation via LCA (Lowest Common Ancestor) intersection
- GO annotation via non-redundant union of closures
- Strong provenance tracking with SHA256 hashes
- Two execution modes: Web server + worker, or CLI-only

## Installation

### Use the Published Docker Image

```bash
docker pull ghcr.io/mriffle/metagomics2:latest

# Or pull a specific release
docker pull ghcr.io/mriffle/metagomics2:0.1.0
```

### Local Python Installation

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install with dev dependencies
pip install -e ".[dev]"
```

> **Note:** The Python package does not include the [DIAMOND](https://github.com/bbuchfink/diamond) sequence aligner. To run the full annotation pipeline (homology search) outside Docker, you must install DIAMOND separately and ensure the `diamond` binary is on your `PATH`. For example:
>
> ```bash
> # Via Conda/Bioconda
> conda install -c bioconda diamond
>
> # Or download a release binary from GitHub
> # https://github.com/bbuchfink/diamond/releases
> ```
>
> Without DIAMOND installed, the pipeline will fail at the homology search stage. The Docker image includes DIAMOND automatically.

### Build the Docker Image Locally

```bash
docker build -t metagomics2 .
```

### Run the Web Interface with Docker Compose

```bash
cp docker-compose.example.yml docker-compose.yml
docker compose up -d
```

## Usage

### Show CLI Help

```bash
# Local Python installation
metagomics2 --help

# Published Docker image
docker run --rm ghcr.io/mriffle/metagomics2:latest --help

# Specific Docker image version
docker run --rm ghcr.io/mriffle/metagomics2:0.1.0 --help
```

### Run the CLI with a Local Python Installation

When running outside Docker, you must have [DIAMOND](https://github.com/bbuchfink/diamond) installed and on your `PATH` (see [Local Python Installation](#local-python-installation)). You must also provide paths to your reference data unless you have arranged bundled reference files at `/app/reference`.

```bash
metagomics2 run \
  --fasta background.fasta \
  --peptides peptides.tsv \
  --outdir results/ \
  --db /path/to/uniprot_sprot.dmnd \
  --annotations-db /path/to/uniprot_sprot.annotations.db \
  --go /path/to/go.obo \
  --taxonomy /path/to/taxonomy/ \
  --max-evalue 1e-5 \
  --min-pident 80
```

### Run the CLI with the Docker Image

The Docker image bundles GO and NCBI taxonomy reference data at `/app/reference`, so you can omit `--go` and `--taxonomy` unless you want to override them.

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/mriffle/metagomics2:latest run \
  --fasta /work/background.fasta \
  --peptides /work/peptides.tsv \
  --outdir /work/results \
  --db /work/databases/uniprot_sprot.dmnd \
  --annotations-db /work/databases/uniprot_sprot.annotations.db \
  --max-evalue 1e-5 \
  --min-pident 80
```

To run a specific released version instead of `latest`, replace the image tag:

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/mriffle/metagomics2:0.1.0 run \
  --fasta /work/background.fasta \
  --peptides /work/peptides.tsv \
  --outdir /work/results \
  --db /work/databases/uniprot_sprot.dmnd \
  --annotations-db /work/databases/uniprot_sprot.annotations.db
```

To override the bundled reference data, mount the files and pass `--go` and `--taxonomy` explicitly:

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/mriffle/metagomics2:latest run \
  --fasta /work/background.fasta \
  --peptides /work/peptides.tsv \
  --outdir /work/results \
  --db /work/databases/uniprot_sprot.dmnd \
  --annotations-db /work/databases/uniprot_sprot.annotations.db \
  --go /work/reference/go.obo \
  --taxonomy /work/reference/taxonomy/ \
  --max-evalue 1e-5 \
  --min-pident 80
```

Use `--peptides` multiple times to process more than one peptide list in a single run.

Both `--db` and `--annotations-db` are required. See [Annotated Databases](#annotated-databases) for how to build these files.

### Run the Web Interface with Docker Compose

The simplest way to run the web interface locally is with Docker Compose. Make sure you have configured your `.env` and `config/databases.json` first (see [Configuration](#configuration)).

```bash
# Copy the example files if you haven't already
cp docker-compose.example.yml docker-compose.yml
cp .env.example .env
cp config/databases.example.json config/databases.json

# Edit .env and config/databases.json for your setup, then start
docker compose up -d
```

The web interface will be available at http://localhost:8000. Docker Compose automatically starts both the web server and the background worker. To view logs:

```bash
docker compose logs -f
```

To stop the service:

```bash
docker compose down
```

### Run the Web Interface without Docker

If you have a local Python installation and want to run without Docker, start the server and worker separately. You must set `METAGOMICS_CONFIG_DIR` (or place a `config/databases.json` in the working directory) and `METAGOMICS_DATA_DIR` so the config loader can find your settings.

```bash
# Start the server
METAGOMICS_CONFIG_DIR=./config METAGOMICS_DATA_DIR=./_data \
  uvicorn metagomics2.server.app:app --host 0.0.0.0 --port 8000

# In another terminal, start the worker
METAGOMICS_CONFIG_DIR=./config METAGOMICS_DATA_DIR=./_data \
  python -m metagomics2.worker.worker
```

Then open http://localhost:8000 in your browser.

## Reference Data

Metagomics 2 bundles Gene Ontology (GO) and NCBI Taxonomy reference data in the Docker image for reproducibility.

### Bundled Versions
- **Gene Ontology**: 2024-01-17 (OBO format)
- **NCBI Taxonomy**: Derived from the upstream `Last-Modified` header for `taxdump.tar.gz`, with a fallback to the UTC fetch date during `docker build`

### Custom Versions
Build with a specific GO version:
```bash
docker build \
  --build-arg GO_VERSION=2024-03-01 \
  -t metagomics2:custom .
```

### Per-Job Snapshots
Each job creates a snapshot of reference data in `/data/jobs/<job_id>/work/ref_snapshot/` with:
- Complete reference files (hardlinked for efficiency)
- SHA256 hashes in manifest
- Version metadata

See [docs/REFERENCE_DATA.md](docs/REFERENCE_DATA.md) for details.

## Annotated Databases

Metagomics 2 requires at least one annotated DIAMOND database for homology search. Each database consists of two files:

| File | Description |
|------|-------------|
| `*.dmnd` | DIAMOND-formatted protein database |
| `*.annotations.db` | Companion SQLite database with taxonomy IDs and GO terms |

Both files must be placed in the directory specified by `METAGOMICS_DATABASES_DIR` in your `.env`.

### Step 1: Download Source Files

Download the UniProt FASTA and GOA annotation file:

```bash
# UniProt SwissProt FASTA (reviewed entries only, ~270 MB compressed)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# GOA annotation file (all UniProt entries, ~10 GB compressed)
wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
```

### Step 2: Build the DIAMOND Database

```bash
diamond makedb --in uniprot_sprot.fasta.gz --db /path/to/databases/uniprot_sprot.dmnd
```

### Step 3: Build the Annotations Database

The annotations database maps UniProt accessions to NCBI taxonomy IDs (from FASTA `OX=` fields) and GO terms (from the GOA GAF file):

```bash
python scripts/build_annotations_db.py \
    --fasta uniprot_sprot.fasta.gz \
    --gaf goa_uniprot_all.gaf.gz \
    --output /path/to/databases/uniprot_sprot.annotations.db
```

Or using the published Docker image:

```bash
docker run --rm \
  -v "$PWD:/work" \
  ghcr.io/mriffle/metagomics2:latest \
  metagomics2-build-annotations \
    --fasta /work/uniprot_sprot.fasta.gz \
    --gaf /work/goa_uniprot_all.gaf.gz \
    --output /work/databases/uniprot_sprot.annotations.db
```

Both `.gz` and plain text inputs are supported. The build script automatically filters the GAF file to only include accessions present in the FASTA, and excludes negative annotations (Qualifier contains "NOT") and ND evidence codes.

### Step 4: Configure the Databases

Create a `databases.json` file inside your config directory (default `./config`). The repository includes an example at `config/databases.example.json`:

```bash
cp config/databases.example.json config/databases.json
```

Edit `config/databases.json` to list your databases:

```json
[
  {
    "name": "UniProt SwissProt 2024_01",
    "description": "Reviewed UniProt entries (SwissProt)",
    "path": "uniprot_sprot.dmnd",
    "annotations": "uniprot_sprot.annotations.db"
  }
]
```

Each entry in the JSON array has the following fields:

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Display name shown in the web UI |
| `description` | Yes | Short description shown as tooltip |
| `path` | Yes | Filename of the `.dmnd` file (relative to `METAGOMICS_DATABASES_DIR`) |
| `annotations` | Yes | Filename of the `.annotations.db` file (relative to `METAGOMICS_DATABASES_DIR`) |

Also set `METAGOMICS_DATABASES_DIR` in your `.env` to point to the directory containing the database files:

```bash
METAGOMICS_DATABASES_DIR=/path/to/databases
```

The service will refuse to start if no databases are configured.

## Input Files

### Background FASTA
Protein sequences (unannotated) that peptides will be matched against.

### Peptide Lists (CSV/TSV)
Required columns:
- `peptide_sequence`: The peptide sequence
- `quantity`: Numeric quantity/count

## Output Files

For each peptide list, the following files are generated:

| File | Description |
|------|-------------|
| `taxonomy_nodes.csv` | Taxonomy nodes with quantities and ratios |
| `go_terms.csv` | GO terms with quantities and ratios |
| `go_taxonomy_combo.csv` | Cross-tabulation of GO terms against taxonomy nodes, used to filter the taxonomy or GO DAG visualizations by the other dimension |
| `coverage.csv` | Coverage statistics |
| `peptide_mapping.parquet` | Per-peptide annotation detail: one row per (peptide, background protein, annotated protein) triple. Contains `peptide_lca_tax_ids` (the LCA taxonomy node and all its ancestors to root, as a list of NCBI tax IDs), `peptide_go_terms` (union of GO annotation closures across all implied subjects), `evalue` and `pident` (e-value and percent identity of the homology search hit between the background protein and the annotated protein). Used by the web UI to display peptide-level details when clicking a node in the taxonomy or GO DAG visualizations |
| `run_manifest.json` | Provenance information |

## Configuration

Metagomics 2 uses a hybrid configuration system:

- **Environment variables / `.env`** — simple scalar values such as paths, thread counts, passwords, and feature flags.
- **JSON config files** — structured data such as database definitions and server settings.

All configuration is loaded and validated once at startup by a central config module (`src/metagomics2/config.py`). The rest of the application works with a clean, validated settings object.

### Quick Start

```bash
# 1. Copy example files
cp .env.example .env
cp config/databases.example.json config/databases.json
cp config/server.example.json config/server.json   # optional

# 2. Edit .env and config/databases.json for your setup

# 3. Start with Docker Compose
cp docker-compose.example.yml docker-compose.yml
docker compose up -d
```

### Environment Variables (`.env`)

The `.env` file holds simple scalar settings. Key variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `METAGOMICS_DATA_DIR` | `/data` | Persistent data directory (jobs, SQLite DB) |
| `METAGOMICS_DATABASES_DIR` | `/databases` | Directory containing `.dmnd` and `.annotations.db` files |
| `METAGOMICS_CONFIG_DIR` | `./config` | Directory containing JSON config files |
| `METAGOMICS_THREADS` | `4` | CPU threads for DIAMOND |
| `METAGOMICS_ADMIN_PASSWORD` | *(empty)* | Password for the admin dashboard |
| `METAGOMICS_MAX_UPLOAD_MB` | `1024` | Max upload size in MB |
| `METAGOMICS_CLEANUP_ON_SUCCESS` | `true` | Delete intermediate files after successful jobs |
| `METAGOMICS_CLEANUP_ON_FAILURE` | `true` | Delete intermediate files after failed jobs |
| `SMTP_HOST` | *(empty)* | SMTP server for email notifications (leave empty to disable) |
| `SMTP_PORT` | `587` | SMTP port |
| `SMTP_USERNAME` | *(empty)* | SMTP username |
| `SMTP_PASSWORD` | *(empty)* | SMTP password |
| `SMTP_FROM` | *(empty)* | Sender address for notifications |
| `SITE_URL` | *(empty)* | Public URL for job links in emails |

See `.env.example` for a fully commented template.

### JSON Config Files

Structured settings live in JSON files inside `METAGOMICS_CONFIG_DIR` (default `./config`).

#### `databases.json` (required)

A JSON array listing the annotated databases available for homology search. See [Annotated Databases](#annotated-databases) for how to build these files and [Step 4](#step-4-configure-the-databases) for the format.

Example: `config/databases.example.json`

#### `server.json` (optional)

Server-specific structured settings such as allowed CORS origins.

Example: `config/server.example.json`

```json
{
  "allowed_origins": ["https://metagomics.example.com"]
}
```

If omitted, CORS allows all origins (`["*"]`).

### Docker Compose

The repository tracks `docker-compose.example.yml` as the default template.
Copy it to `docker-compose.yml` before running Docker locally:

```bash
cp docker-compose.example.yml docker-compose.yml
docker compose up -d
```

The compose file mounts three host directories into the container:

| Host path | Container path | Description |
|-----------|---------------|-------------|
| `METAGOMICS_DATA_DIR` | `/data` | Persistent data |
| `METAGOMICS_DATABASES_DIR` | `/databases` | DIAMOND databases (read-only) |
| `METAGOMICS_CONFIG_DIR` | `/config` | JSON config files (read-only) |

### Filter Parameters

These parameters control how DIAMOND homology hits are filtered before
annotation.  Threshold filters are applied first (as an AND), then the
tie-aware `top_k` ranking selects which hits to keep.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max-evalue` | `1e-10` | Maximum e-value threshold. Also passed to DIAMOND as a pre-filter so low-quality alignments are skipped early. Lower values are more stringent. |
| `--min-pident` | `80` | Minimum percent identity. Applied in post-filtering (not passed to DIAMOND). |
| `--min-qcov` | *(none)* | Minimum query coverage (percent). Applied in post-filtering. |
| `--min-alnlen` | *(none)* | Minimum alignment length (residues). Applied in post-filtering. |
| `--top-k` | `1` | Number of top-scoring hits to keep per query protein, ranked by bitscore. **Tie-aware**: if multiple hits share the same bitscore at the Kth position, all tied hits are retained. For example, with `top_k=1` and five hits tied at the best bitscore, all five are kept. This ensures annotation is not biased by arbitrary tie-breaking. |

### GO Closure Settings

| Parameter | Description |
|-----------|-------------|
| `--go-edge-types` | Edge types for closure (default: `is_a`) |
| `--go-exclude-self` | Exclude terms themselves from closure |

## Project Structure

```
metagomics2/
├── config/                       # Example config files
│   ├── databases.example.json    #   Database definitions
│   └── server.example.json       #   Server settings
├── src/metagomics2/
│   ├── config.py          # Centralized config loader
│   ├── core/              # Core algorithms
│   │   ├── peptides.py        # Peptide parsing
│   │   ├── fasta.py           # FASTA parsing
│   │   ├── matching.py        # Aho-Corasick matching
│   │   ├── filtering.py       # Hit filtering
│   │   ├── taxonomy.py        # Taxonomy + LCA
│   │   ├── go.py              # GO DAG + closure
│   │   ├── annotation.py      # Peptide annotation
│   │   ├── aggregation.py     # Quantity aggregation
│   │   └── reporting.py       # CSV/manifest output
│   ├── pipeline/          # Pipeline orchestration
│   ├── server/            # FastAPI server
│   ├── worker/            # Background worker
│   ├── db/                # SQLite operations
│   ├── models/            # Pydantic models
│   └── cli.py             # CLI entry point
├── tests/
│   ├── unit/              # Unit tests
│   ├── property/          # Property-based tests
│   ├── integration/       # Integration tests
│   └── fixtures/          # Test data
├── frontend/              # React frontend
├── .env.example           # Environment variable template
├── Dockerfile
└── docker-compose.example.yml
```

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=metagomics2 --cov-report=html

# Run only unit tests
pytest tests/unit/

# Run property-based tests
pytest tests/property/

# Skip slow tests
pytest -m "not slow"
```

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/jobs` | Create a new job |
| GET | `/api/jobs/{job_id}` | Get job status |
| GET | `/api/jobs` | List recent jobs |
| GET | `/api/jobs/{job_id}/results/{list_id}/{file}` | Download result file |

## License

Apache 2.0
