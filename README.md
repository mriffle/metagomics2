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
docker pull ghcr.io/mriffle/taxafasta:latest

# Or pull a specific release
docker pull ghcr.io/mriffle/taxafasta:0.1.0
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
docker run --rm ghcr.io/mriffle/taxafasta:latest --help

# Specific Docker image version
docker run --rm ghcr.io/mriffle/taxafasta:0.1.0 --help
```

### Run the CLI with a Local Python Installation

When running outside Docker, provide paths to your reference data unless you have arranged bundled reference files at `/app/reference`.

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
  ghcr.io/mriffle/taxafasta:latest run \
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
  ghcr.io/mriffle/taxafasta:0.1.0 run \
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
  ghcr.io/mriffle/taxafasta:latest run \
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

### Run the Web Interface Locally

```bash
# Start the server
uvicorn metagomics2.server.app:app --host 0.0.0.0 --port 8000

# In another terminal, start the worker
python -m metagomics2.worker.worker
```

Then open http://localhost:8000 in your browser.

## Reference Data

Metagomics 2 bundles Gene Ontology (GO) and NCBI Taxonomy reference data in the Docker image for reproducibility.

### Bundled Versions
- **Gene Ontology**: 2024-01-17 (OBO format)
- **NCBI Taxonomy**: 2024-01-15 (taxonomy dump)

### Custom Versions
Build with specific versions:
```bash
docker build \
  --build-arg GO_VERSION=2024-03-01 \
  --build-arg NCBI_TAXONOMY_DATE=2024-02-28 \
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

Both `.gz` and plain text inputs are supported. The build script automatically filters the GAF file to only include accessions present in the FASTA, and excludes negative annotations (Qualifier contains "NOT") and ND evidence codes.

### Step 4: Configure the Environment

In your `.env` file, set `METAGOMICS_DATABASES_DIR` to the directory containing the database files, and `METAGOMICS_DATABASES` to a JSON array describing each database:

```bash
METAGOMICS_DATABASES_DIR=/path/to/databases

METAGOMICS_DATABASES=[{"name": "UniProt SwissProt", "description": "Reviewed UniProt entries", "path": "uniprot_sprot.dmnd", "annotations": "uniprot_sprot.annotations.db"}]
```

Each entry in the JSON array has the following fields:

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Display name shown in the web UI |
| `description` | Yes | Short description shown as tooltip |
| `path` | Yes | Filename of the `.dmnd` file (relative to `METAGOMICS_DATABASES_DIR`) |
| `annotations` | Yes | Filename of the `.annotations.db` file (relative to `METAGOMICS_DATABASES_DIR`) |

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

### Docker Compose

The repository tracks `docker-compose.example.yml` as the default template.
Copy it to `docker-compose.yml` before running Docker locally so you can keep machine-specific changes out of git:

```bash
cp docker-compose.example.yml docker-compose.yml
docker compose up -d
```

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
├── src/metagomics2/
│   ├── core/           # Core algorithms
│   │   ├── peptides.py     # Peptide parsing
│   │   ├── fasta.py        # FASTA parsing
│   │   ├── matching.py     # Aho-Corasick matching
│   │   ├── filtering.py    # Hit filtering
│   │   ├── taxonomy.py     # Taxonomy + LCA
│   │   ├── go.py           # GO DAG + closure
│   │   ├── annotation.py   # Peptide annotation
│   │   ├── aggregation.py  # Quantity aggregation
│   │   └── reporting.py    # CSV/manifest output
│   ├── pipeline/       # Pipeline orchestration
│   ├── server/         # FastAPI server
│   ├── worker/         # Background worker
│   ├── db/             # SQLite operations
│   ├── models/         # Pydantic models
│   └── cli.py          # CLI entry point
├── tests/
│   ├── unit/           # Unit tests
│   ├── property/       # Property-based tests
│   ├── integration/    # Integration tests
│   └── fixtures/       # Test data
├── frontend/           # React frontend
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
