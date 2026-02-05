# Reference Data System

## Overview

Metagomics 2 bundles Gene Ontology (GO) and NCBI Taxonomy reference data directly in the Docker image to ensure reproducibility and version control. Each job creates a snapshot of these reference files for complete provenance tracking.

## Bundled Reference Data

The Docker image includes:

### Gene Ontology (GO)
- **Location**: `/app/reference/go/go.obo`
- **Format**: OBO (Open Biomedical Ontologies)
- **Version**: Specified via `GO_VERSION` build arg (default: 2024-01-17)
- **Source**: http://purl.obolibrary.org/obo/go/releases/

### NCBI Taxonomy
- **Location**: `/app/reference/taxonomy/`
- **Files**: `nodes.dmp`, `names.dmp`, and other taxonomy dump files
- **Format**: NCBI taxonomy dump
- **Version**: Specified via `NCBI_TAXONOMY_DATE` build arg (default: 2024-01-15)
- **Source**: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

## Version Tracking

Each reference directory contains a `VERSION` file with:
- Version/date identifier
- Source URL

Example `/app/reference/go/VERSION`:
```
2024-01-17
source=http://purl.obolibrary.org/obo/go/releases/2024-01-17/go.obo
```

## Per-Job Snapshots

When a job runs in web mode, the pipeline creates a snapshot of reference data:

1. **Snapshot Location**: `/data/jobs/<job_id>/work/ref_snapshot/`
2. **Method**: Hardlinks (space-efficient) or copies
3. **Structure**:
   ```
   ref_snapshot/
   ├── go/
   │   ├── go.obo
   │   └── VERSION
   └── taxonomy/
       ├── nodes.dmp
       ├── names.dmp
       └── VERSION
   ```

4. **Provenance**: All snapshot files are hashed (SHA256) and recorded in `run_manifest.json`

## Building with Custom Versions

To build the Docker image with specific reference data versions:

```bash
docker build \
  --build-arg GO_VERSION=2024-03-01 \
  --build-arg NCBI_TAXONOMY_DATE=2024-02-28 \
  -t metagomics2:custom .
```

## Supported Formats

### GO Data
- **OBO format** (`.obo`) - Parsed directly
- **JSON format** (`.json`) - Pre-converted format

The pipeline automatically detects and parses OBO files using the built-in parser.

### Taxonomy Data
- **NCBI dump** (directory with `.dmp` files) - Parsed directly
- **JSON format** (`.json`) - Pre-converted format

## CLI Usage

### Using Bundled Reference Data
```bash
# Automatically uses bundled reference data
metagomics2 run --fasta bg.fasta --peptides peptides.tsv --outdir results/
```

### Using Custom Reference Data
```bash
# Override with custom files
metagomics2 run \
  --fasta bg.fasta \
  --peptides peptides.tsv \
  --outdir results/ \
  --go /path/to/custom/go.obo \
  --taxonomy /path/to/custom/taxonomy/
```

## Implementation Details

### OBO Parser
- **Module**: `metagomics2.core.obo_parser`
- **Function**: `parse_obo_file(path) -> GODAG`
- Handles:
  - Term definitions
  - `is_a` relationships
  - `part_of` and other relationships
  - Obsolete term filtering

### NCBI Taxonomy Parser
- **Module**: `metagomics2.core.ncbi_parser`
- **Function**: `parse_ncbi_taxonomy_dump(dir) -> TaxonomyTree`
- Parses:
  - `nodes.dmp` for hierarchy and ranks
  - `names.dmp` for scientific names

### Reference Loader
- **Module**: `metagomics2.core.reference_loader`
- **Functions**:
  - `load_go_data(path)` - Auto-detects OBO or JSON
  - `load_taxonomy_data(path)` - Auto-detects dump or JSON
  - `create_reference_snapshot(src, dest)` - Creates job snapshots
  - `get_bundled_reference_dir()` - Returns `/app/reference`
  - `get_reference_metadata(dir)` - Reads VERSION files

## Manifest Provenance

The `run_manifest.json` includes:

```json
{
  "reference_snapshots": {
    "go_files": {
      "go/go.obo": "sha256:abc123...",
      "go/VERSION": "sha256:def456..."
    },
    "taxonomy_files": {
      "taxonomy/nodes.dmp": "sha256:789abc...",
      "taxonomy/names.dmp": "sha256:012def...",
      "taxonomy/VERSION": "sha256:345678..."
    }
  },
  "parameters": {
    "reference_metadata": {
      "go": "2024-01-17\nsource=http://...",
      "taxonomy": "2024-01-15\nsource=https://..."
    }
  }
}
```

## Reproducibility

To reproduce a job's results:

1. Use the same Docker image tag (contains same reference versions)
2. Or rebuild with the same `GO_VERSION` and `NCBI_TAXONOMY_DATE` build args
3. The manifest's reference file hashes can verify exact versions used

## Updating Reference Data

To update bundled reference data:

1. Update `GO_VERSION` and/or `NCBI_TAXONOMY_DATE` in Dockerfile
2. Rebuild the Docker image
3. Tag with version: `metagomics2:v0.1.0-go2024.03-ncbi2024.02`
4. Document reference versions in release notes
