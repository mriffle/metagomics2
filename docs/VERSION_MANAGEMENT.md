# Version Management

Metagomics 2 uses a comprehensive version management system that ensures version information is consistently available across all interfaces and in data provenance.

## Version Definition

The version is defined in multiple places for different contexts:

### 1. **Package Version** (`pyproject.toml`)
```toml
[project]
version = "0.1.0"
```

### 2. **Docker Build Argument** (`Dockerfile`)
```dockerfile
ARG VERSION=0.1.0
ENV METAGOMICS_VERSION=${VERSION}
```

### 3. **Runtime Detection** (`src/metagomics2/__init__.py`)
```python
__version__ = os.getenv("METAGOMICS_VERSION", "0.1.0")
```

## Building with Custom Version

### Docker Build
```bash
# Build with default version (0.1.0)
docker build -t metagomics2:latest .

# Build with specific version
docker build --build-arg VERSION=1.0.0 -t metagomics2:v1.0.0 .

# Build with version and reference data versions
docker build \
  --build-arg VERSION=1.0.0 \
  --build-arg GO_VERSION=2024-03-01 \
  -t metagomics2:v1.0.0 .
```

### Docker Compose
```yaml
services:
  metagomics2:
    build:
      context: .
      args:
        VERSION: 1.0.0
        GO_VERSION: 2024-03-01
```

## Version Display

### Command Line Interface
```bash
# Show version
metagomics2 --version
# Output: metagomics2 0.1.0

# Version shown at start of run command
metagomics2 run --fasta bg.fasta --peptides pep.tsv --outdir results/
# Output: Metagomics 2 v0.1.0
#         [progress messages...]
```

### Web Interface
- **Header**: Version displayed next to logo (e.g., "Metagomics 2 v0.1.0")
- **Footer**: Version shown in footer text
- **API**: Available via `/api/version` endpoint

### API Endpoints
```bash
# Get version
curl http://localhost:8000/api/version
# Response: {"version": "0.1.0"}

# Health check includes version
curl http://localhost:8000/api/health
# Response: {"status": "healthy", "version": "0.1.0"}
```

## Data Provenance

The version is automatically included in all output manifests:

### `run_manifest.json`
```json
{
  "metagomics2_version": "0.1.0",
  "python_version": "3.11.8",
  "timestamp": "2024-02-04T17:30:00Z",
  "tool_versions": {
    "diamond": "2.1.8"
  },
  ...
}
```

This ensures that:
- Every analysis can be traced to the exact software version
- Results are reproducible with the same version
- Version changes can be tracked in analysis history

## Version Consistency

The system ensures version consistency:

1. **Development**: Uses version from `pyproject.toml` via `__version__`
2. **Docker**: Uses `METAGOMICS_VERSION` environment variable set at build time
3. **Provenance**: Captures actual runtime version in manifests
4. **API**: Exposes version for frontend and external tools

## Best Practices

### Versioning Strategy
Follow [Semantic Versioning](https://semver.org/):
- **MAJOR**: Breaking changes (e.g., 1.0.0 → 2.0.0)
- **MINOR**: New features, backward compatible (e.g., 1.0.0 → 1.1.0)
- **PATCH**: Bug fixes (e.g., 1.0.0 → 1.0.1)

### Release Process
1. Update version in `pyproject.toml`
2. Build Docker image with matching version:
   ```bash
   docker build --build-arg VERSION=1.0.0 -t metagomics2:v1.0.0 .
   ```
3. Tag Docker image:
   ```bash
   docker tag metagomics2:v1.0.0 metagomics2:latest
   ```
4. Tag git repository:
   ```bash
   git tag -a v1.0.0 -m "Release version 1.0.0"
   git push origin v1.0.0
   ```

### Docker Image Tagging
Recommended tagging scheme:
- `metagomics2:latest` - Latest stable release
- `metagomics2:v1.0.0` - Specific version
- `metagomics2:v1.0.0-go2024.03` - Version with reference data versions

Example:
```bash
docker build \
  --build-arg VERSION=1.0.0 \
  --build-arg GO_VERSION=2024-03-01 \
  -t metagomics2:v1.0.0-go2024.03 \
  -t metagomics2:v1.0.0 \
  -t metagomics2:latest \
  .
```

## Verification

### Check Version in Running Container
```bash
# Via environment variable
docker exec <container> printenv METAGOMICS_VERSION

# Via Python
docker exec <container> python -c "from metagomics2 import __version__; print(__version__)"

# Via CLI
docker exec <container> metagomics2 --version

# Via API
curl http://localhost:8000/api/version
```

### Check Version in Image
```bash
# Check label
docker inspect metagomics2:latest | grep version

# Check environment
docker run --rm metagomics2:latest printenv METAGOMICS_VERSION
```
