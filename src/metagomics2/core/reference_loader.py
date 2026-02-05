"""Reference data loading and caching utilities."""

import json
import shutil
from pathlib import Path

from metagomics2.core.go import GODAG, load_go_from_dict
from metagomics2.core.ncbi_parser import convert_ncbi_dump_to_json_dict
from metagomics2.core.obo_parser import convert_obo_to_json_dict
from metagomics2.core.taxonomy import TaxonomyTree, load_taxonomy_from_dict


class ReferenceDataError(Exception):
    """Raised when reference data operations fail."""

    pass


def load_go_data(source_path: Path) -> GODAG:
    """Load GO data from OBO or JSON file.

    Args:
        source_path: Path to GO file (.obo or .json)

    Returns:
        GODAG object

    Raises:
        ReferenceDataError: If loading fails
    """
    if not source_path.exists():
        raise ReferenceDataError(f"GO file not found: {source_path}")

    if source_path.suffix == ".obo":
        # Convert OBO to dict and load
        data = convert_obo_to_json_dict(source_path)
        return load_go_from_dict(data)
    elif source_path.suffix == ".json":
        # Load JSON directly
        with open(source_path) as f:
            data = json.load(f)
        return load_go_from_dict(data)
    else:
        raise ReferenceDataError(f"Unsupported GO file format: {source_path.suffix}")


def load_taxonomy_data(source_path: Path) -> TaxonomyTree:
    """Load taxonomy data from NCBI dump directory or JSON file.

    Args:
        source_path: Path to taxonomy directory (with nodes.dmp/names.dmp) or JSON file

    Returns:
        TaxonomyTree object

    Raises:
        ReferenceDataError: If loading fails
    """
    if not source_path.exists():
        raise ReferenceDataError(f"Taxonomy source not found: {source_path}")

    if source_path.is_dir():
        # NCBI dump directory
        data = convert_ncbi_dump_to_json_dict(source_path)
        return load_taxonomy_from_dict(data)
    elif source_path.suffix == ".json":
        # JSON file
        with open(source_path) as f:
            data = json.load(f)
        return load_taxonomy_from_dict(data)
    else:
        raise ReferenceDataError(f"Unsupported taxonomy format: {source_path}")


def create_reference_snapshot(
    source_dir: Path,
    snapshot_dir: Path,
    use_hardlinks: bool = True,
) -> dict[str, str]:
    """Create a snapshot of reference data files.

    Args:
        source_dir: Source directory containing reference files
        snapshot_dir: Destination directory for snapshot
        use_hardlinks: If True, use hardlinks instead of copying (saves space)

    Returns:
        Dictionary mapping relative file paths to their purposes

    Raises:
        ReferenceDataError: If snapshot creation fails
    """
    if not source_dir.exists():
        raise ReferenceDataError(f"Source directory not found: {source_dir}")

    snapshot_dir.mkdir(parents=True, exist_ok=True)

    snapshot_files = {}

    # Copy/hardlink all files from source
    for source_file in source_dir.rglob("*"):
        if not source_file.is_file():
            continue

        # Calculate relative path
        rel_path = source_file.relative_to(source_dir)
        dest_file = snapshot_dir / rel_path

        # Create parent directories
        dest_file.parent.mkdir(parents=True, exist_ok=True)

        # Copy or hardlink
        if use_hardlinks:
            try:
                # Try hardlink first (same filesystem)
                dest_file.hardlink_to(source_file)
            except (OSError, NotImplementedError):
                # Fall back to copy if hardlink fails
                shutil.copy2(source_file, dest_file)
        else:
            shutil.copy2(source_file, dest_file)

        snapshot_files[str(rel_path)] = str(source_file)

    return snapshot_files


def get_bundled_reference_dir() -> Path:
    """Get the path to bundled reference data in the Docker image.

    Returns:
        Path to /app/reference directory
    """
    return Path("/app/reference")


def get_reference_metadata(ref_dir: Path) -> dict[str, str]:
    """Read reference data metadata from VERSION files.

    Args:
        ref_dir: Reference directory containing subdirectories with VERSION files

    Returns:
        Dictionary with metadata for each reference type
    """
    metadata = {}

    for subdir in ref_dir.iterdir():
        if not subdir.is_dir():
            continue

        version_file = subdir / "VERSION"
        if version_file.exists():
            with open(version_file) as f:
                content = f.read().strip()
                metadata[subdir.name] = content

    return metadata
