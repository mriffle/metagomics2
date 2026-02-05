"""Pytest configuration and shared fixtures."""

import json
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return the path to the fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def small_taxonomy(fixtures_dir: Path) -> dict:
    """Load the small taxonomy fixture."""
    with open(fixtures_dir / "taxonomy" / "small_taxonomy.json") as f:
        return json.load(f)


@pytest.fixture
def small_go(fixtures_dir: Path) -> dict:
    """Load the small GO DAG fixture."""
    with open(fixtures_dir / "go" / "small_go.json") as f:
        return json.load(f)


@pytest.fixture
def subject_annotations(fixtures_dir: Path) -> dict:
    """Load the subject annotations fixture."""
    with open(fixtures_dir / "annotations" / "subjects.json") as f:
        return json.load(f)


@pytest.fixture
def small_background_fasta(fixtures_dir: Path) -> Path:
    """Return path to small background FASTA fixture."""
    return fixtures_dir / "fasta" / "small_background.fasta"


@pytest.fixture
def small_peptides_tsv(fixtures_dir: Path) -> Path:
    """Return path to small peptides TSV fixture."""
    return fixtures_dir / "peptides" / "small_peptides.tsv"


@pytest.fixture
def accepted_hits(fixtures_dir: Path) -> dict:
    """Load the accepted hits mapping fixture."""
    with open(fixtures_dir / "hits" / "accepted_hits.json") as f:
        return json.load(f)
