"""Integration tests for CLI end-to-end execution."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest


class TestCLIEndToEnd:
    """End-to-end CLI tests with mocked homology."""

    def test_cli_run_produces_outputs(self, fixtures_dir: Path, tmp_path: Path):
        """Run CLI and verify output files are created."""
        output_dir = tmp_path / "results"

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "metagomics2.cli",
                "run",
                "--fasta",
                str(fixtures_dir / "fasta" / "small_background.fasta"),
                "--peptides",
                str(fixtures_dir / "peptides" / "small_peptides.tsv"),
                "--outdir",
                str(output_dir),
                "--go",
                str(fixtures_dir / "go" / "small_go.json"),
                "--taxonomy",
                str(fixtures_dir / "taxonomy" / "small_taxonomy.json"),
                "--mock-hits",
                str(fixtures_dir / "hits" / "accepted_hits.json"),
                "--mock-annotations",
                str(fixtures_dir / "annotations" / "subjects.json"),
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Check output files exist
        list_dir = output_dir / "list_000"
        assert list_dir.exists()
        assert (list_dir / "taxonomy_nodes.csv").exists()
        assert (list_dir / "go_terms.csv").exists()
        assert (list_dir / "coverage.csv").exists()
        assert (list_dir / "run_manifest.json").exists()

    def test_cli_version(self):
        """Test version command."""
        result = subprocess.run(
            [sys.executable, "-m", "metagomics2.cli", "version"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "metagomics2" in result.stdout

    def test_cli_help(self):
        """Test help output."""
        result = subprocess.run(
            [sys.executable, "-m", "metagomics2.cli", "--help"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "metagomics2" in result.stdout
        assert "run" in result.stdout

    def test_cli_run_help(self):
        """Test run command help."""
        result = subprocess.run(
            [sys.executable, "-m", "metagomics2.cli", "run", "--help"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "--fasta" in result.stdout
        assert "--peptides" in result.stdout
        assert "--outdir" in result.stdout
        assert "--enrichment-pvalues" in result.stdout

    def test_cli_missing_required_args(self):
        """Test error on missing required arguments."""
        result = subprocess.run(
            [sys.executable, "-m", "metagomics2.cli", "run"],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "error" in result.stderr.lower()

    def test_cli_multiple_peptide_lists(self, fixtures_dir: Path, tmp_path: Path):
        """Test CLI with multiple peptide lists."""
        output_dir = tmp_path / "results"
        peptide_path = str(fixtures_dir / "peptides" / "small_peptides.tsv")

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "metagomics2.cli",
                "run",
                "--fasta",
                str(fixtures_dir / "fasta" / "small_background.fasta"),
                "--peptides",
                peptide_path,
                "--peptides",
                peptide_path,
                "--outdir",
                str(output_dir),
                "--go",
                str(fixtures_dir / "go" / "small_go.json"),
                "--taxonomy",
                str(fixtures_dir / "taxonomy" / "small_taxonomy.json"),
                "--mock-hits",
                str(fixtures_dir / "hits" / "accepted_hits.json"),
                "--mock-annotations",
                str(fixtures_dir / "annotations" / "subjects.json"),
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Should have output for both lists
        assert (output_dir / "list_000").exists()
        assert (output_dir / "list_001").exists()

    def test_cli_with_filter_params(self, fixtures_dir: Path, tmp_path: Path):
        """Test CLI with filter parameters."""
        output_dir = tmp_path / "results"

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "metagomics2.cli",
                "run",
                "--fasta",
                str(fixtures_dir / "fasta" / "small_background.fasta"),
                "--peptides",
                str(fixtures_dir / "peptides" / "small_peptides.tsv"),
                "--outdir",
                str(output_dir),
                "--go",
                str(fixtures_dir / "go" / "small_go.json"),
                "--taxonomy",
                str(fixtures_dir / "taxonomy" / "small_taxonomy.json"),
                "--mock-hits",
                str(fixtures_dir / "hits" / "accepted_hits.json"),
                "--mock-annotations",
                str(fixtures_dir / "annotations" / "subjects.json"),
                "--max-evalue",
                "1e-5",
                "--min-pident",
                "80.0",
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Check manifest contains parameters
        manifest_path = output_dir / "list_000" / "run_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        assert manifest["parameters"]["max_evalue"] == 1e-5
        assert manifest["parameters"]["min_pident"] == 80.0

    def test_cli_with_params_file(self, fixtures_dir: Path, tmp_path: Path):
        """Test CLI with params JSON file."""
        output_dir = tmp_path / "results"

        # Create params file
        params_file = tmp_path / "params.json"
        params_file.write_text(json.dumps({
            "max_evalue": 1e-10,
            "min_pident": 90.0,
            "top_k": 5,
        }))

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "metagomics2.cli",
                "run",
                "--fasta",
                str(fixtures_dir / "fasta" / "small_background.fasta"),
                "--peptides",
                str(fixtures_dir / "peptides" / "small_peptides.tsv"),
                "--outdir",
                str(output_dir),
                "--go",
                str(fixtures_dir / "go" / "small_go.json"),
                "--taxonomy",
                str(fixtures_dir / "taxonomy" / "small_taxonomy.json"),
                "--mock-hits",
                str(fixtures_dir / "hits" / "accepted_hits.json"),
                "--mock-annotations",
                str(fixtures_dir / "annotations" / "subjects.json"),
                "--params",
                str(params_file),
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"CLI failed: {result.stderr}"

        # Check manifest contains parameters from file
        manifest_path = output_dir / "list_000" / "run_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        assert manifest["parameters"]["max_evalue"] == 1e-10
        assert manifest["parameters"]["min_pident"] == 90.0
        assert manifest["parameters"]["top_k"] == 5
