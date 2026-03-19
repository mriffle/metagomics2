"""Unit tests for report generation."""

import csv
import json
from pathlib import Path

import pytest

from metagomics2.core.aggregation import AggregationResult, ComboAggregate, CoverageStats, NodeAggregate
from metagomics2.core.go import load_go_from_dict
from metagomics2.core.reporting import (
    ManifestInfo,
    compute_file_hash,
    write_coverage_csv,
    write_go_taxonomy_combo_csv,
    write_go_terms_csv,
    write_manifest_json,
    write_taxonomy_nodes_csv,
)
from metagomics2.core.taxonomy import load_taxonomy_from_dict


def make_node(node_id: str | int, quantity: float, n_peptides: int) -> NodeAggregate:
    """Helper to create a NodeAggregate."""
    node = NodeAggregate(node_id=node_id)
    node.quantity = quantity
    node.n_peptides = n_peptides
    node.ratio_total = quantity / 100.0  # Assume total=100
    node.ratio_annotated = quantity / 50.0 if quantity <= 50 else 1.0  # Assume annotated=50
    return node


class TestWriteTaxonomyNodesCSV:
    """Tests for taxonomy_nodes.csv generation."""

    def test_columns_present(self, small_taxonomy: dict, tmp_path: Path):
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()
        result.taxonomy_nodes[30] = make_node(30, 10.0, 1)

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert "tax_id" in reader.fieldnames
        assert "name" in reader.fieldnames
        assert "rank" in reader.fieldnames
        assert "parent_tax_id" in reader.fieldnames
        assert "quantity" in reader.fieldnames
        assert "ratio_total" in reader.fieldnames
        assert "ratio_annotated" in reader.fieldnames
        assert "n_peptides" in reader.fieldnames

    def test_node_metadata_correct(self, small_taxonomy: dict, tmp_path: Path):
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()
        result.taxonomy_nodes[30] = make_node(30, 10.0, 1)

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["tax_id"] == "30"
        assert row["name"] == "ClassA"
        assert row["rank"] == "class"
        assert row["parent_tax_id"] == "20"

    def test_stable_sorting(self, small_taxonomy: dict, tmp_path: Path):
        """Output should be sorted by quantity desc, then tax_id."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()
        result.taxonomy_nodes[30] = make_node(30, 10.0, 1)
        result.taxonomy_nodes[20] = make_node(20, 20.0, 2)
        result.taxonomy_nodes[10] = make_node(10, 10.0, 1)  # Same quantity as 30

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # 20 first (highest quantity), then 10, then 30 (same quantity, sorted by id)
        assert rows[0]["tax_id"] == "20"
        assert rows[1]["tax_id"] == "10"
        assert rows[2]["tax_id"] == "30"

    def test_parent_format(self, small_taxonomy: dict, tmp_path: Path):
        """parent_tax_id should be single value."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()
        result.taxonomy_nodes[30] = make_node(30, 10.0, 1)

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        # Should be a single integer string
        assert row["parent_tax_id"] == "20"

    def test_root_has_empty_parent(self, small_taxonomy: dict, tmp_path: Path):
        """Root node should have empty parent_tax_id."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()
        result.taxonomy_nodes[1] = make_node(1, 10.0, 1)

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["parent_tax_id"] == ""


class TestWriteGOTermsCSV:
    """Tests for go_terms.csv generation."""

    def test_columns_present(self, small_go: dict, tmp_path: Path):
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000004"] = make_node("GO:0000004", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert "go_id" in reader.fieldnames
        assert "name" in reader.fieldnames
        assert "namespace" in reader.fieldnames
        assert "parent_go_ids" in reader.fieldnames
        assert "quantity" in reader.fieldnames
        assert "ratio_total" in reader.fieldnames
        assert "ratio_annotated" in reader.fieldnames
        assert "n_peptides" in reader.fieldnames

    def test_multi_parent_format(self, small_go: dict, tmp_path: Path):
        """parent_go_ids should support multiple parents with delimiter."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        # C (GO:0000004) has two parents: A and B
        result.go_terms["GO:0000004"] = make_node("GO:0000004", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, parent_delimiter=";")

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        # Should contain both parents separated by ;
        parents = row["parent_go_ids"].split(";")
        assert "GO:0000002" in parents  # A
        assert "GO:0000003" in parents  # B

    def test_stable_sorting(self, small_go: dict, tmp_path: Path):
        """Output should be sorted by quantity desc, then go_id."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000004"] = make_node("GO:0000004", 10.0, 1)
        result.go_terms["GO:0000002"] = make_node("GO:0000002", 20.0, 2)
        result.go_terms["GO:0000003"] = make_node("GO:0000003", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert rows[0]["go_id"] == "GO:0000002"
        assert rows[1]["go_id"] == "GO:0000003"
        assert rows[2]["go_id"] == "GO:0000004"


    def test_edge_types_is_a_only(self, small_go: dict, tmp_path: Path):
        """edge_types={'is_a'} should only include is_a parents."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        # D (GO:0000005) has is_a parent A (GO:0000002) and part_of parent B (GO:0000003)
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types={"is_a"})

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents = [p for p in row["parent_go_ids"].split(";") if p]
        assert "GO:0000002" in parents  # is_a parent
        assert "GO:0000003" not in parents  # part_of parent excluded

    def test_edge_types_part_of_only(self, small_go: dict, tmp_path: Path):
        """edge_types={'part_of'} should only include part_of parents."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types={"part_of"})

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents = [p for p in row["parent_go_ids"].split(";") if p]
        assert "GO:0000003" in parents  # part_of parent
        assert "GO:0000002" not in parents  # is_a parent excluded

    def test_edge_types_both(self, small_go: dict, tmp_path: Path):
        """edge_types={'is_a', 'part_of'} should include parents from both."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types={"is_a", "part_of"})

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents = [p for p in row["parent_go_ids"].split(";") if p]
        assert "GO:0000002" in parents
        assert "GO:0000003" in parents

    def test_edge_types_none_includes_all(self, small_go: dict, tmp_path: Path):
        """edge_types=None should include parents from all edge types."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types=None)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents = [p for p in row["parent_go_ids"].split(";") if p]
        assert "GO:0000002" in parents
        assert "GO:0000003" in parents

    def test_edge_types_empty_set(self, small_go: dict, tmp_path: Path):
        """edge_types=set() should produce no parents."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types=set())

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["parent_go_ids"] == ""

    def test_edge_types_nonexistent_type(self, small_go: dict, tmp_path: Path):
        """edge_types with a type not in the data should produce no parents for that type."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000005"] = make_node("GO:0000005", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path, edge_types={"regulates"})

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["parent_go_ids"] == ""

    def test_obsolete_go_term_labeled(self, small_go: dict, tmp_path: Path):
        """Obsolete GO terms should show 'OBSOLETE: <name>' in the CSV."""
        # Add an obsolete term to the GO data
        small_go["obsolete_terms"] = {
            "GO:9999999": {
                "name": "positive regulation of something",
                "namespace": "biological_process",
            }
        }
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:9999999"] = make_node("GO:9999999", 5.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["go_id"] == "GO:9999999"
        assert row["name"] == "OBSOLETE: positive regulation of something"
        assert row["namespace"] == "biological_process"
        assert row["parent_go_ids"] == ""

    def test_unknown_go_term_blank_metadata(self, small_go: dict, tmp_path: Path):
        """GO terms not in DAG or obsolete_terms should have blank name/namespace."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:8888888"] = make_node("GO:8888888", 5.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["go_id"] == "GO:8888888"
        assert row["name"] == ""
        assert row["namespace"] == ""


class TestWriteCoverageCSV:
    """Tests for coverage.csv generation."""

    def test_columns_present(self, tmp_path: Path):
        coverage = CoverageStats(
            total_peptide_quantity=18.0,
            annotated_peptide_quantity=10.0,
            unannotated_peptide_quantity=8.0,
            n_peptides_total=3,
            n_peptides_annotated=1,
            n_peptides_unannotated=2,
        )

        output_path = tmp_path / "coverage.csv"
        write_coverage_csv(coverage, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert "total_peptide_quantity" in reader.fieldnames
        assert "annotated_peptide_quantity" in reader.fieldnames
        assert "unannotated_peptide_quantity" in reader.fieldnames
        assert "annotation_coverage_ratio" in reader.fieldnames
        assert "n_peptides_total" in reader.fieldnames
        assert "n_peptides_annotated" in reader.fieldnames
        assert "n_peptides_unannotated" in reader.fieldnames

    def test_values_correct(self, tmp_path: Path):
        coverage = CoverageStats(
            total_peptide_quantity=18.0,
            annotated_peptide_quantity=10.0,
            unannotated_peptide_quantity=8.0,
            n_peptides_total=3,
            n_peptides_annotated=1,
            n_peptides_unannotated=2,
        )

        output_path = tmp_path / "coverage.csv"
        write_coverage_csv(coverage, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert float(row["total_peptide_quantity"]) == pytest.approx(18.0)
        assert float(row["annotated_peptide_quantity"]) == pytest.approx(10.0)
        assert int(row["n_peptides_total"]) == 3


class TestWriteManifestJSON:
    """Tests for run_manifest.json generation."""

    def test_required_keys_present(self, tmp_path: Path):
        manifest = ManifestInfo(
            metagomics2_version="0.1.0",
            python_version="3.11.0",
            search_tool="diamond",
            timestamp_utc="2024-01-01T00:00:00+00:00",
        )

        output_path = tmp_path / "run_manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert "metagomics2_version" in data
        assert "python_version" in data
        assert "search_tool" in data
        assert "timestamp_utc" in data
        assert "inputs" in data
        assert "parameters" in data

    def test_file_hashes_included(self, tmp_path: Path):
        manifest = ManifestInfo(
            metagomics2_version="0.1.0",
            input_fasta_hash="abc123",
            peptide_list_hash="def456",
        )

        output_path = tmp_path / "run_manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["inputs"]["fasta"]["sha256"] == "abc123"
        assert data["inputs"]["peptide_list"]["sha256"] == "def456"

    def test_parameters_included(self, tmp_path: Path):
        manifest = ManifestInfo(
            metagomics2_version="0.1.0",
            parameters={"max_evalue": 1e-5, "min_pident": 80.0},
        )

        output_path = tmp_path / "run_manifest.json"
        write_manifest_json(manifest, output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert data["parameters"]["max_evalue"] == 1e-5
        assert data["parameters"]["min_pident"] == 80.0


class TestComputeFileHash:
    """Tests for file hash computation."""

    def test_hash_deterministic(self, small_background_fasta: Path):
        hash1 = compute_file_hash(small_background_fasta)
        hash2 = compute_file_hash(small_background_fasta)
        assert hash1 == hash2

    def test_hash_is_sha256(self, small_background_fasta: Path):
        file_hash = compute_file_hash(small_background_fasta)
        assert len(file_hash) == 64
        assert all(c in "0123456789abcdef" for c in file_hash)

    def test_known_hash(self, tmp_path: Path):
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello")

        file_hash = compute_file_hash(test_file)
        # Known SHA256 for "hello"
        expected = "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        assert file_hash == expected


class TestRoundTripTreeRebuild:
    """Tests for rebuilding tree structure from CSV."""

    def test_taxonomy_tree_rebuild(self, small_taxonomy: dict, tmp_path: Path):
        """Can rebuild parent-child links from taxonomy_nodes.csv."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        result = AggregationResult()

        # Add several nodes
        for tax_id in [70, 60, 50, 40, 30, 20, 10, 1]:
            result.taxonomy_nodes[tax_id] = make_node(tax_id, 10.0, 1)

        output_path = tmp_path / "taxonomy_nodes.csv"
        write_taxonomy_nodes_csv(result, tree, output_path)

        # Read back and rebuild
        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Build parent map
        parent_map = {}
        for row in rows:
            tax_id = int(row["tax_id"])
            parent = row["parent_tax_id"]
            parent_map[tax_id] = int(parent) if parent else None

        # Verify we can trace from 70 to root
        current = 70
        path = [current]
        while parent_map.get(current) is not None:
            current = parent_map[current]
            path.append(current)

        assert path[-1] == 1  # Ends at root

    def test_go_parents_parseable(self, small_go: dict, tmp_path: Path):
        """Can parse parent_go_ids from go_terms.csv."""
        go_dag = load_go_from_dict(small_go)
        result = AggregationResult()
        result.go_terms["GO:0000004"] = make_node("GO:0000004", 10.0, 1)

        output_path = tmp_path / "go_terms.csv"
        write_go_terms_csv(result, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents_str = row["parent_go_ids"]
        parents = parents_str.split(";") if parents_str else []

        # C has parents A and B
        assert len(parents) == 2
        assert "GO:0000002" in parents
        assert "GO:0000003" in parents


def make_combo(
    tax_id: int,
    go_id: str,
    quantity: float,
    n_peptides: int,
    fraction_of_taxon: float = 0.0,
    fraction_of_go: float = 0.0,
) -> ComboAggregate:
    """Helper to create a ComboAggregate."""
    combo = ComboAggregate(tax_id=tax_id, go_id=go_id)
    combo.quantity = quantity
    combo.n_peptides = n_peptides
    combo.fraction_of_taxon = fraction_of_taxon
    combo.fraction_of_go = fraction_of_go
    return combo


class TestWriteGoTaxonomyComboCSV:
    """Tests for go_taxonomy_combo.csv generation."""

    def test_columns_present(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combos = {
            (30, "GO:0000004"): make_combo(30, "GO:0000004", 10.0, 1, 1.0, 1.0),
        }

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        expected_cols = [
            "tax_id", "tax_name", "tax_rank", "parent_tax_id",
            "go_id", "go_name", "go_namespace", "parent_go_ids",
            "quantity", "fraction_of_taxon", "fraction_of_go", "n_peptides",
            "pvalue_go_for_taxon", "pvalue_taxon_for_go",
            "qvalue_go_for_taxon", "qvalue_taxon_for_go",
        ]
        for col in expected_cols:
            assert col in reader.fieldnames

    def test_metadata_correct(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combos = {
            (30, "GO:0000004"): make_combo(30, "GO:0000004", 10.0, 1, 0.5, 0.8),
        }

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["tax_id"] == "30"
        assert row["tax_name"] == "ClassA"
        assert row["tax_rank"] == "class"
        assert row["parent_tax_id"] == "20"
        assert row["go_id"] == "GO:0000004"
        assert row["go_name"] == "C"
        assert row["go_namespace"] == "biological_process"
        assert float(row["fraction_of_taxon"]) == pytest.approx(0.5)
        assert float(row["fraction_of_go"]) == pytest.approx(0.8)
        assert int(row["n_peptides"]) == 1

    def test_parent_go_ids_present(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        """GO:0000004 (C) has parents GO:0000002 (A) and GO:0000003 (B)."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combos = {
            (30, "GO:0000004"): make_combo(30, "GO:0000004", 10.0, 1),
        }

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        parents = row["parent_go_ids"].split(";")
        assert "GO:0000002" in parents
        assert "GO:0000003" in parents

    def test_stable_sorting(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        """Output sorted by quantity desc, then tax_id, then go_id."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combos = {
            (30, "GO:0000004"): make_combo(30, "GO:0000004", 10.0, 1),
            (20, "GO:0000002"): make_combo(20, "GO:0000002", 20.0, 2),
            (30, "GO:0000002"): make_combo(30, "GO:0000002", 10.0, 1),
        }

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # First: (20, GO:0000002) with qty 20
        assert rows[0]["tax_id"] == "20"
        assert rows[0]["go_id"] == "GO:0000002"
        # Then: (30, GO:0000002) and (30, GO:0000004) both qty 10, same tax_id, sorted by go_id
        assert rows[1]["tax_id"] == "30"
        assert rows[1]["go_id"] == "GO:0000002"
        assert rows[2]["tax_id"] == "30"
        assert rows[2]["go_id"] == "GO:0000004"

    def test_enrichment_values_written(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        """Verify enrichment p-values and q-values are written when populated."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combo = make_combo(30, "GO:0000004", 10.0, 1, 0.5, 0.8)
        combo.pvalue_go_for_taxon = 0.005
        combo.pvalue_taxon_for_go = 0.01
        combo.qvalue_go_for_taxon = 0.05
        combo.qvalue_taxon_for_go = 0.10
        combos = {(30, "GO:0000004"): combo}

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert float(row["pvalue_go_for_taxon"]) == pytest.approx(0.005)
        assert float(row["pvalue_taxon_for_go"]) == pytest.approx(0.01)
        assert float(row["qvalue_go_for_taxon"]) == pytest.approx(0.05)
        assert float(row["qvalue_taxon_for_go"]) == pytest.approx(0.10)

    def test_enrichment_values_empty_when_none(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        """Verify enrichment columns are empty strings when not computed."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        combos = {
            (30, "GO:0000004"): make_combo(30, "GO:0000004", 10.0, 1),
        }

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv(combos, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row["pvalue_go_for_taxon"] == ""
        assert row["pvalue_taxon_for_go"] == ""
        assert row["qvalue_go_for_taxon"] == ""
        assert row["qvalue_taxon_for_go"] == ""

    def test_empty_combos(self, small_taxonomy: dict, small_go: dict, tmp_path: Path):
        """Empty combos produce header-only CSV."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)

        output_path = tmp_path / "go_taxonomy_combo.csv"
        write_go_taxonomy_combo_csv({}, tree, go_dag, output_path)

        with open(output_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 0
        assert "tax_id" in reader.fieldnames
