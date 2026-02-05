"""Unit tests for reference data loader."""

import json
from pathlib import Path

import pytest

from metagomics2.core.go import GODAG
from metagomics2.core.reference_loader import (
    ReferenceDataError,
    create_reference_snapshot,
    get_reference_metadata,
    load_go_data,
    load_taxonomy_data,
)
from metagomics2.core.taxonomy import TaxonomyTree


class TestReferenceLoader:
    """Tests for reference data loading utilities."""

    def test_load_go_from_obo(self, tmp_path: Path):
        """Load GO data from OBO file."""
        obo_content = """[Term]
id: GO:0000001
name: test term
namespace: biological_process
"""
        obo_file = tmp_path / "go.obo"
        obo_file.write_text(obo_content)

        dag = load_go_data(obo_file)

        assert isinstance(dag, GODAG)
        assert "GO:0000001" in dag.terms
        assert dag.terms["GO:0000001"].name == "test term"

    def test_load_go_from_json(self, tmp_path: Path):
        """Load GO data from JSON file."""
        json_data = {
            "terms": {
                "GO:0000001": {
                    "name": "test term",
                    "namespace": "biological_process"
                }
            },
            "edges": {
                "is_a": []
            }
        }
        json_file = tmp_path / "go.json"
        json_file.write_text(json.dumps(json_data))

        dag = load_go_data(json_file)

        assert isinstance(dag, GODAG)
        assert "GO:0000001" in dag.terms

    def test_load_go_file_not_found(self):
        """Loading non-existent GO file should raise error."""
        with pytest.raises(ReferenceDataError, match="GO file not found"):
            load_go_data(Path("/nonexistent/go.obo"))

    def test_load_go_unsupported_format(self, tmp_path: Path):
        """Loading unsupported format should raise error."""
        bad_file = tmp_path / "go.txt"
        bad_file.write_text("not a valid format")

        with pytest.raises(ReferenceDataError, match="Unsupported GO file format"):
            load_go_data(bad_file)

    def test_load_taxonomy_from_ncbi_dump(self, tmp_path: Path):
        """Load taxonomy from NCBI dump directory."""
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        tree = load_taxonomy_data(tmp_path)

        assert isinstance(tree, TaxonomyTree)
        assert 1 in tree.nodes
        assert tree.nodes[1].name == "root"

    def test_load_taxonomy_from_json(self, tmp_path: Path):
        """Load taxonomy from JSON file."""
        json_data = {
            "nodes": {
                "1": {
                    "name": "root",
                    "rank": "no rank",
                    "parent_tax_id": None
                }
            }
        }
        json_file = tmp_path / "taxonomy.json"
        json_file.write_text(json.dumps(json_data))

        tree = load_taxonomy_data(json_file)

        assert isinstance(tree, TaxonomyTree)
        assert 1 in tree.nodes

    def test_load_taxonomy_not_found(self):
        """Loading non-existent taxonomy should raise error."""
        with pytest.raises(ReferenceDataError, match="Taxonomy source not found"):
            load_taxonomy_data(Path("/nonexistent/taxonomy"))

    def test_load_taxonomy_unsupported_format(self, tmp_path: Path):
        """Loading unsupported format should raise error."""
        bad_file = tmp_path / "taxonomy.txt"
        bad_file.write_text("not a valid format")

        with pytest.raises(ReferenceDataError, match="Unsupported taxonomy format"):
            load_taxonomy_data(bad_file)

    def test_create_reference_snapshot(self, tmp_path: Path):
        """Create reference snapshot with hardlinks."""
        # Create source directory with files
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        
        go_dir = source_dir / "go"
        go_dir.mkdir()
        (go_dir / "go.obo").write_text("GO content")
        (go_dir / "VERSION").write_text("2024-01-17")
        
        taxonomy_dir = source_dir / "taxonomy"
        taxonomy_dir.mkdir()
        (taxonomy_dir / "nodes.dmp").write_text("nodes content")
        (taxonomy_dir / "names.dmp").write_text("names content")
        (taxonomy_dir / "VERSION").write_text("2024-01-15")

        # Create snapshot
        snapshot_dir = tmp_path / "snapshot"
        snapshot_files = create_reference_snapshot(source_dir, snapshot_dir)

        # Verify snapshot structure
        assert (snapshot_dir / "go" / "go.obo").exists()
        assert (snapshot_dir / "go" / "VERSION").exists()
        assert (snapshot_dir / "taxonomy" / "nodes.dmp").exists()
        assert (snapshot_dir / "taxonomy" / "names.dmp").exists()
        assert (snapshot_dir / "taxonomy" / "VERSION").exists()

        # Verify content
        assert (snapshot_dir / "go" / "go.obo").read_text() == "GO content"
        assert (snapshot_dir / "taxonomy" / "nodes.dmp").read_text() == "nodes content"

        # Verify snapshot_files mapping
        assert "go/go.obo" in snapshot_files
        assert "taxonomy/nodes.dmp" in snapshot_files

    def test_create_reference_snapshot_source_not_found(self, tmp_path: Path):
        """Creating snapshot from non-existent source should raise error."""
        with pytest.raises(ReferenceDataError, match="Source directory not found"):
            create_reference_snapshot(
                Path("/nonexistent/source"),
                tmp_path / "snapshot"
            )

    def test_create_reference_snapshot_with_copy(self, tmp_path: Path):
        """Create snapshot with copying instead of hardlinks."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        (source_dir / "file.txt").write_text("content")

        snapshot_dir = tmp_path / "snapshot"
        create_reference_snapshot(source_dir, snapshot_dir, use_hardlinks=False)

        assert (snapshot_dir / "file.txt").exists()
        assert (snapshot_dir / "file.txt").read_text() == "content"

    def test_get_reference_metadata(self, tmp_path: Path):
        """Read reference metadata from VERSION files."""
        # Create reference directory structure
        go_dir = tmp_path / "go"
        go_dir.mkdir()
        (go_dir / "VERSION").write_text("2024-01-17\nsource=http://example.com/go.obo")

        taxonomy_dir = tmp_path / "taxonomy"
        taxonomy_dir.mkdir()
        (taxonomy_dir / "VERSION").write_text("2024-01-15\nsource=https://example.com/taxdump.tar.gz")

        metadata = get_reference_metadata(tmp_path)

        assert "go" in metadata
        assert "2024-01-17" in metadata["go"]
        assert "source=http://example.com/go.obo" in metadata["go"]

        assert "taxonomy" in metadata
        assert "2024-01-15" in metadata["taxonomy"]
        assert "source=https://example.com/taxdump.tar.gz" in metadata["taxonomy"]

    def test_get_reference_metadata_no_version_files(self, tmp_path: Path):
        """Metadata should be empty if no VERSION files exist."""
        go_dir = tmp_path / "go"
        go_dir.mkdir()

        metadata = get_reference_metadata(tmp_path)

        assert metadata == {}

    def test_create_snapshot_preserves_directory_structure(self, tmp_path: Path):
        """Snapshot should preserve nested directory structure."""
        source_dir = tmp_path / "source"
        (source_dir / "level1" / "level2" / "level3").mkdir(parents=True)
        (source_dir / "level1" / "level2" / "level3" / "deep.txt").write_text("deep content")

        snapshot_dir = tmp_path / "snapshot"
        create_reference_snapshot(source_dir, snapshot_dir)

        assert (snapshot_dir / "level1" / "level2" / "level3" / "deep.txt").exists()
        assert (snapshot_dir / "level1" / "level2" / "level3" / "deep.txt").read_text() == "deep content"

    def test_create_snapshot_skips_directories(self, tmp_path: Path):
        """Snapshot should only include files, not empty directories."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        (source_dir / "empty_dir").mkdir()
        (source_dir / "file.txt").write_text("content")

        snapshot_dir = tmp_path / "snapshot"
        snapshot_files = create_reference_snapshot(source_dir, snapshot_dir)

        # Should only have the file, not the empty directory
        assert "file.txt" in snapshot_files
        assert len(snapshot_files) == 1
