"""Unit tests for NCBI taxonomy parser."""

from pathlib import Path

import pytest

from metagomics2.core.ncbi_parser import (
    NCBIParsingError,
    convert_ncbi_dump_to_json_dict,
    parse_ncbi_taxonomy_dump,
)


class TestNCBIParser:
    """Tests for NCBI taxonomy dump parsing."""

    def test_parse_minimal_taxonomy_dump(self, tmp_path: Path):
        """Parse a minimal NCBI taxonomy dump."""
        # Create nodes.dmp
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tdomain\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
10\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        # Create names.dmp
        names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\t\t|\tscientific name\t|
10\t|\tProteobacteria\t|\t\t|\tscientific name\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        tree = parse_ncbi_taxonomy_dump(tmp_path)

        assert len(tree.nodes) == 3
        assert 1 in tree.nodes
        assert 2 in tree.nodes
        assert 10 in tree.nodes

        # Check root
        root = tree.nodes[1]
        assert root.name == "root"
        assert root.rank == "no rank"
        assert root.parent_tax_id is None  # Root has no parent

        # Check bacteria
        bacteria = tree.nodes[2]
        assert bacteria.name == "Bacteria"
        assert bacteria.rank == "domain"
        assert bacteria.parent_tax_id == 1

        # Check proteobacteria
        proteo = tree.nodes[10]
        assert proteo.name == "Proteobacteria"
        assert proteo.rank == "phylum"
        assert proteo.parent_tax_id == 2

    def test_parse_taxonomy_with_multiple_name_classes(self, tmp_path: Path):
        """Only scientific names should be used."""
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        # Multiple name types, only scientific name should be used
        names_content = """1\t|\tall\t|\t\t|\tsynonym\t|
1\t|\troot\t|\t\t|\tscientific name\t|
1\t|\troot node\t|\t\t|\tcommon name\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        tree = parse_ncbi_taxonomy_dump(tmp_path)

        assert tree.nodes[1].name == "root"

    def test_parse_taxonomy_missing_nodes_file(self, tmp_path: Path):
        """Missing nodes.dmp should raise error."""
        # Only create names.dmp
        names_file = tmp_path / "names.dmp"
        names_file.write_text("1\t|\troot\t|\t\t|\tscientific name\t|")

        with pytest.raises(NCBIParsingError, match="nodes.dmp not found"):
            parse_ncbi_taxonomy_dump(tmp_path)

    def test_parse_taxonomy_missing_names_file(self, tmp_path: Path):
        """Missing names.dmp should raise error."""
        # Only create nodes.dmp
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text("1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|")

        with pytest.raises(NCBIParsingError, match="names.dmp not found"):
            parse_ncbi_taxonomy_dump(tmp_path)

    def test_parse_taxonomy_directory_not_found(self):
        """Non-existent directory should raise error."""
        with pytest.raises(NCBIParsingError, match="Directory not found"):
            parse_ncbi_taxonomy_dump(Path("/nonexistent/directory"))

    def test_convert_ncbi_dump_to_json_dict(self, tmp_path: Path):
        """Convert NCBI dump to JSON dictionary format."""
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tdomain\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\t\t|\tscientific name\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        result = convert_ncbi_dump_to_json_dict(tmp_path)

        # Check structure
        assert "nodes" in result
        assert "1" in result["nodes"]
        assert "2" in result["nodes"]

        # Check root
        assert result["nodes"]["1"]["name"] == "root"
        assert result["nodes"]["1"]["rank"] == "no rank"
        assert result["nodes"]["1"]["parent_tax_id"] is None

        # Check bacteria
        assert result["nodes"]["2"]["name"] == "Bacteria"
        assert result["nodes"]["2"]["rank"] == "domain"
        assert result["nodes"]["2"]["parent_tax_id"] == 1

    def test_parse_taxonomy_with_various_ranks(self, tmp_path: Path):
        """Parse taxonomy with different rank levels."""
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
2\t|\t1\t|\tdomain\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
10\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
20\t|\t10\t|\tclass\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
30\t|\t20\t|\torder\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
40\t|\t30\t|\tfamily\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
50\t|\t40\t|\tgenus\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
60\t|\t50\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
2\t|\tBacteria\t|\t\t|\tscientific name\t|
10\t|\tProteobacteria\t|\t\t|\tscientific name\t|
20\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|
30\t|\tEnterobacterales\t|\t\t|\tscientific name\t|
40\t|\tEnterobacteriaceae\t|\t\t|\tscientific name\t|
50\t|\tEscherichia\t|\t\t|\tscientific name\t|
60\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        tree = parse_ncbi_taxonomy_dump(tmp_path)

        assert len(tree.nodes) == 8
        assert tree.nodes[1].rank == "no rank"
        assert tree.nodes[2].rank == "domain"
        assert tree.nodes[10].rank == "phylum"
        assert tree.nodes[20].rank == "class"
        assert tree.nodes[30].rank == "order"
        assert tree.nodes[40].rank == "family"
        assert tree.nodes[50].rank == "genus"
        assert tree.nodes[60].rank == "species"

    def test_parse_taxonomy_name_without_scientific_name(self, tmp_path: Path):
        """Taxa without scientific names should use fallback."""
        nodes_content = """1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
999\t|\t1\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
"""
        nodes_file = tmp_path / "nodes.dmp"
        nodes_file.write_text(nodes_content)

        # Only provide scientific name for tax_id 1, not 999
        names_content = """1\t|\troot\t|\t\t|\tscientific name\t|
999\t|\tsome synonym\t|\t\t|\tsynonym\t|
"""
        names_file = tmp_path / "names.dmp"
        names_file.write_text(names_content)

        tree = parse_ncbi_taxonomy_dump(tmp_path)

        # Should use fallback name
        assert tree.nodes[999].name == "taxid_999"
