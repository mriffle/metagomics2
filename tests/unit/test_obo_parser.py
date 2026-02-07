"""Unit tests for OBO parser."""

import tempfile
from pathlib import Path

import pytest

from metagomics2.core.go import GODAG
from metagomics2.core.obo_parser import (
    OBOParsingError,
    convert_obo_to_json_dict,
    parse_obo_file,
)


class TestOBOParser:
    """Tests for OBO format parsing."""

    def test_parse_minimal_obo(self, tmp_path: Path):
        """Parse a minimal valid OBO file."""
        obo_content = """format-version: 1.2

[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process

[Term]
id: GO:0000002
name: mitochondrial genome maintenance
namespace: biological_process
is_a: GO:0000001 ! mitochondrion inheritance
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        assert len(dag.terms) == 2
        assert "GO:0000001" in dag.terms
        assert "GO:0000002" in dag.terms

        term1 = dag.terms["GO:0000001"]
        assert term1.name == "mitochondrion inheritance"
        assert term1.namespace == "biological_process"
        assert len(term1.parents) == 0

        term2 = dag.terms["GO:0000002"]
        assert term2.name == "mitochondrial genome maintenance"
        assert "is_a" in term2.parents
        assert "GO:0000001" in term2.parents["is_a"]

    def test_parse_obo_with_multiple_parents(self, tmp_path: Path):
        """Parse term with multiple is_a parents."""
        obo_content = """[Term]
id: GO:0000003
name: reproduction
namespace: biological_process
is_a: GO:0000001 ! parent A
is_a: GO:0000002 ! parent B
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        term = dag.terms["GO:0000003"]
        assert len(term.parents["is_a"]) == 2
        assert "GO:0000001" in term.parents["is_a"]
        assert "GO:0000002" in term.parents["is_a"]

    def test_parse_obo_with_part_of_relationship(self, tmp_path: Path):
        """Parse term with part_of relationship."""
        obo_content = """[Term]
id: GO:0000004
name: biological regulation
namespace: biological_process
is_a: GO:0000001 ! is_a parent
relationship: part_of GO:0000002 ! part_of parent
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        term = dag.terms["GO:0000004"]
        assert "is_a" in term.parents
        assert "GO:0000001" in term.parents["is_a"]
        assert "part_of" in term.parents
        assert "GO:0000002" in term.parents["part_of"]

    def test_parse_obo_skips_obsolete_terms(self, tmp_path: Path):
        """Obsolete terms should not be in main terms but stored in obsolete_terms."""
        obo_content = """[Term]
id: GO:0000001
name: active term
namespace: biological_process

[Term]
id: GO:0000002
name: obsolete term
namespace: biological_process
is_obsolete: true
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        assert len(dag.terms) == 1
        assert "GO:0000001" in dag.terms
        assert "GO:0000002" not in dag.terms

        # Obsolete term should be stored with metadata
        assert "GO:0000002" in dag.obsolete_terms
        obs = dag.obsolete_terms["GO:0000002"]
        assert obs.name == "obsolete term"
        assert obs.namespace == "biological_process"

    def test_parse_obo_handles_comments(self, tmp_path: Path):
        """Comments and trailing annotations should be handled."""
        obo_content = """! This is a comment

[Term]
id: GO:0000001
name: test term ! with comment
namespace: biological_process ! another comment
is_a: GO:0000002 ! parent term name here
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        term = dag.terms["GO:0000001"]
        assert term.name == "test term"
        assert term.namespace == "biological_process"
        assert "GO:0000002" in term.parents["is_a"]

    def test_parse_obo_skips_typedef_stanzas(self, tmp_path: Path):
        """[Typedef] stanzas should be ignored."""
        obo_content = """[Typedef]
id: part_of
name: part of

[Term]
id: GO:0000001
name: test term
namespace: biological_process
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        # Should only have the term, not the typedef
        assert len(dag.terms) == 1
        assert "GO:0000001" in dag.terms

    def test_parse_obo_file_not_found(self):
        """Parsing non-existent file should raise error."""
        with pytest.raises(OBOParsingError, match="File not found"):
            parse_obo_file(Path("/nonexistent/file.obo"))

    def test_convert_obo_to_json_dict(self, tmp_path: Path):
        """Convert OBO to JSON dictionary format."""
        obo_content = """[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: child A
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000003
name: child B
namespace: biological_process
is_a: GO:0000001
relationship: part_of GO:0000002
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        result = convert_obo_to_json_dict(obo_file)

        # Check structure
        assert "terms" in result
        assert "edges" in result

        # Check terms
        assert "GO:0000001" in result["terms"]
        assert result["terms"]["GO:0000001"]["name"] == "root"
        assert result["terms"]["GO:0000001"]["namespace"] == "biological_process"

        # Check edges
        assert "is_a" in result["edges"]
        assert "part_of" in result["edges"]

        # Verify is_a edges
        is_a_edges = result["edges"]["is_a"]
        assert ["GO:0000002", "GO:0000001"] in is_a_edges
        assert ["GO:0000003", "GO:0000001"] in is_a_edges

        # Verify part_of edges
        part_of_edges = result["edges"]["part_of"]
        assert ["GO:0000003", "GO:0000002"] in part_of_edges

    def test_convert_obo_to_json_dict_includes_obsolete_terms(self, tmp_path: Path):
        """Obsolete terms should be included in JSON dict and survive round-trip."""
        obo_content = """[Term]
id: GO:0000001
name: active term
namespace: biological_process

[Term]
id: GO:0000002
name: obsolete regulation
namespace: molecular_function
is_obsolete: true
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        result = convert_obo_to_json_dict(obo_file)

        # Active term in terms, not in obsolete_terms
        assert "GO:0000001" in result["terms"]
        assert "GO:0000001" not in result.get("obsolete_terms", {})

        # Obsolete term in obsolete_terms, not in terms
        assert "GO:0000002" not in result["terms"]
        assert "GO:0000002" in result["obsolete_terms"]
        assert result["obsolete_terms"]["GO:0000002"]["name"] == "obsolete regulation"
        assert result["obsolete_terms"]["GO:0000002"]["namespace"] == "molecular_function"

        # Round-trip through load_go_from_dict
        from metagomics2.core.go import load_go_from_dict
        dag = load_go_from_dict(result)
        assert "GO:0000001" in dag.terms
        assert "GO:0000002" not in dag.terms
        assert "GO:0000002" in dag.obsolete_terms
        assert dag.obsolete_terms["GO:0000002"].name == "obsolete regulation"

    def test_parse_obo_empty_file(self, tmp_path: Path):
        """Empty OBO file should result in empty DAG."""
        obo_file = tmp_path / "empty.obo"
        obo_file.write_text("")

        dag = parse_obo_file(obo_file)

        assert len(dag.terms) == 0

    def test_parse_obo_only_header(self, tmp_path: Path):
        """OBO file with only header should result in empty DAG."""
        obo_content = """format-version: 1.2
data-version: releases/2024-01-17
"""
        obo_file = tmp_path / "header_only.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        assert len(dag.terms) == 0

    def test_parse_obo_multiple_namespaces(self, tmp_path: Path):
        """Parse terms from different namespaces."""
        obo_content = """[Term]
id: GO:0000001
name: BP term
namespace: biological_process

[Term]
id: GO:0000002
name: MF term
namespace: molecular_function

[Term]
id: GO:0000003
name: CC term
namespace: cellular_component
"""
        obo_file = tmp_path / "test.obo"
        obo_file.write_text(obo_content)

        dag = parse_obo_file(obo_file)

        assert len(dag.terms) == 3
        assert dag.terms["GO:0000001"].namespace == "biological_process"
        assert dag.terms["GO:0000002"].namespace == "molecular_function"
        assert dag.terms["GO:0000003"].namespace == "cellular_component"
