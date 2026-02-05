"""Unit tests for peptide taxonomy annotation semantics."""

import pytest

from metagomics2.core.annotation import (
    SubjectAnnotation,
    annotate_peptide,
    annotate_peptide_taxonomy,
    get_implied_subjects,
    load_subject_annotations_from_dict,
)
from metagomics2.core.go import load_go_from_dict
from metagomics2.core.taxonomy import load_taxonomy_from_dict


class TestGetImpliedSubjects:
    """Tests for getting implied subjects for a peptide."""

    def test_union_across_background_proteins(self):
        peptide_to_proteins = {"PEPTIDE": {"B1", "B2"}}
        protein_to_subjects = {
            "B1": {"U1", "U2"},
            "B2": {"U3"},
        }

        implied = get_implied_subjects("PEPTIDE", peptide_to_proteins, protein_to_subjects)

        assert implied == {"U1", "U2", "U3"}

    def test_no_background_proteins(self):
        peptide_to_proteins = {"PEPTIDE": set()}
        protein_to_subjects = {"B1": {"U1"}}

        implied = get_implied_subjects("PEPTIDE", peptide_to_proteins, protein_to_subjects)

        assert implied == set()

    def test_background_protein_no_subjects(self):
        peptide_to_proteins = {"PEPTIDE": {"B1"}}
        protein_to_subjects = {"B1": set()}

        implied = get_implied_subjects("PEPTIDE", peptide_to_proteins, protein_to_subjects)

        assert implied == set()

    def test_unknown_peptide(self):
        peptide_to_proteins = {}
        protein_to_subjects = {"B1": {"U1"}}

        implied = get_implied_subjects("UNKNOWN", peptide_to_proteins, protein_to_subjects)

        assert implied == set()


class TestAnnotatePeptideTaxonomy:
    """Tests for taxonomy annotation."""

    def test_lca_correctly_computed(self, small_taxonomy: dict, subject_annotations: dict):
        """For PEPTIDE subjects {U1, U2, U3} with tax_ids {70, 71, 72}, LCA should be 30."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U1: tax=70, U2: tax=71, U3: tax=72
        implied_subjects = {"U1", "U2", "U3"}

        lca, nodes = annotate_peptide_taxonomy(implied_subjects, annotations, tree)

        # LCA of 70, 71, 72 is 30 (ClassA)
        assert lca == 30
        # Lineage from 30 to root: 30, 20, 10, 1
        assert nodes == {30, 20, 10, 1}

    def test_lca_sibling_species(self, small_taxonomy: dict, subject_annotations: dict):
        """LCA of U1 (70) and U2 (71) should be 60 (GenusA)."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        implied_subjects = {"U1", "U2"}

        lca, nodes = annotate_peptide_taxonomy(implied_subjects, annotations, tree)

        assert lca == 60
        assert nodes == {60, 50, 40, 30, 20, 10, 1}

    def test_single_subject(self, small_taxonomy: dict, subject_annotations: dict):
        """Single subject returns its own lineage."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        implied_subjects = {"U1"}  # tax=70

        lca, nodes = annotate_peptide_taxonomy(implied_subjects, annotations, tree)

        assert lca == 70
        assert nodes == {70, 60, 50, 40, 30, 20, 10, 1}

    def test_empty_subjects(self, small_taxonomy: dict, subject_annotations: dict):
        """Empty subjects returns None and empty set."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        lca, nodes = annotate_peptide_taxonomy(set(), annotations, tree)

        assert lca is None
        assert nodes == set()

    def test_subjects_without_tax_id(self, small_taxonomy: dict):
        """Subjects without tax_id are ignored."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = {
            "U1": SubjectAnnotation(subject_id="U1", tax_id=None),
        }

        lca, nodes = annotate_peptide_taxonomy({"U1"}, annotations, tree)

        assert lca is None
        assert nodes == set()

    def test_unknown_subject(self, small_taxonomy: dict, subject_annotations: dict):
        """Unknown subjects are ignored."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        lca, nodes = annotate_peptide_taxonomy({"UNKNOWN"}, annotations, tree)

        assert lca is None
        assert nodes == set()


class TestAnnotatePeptideTaxonomyWithFixtures:
    """Tests using the full fixture setup."""

    def test_peptide_annotation_fixture_scenario(
        self, small_taxonomy: dict, small_go: dict, subject_annotations: dict
    ):
        """Test the fixture scenario: PEPTIDE -> B1,B2 -> U1,U2,U3."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        peptide_to_proteins = {
            "PEPTIDE": {"B1", "B2"},
            "ABC": {"B3"},
            "NOMATCH": set(),
        }
        protein_to_subjects = {
            "B1": {"U1", "U2"},
            "B2": {"U3"},
            "B3": set(),
        }

        # PEPTIDE annotation
        result = annotate_peptide(
            peptide="PEPTIDE",
            quantity=10.0,
            peptide_to_proteins=peptide_to_proteins,
            protein_to_subjects=protein_to_subjects,
            subject_annotations=annotations,
            taxonomy_tree=tree,
            go_dag=go_dag,
        )

        assert result.is_annotated is True
        assert result.lca_tax_id == 30
        assert result.taxonomy_nodes == {30, 20, 10, 1}
        assert result.implied_subjects == {"U1", "U2", "U3"}

    def test_unannotated_peptide_no_subjects(
        self, small_taxonomy: dict, small_go: dict, subject_annotations: dict
    ):
        """ABC -> B3 -> {} subjects -> unannotated."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        peptide_to_proteins = {"ABC": {"B3"}}
        protein_to_subjects = {"B3": set()}

        result = annotate_peptide(
            peptide="ABC",
            quantity=5.0,
            peptide_to_proteins=peptide_to_proteins,
            protein_to_subjects=protein_to_subjects,
            subject_annotations=annotations,
            taxonomy_tree=tree,
            go_dag=go_dag,
        )

        assert result.is_annotated is False
        assert result.lca_tax_id is None
        assert result.taxonomy_nodes == set()
        assert result.go_terms == set()

    def test_unannotated_peptide_no_background_hits(
        self, small_taxonomy: dict, small_go: dict, subject_annotations: dict
    ):
        """NOMATCH -> {} background hits -> unannotated."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        peptide_to_proteins = {"NOMATCH": set()}
        protein_to_subjects = {}

        result = annotate_peptide(
            peptide="NOMATCH",
            quantity=3.0,
            peptide_to_proteins=peptide_to_proteins,
            protein_to_subjects=protein_to_subjects,
            subject_annotations=annotations,
            taxonomy_tree=tree,
            go_dag=go_dag,
        )

        assert result.is_annotated is False
        assert result.lca_tax_id is None
        assert result.taxonomy_nodes == set()


class TestLoadSubjectAnnotations:
    """Tests for loading subject annotations."""

    def test_loads_from_dict(self, subject_annotations: dict):
        annotations = load_subject_annotations_from_dict(subject_annotations)

        assert "U1" in annotations
        assert annotations["U1"].tax_id == 70
        assert annotations["U1"].go_terms == {"GO:0000004"}

    def test_handles_empty_go_terms(self, subject_annotations: dict):
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U5 has empty go_terms
        assert annotations["U5"].go_terms == set()
