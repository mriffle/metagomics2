"""Unit tests for peptide GO annotation semantics."""

import pytest

from metagomics2.core.annotation import (
    SubjectAnnotation,
    annotate_peptide,
    annotate_peptide_go,
    load_subject_annotations_from_dict,
)
from metagomics2.core.go import load_go_from_dict
from metagomics2.core.taxonomy import load_taxonomy_from_dict


class TestAnnotatePeptideGO:
    """Tests for GO annotation."""

    def test_go_union_is_non_redundant(self, small_go: dict, subject_annotations: dict):
        """GO(p) must be a set (no duplicates)."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U1: GO={C}, U2: GO={D}, U3: GO={E}
        implied_subjects = {"U1", "U2", "U3"}

        go_terms = annotate_peptide_go(implied_subjects, annotations, go_dag)

        # Should be a set
        assert len(go_terms) == len(set(go_terms))

    def test_go_closure_union_is_a_only(self, small_go: dict, subject_annotations: dict):
        """Test GO closure union with is_a edges only."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U1: GO={C (GO:0000004)}
        # U2: GO={D (GO:0000005)}
        # U3: GO={E (GO:0000006)}
        implied_subjects = {"U1", "U2", "U3"}

        go_terms = annotate_peptide_go(
            implied_subjects, annotations, go_dag, edge_types={"is_a"}
        )

        # Expected closures (is_a only):
        # C (GO:0000004) -> A, B -> root_BP (C has two is_a parents)
        # D (GO:0000005) -> A -> root_BP (is_a only, not part_of)
        # E (GO:0000006) -> B -> root_BP

        # Union should include: C, D, E, A, B, root_BP
        expected = {
            "GO:0000004",  # C
            "GO:0000005",  # D
            "GO:0000006",  # E
            "GO:0000002",  # A
            "GO:0000003",  # B
            "GO:0000001",  # root_BP
        }
        assert go_terms == expected

    def test_go_closure_union_with_part_of(self, small_go: dict, subject_annotations: dict):
        """Test GO closure union with is_a + part_of edges."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U2: GO={D (GO:0000005)} which has part_of -> B
        implied_subjects = {"U2"}

        go_terms_is_a_only = annotate_peptide_go(
            implied_subjects, annotations, go_dag, edge_types={"is_a"}
        )
        go_terms_with_part_of = annotate_peptide_go(
            implied_subjects, annotations, go_dag, edge_types={"is_a", "part_of"}
        )

        # With is_a only: D -> A -> root_BP
        assert go_terms_is_a_only == {"GO:0000005", "GO:0000002", "GO:0000001"}

        # With part_of: D -> A -> root_BP AND D part_of B -> root_BP
        assert "GO:0000003" in go_terms_with_part_of  # B via part_of

    def test_empty_subjects(self, small_go: dict, subject_annotations: dict):
        """Empty subjects returns empty GO set."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        go_terms = annotate_peptide_go(set(), annotations, go_dag)

        assert go_terms == set()

    def test_subjects_without_go_terms(self, small_go: dict, subject_annotations: dict):
        """Subjects without GO terms contribute nothing."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        # U5 has empty go_terms
        go_terms = annotate_peptide_go({"U5"}, annotations, go_dag)

        assert go_terms == set()

    def test_unknown_subject(self, small_go: dict, subject_annotations: dict):
        """Unknown subjects are ignored."""
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        go_terms = annotate_peptide_go({"UNKNOWN"}, annotations, go_dag)

        assert go_terms == set()

    def test_include_self_true(self, small_go: dict):
        """Closure includes the term itself when include_self=True."""
        go_dag = load_go_from_dict(small_go)
        annotations = {
            "U1": SubjectAnnotation(subject_id="U1", go_terms={"GO:0000004"}),
        }

        go_terms = annotate_peptide_go(
            {"U1"}, annotations, go_dag, include_self=True
        )

        assert "GO:0000004" in go_terms

    def test_include_self_false(self, small_go: dict):
        """Closure excludes the term itself when include_self=False."""
        go_dag = load_go_from_dict(small_go)
        annotations = {
            "U1": SubjectAnnotation(subject_id="U1", go_terms={"GO:0000004"}),
        }

        go_terms = annotate_peptide_go(
            {"U1"}, annotations, go_dag, include_self=False
        )

        assert "GO:0000004" not in go_terms
        # But ancestors should still be there
        assert "GO:0000002" in go_terms  # A
        assert "GO:0000001" in go_terms  # root_BP


class TestAnnotatePeptideGOWithFixtures:
    """Tests using the full fixture setup."""

    def test_peptide_go_annotation_fixture_scenario(
        self, small_taxonomy: dict, small_go: dict, subject_annotations: dict
    ):
        """Test the fixture scenario for GO annotation."""
        tree = load_taxonomy_from_dict(small_taxonomy)
        go_dag = load_go_from_dict(small_go)
        annotations = load_subject_annotations_from_dict(subject_annotations)

        peptide_to_proteins = {"PEPTIDE": {"B1", "B2"}}
        protein_to_subjects = {
            "B1": {"U1", "U2"},
            "B2": {"U3"},
        }

        result = annotate_peptide(
            peptide="PEPTIDE",
            quantity=10.0,
            peptide_to_proteins=peptide_to_proteins,
            protein_to_subjects=protein_to_subjects,
            subject_annotations=annotations,
            taxonomy_tree=tree,
            go_dag=go_dag,
            go_edge_types={"is_a"},
        )

        # U1: C, U2: D, U3: E
        # Closures union (is_a only): C, D, E, A, B, root_BP
        expected_go = {
            "GO:0000004",  # C
            "GO:0000005",  # D
            "GO:0000006",  # E
            "GO:0000002",  # A
            "GO:0000003",  # B
            "GO:0000001",  # root_BP
        }
        assert result.go_terms == expected_go

    def test_unannotated_peptide_empty_go(
        self, small_taxonomy: dict, small_go: dict, subject_annotations: dict
    ):
        """Unannotated peptide has empty GO set."""
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

        assert result.go_terms == set()
