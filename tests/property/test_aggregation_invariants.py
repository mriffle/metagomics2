"""Property-based tests for aggregation invariants."""

import pytest
from hypothesis import given, settings, strategies as st

from metagomics2.core.aggregation import (
    aggregate_peptide_annotations,
    validate_aggregation_invariants,
)
from metagomics2.core.annotation import PeptideAnnotation


# Strategy for generating peptide annotations
@st.composite
def peptide_annotation_strategy(draw):
    """Generate a random PeptideAnnotation."""
    peptide = draw(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=3, max_size=20))
    quantity = draw(st.floats(min_value=0.0, max_value=1000.0, allow_nan=False))
    is_annotated = draw(st.booleans())

    if is_annotated:
        # Generate some taxonomy nodes (small set of possible IDs)
        taxonomy_nodes = draw(
            st.sets(st.integers(min_value=1, max_value=100), min_size=1, max_size=10)
        )
        # Generate some GO terms
        go_terms = draw(
            st.sets(
                st.from_regex(r"GO:[0-9]{7}", fullmatch=True),
                min_size=1,
                max_size=10,
            )
        )
    else:
        taxonomy_nodes = set()
        go_terms = set()

    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=is_annotated,
        taxonomy_nodes=taxonomy_nodes,
        go_terms=go_terms,
    )


class TestAggregationInvariants:
    """Property-based tests for aggregation invariants."""

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_quantity_bounds(self, annotations: list[PeptideAnnotation]):
        """For any node: 0 ≤ quantity(n) ≤ total_peptide_quantity."""
        result = aggregate_peptide_annotations(annotations)

        total_qty = result.coverage.total_peptide_quantity

        for node in result.taxonomy_nodes.values():
            assert node.quantity >= 0, f"Negative quantity: {node.quantity}"
            assert node.quantity <= total_qty, (
                f"Quantity {node.quantity} > total {total_qty}"
            )

        for node in result.go_terms.values():
            assert node.quantity >= 0, f"Negative quantity: {node.quantity}"
            assert node.quantity <= total_qty, (
                f"Quantity {node.quantity} > total {total_qty}"
            )

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_ratio_total_bounds(self, annotations: list[PeptideAnnotation]):
        """ratio_total must be in [0, 1]."""
        result = aggregate_peptide_annotations(annotations)

        for node in result.taxonomy_nodes.values():
            assert 0 <= node.ratio_total <= 1, (
                f"ratio_total {node.ratio_total} not in [0,1]"
            )

        for node in result.go_terms.values():
            assert 0 <= node.ratio_total <= 1, (
                f"ratio_total {node.ratio_total} not in [0,1]"
            )

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_ratio_annotated_bounds(self, annotations: list[PeptideAnnotation]):
        """When annotated > 0: ratio_annotated must be in [0, 1]."""
        result = aggregate_peptide_annotations(annotations)

        if result.coverage.annotated_peptide_quantity > 0:
            for node in result.taxonomy_nodes.values():
                if node.ratio_annotated is not None:
                    assert 0 <= node.ratio_annotated <= 1, (
                        f"ratio_annotated {node.ratio_annotated} not in [0,1]"
                    )

            for node in result.go_terms.values():
                if node.ratio_annotated is not None:
                    assert 0 <= node.ratio_annotated <= 1, (
                        f"ratio_annotated {node.ratio_annotated} not in [0,1]"
                    )

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_ratio_ordering(self, annotations: list[PeptideAnnotation]):
        """When annotated > 0: ratio_total ≤ ratio_annotated."""
        result = aggregate_peptide_annotations(annotations)

        if result.coverage.annotated_peptide_quantity > 0:
            for node in result.taxonomy_nodes.values():
                if node.ratio_annotated is not None:
                    assert node.ratio_total <= node.ratio_annotated + 1e-10, (
                        f"ratio_total {node.ratio_total} > ratio_annotated {node.ratio_annotated}"
                    )

            for node in result.go_terms.values():
                if node.ratio_annotated is not None:
                    assert node.ratio_total <= node.ratio_annotated + 1e-10, (
                        f"ratio_total {node.ratio_total} > ratio_annotated {node.ratio_annotated}"
                    )

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_validate_invariants_passes(self, annotations: list[PeptideAnnotation]):
        """validate_aggregation_invariants should return no violations."""
        result = aggregate_peptide_annotations(annotations)
        violations = validate_aggregation_invariants(result)

        assert violations == [], f"Invariant violations: {violations}"

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_coverage_sums(self, annotations: list[PeptideAnnotation]):
        """Coverage quantities should sum correctly."""
        result = aggregate_peptide_annotations(annotations)

        expected_total = sum(a.quantity for a in annotations)
        expected_annotated = sum(a.quantity for a in annotations if a.is_annotated)
        expected_unannotated = sum(a.quantity for a in annotations if not a.is_annotated)

        assert abs(result.coverage.total_peptide_quantity - expected_total) < 1e-10
        assert abs(result.coverage.annotated_peptide_quantity - expected_annotated) < 1e-10
        assert abs(result.coverage.unannotated_peptide_quantity - expected_unannotated) < 1e-10

    @given(st.lists(peptide_annotation_strategy(), min_size=1, max_size=50))
    @settings(max_examples=100, deadline=None)
    def test_n_peptides_bounds(self, annotations: list[PeptideAnnotation]):
        """n_peptides should be positive and bounded by total peptides."""
        result = aggregate_peptide_annotations(annotations)

        n_annotated = sum(1 for a in annotations if a.is_annotated)

        for node in result.taxonomy_nodes.values():
            assert node.n_peptides >= 1
            assert node.n_peptides <= n_annotated

        for node in result.go_terms.values():
            assert node.n_peptides >= 1
            assert node.n_peptides <= n_annotated
