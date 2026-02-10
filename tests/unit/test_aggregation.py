"""Unit tests for aggregation of peptide quantities."""

import pytest

from metagomics2.core.aggregation import (
    AggregationResult,
    ComboAggregate,
    CoverageStats,
    NodeAggregate,
    aggregate_go_taxonomy_combos,
    aggregate_peptide_annotations,
    validate_aggregation_invariants,
)
from metagomics2.core.annotation import PeptideAnnotation


def make_annotation(
    peptide: str,
    quantity: float,
    is_annotated: bool = True,
    taxonomy_nodes: set[int] | None = None,
    go_terms: set[str] | None = None,
) -> PeptideAnnotation:
    """Helper to create a PeptideAnnotation."""
    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=is_annotated,
        taxonomy_nodes=taxonomy_nodes or set(),
        go_terms=go_terms or set(),
    )


class TestCoverageStats:
    """Tests for coverage statistics."""

    def test_coverage_with_fixture_scenario(self):
        """Test coverage with PEPTIDE=10 (annotated), ABC=5, NOMATCH=3 (unannotated)."""
        annotations = [
            make_annotation("PEPTIDE", 10.0, is_annotated=True),
            make_annotation("ABC", 5.0, is_annotated=False),
            make_annotation("NOMATCH", 3.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.coverage.total_peptide_quantity == 18.0
        assert result.coverage.annotated_peptide_quantity == 10.0
        assert result.coverage.unannotated_peptide_quantity == 8.0
        assert result.coverage.n_peptides_total == 3
        assert result.coverage.n_peptides_annotated == 1
        assert result.coverage.n_peptides_unannotated == 2

    def test_coverage_ratio(self):
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True),
            make_annotation("P2", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.coverage.annotation_coverage_ratio == 0.5

    def test_coverage_ratio_all_annotated(self):
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.coverage.annotation_coverage_ratio == 1.0

    def test_coverage_ratio_none_annotated(self):
        annotations = [
            make_annotation("P1", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.coverage.annotation_coverage_ratio == 0.0

    def test_coverage_ratio_empty(self):
        result = aggregate_peptide_annotations([])

        assert result.coverage.annotation_coverage_ratio == 0.0


class TestTaxonomyAggregation:
    """Tests for taxonomy node aggregation."""

    def test_node_quantities_correct(self):
        """Nodes in TAX(p) get quantity added."""
        annotations = [
            make_annotation(
                "PEPTIDE", 10.0, is_annotated=True, taxonomy_nodes={30, 20, 10, 1}
            ),
        ]

        result = aggregate_peptide_annotations(annotations)

        # Each node should have quantity 10
        assert result.taxonomy_nodes[30].quantity == 10.0
        assert result.taxonomy_nodes[20].quantity == 10.0
        assert result.taxonomy_nodes[10].quantity == 10.0
        assert result.taxonomy_nodes[1].quantity == 10.0

    def test_unannotated_peptides_not_counted(self):
        """Unannotated peptides don't contribute to node quantities."""
        annotations = [
            make_annotation("PEPTIDE", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("ABC", 5.0, is_annotated=False, taxonomy_nodes=set()),
        ]

        result = aggregate_peptide_annotations(annotations)

        # Only PEPTIDE contributes
        assert result.taxonomy_nodes[30].quantity == 10.0
        assert 30 in result.taxonomy_nodes
        # No other nodes from ABC

    def test_multiple_peptides_same_node(self):
        """Multiple peptides contributing to same node."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30, 1}),
            make_annotation("P2", 5.0, is_annotated=True, taxonomy_nodes={30, 1}),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.taxonomy_nodes[30].quantity == 15.0
        assert result.taxonomy_nodes[1].quantity == 15.0
        assert result.taxonomy_nodes[30].n_peptides == 2

    def test_n_peptides_counts_distinct(self):
        """n_peptides counts distinct peptides, not occurrences."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("P2", 5.0, is_annotated=True, taxonomy_nodes={30}),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.taxonomy_nodes[30].n_peptides == 2


class TestGOAggregation:
    """Tests for GO term aggregation."""

    def test_go_quantities_correct(self):
        """GO terms in GO(p) get quantity added."""
        annotations = [
            make_annotation(
                "PEPTIDE",
                10.0,
                is_annotated=True,
                go_terms={"GO:0000001", "GO:0000002"},
            ),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.go_terms["GO:0000001"].quantity == 10.0
        assert result.go_terms["GO:0000002"].quantity == 10.0

    def test_multiple_peptides_same_go_term(self):
        """Multiple peptides contributing to same GO term."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, go_terms={"GO:0000001"}),
            make_annotation("P2", 5.0, is_annotated=True, go_terms={"GO:0000001"}),
        ]

        result = aggregate_peptide_annotations(annotations)

        assert result.go_terms["GO:0000001"].quantity == 15.0
        assert result.go_terms["GO:0000001"].n_peptides == 2


class TestRatios:
    """Tests for ratio calculations."""

    def test_ratio_total(self):
        """ratio_total = quantity / total_peptide_quantity."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("P2", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        # total = 20, node 30 quantity = 10
        assert result.taxonomy_nodes[30].ratio_total == 0.5

    def test_ratio_annotated(self):
        """ratio_annotated = quantity / annotated_peptide_quantity."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("P2", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        # annotated = 10, node 30 quantity = 10
        assert result.taxonomy_nodes[30].ratio_annotated == 1.0

    def test_ratio_annotated_none_when_no_annotated(self):
        """ratio_annotated is None when annotated_peptide_quantity = 0."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        # No taxonomy nodes since no annotated peptides
        assert len(result.taxonomy_nodes) == 0


class TestInvariants:
    """Tests for aggregation invariants."""

    def test_quantity_bounds(self):
        """0 ≤ quantity(n) ≤ total_peptide_quantity."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("P2", 8.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        violations = validate_aggregation_invariants(result)
        assert violations == []

        # Check bounds
        assert result.taxonomy_nodes[30].quantity >= 0
        assert result.taxonomy_nodes[30].quantity <= 18.0

    def test_ratio_bounds(self):
        """0 ≤ ratio_total ≤ ratio_annotated ≤ 1."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=True, taxonomy_nodes={30}),
            make_annotation("P2", 10.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        violations = validate_aggregation_invariants(result)
        assert violations == []

        node = result.taxonomy_nodes[30]
        assert 0 <= node.ratio_total <= 1
        assert node.ratio_annotated is not None
        assert 0 <= node.ratio_annotated <= 1
        assert node.ratio_total <= node.ratio_annotated

    def test_validate_catches_violations(self):
        """Validation catches invariant violations."""
        result = AggregationResult()
        result.coverage.total_peptide_quantity = 10.0
        result.coverage.annotated_peptide_quantity = 10.0

        # Create a node with invalid values
        bad_node = NodeAggregate(node_id=1)
        bad_node.quantity = 20.0  # > total
        bad_node.ratio_total = 2.0  # > 1
        bad_node.ratio_annotated = 1.5  # > 1
        result.taxonomy_nodes[1] = bad_node

        violations = validate_aggregation_invariants(result)

        assert len(violations) > 0
        assert any("quantity" in v for v in violations)
        assert any("ratio_total" in v for v in violations)


class TestAggregationWithFixtureScenario:
    """Tests using the fixture scenario."""

    def test_full_fixture_scenario(self):
        """Test with PEPTIDE=10 (annotated), ABC=5, NOMATCH=3 (unannotated)."""
        # PEPTIDE: TAX={30,20,10,1}, GO={GO:0000001, GO:0000002, ...}
        annotations = [
            make_annotation(
                "PEPTIDE",
                10.0,
                is_annotated=True,
                taxonomy_nodes={30, 20, 10, 1},
                go_terms={"GO:0000004", "GO:0000005", "GO:0000006", "GO:0000002", "GO:0000003", "GO:0000001"},
            ),
            make_annotation("ABC", 5.0, is_annotated=False),
            make_annotation("NOMATCH", 3.0, is_annotated=False),
        ]

        result = aggregate_peptide_annotations(annotations)

        # Coverage
        assert result.coverage.total_peptide_quantity == 18.0
        assert result.coverage.annotated_peptide_quantity == 10.0
        assert result.coverage.annotation_coverage_ratio == pytest.approx(10.0 / 18.0)

        # Taxonomy nodes all get 10
        for tax_id in [30, 20, 10, 1]:
            assert result.taxonomy_nodes[tax_id].quantity == 10.0
            assert result.taxonomy_nodes[tax_id].n_peptides == 1
            assert result.taxonomy_nodes[tax_id].ratio_total == pytest.approx(10.0 / 18.0)
            assert result.taxonomy_nodes[tax_id].ratio_annotated == 1.0

        # GO terms all get 10
        for go_id in ["GO:0000004", "GO:0000005", "GO:0000006", "GO:0000002", "GO:0000003", "GO:0000001"]:
            assert result.go_terms[go_id].quantity == 10.0

        # Invariants hold
        violations = validate_aggregation_invariants(result)
        assert violations == []


class TestGoTaxonomyComboAggregation:
    """Tests for GO-taxonomy cross-tabulation."""

    def test_basic_combo(self):
        """Single peptide with one tax node and one GO term produces one combo."""
        annotations = [
            make_annotation(
                "P1", 10.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms={"GO:0000001"},
            ),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        assert len(combos) == 1
        assert (30, "GO:0000001") in combos
        combo = combos[(30, "GO:0000001")]
        assert combo.quantity == 10.0
        assert combo.n_peptides == 1
        assert combo.fraction_of_taxon == pytest.approx(1.0)
        assert combo.fraction_of_go == pytest.approx(1.0)

    def test_cartesian_product(self):
        """Peptide with 2 tax nodes and 2 GO terms produces 4 combos."""
        annotations = [
            make_annotation(
                "P1", 10.0, is_annotated=True,
                taxonomy_nodes={30, 1}, go_terms={"GO:0000001", "GO:0000002"},
            ),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        assert len(combos) == 4
        for tax_id in [30, 1]:
            for go_id in ["GO:0000001", "GO:0000002"]:
                assert (tax_id, go_id) in combos
                assert combos[(tax_id, go_id)].quantity == 10.0

    def test_multiple_peptides_accumulate(self):
        """Two peptides sharing a (tax, GO) pair accumulate quantities."""
        annotations = [
            make_annotation(
                "P1", 10.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms={"GO:0000001"},
            ),
            make_annotation(
                "P2", 5.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms={"GO:0000001"},
            ),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        combo = combos[(30, "GO:0000001")]
        assert combo.quantity == 15.0
        assert combo.n_peptides == 2

    def test_fractions_correct(self):
        """Fractions are computed relative to per-node totals."""
        annotations = [
            make_annotation(
                "P1", 10.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms={"GO:0000001", "GO:0000002"},
            ),
            make_annotation(
                "P2", 10.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms={"GO:0000001"},
            ),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        # tax 30 total = 20, GO:0000001 total = 20, GO:0000002 total = 10
        c1 = combos[(30, "GO:0000001")]
        assert c1.quantity == 20.0
        assert c1.fraction_of_taxon == pytest.approx(1.0)  # 20/20
        assert c1.fraction_of_go == pytest.approx(1.0)  # 20/20

        c2 = combos[(30, "GO:0000002")]
        assert c2.quantity == 10.0
        assert c2.fraction_of_taxon == pytest.approx(0.5)  # 10/20
        assert c2.fraction_of_go == pytest.approx(1.0)  # 10/10

    def test_unannotated_peptides_excluded(self):
        """Unannotated peptides produce no combos."""
        annotations = [
            make_annotation("P1", 10.0, is_annotated=False),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        assert len(combos) == 0

    def test_peptide_with_only_taxonomy_excluded(self):
        """Peptide with taxonomy but no GO terms produces no combos."""
        annotations = [
            make_annotation(
                "P1", 10.0, is_annotated=True,
                taxonomy_nodes={30}, go_terms=set(),
            ),
        ]
        agg = aggregate_peptide_annotations(annotations)
        combos = aggregate_go_taxonomy_combos(annotations, agg)

        assert len(combos) == 0

    def test_empty_annotations(self):
        """Empty annotations produce no combos."""
        agg = aggregate_peptide_annotations([])
        combos = aggregate_go_taxonomy_combos([], agg)

        assert len(combos) == 0
