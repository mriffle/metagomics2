"""Unit tests for single-sample GO x taxonomy enrichment."""

from math import erfc, sqrt

import pytest

from metagomics2.core.annotation import PeptideAnnotation
from metagomics2.core.enrichment import (
    benjamini_hochberg,
    compute_exact_weighted_pvalue,
    compute_go_taxonomy_enrichment,
    compute_weighted_rate_test,
    filter_doubly_annotated_peptides,
)


def make_annotation(
    peptide: str,
    quantity: float,
    *,
    is_annotated: bool = True,
    taxonomy_nodes: set[int] | None = None,
    go_terms: set[str] | None = None,
) -> PeptideAnnotation:
    """Create a test peptide annotation."""
    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=is_annotated,
        taxonomy_nodes=taxonomy_nodes or set(),
        go_terms=go_terms or set(),
    )


class TestFilterDoublyAnnotatedPeptides:
    """Tests for peptide pool filtering."""

    def test_requires_annotation_taxonomy_and_go(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={10}, go_terms={"GO:1"}),
            make_annotation("P2", 5.0, taxonomy_nodes={10}, go_terms=set()),
            make_annotation("P3", 5.0, taxonomy_nodes=set(), go_terms={"GO:1"}),
            make_annotation(
                "P4",
                5.0,
                is_annotated=False,
                taxonomy_nodes={10},
                go_terms={"GO:1"},
            ),
        ]

        pool = filter_doubly_annotated_peptides(annotations)

        assert [ann.peptide for ann in pool] == ["P1"]


class TestBenjaminiHochberg:
    """Tests for BH multiple-testing correction."""

    def test_preserves_original_order_and_monotonic_adjustment(self):
        qvalues = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])

        assert qvalues == pytest.approx([0.04, 0.0533333333, 0.0533333333, 0.2])


class TestWeightedPValues:
    """Tests for exact and approximate p-value helpers."""

    def test_exact_weighted_pvalue_matches_manual_enumeration(self):
        pvalue = compute_exact_weighted_pvalue(
            weights=[2.0, 1.0],
            indicators=[1, 0],
            background_rate=0.25,
        )

        assert pvalue == pytest.approx(0.25)

    def test_large_n_uses_normal_approximation(self):
        weights = [1.0] * 21
        indicators = [1] * 16 + [0] * 5

        pvalue, zscore = compute_weighted_rate_test(weights, indicators, 0.5)

        expected_z = ((16 / 21) - 0.5) / sqrt(0.25 / 21)
        expected_p = erfc(abs(expected_z) / sqrt(2.0))
        assert zscore == pytest.approx(expected_z)
        assert pvalue == pytest.approx(expected_p)


class TestComputeGoTaxonomyEnrichment:
    """Tests for per-pair enrichment statistics."""

    def test_ineligible_pairs_leave_stats_empty(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={1}, go_terms={"GO:1"}),
        ]

        stats = compute_go_taxonomy_enrichment(annotations, [(1, "GO:1")])

        assert stats[(1, "GO:1")].pvalue_go_for_taxon is None
        assert stats[(1, "GO:1")].pvalue_taxon_for_go is None
        assert stats[(1, "GO:1")].qvalue_go_for_taxon is None
        assert stats[(1, "GO:1")].qvalue_taxon_for_go is None

    def test_zero_background_rate_returns_zero_pvalues_and_undefined_zscores(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={10}, go_terms={"GO:1"}),
            make_annotation("P2", 10.0, taxonomy_nodes={20}, go_terms={"GO:2"}),
        ]

        stats = compute_go_taxonomy_enrichment(
            annotations,
            [(10, "GO:1"), (20, "GO:2")],
        )

        first = stats[(10, "GO:1")]
        second = stats[(20, "GO:2")]

        assert first.pvalue_go_for_taxon == 0.0
        assert first.pvalue_taxon_for_go == 0.0
        assert first.qvalue_go_for_taxon == 0.0
        assert first.qvalue_taxon_for_go == 0.0
        assert first.zscore_go_for_taxon == pytest.approx(float("inf"))
        assert first.zscore_taxon_for_go == pytest.approx(float("inf"))

        assert second.pvalue_go_for_taxon == 0.0
        assert second.pvalue_taxon_for_go == 0.0
        assert second.qvalue_go_for_taxon == 0.0
        assert second.qvalue_taxon_for_go == 0.0
        assert second.zscore_go_for_taxon == pytest.approx(float("inf"))
        assert second.zscore_taxon_for_go == pytest.approx(float("inf"))

    def test_zero_background_rate_large_group_stays_significant_without_enumeration(self):
        annotations = [
            make_annotation(f"TAXON_{index}", 1.0, taxonomy_nodes={10}, go_terms={"GO:1"})
            for index in range(21)
        ] + [
            make_annotation(f"OTHER_{index}", 1.0, taxonomy_nodes={20}, go_terms={"GO:2"})
            for index in range(21)
        ]

        stats = compute_go_taxonomy_enrichment(annotations, [(10, "GO:1")])

        result = stats[(10, "GO:1")]
        assert result.pvalue_go_for_taxon == 0.0
        assert result.pvalue_taxon_for_go == 0.0
        assert result.qvalue_go_for_taxon == 0.0
        assert result.qvalue_taxon_for_go == 0.0
        assert result.zscore_go_for_taxon == pytest.approx(float("inf"))
        assert result.zscore_taxon_for_go == pytest.approx(float("inf"))

    def test_pair_specific_eligibility_is_applied_per_direction(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={1, 10}, go_terms={"GO:ROOT", "GO:1"}),
            make_annotation("P2", 10.0, taxonomy_nodes={1, 20}, go_terms={"GO:ROOT", "GO:2"}),
        ]

        stats = compute_go_taxonomy_enrichment(
            annotations,
            [(1, "GO:1"), (10, "GO:ROOT")],
        )

        assert stats[(1, "GO:1")].pvalue_go_for_taxon is None
        assert stats[(1, "GO:1")].qvalue_go_for_taxon is None
        assert stats[(1, "GO:1")].pvalue_taxon_for_go is not None

        assert stats[(10, "GO:ROOT")].pvalue_taxon_for_go is None
        assert stats[(10, "GO:ROOT")].qvalue_taxon_for_go is None
        assert stats[(10, "GO:ROOT")].pvalue_go_for_taxon is not None

    def test_only_doubly_annotated_peptides_contribute_to_background(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={10}, go_terms={"GO:1"}),
            make_annotation("P2", 10.0, taxonomy_nodes={20}, go_terms={"GO:1"}),
            make_annotation("P3", 10.0, taxonomy_nodes={30}, go_terms=set()),
        ]

        stats = compute_go_taxonomy_enrichment(annotations, [(10, "GO:1")])

        assert stats[(10, "GO:1")].pvalue_go_for_taxon == pytest.approx(1.0)
        assert stats[(10, "GO:1")].pvalue_taxon_for_go is None
