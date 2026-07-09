"""Unit tests for single-sample GO x taxonomy enrichment (finite-population null)."""

from itertools import combinations
from math import comb

import pytest

from metagomics2.core.annotation import PeptideAnnotation
from metagomics2.core.enrichment import (
    _exact_fpc_pvalue,
    _saddlepoint_two_sided,
    _subset_sums_by_size,
    benjamini_hochberg,
    collapse_identical_annotations,
    compute_go_taxonomy_enrichment,
    filter_doubly_annotated_peptides,
    format_optional_stat,
    fpc_moments,
)


def make_annotation(
    peptide: str,
    quantity: float,
    *,
    is_annotated: bool = True,
    taxonomy_nodes: set[int] | None = None,
    go_terms: set[str] | None = None,
) -> PeptideAnnotation:
    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=is_annotated,
        taxonomy_nodes=taxonomy_nodes or set(),
        go_terms=go_terms or set(),
    )


# --------------------------------------------------------------------------- #
# Independent brute-force references for the finite-population null
# --------------------------------------------------------------------------- #
def brute_moments(weights, P, n):
    m = len(weights)
    denom = comb(P, n)
    mean = ex2 = 0.0
    for h in range(m + 1):
        if not (0 <= n - h <= P - m):
            continue
        prob = comb(P - m, n - h) / denom
        for subset in combinations(range(m), h):
            s = sum(weights[i] for i in subset)
            mean += prob * s
            ex2 += prob * s * s
    return mean, ex2 - mean * mean


def brute_two_sided_p(weights, P, n, observed):
    m = len(weights)
    denom = comb(P, n)
    mean = (n / P) * sum(weights)
    delta = abs(observed - mean)
    p = 0.0
    for h in range(m + 1):
        if not (0 <= n - h <= P - m):
            continue
        prob = comb(P - m, n - h) / denom
        for subset in combinations(range(m), h):
            s = sum(weights[i] for i in subset)
            if abs(s - mean) >= delta - 1e-9:
                p += prob
    return min(max(p, 0.0), 1.0)


class TestFilterDoublyAnnotatedPeptides:
    def test_requires_annotation_taxonomy_and_go(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={10}, go_terms={"GO:1"}),
            make_annotation("P2", 5.0, taxonomy_nodes={10}, go_terms=set()),
            make_annotation("P3", 5.0, taxonomy_nodes=set(), go_terms={"GO:1"}),
            make_annotation("P4", 5.0, is_annotated=False, taxonomy_nodes={10}, go_terms={"GO:1"}),
        ]
        pool = filter_doubly_annotated_peptides(annotations)
        assert [ann.peptide for ann in pool] == ["P1"]


class TestCollapseIdenticalAnnotations:
    def test_merges_identical_signatures_and_sums_weight(self):
        pool = [
            make_annotation("P1", 10.0, taxonomy_nodes={10, 1}, go_terms={"GO:1"}),
            make_annotation("P2", 5.0, taxonomy_nodes={10, 1}, go_terms={"GO:1"}),
            make_annotation("P3", 7.0, taxonomy_nodes={20, 1}, go_terms={"GO:2"}),
        ]
        collapsed = collapse_identical_annotations(pool)
        assert len(collapsed) == 2
        by_weight = sorted(a.quantity for a in collapsed)
        assert by_weight == pytest.approx([7.0, 15.0])


class TestFpcMoments:
    @pytest.mark.parametrize("P,weights,n", [
        (30, [2.0, 1.0, 3.0], 10),
        (50, [1.0, 1.0, 1.0, 1.0, 5.0], 25),
        (12, [4.0, 1.0], 3),
    ])
    def test_closed_form_matches_exhaustive_enumeration(self, P, weights, n):
        mean, var = fpc_moments(sum(weights), sum(w * w for w in weights), P, n)
        b_mean, b_var = brute_moments(weights, P, n)
        assert mean == pytest.approx(b_mean)
        assert var == pytest.approx(b_var)

    def test_variance_zero_when_feature_universal_or_absent(self):
        _, var0 = fpc_moments(10.0, 40.0, 20, 0)
        _, varP = fpc_moments(10.0, 40.0, 20, 20)
        assert var0 == 0.0 and varP == 0.0


class TestExactFpcPvalue:
    @pytest.mark.parametrize("P,weights,n,obs", [
        (30, [2.0, 1.0, 3.0], 10, 5.0),
        (40, [1.0, 1.0, 1.0, 2.0], 8, 4.0),
        (25, [3.0, 1.5, 2.0, 0.5, 1.0], 12, 7.0),
    ])
    def test_matches_exhaustive_enumeration(self, P, weights, n, obs):
        table = _subset_sums_by_size(tuple(weights))
        p = _exact_fpc_pvalue(table, sum(weights), P, n, obs)
        assert p == pytest.approx(brute_two_sided_p(weights, P, n, obs), abs=1e-9)

    def test_pvalue_is_bounded_and_positive(self):
        table = _subset_sums_by_size((2.0, 1.0, 3.0))
        # observed at the extreme still includes its own mass => p > 0
        p = _exact_fpc_pvalue(table, 6.0, 30, 10, 6.0)
        assert 0.0 < p <= 1.0


class TestSaddlepoint:
    """The double-saddlepoint large-group tail vs the exact finite-population p."""

    def test_matches_exact_fpc_for_continuous_weights(self):
        # A group just above the exact threshold with heavy-tailed (distinct)
        # weights: the double saddlepoint should track the exact p across the tail.
        weights = tuple(round(1.4**i, 4) for i in range(16))
        P, n = 200, 60
        total = sum(weights)
        sq = sum(w * w for w in weights)
        table = _subset_sums_by_size(weights)
        mean, var = fpc_moments(total, sq, P, n)
        for frac in (0.55, 0.65, 0.75):
            obs = frac * total
            exact = _exact_fpc_pvalue(table, total, P, n, obs)
            saddle = _saddlepoint_two_sided(weights, total, P, n, obs, mean, var)
            assert saddle == pytest.approx(exact, rel=0.1)

    def test_pure_group_defers_to_normal_and_is_significant(self):
        # observed == group total (a "pure" group): the saddlepoint defers to the
        # normal, which reports a hugely-significant p-value.
        weights = tuple(1.0 + 0.1 * i for i in range(30))
        P, n = 60, 30
        total = sum(weights)
        mean, var = fpc_moments(total, sum(w * w for w in weights), P, n)
        p = _saddlepoint_two_sided(weights, total, P, n, total, mean, var)
        assert 0.0 < p < 1e-6


class TestBenjaminiHochberg:
    def test_preserves_order_and_monotonic_adjustment(self):
        q = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
        assert q == pytest.approx([0.04, 0.0533333333, 0.0533333333, 0.2])


class TestComputeGoTaxonomyEnrichment:
    def test_ineligible_pairs_leave_stats_empty(self):
        annotations = [make_annotation("P1", 10.0, taxonomy_nodes={1}, go_terms={"GO:1"})]
        stats = compute_go_taxonomy_enrichment(annotations, [(1, "GO:1")])
        # single peptide => taxon group is the whole pool => ineligible both directions
        s = stats[(1, "GO:1")]
        assert s.pvalue_go_for_taxon is None
        assert s.pvalue_taxon_for_go is None
        assert s.qvalue_go_for_taxon is None

    def test_whole_pool_group_is_ineligible(self):
        # taxon 1 is in every peptide (root) => GO-for-taxon ineligible for it
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={1, 10}, go_terms={"GO:1"}),
            make_annotation("P2", 10.0, taxonomy_nodes={1, 20}, go_terms={"GO:1"}),
            make_annotation("P3", 10.0, taxonomy_nodes={1, 30}, go_terms={"GO:2"}),
        ]
        stats = compute_go_taxonomy_enrichment(annotations, [(1, "GO:1"), (10, "GO:1")])
        assert stats[(1, "GO:1")].pvalue_go_for_taxon is None  # root taxon: no complement
        assert stats[(10, "GO:1")].pvalue_go_for_taxon is not None

    def test_detects_strong_enrichment(self):
        # Distinct units (unique second GO term per peptide) so they are NOT
        # pseudo-replicates: taxon 10 carries GO:1, the rest does not.
        annotations = (
            [make_annotation(f"A{i}", 1.0, taxonomy_nodes={10}, go_terms={"GO:1", f"GO:A{i}"}) for i in range(20)]
            + [make_annotation(f"B{i}", 1.0, taxonomy_nodes={20}, go_terms={f"GO:B{i}"}) for i in range(20)]
        )
        stats = compute_go_taxonomy_enrichment(annotations, [(10, "GO:1")])
        result = stats[(10, "GO:1")]
        assert result.pvalue_go_for_taxon < 1e-6
        assert result.zscore_go_for_taxon > 0  # enrichment
        assert result.qvalue_go_for_taxon <= 1.0

    def test_pair_specific_eligibility_is_applied_per_direction(self):
        annotations = [
            make_annotation("P1", 10.0, taxonomy_nodes={1, 10}, go_terms={"GO:R", "GO:1"}),
            make_annotation("P2", 10.0, taxonomy_nodes={1, 20}, go_terms={"GO:R", "GO:2"}),
            make_annotation("P3", 10.0, taxonomy_nodes={1, 30}, go_terms={"GO:R", "GO:3"}),
        ]
        stats = compute_go_taxonomy_enrichment(annotations, [(10, "GO:R"), (10, "GO:1")])
        # GO:R is in every peptide => taxon-for-GO ineligible (GO group == whole pool)
        assert stats[(10, "GO:R")].pvalue_taxon_for_go is None
        assert stats[(10, "GO:R")].pvalue_go_for_taxon is not None

    def test_collapse_makes_result_invariant_to_pseudo_replication(self):
        base = (
            [make_annotation(f"A{i}", 1.0, taxonomy_nodes={10, 1}, go_terms={"GO:1", "GO:R"}) for i in range(8)]
            + [make_annotation(f"B{i}", 1.0, taxonomy_nodes={20, 1}, go_terms={"GO:2", "GO:R"}) for i in range(8)]
            + [make_annotation(f"C{i}", 1.0, taxonomy_nodes={30, 1}, go_terms={"GO:1", "GO:R"}) for i in range(8)]
        )
        triplicated = [
            make_annotation(f"{a.peptide}_{k}", a.quantity,
                            taxonomy_nodes=set(a.taxonomy_nodes), go_terms=set(a.go_terms))
            for a in base for k in range(3)
        ]
        keys = [(10, "GO:1"), (30, "GO:1"), (20, "GO:2")]
        base_stats = compute_go_taxonomy_enrichment(base, keys)
        trip_stats = compute_go_taxonomy_enrichment(triplicated, keys)
        for key in keys:
            assert trip_stats[key].pvalue_go_for_taxon == pytest.approx(
                base_stats[key].pvalue_go_for_taxon
            )

    def test_disabling_collapse_changes_significance(self):
        # identical-annotation pseudo-replicates inflate significance when not collapsed
        pool = (
            [make_annotation(f"A{i}", 1.0, taxonomy_nodes={10, 1}, go_terms={"GO:1", "GO:R"}) for i in range(6)]
            + [make_annotation(f"B{i}", 1.0, taxonomy_nodes={20, 1}, go_terms={"GO:2", "GO:R"}) for i in range(6)]
        )
        keys = [(10, "GO:1")]
        collapsed = compute_go_taxonomy_enrichment(pool, keys, collapse_identical=True)
        raw = compute_go_taxonomy_enrichment(pool, keys, collapse_identical=False)
        # not collapsing treats 6 identical peptides as 6 independent units => smaller p
        assert raw[(10, "GO:1")].pvalue_go_for_taxon <= collapsed[(10, "GO:1")].pvalue_go_for_taxon


class TestFormatOptionalStat:
    def test_scientific_notation_and_specials(self):
        assert format_optional_stat(None) == ""
        assert format_optional_stat(float("inf")) == "+inf"
        assert format_optional_stat(float("-inf")) == "-inf"
        assert format_optional_stat(1e-50) == "1.000000e-50"
        assert float(format_optional_stat(0.05)) == pytest.approx(0.05)
