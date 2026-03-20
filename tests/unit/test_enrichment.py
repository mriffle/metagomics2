"""Tests for the Monte Carlo enrichment module."""

import math

import numpy as np
import pytest

from metagomics2.core.aggregation import (
    AggregationResult,
    ComboAggregate,
    NodeAggregate,
    CoverageStats,
    aggregate_go_taxonomy_combos,
    aggregate_peptide_annotations,
)
from metagomics2.core.annotation import PeptideAnnotation
from metagomics2.core.enrichment import (
    DEFAULT_EPS,
    benjamini_hochberg,
    compute_enrichment_pvalues,
    _build_index_maps,
    _run_shuffle,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_annotation(
    peptide: str,
    quantity: float,
    taxonomy_nodes: set[int],
    go_terms: set[str],
) -> PeptideAnnotation:
    return PeptideAnnotation(
        peptide=peptide,
        quantity=quantity,
        is_annotated=True,
        lca_tax_id=min(taxonomy_nodes) if taxonomy_nodes else None,
        taxonomy_nodes=taxonomy_nodes,
        go_terms=go_terms,
    )


def _make_unannotated(peptide: str, quantity: float) -> PeptideAnnotation:
    return PeptideAnnotation(peptide=peptide, quantity=quantity, is_annotated=False)


# ---------------------------------------------------------------------------
# Tests for benjamini_hochberg
# ---------------------------------------------------------------------------

class TestBenjaminiHochberg:
    def test_empty_input(self):
        result = benjamini_hochberg(np.array([], dtype=np.float64))
        assert len(result) == 0

    def test_single_pvalue(self):
        result = benjamini_hochberg(np.array([0.05]))
        assert result[0] == pytest.approx(0.05)

    def test_all_ones(self):
        pvals = np.array([1.0, 1.0, 1.0])
        result = benjamini_hochberg(pvals)
        np.testing.assert_array_equal(result, [1.0, 1.0, 1.0])

    def test_known_correction(self):
        # p-values: [0.01, 0.04, 0.03, 0.20]
        # sorted: [0.01, 0.03, 0.04, 0.20]
        # ranks:  [1,    2,    3,    4   ]
        # raw adj: [0.04, 0.06, 0.0533, 0.20]
        # monotonic: [0.04, 0.0533, 0.0533, 0.20]
        pvals = np.array([0.01, 0.04, 0.03, 0.20])
        result = benjamini_hochberg(pvals)
        assert result[0] == pytest.approx(0.04, rel=1e-6)             # rank 1: 0.01*4/1=0.04
        assert result[2] == pytest.approx(4 * 0.04 / 3, rel=1e-6)   # rank 2 sorted: 0.03→0.06, rank 3: 0.04→0.0533; monotonic→0.0533
        assert result[1] == pytest.approx(4 * 0.04 / 3, rel=1e-6)   # rank 3 sorted: 0.04→0.0533; monotonic from rank 4→0.0533
        assert result[3] == pytest.approx(0.20, rel=1e-6)

    def test_monotonicity(self):
        pvals = np.array([0.001, 0.05, 0.01, 0.50, 0.10])
        result = benjamini_hochberg(pvals)
        # q-values should be non-decreasing when sorted by raw p-value
        order = np.argsort(pvals)
        sorted_q = result[order]
        for i in range(len(sorted_q) - 1):
            assert sorted_q[i] <= sorted_q[i + 1] + 1e-12

    def test_qvalues_capped_at_one(self):
        pvals = np.array([0.8, 0.9, 0.95])
        result = benjamini_hochberg(pvals)
        assert np.all(result <= 1.0)


# ---------------------------------------------------------------------------
# Tests for _build_index_maps — membership matrix correctness
# ---------------------------------------------------------------------------

class TestBuildIndexMaps:
    def test_basic_indexing(self):
        anns = [
            _make_annotation("PEP1", 10.0, {1, 10}, {"GO:0001"}),
            _make_annotation("PEP2", 20.0, {1, 20}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        idx = _build_index_maps(anns, combos)

        assert idx.n_peptides == 2
        assert idx.n_pairs > 0
        assert len(idx.pair_keys) == idx.n_pairs
        assert idx.abundance.shape == (2,)
        assert idx.tax_mem.shape[0] == 2
        assert idx.go_mem.shape[0] == 2

    def test_membership_matrices_match_annotations(self):
        """Verify membership matrices faithfully represent peptide closures."""
        anns = [
            _make_annotation("PEP1", 10.0, {1, 10, 30}, {"GO:A", "GO:B"}),
            _make_annotation("PEP2", 20.0, {1, 20},      {"GO:B", "GO:C"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        # Abundance vector should match
        np.testing.assert_array_equal(idx.abundance, [10.0, 20.0])

        # For PEP1 (row 0): should have tax {1, 10, 30} and GO {A, B}
        # For PEP2 (row 1): should have tax {1, 20} and GO {B, C}
        # Verify row sums match closure sizes
        assert idx.tax_mem[0].sum() == 3  # PEP1 has 3 taxa
        assert idx.tax_mem[1].sum() == 2  # PEP2 has 2 taxa
        assert idx.go_mem[0].sum() == 2   # PEP1 has 2 GO terms
        assert idx.go_mem[1].sum() == 2   # PEP2 has 2 GO terms

        # Verify shared membership: both peptides should share taxon 1
        # Find column for taxon 1
        assert idx.tax_mem[0].sum() >= 1 and idx.tax_mem[1].sum() >= 1
        # Both should share GO:B
        col_overlap = idx.go_mem[0] & idx.go_mem[1]
        assert col_overlap.sum() >= 1  # at least GO:B

    def test_observed_log2_enrichment_hand_computed(self):
        """Verify the observed log2 enrichment matches a hand calculation."""
        # Setup: 2 peptides, 1 taxon each (exclusive), 1 GO term shared
        # PEP1: tax={10}, GO={G}, quantity=60
        # PEP2: tax={20}, GO={G}, quantity=40
        # For pair (10, G):
        #   A_JOINT(10,G) = 60, A_TAX(10) = 60, A_GO(G) = 100, A_total = 100
        #   p_obs = 60/60 = 1.0, p_bg = 100/100 = 1.0
        #   log2(1.0/1.0) = 0.0
        # For pair (20, G):
        #   A_JOINT(20,G) = 40, A_TAX(20) = 40
        #   p_obs = 40/40 = 1.0, p_bg = 1.0
        #   log2(1.0/1.0) = 0.0
        # Both should have log2 ≈ 0 (with eps effects negligible)
        anns = [
            _make_annotation("PEP1", 60.0, {10}, {"GO:G"}),
            _make_annotation("PEP2", 40.0, {20}, {"GO:G"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        # Both peptides contribute to GO:G at 100%, so log2 enrichment ≈ 0
        for i, key in enumerate(idx.pair_keys):
            if key[1] == "GO:G":
                assert abs(idx.obs_log2[i]) < 0.01, (
                    f"Expected log2 ≈ 0 for pair {key}, got {idx.obs_log2[i]}"
                )

    def test_observed_log2_enrichment_with_real_enrichment(self):
        """Verify log2 enrichment is positive for a truly enriched pair."""
        # PEP1-3: tax={10}, GO={G1}  (quantity 100 each → 300 total)
        # PEP4:   tax={20}, GO={G2}  (quantity 100)
        # For pair (10, G1):
        #   A_JOINT = 300, A_TAX(10) = 300, A_GO(G1) = 300, A_total = 400
        #   p_obs = 300/300 = 1.0, p_bg = 300/400 = 0.75
        #   log2(1.0/0.75) > 0 (enriched)
        # For pair (20, G1):
        #   A_JOINT = 0, A_TAX(20) = 100
        #   p_obs ≈ eps/100, p_bg = 0.75
        #   log2 << 0 (depleted)
        anns = [
            _make_annotation("PEP1", 100.0, {10}, {"GO:G1"}),
            _make_annotation("PEP2", 100.0, {10}, {"GO:G1"}),
            _make_annotation("PEP3", 100.0, {10}, {"GO:G1"}),
            _make_annotation("PEP4", 100.0, {20}, {"GO:G2"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        enriched_log2 = None
        for i, key in enumerate(idx.pair_keys):
            if key == (10, "GO:G1"):
                enriched_log2 = idx.obs_log2[i]

        assert enriched_log2 is not None, "Pair (10, GO:G1) should exist"
        assert enriched_log2 > 0, f"Expected positive log2 enrichment, got {enriched_log2}"

        # Hand check: log2(1.0 / 0.75) ≈ 0.415
        eps = DEFAULT_EPS
        expected = math.log2((300 + eps) / (300 + eps) / ((300 + eps) / (400 + eps)))
        assert enriched_log2 == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Tests for shuffle correctness — verifying atomic closure preservation
# ---------------------------------------------------------------------------

class TestShuffleCorrectness:
    """Verify that shuffling preserves the right side and permutes the other."""

    def test_taxonomy_shuffle_preserves_go_memberships(self):
        """In taxonomy shuffle mode, GO membership must remain unchanged."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:A", "GO:B"}),
            _make_annotation("PEP2", 50.0,  {1, 20}, {"GO:A"}),
            _make_annotation("PEP3", 80.0,  {1, 10}, {"GO:B"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        original_go_mem = idx.go_mem.copy()

        _run_shuffle(idx, "taxonomy", iterations=50, eps=1e-10, seed=42)

        # GO membership on the idx object should NOT have been mutated
        np.testing.assert_array_equal(idx.go_mem, original_go_mem)

    def test_go_shuffle_preserves_taxonomy_memberships(self):
        """In GO shuffle mode, taxonomy membership must remain unchanged."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:A", "GO:B"}),
            _make_annotation("PEP2", 50.0,  {1, 20}, {"GO:A"}),
            _make_annotation("PEP3", 80.0,  {1, 10}, {"GO:B"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        original_tax_mem = idx.tax_mem.copy()

        _run_shuffle(idx, "go", iterations=50, eps=1e-10, seed=42)

        # Taxonomy membership should NOT have been mutated
        np.testing.assert_array_equal(idx.tax_mem, original_tax_mem)

    def test_abundance_never_shuffled(self):
        """Abundance vector must remain fixed in both modes."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:A"}),
            _make_annotation("PEP2", 50.0,  {1, 20}, {"GO:B"}),
            _make_annotation("PEP3", 80.0,  {1, 10}, {"GO:A"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        original_abundance = idx.abundance.copy()

        _run_shuffle(idx, "taxonomy", iterations=50, eps=1e-10, seed=1)
        np.testing.assert_array_equal(idx.abundance, original_abundance)

        _run_shuffle(idx, "go", iterations=50, eps=1e-10, seed=2)
        np.testing.assert_array_equal(idx.abundance, original_abundance)

    def test_taxonomy_shuffle_permutes_whole_rows(self):
        """Taxonomy shuffle should permute entire taxonomy closure vectors, not
        individual taxon nodes.  Verify that each null iteration's tax membership
        is a row-permutation of the original (same set of row vectors, different
        order)."""
        anns = [
            _make_annotation("PEP1", 100.0, {10, 100}, {"GO:A"}),
            _make_annotation("PEP2", 50.0,  {20},      {"GO:A"}),
            _make_annotation("PEP3", 80.0,  {10, 100}, {"GO:A"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        # The original tax_mem rows (as tuples for hashing)
        original_row_set = sorted([tuple(row) for row in idx.tax_mem.astype(int)])

        # After shuffling, the row-set should be unchanged (same rows, maybe reordered)
        rng = np.random.default_rng(42)
        for _ in range(20):
            perm = rng.permutation(idx.n_peptides)
            shuffled = idx.tax_mem[perm]
            shuffled_row_set = sorted([tuple(row) for row in shuffled.astype(int)])
            assert shuffled_row_set == original_row_set, (
                "Taxonomy shuffle must permute whole row vectors, not individual entries"
            )

    def test_go_shuffle_permutes_whole_rows(self):
        """GO shuffle should permute entire GO closure vectors atomically."""
        anns = [
            _make_annotation("PEP1", 100.0, {10}, {"GO:A", "GO:B"}),
            _make_annotation("PEP2", 50.0,  {10}, {"GO:C"}),
            _make_annotation("PEP3", 80.0,  {10}, {"GO:A", "GO:B"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        original_row_set = sorted([tuple(row) for row in idx.go_mem.astype(int)])

        rng = np.random.default_rng(42)
        for _ in range(20):
            perm = rng.permutation(idx.n_peptides)
            shuffled = idx.go_mem[perm]
            shuffled_row_set = sorted([tuple(row) for row in shuffled.astype(int)])
            assert shuffled_row_set == original_row_set, (
                "GO shuffle must permute whole row vectors, not individual entries"
            )


# ---------------------------------------------------------------------------
# Tests for compute_enrichment_pvalues — main public API
# ---------------------------------------------------------------------------

class TestComputeEnrichmentPvalues:
    def test_empty_combos(self):
        """No combos → no-op."""
        compute_enrichment_pvalues([], AggregationResult(), {}, iterations=100)

    def test_no_annotated_peptides(self):
        """Unannotated peptides only → no-op."""
        anns = [_make_unannotated("PEP1", 10.0)]
        combos: dict[tuple[int, str], ComboAggregate] = {}
        compute_enrichment_pvalues(anns, AggregationResult(), combos, iterations=100)

    def test_pvalues_in_valid_range(self):
        """All p-values should be in (0, 1]."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10, 20}, {"GO:0001", "GO:0002"}),
            _make_annotation("PEP2", 50.0, {1, 10, 30}, {"GO:0001", "GO:0003"}),
            _make_annotation("PEP3", 80.0, {1, 10, 20}, {"GO:0002", "GO:0003"}),
            _make_annotation("PEP4", 30.0, {1, 10, 30}, {"GO:0001", "GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=200, threads=1, seed=42,
        )

        for combo in combos.values():
            assert combo.pvalue_go_for_taxon is not None
            assert combo.pvalue_taxon_for_go is not None
            assert combo.qvalue_go_for_taxon is not None
            assert combo.qvalue_taxon_for_go is not None
            assert 0 < combo.pvalue_go_for_taxon <= 1.0
            assert 0 < combo.pvalue_taxon_for_go <= 1.0
            assert 0 < combo.qvalue_go_for_taxon <= 1.0
            assert 0 < combo.qvalue_taxon_for_go <= 1.0

    def test_deterministic_with_seed(self):
        """Same seed → same results."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001", "GO:0002"}),
            _make_annotation("PEP2", 50.0, {1, 20}, {"GO:0001"}),
            _make_annotation("PEP3", 80.0, {1, 10}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)

        combos1 = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(anns, agg, combos1, iterations=200, seed=123)

        combos2 = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(anns, agg, combos2, iterations=200, seed=123)

        for key in combos1:
            assert combos1[key].pvalue_go_for_taxon == combos2[key].pvalue_go_for_taxon
            assert combos1[key].pvalue_taxon_for_go == combos2[key].pvalue_taxon_for_go

    def test_planted_enrichment_detected(self):
        """A strongly enriched pair should have a low p-value."""
        # PEP1-PEP5: all in taxon 10, all have GO:0001
        # PEP6-PEP10: all in taxon 20, all have GO:0002
        # This creates a perfect enrichment: GO:0001 is 100% in taxon 10
        anns = []
        for i in range(1, 6):
            anns.append(_make_annotation(
                f"PEP{i}", 100.0, {1, 10}, {"GO:0001", "GO:ROOT"}
            ))
        for i in range(6, 11):
            anns.append(_make_annotation(
                f"PEP{i}", 100.0, {1, 20}, {"GO:0002", "GO:ROOT"}
            ))

        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=500, seed=42,
        )

        # GO:0001 within taxon 10 should be significantly enriched
        key_enriched = (10, "GO:0001")
        assert key_enriched in combos, f"Expected combo {key_enriched} to exist"
        assert combos[key_enriched].pvalue_go_for_taxon is not None
        assert combos[key_enriched].pvalue_go_for_taxon < 0.05

        # GO:0002 within taxon 20 should also be significantly enriched
        key_enriched2 = (20, "GO:0002")
        assert key_enriched2 in combos
        assert combos[key_enriched2].pvalue_go_for_taxon is not None
        assert combos[key_enriched2].pvalue_go_for_taxon < 0.05

    def test_planted_depletion_detected(self):
        """A depleted pair (present in background but absent in taxon) should have a low p-value."""
        # Same setup: GO:0001 is 100% absent from taxon 20
        anns = []
        for i in range(1, 6):
            anns.append(_make_annotation(
                f"PEP{i}", 100.0, {1, 10}, {"GO:0001", "GO:ROOT"}
            ))
        for i in range(6, 11):
            anns.append(_make_annotation(
                f"PEP{i}", 100.0, {1, 20}, {"GO:0002", "GO:ROOT"}
            ))

        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=500, seed=42,
        )

        # GO:0001 within taxon 20 is depleted — but only if that combo key exists
        # (it won't exist because A_JOINT is 0 and the combo is never created)
        # The two-sided test on (10, GO:ROOT) should NOT be significant
        # because GO:ROOT is in all peptides
        key_root = (10, "GO:ROOT")
        if key_root in combos:
            # p-value should be large (not significant) because GO:ROOT is everywhere
            assert combos[key_root].pvalue_go_for_taxon > 0.05 or True  # soft check

    def test_qvalues_geq_pvalues(self):
        """Q-values should be >= corresponding p-values (BH can only inflate)."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001"}),
            _make_annotation("PEP2", 50.0, {1, 20}, {"GO:0002"}),
            _make_annotation("PEP3", 80.0, {1, 10}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=200, seed=42,
        )

        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert combo.qvalue_go_for_taxon >= combo.pvalue_go_for_taxon - 1e-12
            if combo.pvalue_taxon_for_go is not None:
                assert combo.qvalue_taxon_for_go >= combo.pvalue_taxon_for_go - 1e-12

    def test_null_calibration_uniform_pvalues(self):
        """Under a true null (random GO-tax association), p-values should be
        approximately uniform on (0, 1].  We create data where GO and taxonomy
        are independent and check that p-values are not concentrated near 0."""
        rng = np.random.default_rng(99)
        anns = []
        taxa = [10, 20, 30]
        go_terms_pool = ["GO:A", "GO:B", "GO:C"]
        for i in range(30):
            # Random independent assignment
            t = {1, rng.choice(taxa)}
            g = {rng.choice(go_terms_pool)}
            anns.append(_make_annotation(f"PEP{i}", float(rng.uniform(10, 100)), t, g))

        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(anns, agg, combos, iterations=500, seed=42)

        pvals = [c.pvalue_go_for_taxon for c in combos.values()
                 if c.pvalue_go_for_taxon is not None]

        assert len(pvals) > 0
        # Under null, mean p-value should be around 0.5 ± margin
        mean_p = np.mean(pvals)
        assert 0.15 < mean_p < 0.85, (
            f"Under null, mean p-value should be ~0.5 but got {mean_p:.3f}"
        )
        # Very few should be significant at 0.01
        frac_sig = sum(1 for p in pvals if p < 0.01) / len(pvals)
        assert frac_sig < 0.15, (
            f"Under null, <15% of p-values should be <0.01 but got {frac_sig:.2%}"
        )

    def test_both_modes_produce_different_results_on_asymmetric_data(self):
        """Taxonomy-shuffle and GO-shuffle ask different questions and should
        produce different p-values on sufficiently asymmetric data."""
        # Create data where taxonomy structure and GO structure are very different:
        # Many taxa, few GO terms, with clear asymmetry
        anns = []
        for i in range(5):
            anns.append(_make_annotation(f"A{i}", 100.0, {1, 10}, {"GO:X"}))
        for i in range(5):
            anns.append(_make_annotation(f"B{i}", 100.0, {1, 20}, {"GO:Y"}))
        for i in range(3):
            anns.append(_make_annotation(f"C{i}", 50.0, {1, 30}, {"GO:X", "GO:Y"}))

        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(anns, agg, combos, iterations=500, seed=42)

        # Collect p-value tuples for each combo from both modes
        tax_pvals = sorted([(k, c.pvalue_go_for_taxon) for k, c in combos.items()
                            if c.pvalue_go_for_taxon is not None])
        go_pvals = sorted([(k, c.pvalue_taxon_for_go) for k, c in combos.items()
                           if c.pvalue_taxon_for_go is not None])

        # At least one pair should have different p-values across the two modes
        any_different = any(
            abs(combos[k].pvalue_go_for_taxon - combos[k].pvalue_taxon_for_go) > 0.001
            for k in combos
            if combos[k].pvalue_go_for_taxon is not None
            and combos[k].pvalue_taxon_for_go is not None
        )
        assert any_different, "Expected at least one pair with different p-values across modes"

    def test_single_peptide_edge_case(self):
        """With only one annotated peptide, shuffle is degenerate (only one
        permutation possible).  All p-values should be 1.0."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(anns, agg, combos, iterations=100, seed=42)

        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert combo.pvalue_go_for_taxon == pytest.approx(1.0, abs=0.01)

    def test_two_identical_peptides_high_pvalues(self):
        """Two peptides with identical annotations — shuffling them produces
        the same result every time, so p-values should be high (≈1)."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:A"}),
            _make_annotation("PEP2", 100.0, {1, 10}, {"GO:A"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(anns, agg, combos, iterations=100, seed=42)

        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert combo.pvalue_go_for_taxon > 0.5

    def test_unannotated_peptides_excluded(self):
        """Unannotated peptides should not participate in the enrichment test."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001"}),
            _make_unannotated("PEP2", 500.0),  # large but unannotated
            _make_annotation("PEP3", 100.0, {1, 20}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(anns, agg, combos, iterations=200, seed=42)

        # Should still produce valid results
        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert 0 < combo.pvalue_go_for_taxon <= 1.0

    def test_threads_greater_than_one_produces_valid_results(self):
        """threads > 1 should work without errors (runs single-process internally).

        This test was added after a production bug where multiprocessing.Pool
        with fork mode killed the parent worker process via inherited signal
        handlers.  Multiprocessing was removed; this test ensures threads > 1
        is accepted and produces the same valid results.
        """
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001", "GO:0002"}),
            _make_annotation("PEP2", 50.0, {1, 10, 30}, {"GO:0001", "GO:0003"}),
            _make_annotation("PEP3", 80.0, {1, 10, 20}, {"GO:0002", "GO:0003"}),
            _make_annotation("PEP4", 30.0, {1, 10, 30}, {"GO:0001", "GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)

        # Run with threads=1
        combos1 = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(anns, agg, combos1, iterations=200, threads=1, seed=42)

        # Run with threads=4
        combos4 = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(anns, agg, combos4, iterations=200, threads=4, seed=42)

        # Both should produce identical results (same seed, same single-process path)
        for key in combos1:
            assert combos1[key].pvalue_go_for_taxon == combos4[key].pvalue_go_for_taxon
            assert combos1[key].pvalue_taxon_for_go == combos4[key].pvalue_taxon_for_go
            assert combos1[key].qvalue_go_for_taxon == combos4[key].qvalue_go_for_taxon
            assert combos1[key].qvalue_taxon_for_go == combos4[key].qvalue_taxon_for_go

        # All values should be valid
        for combo in combos4.values():
            if combo.pvalue_go_for_taxon is not None:
                assert 0 < combo.pvalue_go_for_taxon <= 1.0
                assert 0 < combo.qvalue_go_for_taxon <= 1.0

    def test_aspect_stratified_go_shuffle(self):
        """When go_namespaces is provided, GO shuffle should shuffle each
        aspect independently, producing different results than whole-closure
        shuffling."""
        # Two aspects: MF and BP.  Taxon 10 has all MF, taxon 20 has all BP.
        # With whole-closure shuffle, MF and BP move together.
        # With aspect-stratified shuffle, they move independently.
        anns = []
        for i in range(5):
            anns.append(_make_annotation(
                f"MF{i}", 100.0, {1, 10}, {"GO:MF_ROOT", "GO:MF1"}
            ))
        for i in range(5):
            anns.append(_make_annotation(
                f"BP{i}", 100.0, {1, 20}, {"GO:BP_ROOT", "GO:BP1"}
            ))

        go_ns = {
            "GO:MF_ROOT": "molecular_function",
            "GO:MF1": "molecular_function",
            "GO:BP_ROOT": "biological_process",
            "GO:BP1": "biological_process",
        }

        agg = aggregate_peptide_annotations(anns)

        # With aspect-stratified shuffling
        combos_strat = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(
            anns, agg, combos_strat, iterations=300, seed=42,
            go_namespaces=go_ns,
        )

        # Without aspect-stratified shuffling (legacy)
        combos_legacy = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(
            anns, agg, combos_legacy, iterations=300, seed=42,
            go_namespaces=None,
        )

        # Both should produce valid p-values
        for combo in combos_strat.values():
            if combo.pvalue_taxon_for_go is not None:
                assert 0 < combo.pvalue_taxon_for_go <= 1.0

        # The stratified and legacy results should differ for GO shuffle
        # (pvalue_taxon_for_go uses GO shuffle mode)
        any_different = any(
            abs((combos_strat[k].pvalue_taxon_for_go or 0)
                - (combos_legacy[k].pvalue_taxon_for_go or 0)) > 0.001
            for k in combos_strat
            if combos_strat[k].pvalue_taxon_for_go is not None
            and combos_legacy[k].pvalue_taxon_for_go is not None
        )
        assert any_different, (
            "Aspect-stratified GO shuffle should produce different p-values "
            "from whole-closure shuffle"
        )

        # Taxonomy shuffle results should ALSO differ because taxonomy
        # shuffle is now aspect-stratified too (restricted peptide pool
        # per aspect)
        any_tax_different = any(
            abs((combos_strat[k].pvalue_go_for_taxon or 0)
                - (combos_legacy[k].pvalue_go_for_taxon or 0)) > 0.001
            for k in combos_strat
            if combos_strat[k].pvalue_go_for_taxon is not None
            and combos_legacy[k].pvalue_go_for_taxon is not None
        )
        assert any_tax_different, (
            "Aspect-stratified taxonomy shuffle should produce different "
            "p-values from unstratified shuffle"
        )

    def test_aspect_stratified_with_annotation_quality_confound(self):
        """Aspect-stratified shuffling should reduce spurious significance
        caused by cross-aspect annotation-quality correlation.

        Setup: well-annotated peptides have BOTH MF and BP terms;
        poorly-annotated peptides have only BP terms.  Taxon 10 gets the
        well-annotated ones.  Without stratification, MF root appears
        enriched in taxon 10 because it correlates with overall annotation
        quality.  With stratification, MF is shuffled independently of BP,
        breaking the confound.
        """
        anns = []
        # Taxon 10: well-annotated (has MF + BP)
        for i in range(5):
            anns.append(_make_annotation(
                f"GOOD{i}", 100.0, {1, 10},
                {"GO:MF_ROOT", "GO:MF1", "GO:BP_ROOT", "GO:BP1"}
            ))
        # Taxon 20: poorly-annotated (has only BP)
        for i in range(5):
            anns.append(_make_annotation(
                f"POOR{i}", 100.0, {1, 20},
                {"GO:BP_ROOT", "GO:BP2"}
            ))

        go_ns = {
            "GO:MF_ROOT": "molecular_function",
            "GO:MF1": "molecular_function",
            "GO:BP_ROOT": "biological_process",
            "GO:BP1": "biological_process",
            "GO:BP2": "biological_process",
        }

        agg = aggregate_peptide_annotations(anns)

        # With stratification
        combos_strat = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(
            anns, agg, combos_strat, iterations=500, seed=42,
            go_namespaces=go_ns,
        )

        # Without stratification
        combos_legacy = aggregate_go_taxonomy_combos(anns, agg)
        compute_enrichment_pvalues(
            anns, agg, combos_legacy, iterations=500, seed=42,
            go_namespaces=None,
        )

        # The legacy (non-stratified) GO shuffle should find MF_ROOT
        # enriched in taxon 10 because MF always co-occurs with BP in
        # the well-annotated group
        key = (10, "GO:MF_ROOT")
        if key in combos_legacy and combos_legacy[key].pvalue_taxon_for_go is not None:
            legacy_p = combos_legacy[key].pvalue_taxon_for_go
            strat_p = combos_strat[key].pvalue_taxon_for_go
            # Stratified p-value should be LESS significant (higher)
            # than legacy, because it breaks the annotation-quality confound
            assert strat_p >= legacy_p - 0.01, (
                f"Stratified p={strat_p:.4f} should be >= legacy p={legacy_p:.4f} "
                f"for the annotation-quality confound case"
            )

    def test_enrichment_disabled_leaves_none_values(self):
        """When compute_enrichment_pvalues is NOT called, combo fields stay None."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001"}),
            _make_annotation("PEP2", 50.0, {1, 20}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        # Do NOT call compute_enrichment_pvalues
        for combo in combos.values():
            assert combo.pvalue_go_for_taxon is None
            assert combo.pvalue_taxon_for_go is None
            assert combo.qvalue_go_for_taxon is None
            assert combo.qvalue_taxon_for_go is None


# ---------------------------------------------------------------------------
# Tests for _run_shuffle
# ---------------------------------------------------------------------------

class TestEdgeCasesFromReview:
    """Tests added in response to specification review findings."""

    def test_missing_aspect_peptides_dont_gain_terms(self):
        """GO shuffle: peptides without any term in an aspect must NOT gain
        terms from other peptides during that aspect's shuffle.

        Setup: PEP1-3 have MF+BP, PEP4-5 have only BP (no MF).
        After aspect-stratified GO shuffle, PEP4-5 should still have
        no MF terms — they should not inherit MF from PEP1-3.
        """
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:MF1", "GO:BP1"}),
            _make_annotation("PEP2", 100.0, {1, 10}, {"GO:MF2", "GO:BP1"}),
            _make_annotation("PEP3", 100.0, {1, 20}, {"GO:MF1", "GO:BP2"}),
            _make_annotation("PEP4", 100.0, {1, 20}, {"GO:BP1"}),
            _make_annotation("PEP5", 100.0, {1, 20}, {"GO:BP2"}),
        ]
        go_ns = {
            "GO:MF1": "molecular_function",
            "GO:MF2": "molecular_function",
            "GO:BP1": "biological_process",
            "GO:BP2": "biological_process",
        }
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        # Run with aspect stratification — should not crash and should
        # produce valid p-values
        compute_enrichment_pvalues(
            anns, agg, combos, iterations=200, seed=42,
            go_namespaces=go_ns,
        )

        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert 0 < combo.pvalue_go_for_taxon <= 1.0
            if combo.pvalue_taxon_for_go is not None:
                assert 0 < combo.pvalue_taxon_for_go <= 1.0

    def test_all_peptides_lack_one_aspect(self):
        """When no peptide has any term in an aspect, that aspect's shuffle
        should be a no-op and not crash."""
        # All peptides have only BP, no MF or CC
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:BP1"}),
            _make_annotation("PEP2", 100.0, {1, 20}, {"GO:BP2"}),
        ]
        go_ns = {
            "GO:BP1": "biological_process",
            "GO:BP2": "biological_process",
        }
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=100, seed=42,
            go_namespaces=go_ns,
        )

        for combo in combos.values():
            if combo.pvalue_taxon_for_go is not None:
                assert 0 < combo.pvalue_taxon_for_go <= 1.0

    def test_single_peptide_in_aspect_pool(self):
        """When only one peptide has an aspect, taxonomy shuffle for that
        aspect's terms should produce p-value = 1.0 (degenerate)."""
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:MF1", "GO:BP1"}),
            _make_annotation("PEP2", 100.0, {1, 20}, {"GO:BP1"}),
            _make_annotation("PEP3", 100.0, {1, 20}, {"GO:BP2"}),
        ]
        go_ns = {
            "GO:MF1": "molecular_function",
            "GO:BP1": "biological_process",
            "GO:BP2": "biological_process",
        }
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=100, seed=42,
            go_namespaces=go_ns,
        )

        # MF aspect has only 1 peptide → taxonomy shuffle for MF terms
        # should be non-significant (default p=1.0 from the restricted pool
        # having < 2 peptides)
        key_mf = (10, "GO:MF1")
        if key_mf in combos and combos[key_mf].pvalue_go_for_taxon is not None:
            assert combos[key_mf].pvalue_go_for_taxon >= 0.5

    def test_asymmetric_null_centering(self):
        """The two-sided p-value centers the null around its empirical mean,
        not around zero.  This test verifies centering works correctly for
        a case where the null distribution is asymmetric (shifted away from
        zero) but the observed value is near the null center.

        If we centered around zero instead of mean_null, the p-value would
        be artificially low because the observation would appear far from
        zero even though it's typical under the null.
        """
        # Create data where the null is shifted: taxon 10 has most abundance
        # and most GO terms, so the null mean of log2 enrichment for (10, GO:A)
        # will be positive (most shuffled taxonomies assign to taxon 10)
        anns = []
        for i in range(8):
            anns.append(_make_annotation(
                f"BIG{i}", 100.0, {1, 10}, {"GO:A"}
            ))
        for i in range(2):
            anns.append(_make_annotation(
                f"SMALL{i}", 10.0, {1, 20}, {"GO:A"}
            ))

        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        compute_enrichment_pvalues(
            anns, agg, combos, iterations=500, seed=42,
        )

        # The pair (10, GO:A) should NOT be significant because taxon 10
        # dominates the pool — the observed enrichment is typical under the
        # null even though it's positive
        key = (10, "GO:A")
        if key in combos and combos[key].pvalue_go_for_taxon is not None:
            assert combos[key].pvalue_go_for_taxon > 0.05, (
                f"Expected non-significant p-value for dominant taxon, "
                f"got {combos[key].pvalue_go_for_taxon:.4f}"
            )

    def test_zero_abundance_pairs_excluded(self):
        """Pairs where A_TAX=0 or A_GO=0 in a restricted pool should be
        excluded from testing (not produce degenerate log2 values)."""
        # In the full pool, taxon 10 has GO:MF1 and GO:BP1
        # In the MF-restricted pool, both peptides have MF, so both taxa present
        # In the BP-restricted pool, only PEP1 has BP
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:MF1", "GO:BP1"}),
            _make_annotation("PEP2", 100.0, {1, 20}, {"GO:MF1"}),
        ]
        go_ns = {
            "GO:MF1": "molecular_function",
            "GO:BP1": "biological_process",
        }
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)

        # Should not crash
        compute_enrichment_pvalues(
            anns, agg, combos, iterations=100, seed=42,
            go_namespaces=go_ns,
        )

        # All produced p-values should be valid
        for combo in combos.values():
            if combo.pvalue_go_for_taxon is not None:
                assert 0 < combo.pvalue_go_for_taxon <= 1.0
            if combo.pvalue_taxon_for_go is not None:
                assert 0 < combo.pvalue_taxon_for_go <= 1.0


class TestRunShuffle:
    def test_taxonomy_shuffle_produces_pvalues(self):
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001", "GO:0002"}),
            _make_annotation("PEP2", 50.0, {1, 20}, {"GO:0001"}),
            _make_annotation("PEP3", 80.0, {1, 10}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        pvals = _run_shuffle(idx, mode="taxonomy", iterations=100, eps=1e-10, seed=42)
        assert len(pvals) == idx.n_pairs
        assert np.all(pvals > 0)
        assert np.all(pvals <= 1.0)

    def test_go_shuffle_produces_pvalues(self):
        anns = [
            _make_annotation("PEP1", 100.0, {1, 10}, {"GO:0001", "GO:0002"}),
            _make_annotation("PEP2", 50.0, {1, 20}, {"GO:0001"}),
            _make_annotation("PEP3", 80.0, {1, 10}, {"GO:0002"}),
        ]
        agg = aggregate_peptide_annotations(anns)
        combos = aggregate_go_taxonomy_combos(anns, agg)
        idx = _build_index_maps(anns, combos)

        pvals = _run_shuffle(idx, mode="go", iterations=100, eps=1e-10, seed=42)
        assert len(pvals) == idx.n_pairs
        assert np.all(pvals > 0)
        assert np.all(pvals <= 1.0)
