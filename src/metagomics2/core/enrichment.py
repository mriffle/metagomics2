"""Single-sample GO x taxonomy enrichment statistics."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import erfc, isinf, sqrt

from metagomics2.core.annotation import PeptideAnnotation


DEFAULT_EXACT_ENUMERATION_THRESHOLD = 14
_EXACT_PVALUE_TOLERANCE = 1e-12


@dataclass
class ComboEnrichmentStats:
    """Enrichment statistics for a single (taxon, GO) pair."""

    pvalue_go_for_taxon: float | None = None
    pvalue_taxon_for_go: float | None = None
    qvalue_go_for_taxon: float | None = None
    qvalue_taxon_for_go: float | None = None
    zscore_go_for_taxon: float | None = None
    zscore_taxon_for_go: float | None = None


@dataclass
class GroupSummary:
    """Cached summary statistics for one tested group."""

    count: int = 0
    total_weight: float = 0.0
    sum_weight_sq: float = 0.0
    exact_peptides: tuple[PeptideAnnotation, ...] | None = None


def summarize_group(
    peptides: list[PeptideAnnotation],
    exact_enumeration_threshold: int = DEFAULT_EXACT_ENUMERATION_THRESHOLD,
) -> GroupSummary:
    """Summarize a peptide group once for reuse across many pair tests."""
    total_weight = 0.0
    sum_weight_sq = 0.0
    for peptide in peptides:
        total_weight += peptide.quantity
        sum_weight_sq += peptide.quantity * peptide.quantity

    exact_peptides = (
        tuple(peptides)
        if len(peptides) <= exact_enumeration_threshold
        else None
    )
    return GroupSummary(
        count=len(peptides),
        total_weight=total_weight,
        sum_weight_sq=sum_weight_sq,
        exact_peptides=exact_peptides,
    )


def filter_doubly_annotated_peptides(
    annotations: list[PeptideAnnotation],
) -> list[PeptideAnnotation]:
    """Keep only peptides eligible for enrichment testing."""
    return [
        ann
        for ann in annotations
        if ann.is_annotated and ann.taxonomy_nodes and ann.go_terms
    ]


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction."""
    if not pvalues:
        return []

    ranked = sorted(enumerate(pvalues), key=lambda item: (item[1], item[0]))
    n = len(ranked)
    adjusted = [0.0] * n

    for rank, (_, pvalue) in enumerate(ranked, start=1):
        adjusted[rank - 1] = min((pvalue * n) / rank, 1.0)

    running_min = 1.0
    for i in range(n - 1, -1, -1):
        running_min = min(running_min, adjusted[i])
        adjusted[i] = running_min

    restored = [0.0] * n
    for adjusted_value, (original_index, _) in zip(adjusted, ranked):
        restored[original_index] = adjusted_value
    return restored


def compute_signed_z_score_from_summary(
    total_weight: float,
    sum_weight_sq: float,
    observed_rate: float,
    background_rate: float,
) -> float | None:
    """Compute the signed z-score from precomputed weight moments."""
    if total_weight <= 0:
        return None

    variance_numerator = background_rate * (1.0 - background_rate) * sum_weight_sq
    if variance_numerator <= 0:
        return None

    variance = variance_numerator / (total_weight * total_weight)
    if variance <= 0:
        return None

    return (observed_rate - background_rate) / sqrt(variance)


def compute_boundary_rate_pvalue(
    observed_rate: float,
    background_rate: float,
) -> float | None:
    """Handle degenerate background rates without exact enumeration."""
    if background_rate <= 0.0:
        return 1.0 if observed_rate <= 0.0 else 0.0
    if background_rate >= 1.0:
        return 1.0 if observed_rate >= 1.0 else 0.0
    return None


def compute_boundary_zscore(
    observed_rate: float,
    background_rate: float,
) -> float | None:
    """Return a signed infinite z-score when the boundary direction is known."""
    if observed_rate > background_rate:
        return float("inf")
    if observed_rate < background_rate:
        return float("-inf")
    return None


def compute_exact_weighted_pvalue(
    weights: list[float],
    indicators: list[int],
    background_rate: float,
) -> float:
    """Compute an exact two-sided p-value by enumerating all assignments."""
    if not weights:
        return 1.0

    total_weight = sum(weights)
    if total_weight <= 0:
        return 1.0

    if background_rate <= 0.0:
        observed_rate = (
            sum(weight for weight, indicator in zip(weights, indicators) if indicator)
            / total_weight
        )
        return 1.0 if observed_rate <= 0.0 else 0.0

    if background_rate >= 1.0:
        observed_rate = (
            sum(weight for weight, indicator in zip(weights, indicators) if indicator)
            / total_weight
        )
        return 1.0 if observed_rate >= 1.0 else 0.0

    observed_rate = (
        sum(weight for weight, indicator in zip(weights, indicators) if indicator)
        / total_weight
    )
    observed_delta = abs(observed_rate - background_rate)

    n = len(weights)
    pvalue = 0.0
    for mask in range(1 << n):
        weighted_sum = 0.0
        probability = 1.0
        for index, weight in enumerate(weights):
            has_feature = (mask >> index) & 1
            if has_feature:
                weighted_sum += weight
                probability *= background_rate
            else:
                probability *= 1.0 - background_rate

            if probability == 0.0:
                break

        if probability == 0.0:
            continue

        rate = weighted_sum / total_weight
        if abs(rate - background_rate) >= observed_delta - _EXACT_PVALUE_TOLERANCE:
            pvalue += probability

    return min(max(pvalue, 0.0), 1.0)


def _resolve_pvalue_and_zscore(
    total_weight: float,
    sum_weight_sq: float,
    observed_rate: float,
    background_rate: float,
    exact_weights: list[float],
    exact_indicators: list[int],
    use_exact: bool,
) -> tuple[float, float | None]:
    """Select the boundary, exact, or normal-approximation result.

    This is the single decision path shared by the batch enrichment engine
    (via ``_test_pair_direction``) and the standalone
    ``compute_weighted_rate_test`` helper, so the tested logic and the logic the
    pipeline actually runs cannot drift apart.
    """
    zscore = compute_signed_z_score_from_summary(
        total_weight, sum_weight_sq, observed_rate, background_rate
    )

    boundary_pvalue = compute_boundary_rate_pvalue(observed_rate, background_rate)
    if boundary_pvalue is not None:
        return boundary_pvalue, compute_boundary_zscore(observed_rate, background_rate)

    if use_exact:
        pvalue = compute_exact_weighted_pvalue(
            exact_weights, exact_indicators, background_rate
        )
        return pvalue, zscore

    if zscore is None:
        # Interior background rate whose weighted variance underflowed to zero:
        # the normal approximation is undefined and the group is too large to
        # enumerate, so report a non-significant result with no effect size.
        return 1.0, None

    pvalue = min(max(erfc(abs(zscore) / sqrt(2.0)), 0.0), 1.0)
    return pvalue, zscore


def compute_weighted_rate_test(
    weights: list[float],
    indicators: list[int],
    background_rate: float,
    exact_enumeration_threshold: int = DEFAULT_EXACT_ENUMERATION_THRESHOLD,
) -> tuple[float, float | None]:
    """Test whether a weighted Bernoulli rate differs from background."""
    if not weights:
        return 1.0, None

    total_weight = sum(weights)
    if total_weight <= 0:
        return 1.0, None

    sum_weight_sq = sum(weight * weight for weight in weights)
    observed_rate = (
        sum(weight for weight, indicator in zip(weights, indicators) if indicator)
        / total_weight
    )
    return _resolve_pvalue_and_zscore(
        total_weight,
        sum_weight_sq,
        observed_rate,
        background_rate,
        weights,
        indicators,
        use_exact=len(weights) <= exact_enumeration_threshold,
    )


def _test_pair_direction(
    summary: GroupSummary,
    joint: float,
    other_group_total: float,
    background_denominator: float,
    feature_id: int | str,
    feature_collection: Callable[[PeptideAnnotation], set[int] | set[str]],
) -> tuple[float, float | None]:
    """Run one directional leave-one-out test for a single (taxon, GO) pair."""
    background_rate = min(
        max((other_group_total - joint) / background_denominator, 0.0),
        1.0,
    )
    observed_rate = joint / summary.total_weight

    if summary.exact_peptides is not None:
        exact_weights = [ann.quantity for ann in summary.exact_peptides]
        exact_indicators = [
            1 if feature_id in feature_collection(ann) else 0
            for ann in summary.exact_peptides
        ]
        use_exact = True
    else:
        exact_weights = []
        exact_indicators = []
        use_exact = False

    return _resolve_pvalue_and_zscore(
        summary.total_weight,
        summary.sum_weight_sq,
        observed_rate,
        background_rate,
        exact_weights,
        exact_indicators,
        use_exact=use_exact,
    )


def compute_go_taxonomy_enrichment(
    annotations: list[PeptideAnnotation],
    combo_keys: list[tuple[int, str]],
    exact_enumeration_threshold: int = DEFAULT_EXACT_ENUMERATION_THRESHOLD,
) -> dict[tuple[int, str], ComboEnrichmentStats]:
    """Compute enrichment statistics for observed GO x taxonomy pairs."""
    stats_by_pair = {
        pair: ComboEnrichmentStats() for pair in combo_keys
    }
    pool = filter_doubly_annotated_peptides(annotations)
    if not pool or not combo_keys:
        return stats_by_pair

    total_abundance = sum(ann.quantity for ann in pool)
    if total_abundance <= 0:
        return stats_by_pair

    tax_totals: dict[int, float] = {}
    go_totals: dict[str, float] = {}
    joint_totals: dict[tuple[int, str], float] = {}
    tax_groups: dict[int, list[PeptideAnnotation]] = {}
    go_groups: dict[str, list[PeptideAnnotation]] = {}

    for ann in pool:
        for tax_id in ann.taxonomy_nodes:
            tax_totals[tax_id] = tax_totals.get(tax_id, 0.0) + ann.quantity
            tax_groups.setdefault(tax_id, []).append(ann)

        for go_id in ann.go_terms:
            go_totals[go_id] = go_totals.get(go_id, 0.0) + ann.quantity
            go_groups.setdefault(go_id, []).append(ann)

        for tax_id in ann.taxonomy_nodes:
            for go_id in ann.go_terms:
                key = (tax_id, go_id)
                if key in stats_by_pair:
                    joint_totals[key] = joint_totals.get(key, 0.0) + ann.quantity

    tax_summaries = {
        tax_id: summarize_group(peptides, exact_enumeration_threshold)
        for tax_id, peptides in tax_groups.items()
    }
    go_summaries = {
        go_id: summarize_group(peptides, exact_enumeration_threshold)
        for go_id, peptides in go_groups.items()
    }

    go_for_taxon_pairs: list[tuple[tuple[int, str], float]] = []
    taxon_for_go_pairs: list[tuple[tuple[int, str], float]] = []

    for pair, stats in stats_by_pair.items():
        tax_id, go_id = pair
        joint = joint_totals.get(pair, 0.0)
        if joint <= 0.0:
            continue

        tax_total = tax_totals.get(tax_id, 0.0)
        go_total = go_totals.get(go_id, 0.0)

        tax_summary = tax_summaries.get(tax_id)
        tax_background_denominator = total_abundance - tax_total
        if tax_summary and tax_total > 0.0 and tax_background_denominator > 0.0:
            pvalue, zscore = _test_pair_direction(
                tax_summary,
                joint,
                go_total,
                tax_background_denominator,
                go_id,
                lambda ann: ann.go_terms,
            )
            stats.pvalue_go_for_taxon = pvalue
            stats.zscore_go_for_taxon = zscore
            go_for_taxon_pairs.append((pair, pvalue))

        go_summary = go_summaries.get(go_id)
        go_background_denominator = total_abundance - go_total
        if go_summary and go_total > 0.0 and go_background_denominator > 0.0:
            pvalue, zscore = _test_pair_direction(
                go_summary,
                joint,
                tax_total,
                go_background_denominator,
                tax_id,
                lambda ann: ann.taxonomy_nodes,
            )
            stats.pvalue_taxon_for_go = pvalue
            stats.zscore_taxon_for_go = zscore
            taxon_for_go_pairs.append((pair, pvalue))

    go_qvalues = benjamini_hochberg([pvalue for _, pvalue in go_for_taxon_pairs])
    for (pair, _), qvalue in zip(go_for_taxon_pairs, go_qvalues):
        stats_by_pair[pair].qvalue_go_for_taxon = qvalue

    taxon_qvalues = benjamini_hochberg([pvalue for _, pvalue in taxon_for_go_pairs])
    for (pair, _), qvalue in zip(taxon_for_go_pairs, taxon_qvalues):
        stats_by_pair[pair].qvalue_taxon_for_go = qvalue

    return stats_by_pair


def format_optional_stat(value: float | None) -> str:
    """Format an optional enrichment statistic for CSV output.

    Uses scientific notation so that very small p-values and q-values (which
    are exactly the most significant results) keep their magnitude instead of
    collapsing to a string of zeros.
    """
    if value is None:
        return ""
    if isinf(value):
        return "+inf" if value > 0 else "-inf"
    return f"{value:.6e}"
