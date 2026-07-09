"""Single-sample GO x taxonomy enrichment statistics.

The null is a weighted finite-population (sampling-without-replacement) permutation
null, equivalent to a weighted Fisher/hypergeometric test that conditions on both
margins. For an observed ``(taxon t, GO g)`` pair the two directional tests are:

- GO-for-taxon: hold the taxon-``t`` peptide group (weights) fixed and place the
  ``n_g`` "carries-``g``" labels uniformly without replacement among all ``P``
  peptides; ask whether the group's labelled weight (``A_JOINT``) is extreme.
- taxon-for-GO: symmetric, with the GO-``g`` group fixed and the ``n_t`` "in-``t``"
  labels placed.

Under this null the statistic has a closed-form mean/variance, so large groups use a
normal-tail p-value while small groups use an exact enumeration of the group's
weighted subsets (computed once per group and reused across every pair that shares
it). Identical-annotation peptides are collapsed to one weighted unit first, so
statistically-redundant pseudo-replicates are not double counted.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from dataclasses import dataclass
from math import copysign, erfc, exp, isinf, lgamma, log, log1p, pi, sqrt, tanh

from metagomics2.core.annotation import PeptideAnnotation

# Group sizes at or below this use exact enumeration; larger groups use the
# double-saddlepoint tail approximation. Not user-configurable.
EXACT_GROUP_SIZE_MAX = 14

# Below this |z| the observed value is near the null mean, where the normal
# approximation is accurate, the Lugannani-Rice formula is singular, and the
# p-value is far from significant; fall back to the normal there.
_SADDLEPOINT_NEAR_MEAN_Z = 0.5

_SQRT2 = sqrt(2.0)
_SQRT2PI = sqrt(2.0 * pi)


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
    """Cached summary statistics for one tested group.

    ``exact`` marks groups small enough for exact enumeration; larger groups use
    the saddlepoint tail. ``weights`` is retained for both paths.
    """

    count: int = 0
    total_weight: float = 0.0
    sum_weight_sq: float = 0.0
    weights: tuple[float, ...] = ()
    exact: bool = False


def filter_doubly_annotated_peptides(
    annotations: list[PeptideAnnotation],
) -> list[PeptideAnnotation]:
    """Keep only peptides eligible for enrichment testing."""
    return [
        ann
        for ann in annotations
        if ann.is_annotated and ann.taxonomy_nodes and ann.go_terms
    ]


def collapse_identical_annotations(
    pool: list[PeptideAnnotation],
) -> list[PeptideAnnotation]:
    """Collapse peptides sharing an identical (taxonomy, GO) annotation into one
    weighted unit (summed quantity).

    Peptides with identical annotations are statistically indistinguishable to the
    enrichment test -- they always co-occur or co-absent for every ``(t, g)`` pair --
    so treating them as separate independent observations inflates significance.
    Collapsing them yields the correct independent-unit count for the null.
    """
    aggregated: dict[tuple[frozenset, frozenset], float] = defaultdict(float)
    representative: dict[tuple[frozenset, frozenset], PeptideAnnotation] = {}
    for ann in pool:
        key = (frozenset(ann.taxonomy_nodes), frozenset(ann.go_terms))
        aggregated[key] += ann.quantity
        representative.setdefault(key, ann)

    collapsed: list[PeptideAnnotation] = []
    for key, weight in aggregated.items():
        rep = representative[key]
        collapsed.append(
            PeptideAnnotation(
                peptide=rep.peptide,
                quantity=weight,
                is_annotated=True,
                taxonomy_nodes=set(rep.taxonomy_nodes),
                go_terms=set(rep.go_terms),
            )
        )
    return collapsed


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


# --------------------------------------------------------------------------- #
# Finite-population null: moments (closed form) and exact enumeration
# --------------------------------------------------------------------------- #
def fpc_moments(
    group_total: float,
    group_sum_sq: float,
    pool_size: int,
    label_count: int,
) -> tuple[float, float]:
    """Mean and variance of the group's labelled weight under the null.

    ``label_count`` labels are placed uniformly without replacement among
    ``pool_size`` peptides; the returned moments are for the summed weight of the
    tested group's peptides that get labelled.
    """
    f = label_count / pool_size
    mean = f * group_total
    if pool_size <= 1:
        return mean, 0.0
    variance = (
        label_count
        * (pool_size - label_count)
        / (pool_size * pool_size * (pool_size - 1))
    ) * (pool_size * group_sum_sq - group_total * group_total)
    return mean, max(variance, 0.0)


def _log_comb(n: int, k: int) -> float:
    """Log of the binomial coefficient C(n, k); -inf when out of range."""
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def _subset_sums_by_size(weights: tuple[float, ...]) -> list[list[float]]:
    """Sorted subset sums grouped by subset size, via O(2^m) DP over the weights."""
    sums_by_size: list[list[float]] = [[0.0]]
    for weight in weights:
        extended: list[list[float]] = [[] for _ in range(len(sums_by_size) + 1)]
        for size, sums in enumerate(sums_by_size):
            extended[size].extend(sums)  # exclude this weight
            extended[size + 1].extend(value + weight for value in sums)  # include it
        sums_by_size = extended
    for sums in sums_by_size:
        sums.sort()
    return sums_by_size


def _exact_fpc_pvalue(
    sums_by_size: list[list[float]],
    group_total: float,
    pool_size: int,
    label_count: int,
    observed: float,
) -> float:
    """Exact two-sided p-value under the finite-population null for a small group.

    A specific size-``h`` subset of the group is the labelled set with probability
    ``C(P-m, n-h) / C(P, n)`` (only the size matters), so the tail probability is a
    binary-search scan over the cached subset sums.
    """
    group_size = len(sums_by_size) - 1
    mean = (label_count / pool_size) * group_total
    delta = abs(observed - mean)
    tol = 1e-9 * (abs(mean) + delta + 1.0)
    low = mean - delta + tol
    high = mean + delta - tol

    log_denominator = _log_comb(pool_size, label_count)
    pvalue = 0.0
    for size, sums in enumerate(sums_by_size):
        log_weight = _log_comb(pool_size - group_size, label_count - size)
        if isinf(log_weight):
            continue
        per_subset_prob = exp(log_weight - log_denominator)
        if per_subset_prob == 0.0:
            continue
        n_low = bisect.bisect_right(sums, low)
        n_high = len(sums) - bisect.bisect_left(sums, high)
        in_tail = n_low + n_high
        if in_tail:
            pvalue += per_subset_prob * in_tail
    return min(max(pvalue, 0.0), 1.0)


def _std_normal_sf(x: float) -> float:
    """Upper tail of the standard normal, 1 - Phi(x)."""
    return 0.5 * erfc(x / _SQRT2)


def _std_normal_pdf(x: float) -> float:
    return exp(-0.5 * x * x) / _SQRT2PI


def _dsp_moments(
    t: float, theta: float, weights: tuple[float, ...], non_group: int
) -> tuple[float, float, float, float, float, float]:
    """Cumulants of the joint (group weight ``S``, total label count) under the
    base Bernoulli(1/2) model tilted by ``t`` on ``S`` and ``theta`` on the count.

    Returns ``(K, K_t, K_theta, K_tt, K_ttheta, K_thetatheta)``. The constant
    ``-P*log2`` is omitted because it cancels in the Skovgaard difference.
    """
    cgf = k_t = k_th = k_tt = k_tth = k_thth = 0.0
    for w in weights:
        a = t * w + theta
        s = 0.5 * (1.0 + tanh(0.5 * a))                            # sigmoid(a)
        sp = s * (1.0 - s)
        cgf += (a if a > 0.0 else 0.0) + log1p(exp(-abs(a)))       # softplus(a)
        k_t += w * s
        k_th += s
        k_tt += w * w * sp
        k_tth += w * sp
        k_thth += sp
    s0 = 0.5 * (1.0 + tanh(0.5 * theta))
    sp0 = s0 * (1.0 - s0)
    cgf += non_group * ((theta if theta > 0.0 else 0.0) + log1p(exp(-abs(theta))))
    k_th += non_group * s0
    k_thth += non_group * sp0
    return cgf, k_t, k_th, k_tt, k_tth, k_thth


def _dsp_upper(
    target: float, weights: tuple[float, ...], pool_size: int, label_count: int
) -> float | None:
    """Double-saddlepoint (Skovgaard) approximation to P(S >= target | total = n).

    Conditions on both margins, so it is accurate across all group/pool ratios.
    Returns None if the 2-D saddlepoint solve fails to converge (caller falls back
    to the normal).
    """
    non_group = pool_size - len(weights)
    theta0 = log(label_count / (pool_size - label_count))
    k0 = pool_size * ((theta0 if theta0 > 0.0 else 0.0) + log1p(exp(-abs(theta0))))
    k_thth0 = label_count * (pool_size - label_count) / pool_size

    t, theta = 0.0, theta0
    converged = False
    for _ in range(100):
        _, k_t, k_th, k_tt, k_tth, k_thth = _dsp_moments(t, theta, weights, non_group)
        grad_t = target - k_t
        grad_th = label_count - k_th
        tol_t = 1e-11 * (abs(target) + 1.0)
        tol_th = 1e-11 * (label_count + 1.0)
        if abs(grad_t) <= tol_t and abs(grad_th) <= tol_th:
            converged = True
            break
        det = k_tt * k_thth - k_tth * k_tth
        if det <= 0.0:
            break
        t += (k_thth * grad_t - k_tth * grad_th) / det
        theta += (-k_tth * grad_t + k_tt * grad_th) / det
    if not converged:
        return None

    cgf, _, _, k_tt, k_tth, k_thth = _dsp_moments(t, theta, weights, non_group)
    det = k_tt * k_thth - k_tth * k_tth
    if det <= 0.0:
        return None
    inside = 2.0 * ((t * target + theta * label_count - cgf) - (theta0 * label_count - k0))
    w_root = copysign(sqrt(max(inside, 0.0)), t)
    u = t * sqrt(max(det / k_thth0, 1e-300))
    if w_root == 0.0 or u == 0.0:
        return 0.5
    value = _std_normal_sf(w_root) + _std_normal_pdf(w_root) * (1.0 / u - 1.0 / w_root)
    return min(max(value, 0.0), 1.0)


def _saddlepoint_two_sided(
    weights: tuple[float, ...],
    total_weight: float,
    pool_size: int,
    label_count: int,
    observed: float,
    mean: float,
    variance: float,
) -> float:
    """Two-sided (distance-from-mean) p-value via the double saddlepoint.

    Conditions on both margins (the finite-population null), so it is accurate
    across all group/pool ratios and captures the skewness the normal misses. It
    defers to the normal near the mean (normal accurate, Lugannani-Rice singular)
    and at the support edge (a "pure" group is an atom the continuous approximation
    mishandles, but is hugely significant regardless).
    """
    if variance <= 0.0:
        return 1.0
    delta = abs(observed - mean)
    z = delta / sqrt(variance)
    if z < _SADDLEPOINT_NEAR_MEAN_Z:
        return min(max(erfc(z / _SQRT2), 0.0), 1.0)

    boundary = 1e-9 * (total_weight + 1.0)
    if observed >= total_weight - boundary or observed <= boundary:
        return min(max(erfc(z / _SQRT2), 0.0), 1.0)

    high = mean + delta
    low = mean - delta
    above_support = high >= total_weight - boundary
    upper = 0.0 if above_support else _dsp_upper(high, weights, pool_size, label_count)
    lower_sf = 1.0 if low <= boundary else _dsp_upper(low, weights, pool_size, label_count)
    if upper is None or lower_sf is None:
        return min(max(erfc(z / _SQRT2), 0.0), 1.0)
    return min(max(upper + (1.0 - lower_sf), 0.0), 1.0)


def _test_direction(
    summary: GroupSummary,
    cache_key: tuple[str, object],
    subset_cache: dict[tuple[str, object], list[list[float]]],
    pool_size: int,
    label_count: int,
    observed: float,
) -> tuple[float, float | None]:
    """Run one directional finite-population test; returns (p-value, signed z)."""
    mean, variance = fpc_moments(
        summary.total_weight, summary.sum_weight_sq, pool_size, label_count
    )
    zscore = (observed - mean) / sqrt(variance) if variance > 0.0 else None

    if label_count <= 0 or label_count >= pool_size:
        # Feature absent or universal: the statistic is deterministic.
        return 1.0, zscore

    if summary.exact:
        table = subset_cache.get(cache_key)
        if table is None:
            table = _subset_sums_by_size(summary.weights)
            subset_cache[cache_key] = table
        pvalue = _exact_fpc_pvalue(
            table, summary.total_weight, pool_size, label_count, observed
        )
    else:
        pvalue = _saddlepoint_two_sided(
            summary.weights, summary.total_weight, pool_size, label_count, observed, mean, variance
        )
    return pvalue, zscore


def _summarize_group(
    weights: list[float],
    indices: list[int],
    exact_group_size_max: int,
) -> GroupSummary:
    total = 0.0
    sum_sq = 0.0
    for i in indices:
        w = weights[i]
        total += w
        sum_sq += w * w
    m = len(indices)
    return GroupSummary(
        count=m,
        total_weight=total,
        sum_weight_sq=sum_sq,
        weights=tuple(weights[i] for i in indices),
        exact=m <= exact_group_size_max,
    )


def compute_go_taxonomy_enrichment(
    annotations: list[PeptideAnnotation],
    combo_keys: list[tuple[int, str]],
    *,
    collapse_identical: bool = True,
    exact_group_size_max: int = EXACT_GROUP_SIZE_MAX,
) -> dict[tuple[int, str], ComboEnrichmentStats]:
    """Compute finite-population enrichment statistics for observed GO x taxonomy pairs."""
    stats_by_pair = {pair: ComboEnrichmentStats() for pair in combo_keys}

    pool = filter_doubly_annotated_peptides(annotations)
    if collapse_identical:
        pool = collapse_identical_annotations(pool)
    if not pool or not combo_keys:
        return stats_by_pair

    pool_size = len(pool)
    weights = [ann.quantity for ann in pool]

    tax_groups: dict[int, list[int]] = defaultdict(list)
    go_groups: dict[str, list[int]] = defaultdict(list)
    joint_totals: dict[tuple[int, str], float] = defaultdict(float)
    for i, ann in enumerate(pool):
        for tax_id in ann.taxonomy_nodes:
            tax_groups[tax_id].append(i)
        for go_id in ann.go_terms:
            go_groups[go_id].append(i)
        for tax_id in ann.taxonomy_nodes:
            for go_id in ann.go_terms:
                key = (tax_id, go_id)
                if key in stats_by_pair:
                    joint_totals[key] += ann.quantity

    tax_summaries = {
        tax_id: _summarize_group(weights, indices, exact_group_size_max)
        for tax_id, indices in tax_groups.items()
    }
    go_summaries = {
        go_id: _summarize_group(weights, indices, exact_group_size_max)
        for go_id, indices in go_groups.items()
    }
    tax_count = {tax_id: len(indices) for tax_id, indices in tax_groups.items()}
    go_count = {go_id: len(indices) for go_id, indices in go_groups.items()}

    subset_cache: dict[tuple[str, object], list[list[float]]] = {}
    go_for_taxon_pairs: list[tuple[tuple[int, str], float]] = []
    taxon_for_go_pairs: list[tuple[tuple[int, str], float]] = []

    for pair, stats in stats_by_pair.items():
        tax_id, go_id = pair
        joint = joint_totals.get(pair, 0.0)
        if joint <= 0.0:
            continue

        # GO-for-taxon: taxon group fixed, place the g-carrier labels.
        tax_summary = tax_summaries.get(tax_id)
        if (
            tax_summary is not None
            and tax_summary.total_weight > 0.0
            and tax_summary.count < pool_size
        ):
            pvalue, zscore = _test_direction(
                tax_summary, ("tax", tax_id), subset_cache,
                pool_size, go_count[go_id], joint,
            )
            stats.pvalue_go_for_taxon = pvalue
            stats.zscore_go_for_taxon = zscore
            go_for_taxon_pairs.append((pair, pvalue))

        # taxon-for-GO: GO group fixed, place the taxon-t member labels.
        go_summary = go_summaries.get(go_id)
        if (
            go_summary is not None
            and go_summary.total_weight > 0.0
            and go_summary.count < pool_size
        ):
            pvalue, zscore = _test_direction(
                go_summary, ("go", go_id), subset_cache,
                pool_size, tax_count[tax_id], joint,
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

    Uses scientific notation so that very small p-values and q-values (which are
    exactly the most significant results) keep their magnitude instead of
    collapsing to a string of zeros.
    """
    if value is None:
        return ""
    if isinf(value):
        return "+inf" if value > 0 else "-inf"
    return f"{value:.6e}"
