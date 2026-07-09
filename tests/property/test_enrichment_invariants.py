"""Property-based tests for the finite-population enrichment null.

Verifies the closed-form moments and the analytic-vs-exact machinery against
exhaustive enumeration of the exact weighted-hypergeometric placement.
"""

from itertools import combinations
from math import comb

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from metagomics2.core.enrichment import (
    _exact_fpc_pvalue,
    _subset_sums_by_size,
    fpc_moments,
)

weights_strategy = st.lists(
    st.floats(min_value=0.01, max_value=100.0, allow_nan=False, allow_infinity=False),
    min_size=1,
    max_size=6,
)


def _brute(weights, P, n):
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


def _brute_p(weights, P, n, observed):
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
            if abs(sum(weights[i] for i in subset) - mean) >= delta - 1e-9:
                p += prob
    return min(max(p, 0.0), 1.0)


@given(
    weights=weights_strategy,
    extra=st.integers(min_value=1, max_value=34),
    n_frac=st.floats(min_value=0.05, max_value=0.95),
)
@settings(max_examples=200, deadline=None)
def test_closed_form_moments_match_enumeration(weights, extra, n_frac):
    m = len(weights)
    P = m + extra
    n = max(1, min(P - 1, round(n_frac * P)))
    mean, var = fpc_moments(sum(weights), sum(w * w for w in weights), P, n)
    b_mean, b_var = _brute(weights, P, n)
    assert mean == pytest.approx(b_mean, rel=1e-9, abs=1e-9)
    assert var == pytest.approx(b_var, rel=1e-6, abs=1e-6)


@given(
    weights=weights_strategy,
    extra=st.integers(min_value=1, max_value=34),
    n_frac=st.floats(min_value=0.05, max_value=0.95),
    obs_frac=st.floats(min_value=0.0, max_value=1.0),
)
@settings(max_examples=200, deadline=None)
def test_exact_pvalue_matches_enumeration_and_is_bounded(weights, extra, n_frac, obs_frac):
    m = len(weights)
    P = m + extra
    n = max(1, min(P - 1, round(n_frac * P)))
    observed = obs_frac * sum(weights)
    table = _subset_sums_by_size(tuple(weights))
    p = _exact_fpc_pvalue(table, sum(weights), P, n, observed)
    assert 0.0 <= p <= 1.0
    assert p == pytest.approx(_brute_p(weights, P, n, observed), abs=1e-9)
