# Implementation Plan — Recalibrate Single-Sample Enrichment (FPC / weighted hypergeometric null)

Status: **implemented** on `single-file-p-values` (D1 replace-in-place and D4 collapse-by-default
both confirmed; the enrichment feature is not yet released to `main`, so this was a pre-release
change, not a migration of shipped behavior). The identical-annotation collapse ships as an
internal `collapse_identical` parameter (default on); no user-facing disable flag was added.

## 1. Motivation

The current enrichment method (`core/enrichment.py`) models each peptide in a tested group
as an independent weighted Bernoulli trial against a leave-one-out background rate. Empirical
evaluation showed this null is **anti-conservative**:

- Clean synthetic H0, large-fraction groups: ~2.7× too many hits at p<.05, ~3.7× at p<.01.
- **Real data** (`test-data/list_015`, reconstruction validated to exact agreement with the
  pipeline's own `taxonomy_nodes.csv`/`go_terms.csv`): under a real-structure permutation null,
  the current method produced **frac<.05 = .134 (2.7×)** and **frac<.01 = .076 (7.6×)**.
- Consequence on that sample: the current method calls 21,981 pairs significant at q<.05 vs
  11,206 for the recalibrated null; **~59% of its q<.05 hits do not survive** a calibrated test.
- It is also slow: the exact `O(2^N)` enumeration is the runtime hotspot (9.4 s on 390 peptides;
  the recalibrated analytic form ran the same comparison in 0.3 s, ~30×).

The replacement null is a **weighted finite-population (sampling-without-replacement)
permutation null**, equivalent to a weighted Fisher/hypergeometric test that conditions on both
margins. It is the exact null for "taxonomy ⟂ function given the observed margins," was validated
to match the honest whole-object permutation (closed-form mean/var matched permutation to 3–4
sig figs), is near-calibrated under real-structure H0 (frac<.05 = .054), retains full power
(planted signals detected 100%), and is dramatically cheaper.

See the prototype scripts (`scratchpad/enrichment_fpc.py`, `study.py`, `real_run.py`,
`real_calib.py`) for the validation harness that produced these numbers.

## 2. The method

For each observed `(taxon t, GO g)` pair in the doubly-annotated pool of size `P`, two
directional tests (unchanged from today):

- **GO-for-taxon**: tested group = peptides in taxon `t` (weights `w_i`); the "labels" are the
  `n_g` peptides carrying `g`, placed uniformly without replacement among all `P` slots.
- **taxon-for-GO**: symmetric — tested group = peptides carrying `g`; labels are the `n_t`
  peptides in taxon `t`.

The test statistic is the joint abundance `A_JOINT(t,g)` (equivalently the weighted rate).
Under the null its moments are closed-form:

```
mean = (n / P) * W_S
var  = [ n (P - n) / (P^2 (P - 1)) ] * ( P * Σ_S w_i^2  -  W_S^2 )
```

where `W_S = Σ_S w_i` is the tested group's total weight, `Σ_S w_i^2` its sum of squares, and
`n` is the *count* of labels (n_g or n_t). Two-sided p-value:

- **Large groups**: normal tail, `p = erfc(|z|/√2)`, `z = (A_JOINT - mean)/√var`.
- **Small groups**: exact enumeration of the group's `2^m` weighted subsets. Each subset of size
  `h` contributes probability `C(P-m, n-h)/C(P,n)` (the specific members don't matter, only the
  size), so the exact two-sided p is a scan over cached `(size, sum)` pairs. This enumeration is
  done **once per group** and reused across every pair that shares that group — unlike the current
  code, which re-enumerates `2^N` per pair.

The leave-one-out background disappears (the conditioning handles self-inclusion, like Fisher's
exact test). `z` remains a signed effect size (>0 enrichment). Degenerate cases (`n ∈ {0, P}`) are
deterministic: `p = 1`, `z` empty — there are no ±inf boundary cases anymore.

## 3. Design decisions

| # | Decision | Recommendation | Rationale |
|---|----------|----------------|-----------|
| D1 | Replace current null, or keep it as a selectable legacy option? | **Replace in place** | Current null is anti-conservative (a correctness issue), and the feature is unreleased. Keeping two nulls doubles the test surface for no user benefit. |
| D2 | Keep the CSV/column contract? | **Yes, unchanged** | Same six columns (`pvalue/qvalue/zscore_*`), same scientific-notation format, same `compute_enrichment_pvalues` flag. Frontend and any consumers need no changes. Only the numeric values change. |
| D3 | Small-group accuracy strategy | **Exact enumeration (cached per group) + analytic for large; no Monte Carlo in the shipped path** | Deterministic (preserves reproducibility/provenance, a project value). MC stays a test-only validation tool. |
| D4 | Collapse pseudo-replicate peptides before testing? | **Default on for identical `(taxonomy, go)` signatures**, flag to disable | Peptides with identical annotations are statistically indistinguishable to the test; collapsing (summing weights) is unambiguously correct and is the dominant calibration lever on rich samples. Protein-level collapse is deferred (D-future). |
| D5 | Residual extreme-tail error on large groups | **Implemented: double-saddlepoint (Skovgaard) tail** | Single-saddlepoint was found inadequate for large m/P; the double-saddlepoint conditions on both margins and matches the exact tail to a few percent across all m/P. Residual mild `p<.01` anti-conservatism (~2-3x, vs ~8x for the old method) now traces to the two-sided *definition* (distance-from-mean), not the tail approximation. |

D1, D4 are the two that change scientific output most — please confirm before implementation.

## 4. Performance

- Preprocessing: `O(P · t · g)` — reuse the totals/groups already built by
  `aggregate_go_taxonomy_combos` instead of recomputing them (removes today's duplicate pass).
- Per pair: `O(1)` analytic for large groups; for small groups `O(m)` per pair after a one-time
  `O(2^m)` per-group enumeration that is **cached and reused across all pairs sharing the group**.
- BH: `O(K log K)` (unchanged).

Net: strictly faster than the current per-pair `2^N` (measured ~30× on real data), while being
correct. Memory: process small groups group-major so only one group's `2^m` table is resident at
a time.

## 5. File-by-file changes

### Backend
- **`src/metagomics2/core/enrichment.py`** (core rewrite):
  - Add `fpc_moments(group_total, group_sq, P, n)` and `fpc_analytic_pvalue(...)`.
  - Add exact small-group null: per-group cached subset enumeration + two-sided tail query,
    keyed by group; reuse `GroupSummary.exact_peptides` (already caches small-group weights).
  - Rewrite `_test_pair_direction` / the two direction blocks in `compute_go_taxonomy_enrichment`
    to use FPC. Keep the existing scaffolding (`tax_groups`, `go_groups`, `summarize_group`,
    `joint_totals`) — it already provides `W`, `Σw²`, group weights, and counts (`len(group)`).
  - Restructure the small-group path to be group-major (enumerate once per group, apply to all its
    pairs) to realize the caching win and bound memory.
  - Remove now-dead Bernoulli-LOO helpers: `compute_boundary_rate_pvalue`,
    `compute_exact_weighted_pvalue`, `compute_weighted_rate_test`, `_resolve_pvalue_and_zscore`.
    Keep `benjamini_hochberg`, `format_optional_stat`, `ComboEnrichmentStats`, `GroupSummary`.
  - Add `collapse_identical_annotations(pool)` (D4), applied at the top of
    `compute_go_taxonomy_enrichment` when enabled.
- **`src/metagomics2/pipeline/runner.py`**: no change to wiring; still calls
  `compute_go_taxonomy_enrichment(annotations, combo_keys)`. Manifest still records the flag.
  (If D4 exposes a disable flag, thread it like `compute_enrichment_pvalues`.)
- **`reporting.py`, `models/job.py`, `cli.py`, `worker/worker.py`**: unchanged (contract + flag
  unchanged), unless D4 adds a collapse toggle — then thread it through the same four places.

### Frontend
- **No code changes.** Same columns/parsing. Verify `tsc --noEmit` + `vitest` still pass; the GO
  DAG q-value coloring and tooltips work unchanged (values are just better calibrated).

### Docs
- **`docs/SPECIFICATION_SINGLE_FILE_ENRICHMENT.md`**: rewrite §5 (statistical tests),
  §6 (p-values/z), §8 (performance), §14 (limitations) to describe the FPC null, exact/analytic
  split, per-group caching, and the identical-annotation collapse. Replace the leave-one-out
  framing with both-margins conditioning; add a short "Calibration evidence" subsection.
- **`docs/SPECIFICATION_BACKEND.md`**: update §5 Stage 6 and §7.6 to match.
- Delete this plan doc (or mark it "implemented") once merged.

## 6. Testing plan

- **`tests/unit/test_enrichment.py`** (rewrite the null-specific tests; keep BH/eligibility/format):
  - `fpc_moments` vs brute-force permutation moments (small deterministic case).
  - Exact small-group p vs brute-force enumeration of the hypergeometric placement.
  - Analytic vs exact agreement in an overlap regime (medium group) within tolerance.
  - Directional eligibility (both directions), degenerate `n ∈ {0, P}` → `p=1`, `z` empty.
  - `collapse_identical_annotations` correctness (weights summed, signatures deduped).
- **`tests/property/`**: add a Hypothesis test that `fpc_moments` matches empirical
  permutation moments over random weighted pools (fits the project's property-testing philosophy).
- **Regression**: `test_reporting.py`, `test_worker.py`, `test_server_api.py`, integration tests
  assert the column contract / flag plumbing, not method values — they should pass unchanged.
- **Calibration acceptance** (checked-in as a `@pytest.mark.slow` test or a `scripts/` harness):
  under a permutation H0 on synthetic data, require FPC `frac<.05 ∈ [.04, .06]` and its `.01`
  tail materially below the current method's.

## 7. Validation / acceptance criteria (definition of done)

1. All existing + new unit/property tests pass (`./venv/bin/python -m pytest tests/`).
2. Frontend `tsc --noEmit` and `vitest run` pass in Docker.
3. Re-running `real_run.py`/`real_calib.py` on `test-data/`:
   - pool reconstruction still validates exactly against the pipeline CSVs;
   - FPC calibration under real-structure H0: `frac<.05 ≈ .05`, `.01` tail ≪ current;
   - a summary of how the significant-set changes vs the current method is recorded.
4. A full pipeline run with `--enrichment-pvalues` produces a well-formed `go_taxonomy_combo.csv`
   (six columns, scientific notation, empties where ineligible) and the manifest records the flag.
5. Docs updated to describe the shipped method.

## 8. Rollout sequence

1. **Core** — implement FPC (analytic + exact-cached), rewrite unit/property tests, validate
   against permutation MC. (D1, D3)
2. **Collapse** — add identical-annotation collapse (D4), re-validate calibration on real data.
3. **Docs + frontend verification** — update specs; confirm no frontend change needed.
4. **Future (separate PR)** — saddlepoint tail for large-group extreme tails (D5); optional
   protein-level collapse; expose method selection only if a legacy null is ever required.

## 9. Risks & backward compatibility

- **Numeric values change** — this is the intended effect. Because the feature is unreleased,
  no released behavior breaks. Anyone comparing to prior branch runs will see fewer, better-
  calibrated hits; call this out in the docs and PR description.
- **CSV/format/flag unchanged** → frontend, `comboParser`, and any external consumers are
  unaffected; old combo CSVs still parse.
- **Determinism preserved** (no MC in shipped path) → provenance/reproducibility intact.
- **Small residual anti-conservatism at the extreme tail** for large groups (normal
  approximation) — documented as a known limitation; saddlepoint deferred.
- **Interpretation** — outputs remain within-sample, single-replicate enrichment; keep the
  "scores, not population inference" framing in the docs (no biological replication is modeled).
