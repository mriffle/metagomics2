# Metagomics 2 — Single-Sample GO x Taxonomy Enrichment

## 1. Purpose

This document describes the single-sample GO x taxonomy enrichment feature as it is actually implemented in Metagomics 2.

The feature adds within-sample enrichment statistics to `go_taxonomy_combo.csv` so users can ask two directional questions for each observed `(taxon, GO term)` pair:

1. Is GO term `g` enriched or depleted within taxon `t` relative to the rest of the sample?
2. Is taxon `t` enriched or depleted within GO term `g` relative to the rest of the sample?

This is a within-sample weighted enrichment analysis. It is not a between-sample differential test.

---

## 2. Scope

The analysis operates on a single peptide list after peptide annotation is complete.

It uses:

- peptide abundance as a weight
- taxonomy lineage closure from LCA to root
- GO closure union from direct terms to ancestors

It does not model uncertainty in:

- peptide identification
- protein inference
- GO annotation quality
- taxonomy assignment quality
- biological replication

---

## 3. Inputs and Eligibility

The enrichment engine consumes `PeptideAnnotation` objects from the annotation pipeline.

Relevant fields:

| Field | Type | Meaning |
|-------|------|---------|
| `peptide` | `str` | Peptide sequence |
| `quantity` | `float` | Peptide abundance |
| `is_annotated` | `bool` | Whether the peptide received annotation |
| `taxonomy_nodes` | `set[int]` | Taxonomy lineage from LCA to root |
| `go_terms` | `set[str]` | GO closure union |

Only peptides satisfying all of the following are included in enrichment:

- `is_annotated == True`
- `taxonomy_nodes` is non-empty
- `go_terms` is non-empty

This means enrichment uses the doubly annotated peptide pool, not all sample abundance.

### 3.1 Collapsing pseudo-replicates

Before testing, peptides that share an identical `(taxonomy_nodes, go_terms)` annotation are collapsed into a single weighted unit whose quantity is the sum of the merged peptides. Such peptides are statistically indistinguishable to the test — they co-occur in, and are co-absent from, every `(taxon, GO)` pair — so counting them as separate independent observations would inflate significance. The collapsed pool defines the population size `P` used by the null. Collapsing is on by default (`collapse_identical=True`).

### 3.2 Eligibility

Only observed `(taxon, GO)` pairs are tested. A pair must have positive joint abundance in the collapsed pool to receive enrichment statistics.

A test direction is ineligible when its tested group is the entire pool (there is no complement to compare against):

- GO-for-taxon is ineligible when taxon `t` is present in every unit
- taxon-for-GO is ineligible when GO term `g` is present in every unit

If a direction is ineligible, its p-value, q-value, and z-score fields are left empty in the CSV.

---

## 4. Observed Quantities

For the filtered peptide pool:

| Symbol | Definition |
|--------|-----------|
| `P` | Number of units in the collapsed pool |
| `A_total` | Total abundance across all eligible units |
| `A_TAX(t)` | Total abundance of units whose taxonomy closure contains `t` |
| `A_GO(g)` | Total abundance of units whose GO closure contains `g` |
| `A_JOINT(t, g)` | Total abundance of units containing both `t` and `g` |
| `n_TAX(t)` | Number of units whose taxonomy closure contains `t` |
| `n_GO(g)` | Number of units whose GO closure contains `g` |

These quantities are computed over propagated closures, so parent and child nodes in either hierarchy share signal and therefore produce correlated results. That correlation is expected.

---

## 5. Statistical Tests

The null hypothesis is that taxonomy annotation and GO annotation are independent across units, conditional on the observed margins. It is realized as a weighted finite-population (sampling-without-replacement) permutation: the labels of one axis are placed uniformly at random among the `P` units, holding the other axis and the weights fixed. This is equivalent to a weighted Fisher/hypergeometric test and conditions on both margins, so no leave-one-out background is required.

### 5.1 GO for Taxon

Holds the taxon-`t` group of units (with their weights) fixed and places the `n_GO(g)` "carries `g`" labels uniformly without replacement among all `P` units. The statistic is the group's labelled weight:

```text
A_JOINT(t, g) = Σ_{u ∈ TAX(t)} w_u · 1[u carries g]
```

compared against its null distribution. Large `A_JOINT` relative to the null means `g` is enriched within `t`.

### 5.2 Taxon for GO

Symmetric: holds the GO-`g` group fixed and places the `n_TAX(t)` "in taxon `t`" labels. The statistic is the same `A_JOINT(t, g)`, compared against the null defined by the GO-`g` group and `n_TAX(t)`.

### 5.3 Null Moments

Under placement of `n` labels among `P` units, the tested group's labelled weight has closed-form moments:

```text
mean = (n / P) · W_S
var  = [ n (P − n) / (P² (P − 1)) ] · ( P · Σ_S w_u² − W_S² )
```

where `W_S` is the group's total weight and `Σ_S w_u²` its sum of squared weights. The variance carries the finite-population correction that the earlier independent-Bernoulli model lacked, and conditioning on both margins removes the need for a leave-one-out background.

---

## 6. P-Values and Z-Scores

### 6.1 Exact vs Approximate Computation

Two paths, selected by the size `m` of the tested group:

- Exact enumeration for groups with `m <= 14`
- Normal approximation for groups with `m > 14`

The threshold is implemented as:

```python
EXACT_GROUP_SIZE_MAX = 14
```

This threshold is not user-configurable in the UI or CLI.

### 6.2 Exact Two-Sided P-Value

For small groups the code enumerates the group's `2^m` weighted subsets once (grouped by subset size and sorted) and caches the result per group. A specific size-`h` subset is the labelled set with probability `C(P − m, n − h) / C(P, n)` (only its size matters), so the two-sided p-value is the total probability mass of subsets whose sum is at least as far from the null mean as the observed `A_JOINT`. Because the enumeration is cached per group, every `(taxon, GO)` pair that shares the group reuses it rather than re-enumerating.

### 6.3 Large-Group Tail (Double Saddlepoint)

For larger groups the two-sided p-value comes from a **double-saddlepoint (Skovgaard) approximation** to the finite-population null. It conditions on both margins, so it captures the skewness a normal approximation misses and stays accurate across all group/pool ratios (validated against the exact enumeration of §6.2 to within a few percent in the tail). The signed z-score (§6.4) is still taken from the closed-form moments:

```text
z = (A_JOINT − mean) / sqrt(var)
```

Two regimes defer to the plain normal tail `erfc(|z| / sqrt(2))`:

- **near the mean** (`|z| < 0.5`), where the normal is accurate and the saddlepoint formula is singular; and
- **at the support edge** — a "pure" group whose entire weight carries (or lacks) the feature — an atom the continuous approximation mishandles but which is hugely significant regardless.

### 6.4 Signed Z-Score Semantics

The z-score is a signed directional effect measure:

- `z > 0` means enrichment
- `z < 0` means depletion

It is computed from the closed-form moments and reported for both paths (exact and approximate) as a descriptive effect size whenever the variance is defined.

### 6.5 Degenerate Cases

When a feature is absent (`n = 0`) or universal (`n = P`) the statistic is deterministic: the variance is zero, so the p-value is `1` and the z-score is left empty. This null does not produce infinite z-scores.

---

## 7. Multiple Testing Correction

Benjamini-Hochberg FDR correction is applied separately for the two test directions:

- `qvalue_go_for_taxon` is computed across all tested GO-for-taxon p-values
- `qvalue_taxon_for_go` is computed across all tested taxon-for-GO p-values

The implementation restores q-values to the original pair order after ranking.

---

## 8. Performance Characteristics

### 8.1 Cached Group Summaries

For each taxon group and GO group, the engine caches total group weight, sum of squared weights, and the group's weights. Small groups (at or below the exact threshold) additionally get a lazily-built, sorted subset-sum table.

### 8.2 Complexity

- preprocessing over the collapsed pool: proportional to propagated taxonomy and GO memberships
- large-group pair testing: `O(1)` near the mean (normal); `O(m)` per tail in the tail region (the double-saddlepoint's Newton solve)
- small-group pair testing: `O(2^m)` **once per group** (cached), then `O(m)` per pair that shares the group

The exact enumeration is amortized per group rather than repeated per pair, and most large-group pairs fall near the mean (the fast normal path), so the method is much faster than a per-pair exact approach.

---

## 9. Output Columns

`go_taxonomy_combo.csv` always includes these six enrichment columns:

| Column | Meaning |
|--------|---------|
| `pvalue_go_for_taxon` | Raw p-value for GO-within-taxon |
| `pvalue_taxon_for_go` | Raw p-value for taxon-within-GO |
| `qvalue_go_for_taxon` | BH-adjusted q-value for GO-within-taxon |
| `qvalue_taxon_for_go` | BH-adjusted q-value for taxon-within-GO |
| `zscore_go_for_taxon` | Signed z-score for GO-within-taxon |
| `zscore_taxon_for_go` | Signed z-score for taxon-within-GO |

Formatting rules:

- disabled enrichment: empty strings
- ineligible direction: empty strings
- finite values: scientific-notation string (e.g. `1.234560e-02`), which preserves the magnitude of very small p-values and q-values

The finite-population null does not produce infinite z-scores; for backward compatibility the formatter still writes `+inf` / `-inf` if an infinite value is ever encountered. This behavior is part of the CSV contract used by the frontend parser.

---

## 10. User-Facing Configuration

### 10.1 Web UI

`NewJobPage.tsx` includes a checkbox:

- `Calculate enrichment p-values`

When checked, the frontend submits:

```json
{
  "compute_enrichment_pvalues": true
}
```

### 10.2 CLI

The CLI exposes:

```text
--enrichment-pvalues
```

### 10.3 Pipeline and Job Models

The flag is threaded through:

- `PipelineConfig.compute_enrichment_pvalues`
- `JobParams.compute_enrichment_pvalues`
- worker job execution
- run manifest output

The exact enumeration threshold is not exposed as a user-facing parameter.

---

## 11. Frontend Integration

### 11.1 Combo CSV Parsing

The frontend combo parser:

- parses old CSVs without enrichment columns
- parses enriched CSVs with all six enrichment fields
- parses `+inf` and `-inf` z-score strings into numeric `Infinity` and `-Infinity`

### 11.2 GO DAG Page

When a taxonomy filter is active and combo enrichment data is present:

- the GO DAG page can color by `Q-value (GO for Taxon)`
- q-values are visualized as `-log10(q + eps)`
- the GO DAG page can filter by a **max q-value (GO for taxon)** threshold. This iteratively prunes leaf GO terms whose `qvalue_go_for_taxon` exceeds the threshold until every remaining leaf passes. Only leaves are pruned: internal GO terms are preserved even when their own q-value is worse than the threshold, so the DAG spine from significant leaves up to the root(s) stays intact. Because q-value eligibility for the go-for-taxon direction is an all-or-nothing property of the selected taxon (a q-value is empty only when the taxon is present in every unit, which also hides this control), every candidate leaf has a q-value whenever the filter is usable.
- GO tooltips show:
  - fraction of taxon
  - fraction of GO
  - q-value
  - z-score

### 11.3 Taxonomy Page

When a GO filter is active and combo enrichment data is present:

- taxonomy tooltips show:
  - fraction of taxon
  - fraction of GO
  - q-value for taxon-within-GO
  - z-score for taxon-within-GO

The taxonomy page does not currently add a q-value coloring metric analogous to the GO DAG page.

---

## 12. Implementation Locations

Primary backend code:

- `src/metagomics2/core/enrichment.py`
- `src/metagomics2/core/aggregation.py`
- `src/metagomics2/core/reporting.py`
- `src/metagomics2/pipeline/runner.py`

User-facing plumbing:

- `src/metagomics2/cli.py`
- `src/metagomics2/models/job.py`
- `src/metagomics2/worker/worker.py`
- `frontend/src/pages/NewJobPage.tsx`
- `frontend/src/utils/comboParser.ts`
- `frontend/src/pages/GoDagPage.tsx`
- `frontend/src/components/GoDagViewer.tsx`
- `frontend/src/components/TaxonomyChart.tsx`

---

## 13. Test Coverage

The implemented feature is covered by:

- backend unit tests for enrichment math and edge cases
- reporting tests for CSV column formatting, including `+inf` and `-inf`
- API, worker, and pipeline integration tests for flag plumbing and output generation
- frontend combo parser tests for backward compatibility and infinity parsing

Backend verification used during implementation:

```bash
./venv/bin/python -m pytest tests/unit/test_enrichment.py \
  tests/unit/test_reporting.py \
  tests/unit/test_worker.py \
  tests/unit/test_server_api.py \
  tests/integration/test_pipeline_mocked_homology.py \
  tests/integration/test_pipeline_cli_end_to_end.py
```

Frontend verification follows the project frontend spec and runs in Docker:

```bash
docker build --target frontend-builder -t metagomics2-frontend-test .
docker run --rm metagomics2-frontend-test npx vitest run
docker run --rm metagomics2-frontend-test npx tsc --noEmit
```

---

## 14. Known Limitations

- The null treats units as exchangeable. Collapsing identical annotations removes exact pseudo-replication, but distinct peptides with correlated (non-identical) annotations from homologous proteins are still treated as independent.
- Parent and child nodes in GO and taxonomy remain highly correlated because closures are propagated.
- Very broad GO terms and very high-level taxonomy nodes can produce biologically uninformative but statistically valid results.
- The two-sided p-value uses distance from the null mean. For strongly skewed groups this is mildly anti-conservative in the extreme tail — empirically ~2-3x nominal at `p < 0.01` under a real-structure null, versus ~8x for the previous method — while the `p < 0.05` level is well calibrated. This residual is a property of the two-sided definition, not the tail approximation (the double-saddlepoint tracks the exact tail to a few percent, and the "twice the smaller tail" alternative is worse). Large groups use the double-saddlepoint tail; small groups use the exact enumeration.
- The analysis is within-sample and single-replicate; it does not model biological variability. Outputs are best read as enrichment scores rather than population-level inference.
