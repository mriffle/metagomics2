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

Only observed `(taxon, GO)` pairs are tested. A pair must have positive joint abundance in the filtered pool to receive enrichment statistics.

Each test direction also requires a positive leave-one-out denominator:

- GO-for-taxon requires `A_total - A_TAX(t) > 0`
- taxon-for-GO requires `A_total - A_GO(g) > 0`

If a direction is ineligible, its p-value, q-value, and z-score fields are left empty in the CSV.

---

## 4. Observed Quantities

For the filtered peptide pool:

| Symbol | Definition |
|--------|-----------|
| `A_total` | Total abundance across all eligible peptides |
| `A_TAX(t)` | Total abundance of peptides whose taxonomy closure contains `t` |
| `A_GO(g)` | Total abundance of peptides whose GO closure contains `g` |
| `A_JOINT(t, g)` | Total abundance of peptides containing both `t` and `g` |

These quantities are computed over propagated closures, so parent and child nodes in either hierarchy share signal and therefore produce correlated results. That correlation is expected.

---

## 5. Statistical Tests

### 5.1 GO for Taxon

This asks whether GO term `g` is unusually concentrated within taxon `t`.

Observed rate:

```text
p_obs = A_JOINT(t, g) / A_TAX(t)
```

Leave-one-out background rate:

```text
p_bg = (A_GO(g) - A_JOINT(t, g)) / (A_total - A_TAX(t))
```

### 5.2 Taxon for GO

This asks whether taxon `t` claims an unusual share of GO term `g`.

Observed rate:

```text
p_obs = A_JOINT(t, g) / A_GO(g)
```

Leave-one-out background rate:

```text
p_bg = (A_TAX(t) - A_JOINT(t, g)) / (A_total - A_GO(g))
```

### 5.3 Why Leave-One-Out

The implementation uses leave-one-out backgrounds so the tested slice is compared against "the rest of the sample" rather than against a background partly containing itself. This avoids self-dilution for dominant taxa or broad GO groups.

---

## 6. P-Values and Z-Scores

### 6.1 Exact vs Approximate Computation

The implementation uses two paths:

- Exact weighted enumeration for groups with `N <= 14`
- Normal approximation for groups with `N > 14`

The threshold is implemented as:

```python
DEFAULT_EXACT_ENUMERATION_THRESHOLD = 14
```

This threshold is not user-configurable in the UI or CLI.

### 6.2 Exact Weighted P-Value

For small groups, the code enumerates all `2^N` Bernoulli assignments for the tested group and computes an exact two-sided p-value using peptide abundances as weights.

For a group with weights `w_i` and binary indicators `x_i`:

```text
p_obs = sum(w_i * x_i) / sum(w_i)
```

The exact p-value is the probability, under the null Bernoulli rate `p_bg`, of observing a weighted rate at least as far from `p_bg` as the observed one.

### 6.3 Normal Approximation

For larger groups, the weighted rate is approximated as normal:

```text
Var[S] = p_bg * (1 - p_bg) * sum(w_i^2) / (sum(w_i)^2)
z = (p_obs - p_bg) / sqrt(Var[S])
```

The two-sided p-value is computed as:

```text
erfc(|z| / sqrt(2))
```

### 6.4 Signed Z-Score Semantics

The z-score is a signed directional effect measure:

- `z > 0` means enrichment
- `z < 0` means depletion

For large groups, the z-score is the same value used for the approximate p-value.

For small groups, the p-value comes from exact enumeration, but the z-score is still reported as a descriptive effect size when the variance is defined.

### 6.5 Boundary Cases

When `p_bg` is exactly `0` or `1`, the usual variance term is zero and a finite z-score does not exist.

The shipped implementation handles these cases explicitly:

- if direction is enriched, z-score is `+inf`
- if direction is depleted, z-score is `-inf`
- if `p_obs == p_bg`, z-score remains empty

The CSV formatter writes these values as `+inf` and `-inf`.

The p-value behavior at the boundary is:

- `pvalue = 0` when the observed rate contradicts a boundary null
- `pvalue = 1` when the observed rate matches the boundary null

So a row can legitimately contain `pvalue = 0` and `zscore = +inf` or `-inf`.

---

## 7. Multiple Testing Correction

Benjamini-Hochberg FDR correction is applied separately for the two test directions:

- `qvalue_go_for_taxon` is computed across all tested GO-for-taxon p-values
- `qvalue_taxon_for_go` is computed across all tested taxon-for-GO p-values

The implementation restores q-values to the original pair order after ranking.

---

## 8. Performance Characteristics

The current implementation is optimized for the large-group path.

### 8.1 Cached Group Summaries

For each taxon group and GO group, the engine caches:

- total group weight
- sum of squared weights
- the exact peptide tuple only when the group size is at or below the exact threshold

This means the normal-approximation path no longer rebuilds per-pair weight vectors for large groups.

### 8.2 Current Complexity

Approximate high-level behavior:

- preprocessing over the peptide pool: proportional to propagated taxonomy and GO memberships
- large-group pair testing: approximately `O(1)` per tested pair
- small-group exact testing: `O(2^N)` per tested pair, but only for groups with `N <= 14`

In practice, the exact branch is still the dominant hotspot on workloads containing many small tested groups with many observed pairs.

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
- finite values: fixed decimal string
- infinite z-scores: `+inf` or `-inf`

This behavior is part of the backward-compatible CSV contract used by the frontend parser.

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

- The method assumes independent Bernoulli behavior within a tested group; peptides from the same protein can violate that assumption.
- Parent and child nodes in GO and taxonomy remain highly correlated because closures are propagated.
- Very broad GO terms and very high-level taxonomy nodes can produce biologically uninformative but statistically valid results.
- Exact enumeration can still be slow on pathological workloads with many small tested groups and many observed pairs, though lowering the threshold to `14` substantially reduces this cost.
