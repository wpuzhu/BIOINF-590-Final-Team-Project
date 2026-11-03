# U1 · Bootstrap

## Cluster (unit) bootstrap: the default
Resample **units** (patients, slides, blocks) with replacement; recompute the statistic; take percentiles. For a **contrast** between two input conditions, bootstrap the **per-unit paired differences**, not the two conditions separately. Pairing removes between-unit variance.

```python
from uq import fold_scores, unit_folds, cluster_bootstrap_ci
fs = fold_scores(d, make_X, unit_folds(d.sample), "cond")          # one score per (protein, unit)
est, lo, hi = cluster_bootstrap_ci(per_unit_diffs, B=2000)
```
- **Few units (< 8):** percentile CIs under-cover. Prefer BCa, or report the sign-flip p-value alongside, and say n.
- **Two-level designs** (blocks within patients): resample patients, then blocks within each resampled patient (hierarchical bootstrap).

## Spatial block bootstrap (single slide)
Tile the slide into blocks larger than the spatial correlation range (estimate it from a variogram or Moran's I correlogram of the residuals), and resample blocks. Too-small blocks under-estimate variance; too-large blocks leave few units. Report the CI across 2–3 block sizes.

## Bootstrap of curves and argmax (E3)
Pointwise percentile bands for R²(r); for r*, bootstrap units → recompute the curve → take the argmax → percentile CI on r*. A flat curve gives a wide r* CI: say "no preferred scale".

## Stability selection (N3, E2 gene lists)
Repeat selection (LASSO / forward selection) on 100 unit-bootstrap resamples; report each gene's **selection frequency**; keep genes selected in > 60–80% of resamples. This controls the expected number of false selections (Meinshausen & Bühlmann 2010).

## What not to do
- Bootstrapping cells inside one slide: autocorrelation makes the CIs far too narrow.
- Bootstrapping after model selection on all data: the selection leaks (A8). Bootstrap the whole pipeline.
