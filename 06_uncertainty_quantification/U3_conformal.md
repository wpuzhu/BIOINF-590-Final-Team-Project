# U3 · Conformal prediction: uncertainty that tells you whether information transferred

## Split conformal (Lei et al. 2018; Angelopoulos & Bates 2021 tutorial)
Train on units A; compute absolute residuals on calibration units B; q = the ⌈(n+1)(1−α)⌉-th smallest residual. The interval ŷ ± q covers new points with probability ≥ 1−α **if calibration and test points are exchangeable**.

## Mondrian (class-conditional) conformal: per cell type
Compute q separately in each cell type, which gives coverage **within** each type (`uq.mondrian_conformal`). Essential here because RNA explains protein very differently across types.

## Coverage as a transferability test
Calibrate on held-out units **from the source** population, then measure empirical coverage on the **target** (new region, dataset, platform):
- coverage ≈ 1−α: the residual distribution transferred (the mapping did too, at least in its error profile);
- coverage < 1−α: shift. Its size per cell type localises *where* transfer fails.

Report coverage with a **binomial CI** (n_test cells in type k; but cells are dependent, so treat the CI as optimistic, or compute coverage per block and bootstrap blocks).

## Weighted conformal under covariate shift (Tibshirani et al., NeurIPS 2019)
Re-weight calibration residuals by w(x) = p_target(x) / p_source(x), estimated with a classifier (`uq.weighted_conformal_coverage`).
- **Estimate w on low-dimensional context features** (own cell type + neighbourhood composition), not on 100+ RNA dimensions; high-dimensional density ratios are unstable (we observed this in development).
- If weighting **restores** coverage, the failure was covariate shift: transferable, conditional on context. If not, it was concept shift.
- In simulation: PD-L1 core coverage went from 0.72 to 0.95 after weighting (covariate shift); pS6 went from 0.71 to 0.85 (only partly fixed).

## Spatial variants
- **Block conformal:** calibration units are spatial blocks disjoint from training blocks (exchangeability between blocks, not cells).
- **GeoConformal** (Lou, Luo & Meng, *Ann Am Assoc Geogr* 2025): geographically weighted calibration that gives *spatially varying* interval widths. Use it to map where predictions are uncertain on a slide.
- **Conformalised quantile regression** (Romano et al. 2019): adaptive widths when heteroscedastic (e.g. protein variance grows with expression).

## For classification (H1/H2)
Adaptive prediction sets (APS) per cell type; set size is an uncertainty measure for cell typing on new slides. Report coverage and mean set size per arm (E vs F vs G).
