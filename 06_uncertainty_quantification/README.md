# 06 · Uncertainty quantification

**One rule above all:** the unit you resample or permute must be the unit of independence. That means patients or slides, or spatial blocks larger than the spatial correlation range. Never individual cells. Our simulation shows cell-level nulls can be **13× too narrow** (`08_results/simulation/uq/`).

| Page | Methods | Answers |
|---|---|---|
| [U1 · bootstrap](U1_bootstrap.md) | cluster bootstrap, spatial block bootstrap, stability selection | how precise is this gain / curve / gene list? |
| [U2 · permutation tests](U2_permutation_tests.md) | sign-flip, within-cell-type, toroidal shift, block permutation, random-gene-set null | is this gain / association real? |
| [U3 · conformal prediction](U3_conformal.md) | split, Mondrian (per cell type), weighted (covariate shift), spatial-block conformal | does uncertainty transfer to new regions/datasets? |
| [U4 · model uncertainty](U4_model_uncertainty.md) | seeds, deep ensembles, MC dropout, Bayesian ridge | how much of the result is the model's randomness? |
| [U5 · multiplicity and power](U5_multiplicity_power.md) | Holm, BH, hierarchical testing, knockoffs, power | how many claims can we make with n patients? |

## Which UQ goes with which experiment

| Experiment | Estimate | Interval / test |
|---|---|---|
| E1 variance ladder | ΔR² per step | cluster bootstrap CI on per-unit Δ; sign-flip p |
| E2 gene ablation | ΔR² / importance | bootstrap CI; empirical p vs random-gene-set null; BH |
| E3 radius sweep | R²(r), r* | pointwise bootstrap bands; bootstrap CI of argmax r* |
| E4 context controls | arrangement, composition contrasts | bootstrap CI; sign-flip p (average ≥5 shuffles per unit first) |
| E5 region hold-out | coverage, shift parts | Mondrian / weighted conformal; binomial CI on coverage; bootstrap the shift parts |
| E6 residual maps | Moran's I, residual–covariate r | within-type permutation; **toroidal shift** for associations |
| E7 cross-dataset | transfer ratio, rank agreement | bootstrap over target units; permutation of protein labels |
| E8 pairing | δ½ | bootstrap CI over units; repeat jitter draws |
| H1/H2 | macro-F1 differences | patient bootstrap, paired Wilcoxon/sign-flip over patients × seeds, Holm |

## Code

`07_code/ablation/uq.py`: `spatial_blocks`, `unit_folds`, `fold_scores`, `cluster_bootstrap_ci`, `signflip_test`, `group_permutation_importance`, `toroidal_shift_test`, `mondrian_conformal`, `weighted_conformal_coverage`, `shift_decomposition`, `holm`, `bh`. Demo: `run_uq_demo.py` → `08_results/simulation/uq/` (with `uq_overview.png`).
