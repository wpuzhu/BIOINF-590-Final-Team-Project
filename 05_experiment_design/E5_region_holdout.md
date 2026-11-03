# E5 · Region hold-out (spatial out-of-distribution transfer)

**Question.** Does the RNA → protein mapping learned in some tissue regions hold in a region never seen in training?

**Code.** `region_holdout(d, F)`; `uq.mondrian_conformal`, `uq.weighted_conformal_coverage`, `uq.shift_decomposition`

## Manipulation (input level only)
Train on all regions but one (across all training samples); test on the held-out region in held-out samples. Regions from pathology annotation, RNA-only domains (Novae/BANKSY on RNA), or protein niches. Report which, since protein niches make the split partly circular.

## Controls
In-distribution baseline: same regions, held-out samples. Covariate-shift reweighting.

## Read-out
Pearson r and MSE (R² is unstable when held-out variance is small); per-cell-type conformal coverage; shift decomposition (covariate part, concept part, target-specific gain).

## Uncertainty quantification
Bootstrap over held-out samples; conformal coverage with a binomial CI; domain-classifier AUC for the size of the shift.

## Interpretation
See the decision tree in `03_assumptions/why_information_does_not_transfer.md`. If weighting restores coverage → covariate shift (transferable given type/context). If not → concept shift (biological or technical).

## Failure modes
Unseen cell types in the held-out region (support shift); too few cells of a type for conformal (use ≥20 calibration cells).

## Simulation self-check
PD-L1 core coverage 0.72 → 0.95 after covariate weighting (covariate shift); pS6 0.71 → 0.85 (partly concept shift); lineage markers stay between 0.84 and 0.93.
