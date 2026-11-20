# E6 · Residual spatial maps

**Question.** Which protein spatial structure remains after RNA (own + neighbourhood) is accounted for?

**Code.** `residual_map(...)` (within-type Moran's I); `uq.toroidal_shift_test` for associations

## Manipulation (input level only)
Residual = protein − **out-of-fold** prediction (cross-fitting). Then within-cell-type spatial autocorrelation, and association with candidate drivers (hypoxia score, vessel distance, immune density).

## Controls
Lineage markers as negative controls (residual should be flat); the predictor with all RNA shuffled (residual = protein).

## Read-out
Moran's I of protein vs residual per type; association r with covariates.

## Uncertainty quantification
Within-type permutation for Moran's I; **toroidal shift** (not cell shuffling) for associations between two spatial fields; BH across protein × type × covariate.

## Interpretation
Residual structure that is associated with a biological covariate is a candidate post-transcriptional or microenvironmental mechanism (B1–B5).

## Failure modes
Model misspecification leaves RNA signal in residuals, so confirm with a non-linear predictor.

## Simulation self-check
pS6 residual Moran's I = 0.84 (p = 0.005); the residual–y-coordinate association is significant under the toroidal null, while the toroidal null SD is 13× the naive one.
