# E7 · Cross-dataset / cross-platform transfer

**Question.** Is the information map tissue- and platform-general?

**Code.** `fold_scores` with folds = datasets; shared-feature restriction helper (to add)

## Manipulation (input level only)
Fit on dataset A (e.g. Xenium Protein RCC), predict dataset B (DBiTplus lymph node, spatial-CITE-seq tonsil, CosMx breast) on the **shared** genes × proteins; harmonise protein scales per dataset (z-score) and compare rank metrics.

## Controls
Within-B cross-validated model (oracle); A → A held-out (in-distribution).

## Read-out
Transfer ratio = R²(A→B) / R²(B→B); Spearman per protein; information-map rank agreement between datasets (Kendall τ of per-protein gains).

## Uncertainty quantification
Bootstrap over B's units; permutation of protein labels for a rank-agreement null.

## Interpretation
High transfer for lineage proteins and low for state proteins is the expected pattern; deviations point to platform effects (T4) or tissue-specific regulation.

## Failure modes
Resolution differences (spot vs cell) and panel mismatch dominate. Run at matched resolution (aggregate cells to pseudo-spots).

## Simulation self-check
Not simulated (needs real datasets).
