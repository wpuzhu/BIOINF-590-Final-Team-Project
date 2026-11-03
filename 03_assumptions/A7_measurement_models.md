# A7. Measurement noise is modelled appropriately for each modality

**Statement.** RNA counts (sparse, Poisson/NB-like, dropout) and protein intensities (dense, continuous, with background, spillover and saturation) are each transformed so that a linear or Gaussian model is reasonable.

**Violations and fixes.**

| Issue | Modality | Fix |
|---|---|---|
| Library size, dropout | RNA | size-factor normalisation + log1p; PCA; optionally neighbourhood smoothing reported separately |
| Antibody background / non-specific binding | protein (ADT, imaging) | CLR (ADT), totalVI background model, or arcsinh with per-marker cofactor; isotype controls as negative "proteins" |
| **Lateral spillover** (signal leaks into adjacent cells) | protein imaging | spillover makes a neighbour's protein look like the cell's own. Check E4: a "neighbourhood RNA helps" effect for a lineage marker of the *neighbour's* type is a spillover warning. Use REDSEA-style compensation, or nucleus-restricted intensities |
| Batch/slide scale differences | both | z-score per sample before pooling; never fit scalers on test units |

**If violated.** Spurious information. Spillover is the most dangerous case, because it fakes exactly the neighbourhood effects we want to find.
