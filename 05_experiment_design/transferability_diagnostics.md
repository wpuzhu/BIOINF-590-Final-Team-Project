# Informative vs transferable, and why things fail to transfer

## Definitions used in this project

- **Informative** (for protein *p*, cell type *k*, input set *S*): adding *S* to the inputs raises held-out R² for *p* in *k*, on **held-out units drawn from the same population** (same tissue, same platform, same regions). Measured by E1–E4. The CI must exclude 0.
- **Transferable:** the gain from *S* persists when the held-out units come from a **different population**: a new region (E5), a new dataset/platform/tissue (E7), or pairing with error (E8). Measured as the transfer gap between the in-distribution and out-of-distribution gain, with a CI.
- **Non-transferable but informative:** informative in-distribution, gap ≫ 0 out of distribution. This is the interesting category. Diagnose it with the steps below.

## Diagnostic battery

| Step | Tool (code) | Question |
|---|---|---|
| 1. How different is the target? | domain classifier AUC on context features: cell type + neighbourhood composition (`uq._density_ratio`) | size of covariate shift; AUC ≈ 1 for a subset = support shift |
| 2. Is there support? | per-cell-type counts in source vs target; nearest-neighbour distance of target cells to the source cloud | extrapolation flag |
| 3. Does coverage hold? | Mondrian conformal per cell type (`uq.mondrian_conformal`) | distribution-free transfer check |
| 4. Does reweighting fix it? | weighted conformal (`uq.weighted_conformal_coverage`); shift decomposition (`uq.shift_decomposition`) | covariate shift (fixable) vs concept shift |
| 5. Could a target-specific model do better? | oracle refit in target (`shift_decomposition.mse_target_oracle`) | is the failure in the *mapping* (non-transferable) or the *protein* (unpredictable anywhere)? |
| 6. What is in the residual? | E6 residual Moran's I + toroidal association tests | biological mechanism (B1–B8) vs noise |

## Reading the shift decomposition (MSE units)

| Pattern | Meaning | Example in our simulation |
|---|---|---|
| covariate part large, concept part ≈ 0, oracle ≈ transfer | transferable once you condition on cell type/context | — |
| covariate part large, oracle ≪ transfer | the context mix changed **and** the model class can't extrapolate; a region-specific model does better → "transferable only with region-aware models" | PD-L1 into the T-cell-poor core (covariate 0.73, target-specific gain 0.37) |
| oracle ≈ transfer ≈ target variance | protein unpredictable from RNA anywhere; the shift parts only reflect variance differences between regions, not a mapping failure | pS6 (oracle 1.59 ≈ transfer 1.58 ≈ target variance 1.58) |
| all parts ≈ 0 | transfers | CD3, Ki-67, CD68 |

## Why some regions are non-transferable (summary)

1. **Composition outside the training range:** e.g. tumour core vs margin. Covariate/support shift.
2. **Regional regulation of translation or stability:** hypoxia, necrosis, nutrient gradients. Concept shift with a biological covariate (B1, B2).
3. **Region-specific signalling:** tertiary lymphoid structures, invasive front. Context the model never saw (B8).
4. **Technical regional effects:** edge artefacts, tissue folding, uneven staining, registration error that grows towards section edges (T3–T5).

The full mechanism table and decision tree are in [`../03_assumptions/why_information_does_not_transfer.md`](../03_assumptions/why_information_does_not_transfer.md).
