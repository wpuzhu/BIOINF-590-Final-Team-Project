# N2 · Residual ctSVPs: protein spatial structure that RNA cannot explain

| | |
|---|---|
| **Status** | E6 implemented (Moran's I); full detectors via `ctsvp_baseline` |
| **Question** | After removing what RNA (own + neighbourhood) predicts, is the protein still spatially variable within a cell type, and what is it associated with? |
| **Assumptions** | A5, A7, and correct model specification |
| **Data** | same-cell data (P1/P2/P6) |
| **Methods** | cross-fitted RNA → protein predictor (ridge, then DGAT); ctSV detectors on residuals |
| **Experiments** | E6 |
| **UQ** | **cross-fitting** (residuals from out-of-fold predictions only); within-type permutation for Moran's I; **toroidal shift** for association with covariates (hypoxia score, distance to vessels/T cells); BH across proteins × types |
| **Cost** | minutes–hours |

## Motivation
The non-transferable part is not noise to throw away: it is where protein biology departs from transcription.

## Design
Residual = protein − out-of-fold prediction. Test spatial structure, then test association with candidate drivers using spatial nulls.

## Expected outcomes and how to read them
Phospho/stability proteins should show residual structure; lineage markers should not (negative controls).

## Risks
Misspecification leaves RNA signal in the residual; confirm with a non-linear predictor.
