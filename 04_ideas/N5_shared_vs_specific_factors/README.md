# N5 · Shared vs modality-specific latent factors

| | |
|---|---|
| **Status** | design only |
| **Question** | Which variance components are shared by RNA and protein, which are specific, and which are spatially smooth? |
| **Assumptions** | A1, A7 |
| **Data** | same-cell data |
| **Methods** | MOFA+ / MEFISTO (spatial smoothness), SpatialGlue attention as a descriptive check |
| **Experiments** | fit on training units; project held-out units; check factor stability |
| **UQ** | factor stability across bootstrap resamples of units (match factors by correlation; report the reproducibility rate) |
| **Cost** | minutes–hours |

## Motivation
A model-based complement to E1: instead of predicting protein, decompose both modalities jointly.

## Design
Variance explained per factor per modality; spatial length-scale per factor (MEFISTO).

## Expected outcomes and how to read them
Shared smooth factors = transferable tissue programmes; protein-specific factors = candidates for non-transferable biology.

## Risks
Factor rotation ambiguity; choose K by stability.
