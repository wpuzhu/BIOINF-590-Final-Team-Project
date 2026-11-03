# N1 · Spatial-variability transfer: do RNA ctSVGs become protein ctSVPs?

| | |
|---|---|
| **Status** | ready to run with v1 detectors (`07_code/ctsvp_baseline/`) |
| **Question** | For a gene/protein pair within one cell type, is spatial variability concordant between modalities? |
| **Assumptions** | A1, A2, A5 |
| **Data** | P1/P2/P6 (same cell), D13 (spot) |
| **Methods** | CELINA, SPARK-X, within-type Moran's I ([G_spatial_variability](../../01_literature_review/algorithms/G_spatial_variability.md)) |
| **Experiments** | run detectors on cognate gene and protein separately in each cell type; 2×2 concordance table |
| **UQ** | within-type permutation calibration of each detector; bootstrap over spatial blocks for the concordance (Cohen's κ with CI); BH across pairs |
| **Cost** | hours (R methods) |

## Motivation
Prediction accuracy says whether levels transfer. This asks whether **spatial patterns** transfer, which is the more biological question.

## Design
Four boxes: concordant SV (transfers); RNA-only SV (buffered at protein level); protein-only SV (post-transcriptional spatial regulation); neither.

## Expected outcomes and how to read them
A table of proteins per box, per cell type. Protein-only SV pairs are prime candidates for B1/B2/B4 mechanisms ([why_information_does_not_transfer](../../03_assumptions/why_information_does_not_transfer.md)).

## Risks
Different detection power in each modality (RNA sparser), so compare effect sizes, not only calls; detector miscalibration (STANCE).
