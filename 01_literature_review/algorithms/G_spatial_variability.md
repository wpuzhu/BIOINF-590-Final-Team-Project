# G. Spatial-variability detectors (from v1, reused in v2)

The v1 project benchmarked cell-type-specific spatially variable gene (ctSVG) detectors for **proteins** (ctSVPs). The template notebooks are in [`../notebooks_v1/`](../notebooks_v1/), the runnable pipeline is in [`../../07_code/ctsvp_baseline/`](../../07_code/ctsvp_baseline/), and the full review is in [`../archive_v1/`](../archive_v1/).

## Why they matter in v2

| v2 use | What is tested | Detector input |
|---|---|---|
| **E6 / idea N2: residual ctSVPs** | is *protein minus RNA-predicted protein* still spatially variable within a cell type? | residual matrix (cells × proteins) + coordinates + cell types |
| **Idea N1: spatial-variability transfer** | is a gene ctSV *and* its protein ctSV in the same cell type? | gene matrix and protein matrix separately, same cells |
| **Calibration check** | p-values under within-type coordinate permutation | as above, permuted |

## Detectors and their inputs (summary)

| Method | Type | Input | Notes from v1 runs |
|---|---|---|---|
| **CELINA** | ctSVG, spatially varying coefficient, random effects, 11 kernels | counts + coords + cell-type matrix | best calibrated in v1 (type-I ≈ 0.02–0.04) |
| **STANCE** | ctSVG + SVG, LMM, two-stage | counts + coords + proportions | **badly anti-conservative** in v1 (type-I 0.33–0.47); calls marker genes ctSV |
| **SPARK-X** | global SVG, kernel covariance test | counts/intensities + coords | well calibrated; run within cell type for ctSV |
| **Within-type Moran's I** | autocorrelation on same-type edges | values + coords + types | in `07_code/ablation/xmodal_ablation.py::morans_i_within_type` with a within-type permutation null |
| **C-SIDE, CTSV, spVC, ctSVG** | spot- or cell-level ctSVG | see notebooks | spot-based tools need deconvolution weights |

## Statistical cautions (carry over to every spatial test in v2)

- **Composition confounding:** a global SVG test mostly rediscovers lineage markers (v1: Spearman ρ up to 0.55 between global Moran's I and cell-type η²). Always test **within** cell type.
- **Nulls must respect spatial structure** when testing *association* between two spatial fields (e.g. residual vs T-cell density). Use toroidal shifts or block permutations ([`../../06_uncertainty_quantification/`](../../06_uncertainty_quantification/)), not cell-level shuffles.
- **Residuals inherit the predictor's misspecification.** A linear predictor leaves non-linear RNA signal in the residual, which can look like "post-transcriptional" structure. Confirm residual ctSVPs with a non-linear predictor (DGAT) before interpreting them.
