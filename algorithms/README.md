# Algorithms

> **v2 note.** These notebooks are from v1 (ctSVG/ctSVP biomarker detection). In v2 they still serve three roles: ctSV detectors for residual maps (E6) and RNA-vs-protein spatial-variability transfer (N1); protein-native baselines (mistyR, Statial, scimap, CytoCommunity); and protein foundation models (KRONOS, VirTues). Cross-modal RNA ↔ protein methods (SpatialGlue, MultiGATE, MaxFuse, DGAT, Nicheformer, Spatium, HEIST, ...) are listed in [`../docs/03_methods_catalog.md`](../docs/03_methods_catalog.md).

Every notebook is a **template** (not executed) that shows installation, the input format, and the core calls.
R notebooks are `.Rmd` (open in RStudio; chunks have `eval = FALSE` until you install the packages). Python notebooks are `.ipynb`.

| Folder | Notebook | Methods |
|---|---|---|
| `01_ctSVG_methods/` | `ctSVG_methods_tutorial.Rmd` | CELINA, STANCE, C-SIDE, CTSV, spVC, MMM (MCube), ctSVG |
| `02_global_SVG_methods/` | `global_SVG_methods_R.Rmd` | SPARK-X, nnSVG |
|  | `global_SVG_methods_python.ipynb` | Moran's I (Squidpy), within-cell-type Moran's I, SpatialDE, Hotspot |
| `03_spatial_proteomics_biomarker_methods/` | `spatial_proteomics_biomarkers_R.Rmd` | Statial (SpatioMark, Kontextual), spicyR, lisaClust, mistyR, imcRtools |
|  | `scimap_cytocommunity_python.ipynb` | scimap neighbourhoods/interactions, CytoCommunity |
| `04_pretrained_foundation_models/` | `pretrained_models.ipynb` | KRONOS, VirTues, Mesmer, Cellpose |
| `05_evaluation_framework/` | `metrics.py`, `evaluation_and_simulation.ipynb` | EP, AUPRC, FDP, power, type I error, KS; simulator; permutation/rotation tests |

**Common output format.** Every method should write a long table `method, cell_type, feature, p_value` (use `-1` for "not tested") so that `05_evaluation_framework/metrics.py` can score all of them the same way. `metrics.py` has been self-tested. The simulator and metrics were smoke-tested together on synthetic coordinates.

---

## 1. Cell-type-specific spatially variable (ctSV) methods

These are the methods in the anchor benchmark. Install sources match the pinned `renv.lock` of the official benchmark repo [weiqi-grace-li/ctsvg-benchmark](https://github.com/weiqi-grace-li/ctsvg-benchmark).

| Method | Model | Input | Strengths (benchmark) | Weaknesses (benchmark) | Install |
|---|---|---|---|---|---|
| **CELINA** | spatially varying coefficient, random effects, 11 kernels, Cauchy combination | counts, coords, cell-type matrix, scRNA ref (for QC) | lowest FDP (41%) in realistic sims; robust for minority types | QC filter removes true ctSVGs; scalability | `devtools::install_github("pekjoonwu/CELINA")` |
| **STANCE** | LMM, single kernel, two-stage (SVG → ctSVG) | counts, coords, proportions | rotation-invariant; robust for rare types | calls marker genes ctSVGs; FDP 98%; slow growth | `devtools::install_github("Cui-STT-Lab/STANCE")` · [tutorial](https://haroldsu.github.io/STANCE/tutorial.html) |
| **C-SIDE** | fixed effects, spline basis (`run.CSIDE.nonparam`), on RCTD | counts, coords, RCTD weights | flexible covariates, mature package | Bonferroni; QC removes true signal; FDP 81% | `devtools::install_github("dmcable/spacexr")` |
| **CTSV** | ZINB, x/y kernels, Cauchy | `SpatialExperiment` + W (n × K) | simple Bioconductor API | not rotation-invariant; scalability; FDP 73% | `BiocManager::install("CTSV")` |
| **spVC** | generalized Poisson, constant + spatially varying coefficients on a triangulation | counts, coords, covariates, boundary mesh | interpretable spatial coefficient maps | needs tissue boundary; FDP 72% | `devtools::install_github("shanyu-stat/spVC")` (+ BPST, Triangulation) · [tutorial](https://github.com/shanyu-stat/spVC/wiki/Tutorial-on-the-spVC) |
| **MMM (MCube)** | mixed model on RCTD reference | counts, coords, proportions, RCTD reference | sixth method in the official repo | — | `devtools::install_github("YangLabHKUST/MCube")` |
| **ctSVG** (Zhuang et al.) | per-cluster B-spline + permutation | single-cell-resolution data (Visium HD / Xenium) | made for single-cell resolution, like multiplexed imaging | requires clustering choices | `devtools::install_github("haotian-zhuang/ctSVG")` · [vignettes](https://haotian-zhuang.github.io/ctSVG/) |

**Adapting to protein data.** Single-cell imaging gives a one-hot cell-type matrix, so no deconvolution is needed. Count-based models (C-SIDE, CTSV, spVC, MMM) expect integers. Use spot-level protein UMIs (Visium CytAssist, D13), or treat pseudo-counts as an explicit experimental factor. Exclude lineage markers used for phenotyping when testing within that type.

## 2. Global SVG methods (baselines; also usable within one cell type)

| Method | Model | Language | Install | Notes |
|---|---|---|---|---|
| **SPARK-X** | non-parametric kernel covariance test | R | `devtools::install_github("xzhoulab/SPARK")` | best overall in the 2025 Genome Biology SVG benchmark; fast and well calibrated |
| **nnSVG** | NNGP, feature-specific length scale | R | `BiocManager::install("nnSVG")` | gives effect size `prop_sv` |
| **Moran's I** | spatial autocorrelation | Python (Squidpy) | `pip install squidpy` | strong cheap baseline |
| **SpatialDE** | Gaussian process | Python | `pip install spatialde` | O(n³), so subsample |
| **Hotspot** | kNN autocorrelation, `normal` model | Python | `pip install hotspotsc` | also finds marker modules |

## 3. Spatial-proteomics-native biomarker methods

| Method | Biomarker question | Language | Install |
|---|---|---|---|
| **Statial: SpatioMark** | does marker *m* in type A change with distance to / abundance of type B? (ctSV analogue) | R | `BiocManager::install("Statial")` |
| **Statial: Kontextual** | context-aware co-localisation robust to tissue structure | R | same |
| **spicyR** | does co-localisation differ by condition or predict survival? | R | `BiocManager::install("spicyR")` |
| **lisaClust** | spatial regions (local L-function clusters) | R | `BiocManager::install("lisaClust")` |
| **mistyR** | intra- vs. juxta- vs. paracrine marker dependencies | R | `BiocManager::install("mistyR")` |
| **imcRtools** | spatial graphs, neighbourhood aggregation, permutation interaction tests | R | `BiocManager::install("imcRtools")` |
| **scimap** | phenotyping, spatial distance/interaction, neighbourhoods, spatial LDA | Python | `pip install scimap` |
| **CytoCommunity** | condition-specific niches via supervised GNN | Python | [huBioinfo/CytoCommunity](https://github.com/huBioinfo/CytoCommunity) |
| **SPIAT** | infiltration, margins, spatial heterogeneity metrics | R | `BiocManager::install("SPIAT")` |
| **NaroNet** | weakly supervised patient-level biomarker discovery from images | Python | [djimenezsanchez/NaroNet](https://github.com/djimenezsanchez/NaroNet) |

## 4. Pretrained / foundation models

| Model | Training data | Weights | How to connect |
|---|---|---|---|
| **KRONOS** | more than 47M patches, 175 markers, 16 tissues, 8 platforms | HF `MahmoodLab/kronos` (**gated**: request access) | `create_model_from_pretrained("hf_hub:MahmoodLab/kronos")`. See notebook and [repo tutorials](https://github.com/mahmoodlab/KRONOS/tree/main/tutorials) |
| **VirTues** | 31–32 spatial proteomics datasets (spora) | HF `bunnelab/virtues`: `virtues-sp31` (MIT), `virtues-sp32` / `virtues-imc14` (CC BY-NC) | `hf_hub_download("bunnelab/virtues", "virtues-sp31/model.safetensors")` + [demo notebooks](https://github.com/bunnelab/virtues/tree/main/notebooks) |
| **Mesmer** | TissueNet | auto-download (DeepCell access token) | `deepcell.applications.Mesmer()` |
| **Cellpose** | generalist | auto-download | `cellpose.models.CellposeModel(model_type="cyto3")` |

**How these feed biomarker discovery:** segmentation → cell × marker matrix (input to everything else). Foundation-model embeddings → finer cell states, used as cell-type covariates in ctSV tests, or MIL-based patient stratification with attention maps as spatial biomarker candidates.

## 5. Evaluation framework

`05_evaluation_framework/metrics.py` implements the benchmark's metrics:
early precision, AUPRC, FDP and power at BH-q, type I error, KS distance of null p-values, and a count of **true pairs a method's QC dropped**. It also provides `rotate_coords`, `permute_within_type`, and `jaccard_topk` for ctSVGbench-style robustness tests.
The notebook includes a lightweight protein simulator (real coordinates and types, planted within-type patterns, optional spillover).
