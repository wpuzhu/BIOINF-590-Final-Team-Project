> **v1 component (kept in v2).** This was the original project's main pipeline. In v2 it is reused to test whether the *residual* protein signal left after RNA prediction is still spatially variable within cell types (E6 / idea N2), and to compare RNA vs protein spatial variability (idea N1). Paths are relative, so it runs from this folder as before.

# Simple pipeline: small-batch, zero-shot ctSVP detection in spatial proteomics

The goal is to find **cell-type-specific spatially variable proteins (ctSVPs)** on public multiplexed-imaging data, using the methods and evaluation logic of the baseline benchmark:

> Li W, Ge X, Jiang Y. *Benchmarking cell-type-specific spatially variable gene detection methods* (bioRxiv 2025, doi:10.1101/2025.11.26.690782; code: [weiqi-grace-li/ctsvg-benchmark](https://github.com/weiqi-grace-li/ctsvg-benchmark)).

The pipeline aggregates the methods' calls to find what they agree on, then checks the agreed calls against the **original dataset annotations**.

* **Small batch:** 10 images in total, from 3 public datasets on 3 imaging platforms.
* **Zero-shot:** every method runs with fixed package defaults. Nothing is tuned or trained, and the validation annotations are never used for detection; the detectors see only marker values, x/y coordinates and cell-type labels.

**Results:** [`results/SUMMARY.md`](results/SUMMARY.md), auto-generated with tables and figures. The headline findings are at the bottom of this page.

---

## Pipeline

```
download ─► prep ─► detect (6 ct methods + 2 global baselines) ─► aggregate ─► validate ─► report
            │        │                                             │             │
            │        ├ observed coordinates                        │             └ original annotations
            │        └ coordinates permuted within cell type (null)│               + permutation null
            └ one table per image: cells × markers, x, y, cell_type (+ niche/group for validation only)
```

| Step | Script | What it does |
|---|---|---|
| 0 | `scripts/install_r_methods.sh` | Installs SPARK, STANCE and CELINA into a project-local `.Rlib` (your global R library is untouched) |
| 1 | `scripts/download_data.sh`, `src/prep.py` | Downloads and checksum-verifies the data. Applies the same arcsinh transform to all raw intensities (cofactor = 20% of each marker's 99th percentile) and writes standard per-image tables |
| 2a | `src/methods_py.py` | Python methods (below), on observed coordinates and on coordinates permuted within each cell type |
| 2b | `src/methods_r.R`, `src/run_r_parallel.sh` | The benchmark's own R methods, one image per process |
| 3 | `src/aggregate.py` | Per-method BH-FDR, then a consensus call, method agreement, cross-image replication and cross-dataset themes |
| 4 | `src/validate.py` | Validation against the annotations and the null (V1–V6 below) |
| 5 | `src/report.py` | Figures plus `results/SUMMARY.md` |

Run everything with `bash run_pipeline.sh` (about 30–40 min on an 8-core laptop), or the Python-only part with `SKIP_R=1 bash run_pipeline.sh` (about 1 min).
All settings are in [`src/config.py`](src/config.py).

## Datasets (small batch, all annotated, all public)

| Dataset | Platform | Images used | Annotations used for validation |
|---|---|---|---|
| Schürch et al. 2020, *Cell*: colorectal cancer | CODEX, 56 markers | 6 TMA spots: 3 CLR + 3 DII patients, the spots closest to the median size | 29 cell types; **authors' 9 cellular neighbourhoods**; **authors' per-channel descriptions** ("CD5 - T cells"); CLR/DII group |
| Jackson et al. 2020, *Nature*: breast cancer (Squidpy subset) | IMC, 31 markers | 1 image | 11 cell types |
| Hartmann et al. 2021, *Nat Biotech*: colorectal cancer (Squidpy subset) | MIBI-TOF, 35 markers | 3 points, 2 donors | 8 cell types |

## Methods (one per model class from the benchmark)

Single-cell data means the cell-type "composition" is one-hot. A ct test of marker *m* in type *k* asks whether *m* varies in space **among the cells of type k**.

| Method | Class (as in the benchmark) | Implementation |
|---|---|---|
| `ct_STANCE` | random-effect kernel, two-stage | official R package (Cui-STT-Lab/STANCE), `normalized = TRUE`, on a stratified ~600-cell subsample per image (runtime) |
| `ct_CELINA` | random-effect spatially varying coefficient, 11 kernels | official R package (pekjoonwu/CELINA). Its count-based gene QC is **bypassed** so every marker is tested; the benchmark showed that this QC drops true signals |
| `ct_sparkx_R` | fixed-effect projection kernels, Cauchy combination | official SPARK-X (xzhoulab/SPARK), run within each type |
| `ct_sparkx_py` | same as above | Python port of the SPARK-X statistic (cross-checked against R) |
| `ct_kernel_score` | random-effect Gaussian-kernel score test (3 bandwidths) | Python |
| `ct_moran` | spatial autocorrelation | Python, kNN (k = 10) Moran's I |
| `global_moran`, `global_sparkx_py` | ignore cell type: the confounded analysis | Python baselines |

C-SIDE, CTSV, spVC and MMM are not included. They assume RNA counts plus RCTD deconvolution, which doesn't fit continuous single-cell protein intensities without extra modelling; see `../algorithms/01_ctSVG_methods/` to add them.

**Consensus:** ACAT-combined q < 0.05 **and** a strict majority of the methods that tested the pair call it at BH q < 0.05.

## Validation (against original annotations + nulls)

| Check | Question | Uses |
|---|---|---|
| V1 / V1b | Do global SVPs just rediscover cell-type identity (the benchmark's central warning)? | annotated cell types (η² per marker) |
| V2 | Is a ct call associated with tissue context? | authors' neighbourhood labels (Schürch); annotated neighbour composition |
| V3 | **Spillover**: does the marker rise where neighbours that express it much more are abundant? | annotated cell types |
| V3b | Is the marker native to the cell type (within 1.5 SD of the best-expressing type)? | annotated cell types |
| V3c | Does the pattern survive removing each cell's overall brightness (shared regional or tile background)? | — |
| V3d | Is the marker documented for a different lineage than the tested type? | authors' channel descriptions (Schürch) |
| V4 | Calibration: type-I error and KS under within-type permutation | null |
| V5 | Replication across images, vs. the permuted null | multi-image datasets |
| V6 | Python SPARK-X port vs. official R | — |

A call is **validated** if it is a consensus call that is annotation-supported, passes V3/V3b/V3c/V3d, and replicates in ≥ 2 images. Single-image datasets count as replicated by default, because replication can't be tested there.

> "Validated" means *it survived every automated check against the annotations*. It is a short-list for follow-up, not proof of biology. V2 is partly circular: any spatial signal tends to correlate with a spatial annotation. That is why V3–V3d exist.

## Headline findings (run of 21 Oct 2026)

See [`results/SUMMARY.md`](results/SUMMARY.md) for every table.

1. **The benchmark's central warning holds on protein data.** Global SVG scores track how much of a marker's variance cell-type identity explains (mean Spearman ρ with η² = 0.19–0.55). Within-type methods mostly remove this (ρ ≈ −0.05 to 0.22), except STANCE (ρ = 0.22–0.52). This matches the benchmark's finding that STANCE mistakes cell-type markers for ctSVGs, because it does not model a constant cross-cell-type effect. [fig 1](results/figures/fig1_confounding_eta2.png)
2. **Calibration under within-type permutation:** CELINA 2–4%, kernel score / Moran / SPARK-X 3–7%, **STANCE 33–47%**, with 125–474 false BH calls per dataset. This mirrors the benchmark's STANCE FDP of 98%. [fig 2](results/figures/fig2_calibration.png)
3. **What the methods have in common:** the five calibrated methods agree on most calls (Jaccard 0.65–0.80). STANCE is the outlier (0.52–0.59). The Python SPARK-X port reproduces official R SPARK-X (Spearman 0.96–0.99, Jaccard 0.95). [fig 3](results/figures/fig3_method_agreement.png)
4. **Raw consensus is not enough.** 1,546 consensus calls shrink to **412 validated candidates** (156 MIBI, 47 IMC, 209 CODEX). Removed calls were not annotation-supported, or were spillover (e.g. SMA in CD4 T cells next to fibroblasts; CK in T cells next to epithelium), not native to the cell type (e.g. FOXP3 in tumour cells), background-driven, or documented for another lineage. Before validation, immune markers (CD3, CD45, CD68) inside tumour cells were among the most common cross-dataset hits. [fig 4](results/figures/fig4_validation_funnel.png)
5. **Validated themes shared across datasets:** Ki67, pan-cytokeratin and EGFR/E-cadherin within tumour cells, and SMA/vimentin within stromal or endothelial cells. MIBI-TOF epithelial cells show spatially variable **metabolic** markers (GLUT1, LDHA, VDAC1, ATP5A, NRF2p, CD36) in all 3 images, consistent with Hartmann et al.'s report of spatially organised metabolic states. 65% of calls replicate in ≥ 2 images, and none do under the permutation null.
6. **Caveats:**
   - The Jackson IMC data are a single image, so replication can't be tested there; its candidates (e.g. Slug or fibronectin in T cells) are weaker.
   - A few implausible calls survive, e.g. PD-1 within tumour cells. The checkpoint channel is "function", so the lineage check allows it.
   - STANCE failed on CODEX spot16_A ("sample is too sparse").
   - With a small batch, ACAT q-values saturate at about 1e-14, so the ranks within the top lists are ties.

## Next steps for the team
- Scale the batch up: all 140 CODEX spots, the full Jackson cohort via `imcdatasets`, and Keren MIBI. Then test CLR vs DII associations of validated ctSVPs (the `group` column is already in the processed metadata).
- Add a protein version of the benchmark's realistic simulation (see `../algorithms/05_evaluation_framework/`), so that precision/FDP can be measured against known truth instead of annotation proxies.
- Add C-SIDE / CTSV / spVC on spot-level protein counts (10x Visium CytAssist), where their count models apply.
