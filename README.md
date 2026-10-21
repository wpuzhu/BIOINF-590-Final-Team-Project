# BIOINF 590 Final Team Project

**Biomarker identification in spatial proteomics data: a literature review of algorithms and public validation datasets**

This repo is our team's starting point. It collects:

1. a literature review of algorithms that find **spatially variable** and **cell-type-specific spatially variable** features (genes or proteins). In spatial proteomics we treat these features as candidate biomarkers.
2. a list of **public datasets** (links plus download scripts) for validating those algorithms.
3. **starter notebooks** that show how to call each algorithm, or how to load each pretrained model.

The conceptual anchor is the benchmark paper
> Li W., Ge X., Jiang Y. *Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework.* bioRxiv (2025). doi:[10.1101/2025.11.26.690782](https://www.biorxiv.org/content/10.1101/2025.11.26.690782v1)

- **Official code:** [weiqi-grace-li/ctsvg-benchmark](https://github.com/weiqi-grace-li/ctsvg-benchmark) (MIT). It includes method wrappers, scDesign3-based realistic simulators, and a pinned `renv.lock`. Intermediate data are on [Zenodo 21405877](https://zenodo.org/records/21405877). The repo uses a newer title, *"Realistic and Decomposable Benchmarking Reveals Limitations of Cell-Type-Specific Spatially Variable Gene Detection Methods"*, and adds a sixth method, MMM (MCube).
- **Companion benchmark:** [FangZY-Lab/ctSVGbench](https://github.com/FangZY-Lab/ctSVGbench), an independent ctSVG benchmark with rotation tests, deconvolution sensitivity, permutation nulls, and scalability. We borrow its stress tests too. Our own takeaways from both benchmarks are in [`docs/literature_review.md`](docs/literature_review.md).

> **Nothing here has been run yet.** The notebooks are templates. Function calls follow each package's documentation at the time of writing (Oct 2026). Check the linked docs before running anything.

---

## Repository layout

```
BIOINF-590-Final-Team-Project/
├── README.md                     ← you are here (brief intro to all datasets and algorithms)
├── simple_pipeline/              ← RUNNABLE small-batch zero-shot ctSVP pipeline + results (start here)
├── docs/
│   ├── literature_review.md      ← full review: ctSVG benchmark ideas, transfer to proteomics
│   └── project_plan.md           ← proposed validation design and team task split
├── datasets/
│   ├── README.md                 ← dataset descriptions, access notes, and recommended use
│   ├── datasets_catalog.csv      ← machine-readable table of every dataset and link
│   └── download_scripts/         ← R / Python / shell helpers that fetch the data
└── algorithms/
    ├── README.md                 ← algorithm table: model type, input, language, link
    ├── 01_ctSVG_methods/         ← CELINA, STANCE, C-SIDE, CTSV, spVC, MMM, ctSVG (R)
    ├── 02_global_SVG_methods/    ← SPARK-X, nnSVG (R); Moran's I/Squidpy, SpatialDE, Hotspot (Py)
    ├── 03_spatial_proteomics_biomarker_methods/ ← Statial (Kontextual/SpatioMark), spicyR,
    │                                              lisaClust, mistyR, imcRtools, scimap, CytoCommunity
    ├── 04_pretrained_foundation_models/         ← KRONOS, VirTues, Mesmer, Cellpose
    └── 05_evaluation_framework/  ← metrics from the benchmark (EP, AUPRC, FDP, type I error, KS)
```

> **New: [`simple_pipeline/`](simple_pipeline/)** runs 6 ctSVG/ctSVP methods (including the benchmark's own STANCE, CELINA and SPARK-X) zero-shot on 10 images from 3 public protein-imaging datasets. It aggregates their calls and validates them against the original annotations. Results: [`simple_pipeline/results/SUMMARY.md`](simple_pipeline/results/SUMMARY.md).

---

## Why the ctSVG benchmark matters for spatial proteomics

Multiplexed protein imaging (IMC, MIBI-TOF, CODEX/PhenoCycler, CyCIF) gives **single-cell resolution** data, typically 30–60 markers, with cell-type labels from phenotyping. A plain "spatially variable protein" test mostly rediscovers **lineage markers** (CD3, CD20, pan-CK). Those vary in space only because the cell types themselves are unevenly placed. Li et al. identify the same confounding in transcriptomics. The more useful biomarker question is:

> *Within one cell type (e.g., tumour cells or CD8 T cells), does a functional marker (PD-L1, Ki-67, HLA-DR, GZMB) vary in space, and does that pattern relate to the niche or to the outcome?*

That question is exactly the **cell-type-specific spatially variable (ct-SV)** problem. In this repo we call the protein version **ctSVPs** (cell-type-specific spatially variable proteins).

The benchmark's main lessons carry over directly:
- Most methods look fine on idealized simulations but **fail on realistic ones**. Early precision falls from above 80% to below 20%, and FDP at nominal 5% FDR ranges from 41% to 98%.
- **Rare cell types** get little power. Random-effect models (CELINA, STANCE) degrade the least.
- **QC filters** (CELINA, C-SIDE) can quietly remove true signals.
- Methods that don't model a constant cell-type effect (STANCE) **mistake marker genes for ctSVGs**.
- **No single method wins** on accuracy, calibration, scalability, and usability together.

---

## Algorithms at a glance

Details, install commands, and usage are in [`algorithms/README.md`](algorithms/README.md) and in each notebook.

| # | Method | Category | Model | Language | Code |
|---|---|---|---|---|---|
| 1 | **CELINA** (Shang, Wu, Zhou 2025, *Nat Commun*) | ctSVG | spatially varying coefficient, random effects, 11 kernels | R | [pekjoonwu/CELINA](https://github.com/pekjoonwu/CELINA) |
| 2 | **STANCE** (Su et al. 2025, *Nat Commun*) | ctSVG + SVG | linear mixed model, two-stage test, rotation-invariant | R | [Cui-STT-Lab/STANCE](https://github.com/Cui-STT-Lab/STANCE) |
| 3 | **C-SIDE** (Cable et al. 2022, *Nat Methods*) | ctSVG / ctDE | fixed effects, spline basis, on RCTD | R | [dmcable/spacexr](https://github.com/dmcable/spacexr) |
| 4 | **CTSV** (Yu & Luo 2022, *Bioinformatics*) | ctSVG | ZINB, 10 kernels, Cauchy combination | R (Bioc) | [Bioconductor: CTSV](https://bioconductor.org/packages/CTSV) |
| 5 | **spVC** (Yu & Li 2024, *Genome Biol*) | ctSVG / covariate-SV | generalized Poisson, spatially varying coefficients (bivariate splines over triangulation) | R | [shanyu-stat/spVC](https://github.com/shanyu-stat/spVC) |
| 5b | **MMM / MCube** (in the official benchmark repo) | ctSVG | mixed model on RCTD reference | R | [YangLabHKUST/MCube](https://github.com/YangLabHKUST/MCube) |
| 6 | **ctSVG** (Zhuang, Shang, Hou, Ji 2025, *Genome Biol*) | ctSVG at single-cell resolution | B-spline regression, permutation null | R | [haotian-zhuang/ctSVG](https://github.com/haotian-zhuang/ctSVG) |
| 7 | **SPARK-X** (Zhu, Sun, Zhou 2021) | global SVG | non-parametric kernel covariance test | R | [xzhoulab/SPARK](https://github.com/xzhoulab/SPARK) |
| 8 | **nnSVG** (Weber et al. 2023) | global SVG | nearest-neighbour Gaussian process | R (Bioc) | [Bioconductor: nnSVG](https://bioconductor.org/packages/nnSVG) |
| 9 | **Moran's I / Geary's C** (Squidpy) | global SVG baseline | spatial autocorrelation, permutation | Python | [scverse/squidpy](https://github.com/scverse/squidpy) |
| 10 | **SpatialDE / SpatialDE2** | global SVG | Gaussian process | Python | [Teichlab/SpatialDE](https://github.com/Teichlab/SpatialDE) |
| 11 | **Hotspot** | global SVG + modules | kNN autocorrelation, `normal` model suits protein data | Python | [YosefLab/Hotspot](https://github.com/YosefLab/Hotspot) |
| 12 | **Statial: Kontextual + SpatioMark** | proteomics, marker changes vs. proximity | context-aware L-function; marker ~ distance/abundance regression | R (Bioc) | [Bioconductor: Statial](https://bioconductor.org/packages/Statial) |
| 13 | **spicyR** | proteomics, colocalization ↔ outcome | mixed-effects / survival on L-function | R (Bioc) | [Bioconductor: spicyR](https://bioconductor.org/packages/spicyR) |
| 14 | **lisaClust** | proteomics, spatial regions | local L-function clustering | R (Bioc) | [Bioconductor: lisaClust](https://bioconductor.org/packages/lisaClust) |
| 15 | **mistyR** | proteomics, intra- and intercellular marker dependencies | multi-view random forest | R (Bioc) | [Bioconductor: mistyR](https://bioconductor.org/packages/mistyR) |
| 16 | **imcRtools** | proteomics workflow | neighbourhoods, interaction tests (histoCAT-style) | R (Bioc) | [Bioconductor: imcRtools](https://bioconductor.org/packages/imcRtools) |
| 17 | **scimap** | proteomics workflow | phenotyping, neighbourhoods, spatial LDA | Python | [labsyspharm/scimap](https://github.com/labsyspharm/scimap) |
| 18 | **CytoCommunity** | proteomics, condition-specific niches (supervised GNN) | graph neural network | Python | [huBioinfo/CytoCommunity](https://github.com/huBioinfo/CytoCommunity) |
| 19 | **KRONOS** (Mahmood lab 2025) | pretrained foundation model, multiplex images | ViT, 47M patches, 175 markers | Python | [mahmoodlab/KRONOS](https://github.com/mahmoodlab/KRONOS) · HF `MahmoodLab/kronos` (gated) |
| 20 | **VirTues** (Bunne lab, 2026 *Nature*) | pretrained foundation model, virtual biomarkers | multi-scale transformer + ESM-2 marker embeddings | Python | [bunnelab/virtues](https://github.com/bunnelab/virtues) · HF `bunnelab/virtues` |
| 21 | **Mesmer / DeepCell**, **Cellpose** | pretrained segmentation (upstream) | CNN | Python | [vanvalenlab/deepcell-tf](https://github.com/vanvalenlab/deepcell-tf), [MouseLand/cellpose](https://github.com/MouseLand/cellpose) |

---

## Datasets at a glance

Full descriptions are in [`datasets/README.md`](datasets/README.md). Every link is also in [`datasets/datasets_catalog.csv`](datasets/datasets_catalog.csv).

| # | Dataset | Technology | Tissue / disease | Size | Why useful | Access |
|---|---|---|---|---|---|---|
| D1 | Keren et al. 2018 | MIBI-TOF | TNBC | 41 pts, 36 markers | survival, tumour–immune architecture; `kerenSCE` in Statial | [Angelo lab](https://www.angelolab.com/mibi-data) |
| D2 | Jackson, Fischer et al. 2020 | IMC | breast cancer (Basel + Zurich) | ~350 pts, 35 markers | two cohorts for **replication**; clinical outcome | [Zenodo 3518284](https://zenodo.org/records/3518284), `imcdatasets` |
| D3 | Danenberg et al. 2022 (METABRIC) | IMC | breast cancer | ~700 pts, 37 markers | genomics and survival linkage | [Zenodo 5850952](https://zenodo.org/records/5850952) |
| D4 | Ali et al. 2020 (METABRIC) | IMC | breast cancer | 483 tumours | independent METABRIC IMC cohort | [IDR idr0076](https://idr.openmicroscopy.org/study/idr0076) |
| D5 | Damond et al. 2019 | IMC | pancreas, type 1 diabetes | 12 donors, 38 markers | disease stage contrast | `imcdatasets`, [Mendeley](https://data.mendeley.com/datasets/cydmwsfztj) |
| D6 | Hoch, Schulz et al. 2022 | IMC (+ RNAscope) | metastatic melanoma | 2 panels | protein + chemokine RNA on same tissue | `imcdatasets` |
| D7 | IMMUcan 2022 example | IMC | NSCLC and others | small | tutorial-sized, quick tests | `imcdatasets` |
| D8 | Schürch et al. 2020 | CODEX | colorectal cancer | 35 pts, 56 markers | cellular neighbourhoods ↔ survival (CLR vs DII) | [Mendeley mpjzbtfgfr](https://data.mendeley.com/datasets/mpjzbtfgfr/1) |
| D9 | Hickey et al. 2023 / HuBMAP | CODEX | healthy intestine | 64 sections, 57 markers | large, healthy reference; used by STELLAR | [HuBMAP portal](https://portal.hubmapconsortium.org/) |
| D10 | HTAN CyCIF atlases (e.g., Lin et al. 2023 CRC) | CyCIF / Orion | CRC, many cancers | very large | whole-slide multiplexed IF | [HTAN portal](https://humantumoratlas.org/explore) |
| D11 | **spora [data]** (Bunne lab) | mixed (IMC, MIBI, CODEX…) | 31+ harmonized datasets | large | one-stop harmonized collection | [spora.epfl.ch](https://spora.epfl.ch/datasets.html) |
| D12 | Squidpy built-ins (`imc`, `mibitof`) | IMC / MIBI | breast / colorectal | small | one-line loaders for notebooks | `squidpy.datasets` |
| D13 | **10x Visium CytAssist Gene + Protein**, human tonsil | spot RNA + 35 antibody UMIs | tonsil | 1–2 sections | **protein counts at spot level**: runs spot-based ctSVG tools directly, with RNA for deconvolution | [10x dataset](https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard) |
| D14 | 10x Xenium breast / ovarian / lymph node | Xenium (RNA, subcellular) | the benchmark's source data | — | reproduces the benchmark's realistic simulation | [10x datasets](https://www.10xgenomics.com/datasets) |
| D15 | STANCE / CELINA reference data (HER2+ breast ST, kidney cancer, mouse OB) | ST / Visium | — | — | positive controls for ctSVG pipelines | [Zenodo 4751624](https://zenodo.org/records/4751624), [Mendeley g67bkbnhhg](https://data.mendeley.com/datasets/g67bkbnhhg/1) |
| D16 | ctSVG paper datasets (7 Visium HD + 3 Xenium) | Visium HD / Xenium | CRC, lung, pancreas, mouse | — | single-cell-resolution ctSVG examples | see `datasets/README.md` |
| D17 | Spatial-CITE-seq / SPOTS | spatial RNA + ADT protein | human tonsil, skin; mouse spleen | — | spatial protein UMIs + transcriptome | GEO (see catalog) |

---

## Suggested validation strategy (summary)

The full plan is in [`docs/project_plan.md`](docs/project_plan.md).

1. **Simulation with known ground truth.** Follow Li et al. Fit a generator to real IMC/CODEX single-cell data (scDesign3 with a Gaussian family on arcsinh intensities, or a custom simulator). Plant known ctSVPs inside chosen cell types, then score with EP, AUPRC, FDP@5%, type I error, and KS distance (code in `algorithms/05_evaluation_framework`).
2. **Calibration on real data.** Permute coordinates within each cell type and check that p-values are uniform. Also rotate the tissue to check rotation invariance, as ctSVGbench does.
3. **Biological positive controls.** Use known spatially organised functional markers, e.g., Ki-67 at the invasive margin, PD-L1/IDO on tumour cells next to CD8 T cells, HLA-DR on macrophages.
4. **Cross-cohort replication.** Jackson Basel ↔ Zurich, Danenberg ↔ Ali (METABRIC), Keren MIBI ↔ IMC breast.
5. **Clinical relevance.** Associate the ctSVP scores per patient with survival (spicyR/Statial style).

---

## Team

BIOINF 590, University of Michigan. Add names and roles here.

## Citation note

Each dataset and algorithm keeps its original license. Please cite the original papers (listed in `docs/literature_review.md`) in the final report.
