# Literature Review: Biomarker Identification Algorithms for Spatial Proteomics

*BIOINF 590 Final Team Project, draft (Oct 2026)*

## 1. Scope and framing

**Spatial proteomics** here means multiplexed antibody-based imaging at single-cell resolution: Imaging Mass Cytometry (IMC), MIBI-TOF, CODEX/PhenoCycler, CyCIF/Orion, and sequencing-based spatial protein assays (10x Visium CytAssist Gene + Protein, Spatial-CITE-seq, SPOTS). A typical processed dataset is a cell × marker intensity matrix with 30–60 markers, plus x/y centroids, a phenotype label for each cell, and patient/clinical metadata.

"Biomarker identification" in this setting covers three related tasks:

| Task | Unit | Typical question | Representative methods |
|---|---|---|---|
| **A. Spatially variable features** | marker | Does marker *m* show spatial structure? | SPARK-X, nnSVG, Moran's I, SpatialDE, Hotspot |
| **B. Cell-type-specific spatially variable features (ctSV)** | marker × cell type | Does *m* vary in space **within** cell type *k*, beyond composition? | CELINA, STANCE, C-SIDE, CTSV, spVC, ctSVG; Statial-SpatioMark (proteomics analogue) |
| **C. Spatial features ↔ outcome** | patient | Does a spatial pattern (co-localization, niche, ctSV score) predict survival or response? | spicyR, Statial-Kontextual, lisaClust, CytoCommunity, KRONOS / VirTues + MIL |

Our project centres on **Task B**. It is the focus of the anchor benchmark, and it is the most biologically specific route to functional biomarkers. Tasks A and C are the baseline and the downstream readout.

## 2. The anchor paper: Li, Ge & Jiang (2025)

*Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework*, bioRxiv 10.1101/2025.11.26.690782 (Oregon State University).

### 2.1 Problem statement
Cell types are unevenly distributed across tissue. A gene that is simply a marker of a spatially clustered cell type looks spatially variable even when its expression **within** that cell type is constant. ctSVG methods try to separate (i) composition-driven variation from (ii) true within-cell-type spatial variation.

**Code:** https://github.com/weiqi-grace-li/ctsvg-benchmark (MIT). Layout: `method_wrappers/`, `simulators/realistic/`, `run/{idealized,realistic,decomposed,real_data,scalability}.R`, pinned `renv.lock`. Data: Zenodo 21405877. The repo's current title is *Realistic and Decomposable Benchmarking Reveals Limitations of Cell-Type-Specific Spatially Variable Gene Detection Methods*. It benchmarks **six** methods: the five below plus **MMM** (via the MCube package, YangLabHKUST/MCube).

### 2.2 Methods benchmarked
| Method | Effect type | Spatial basis | P-value combination | Notes from benchmark |
|---|---|---|---|---|
| CELINA | random | 11 kernels | Cauchy | best FDP (41%) in realistic settings; QC removes ~8% of true ctSVGs even in idealized settings and most under realistic patterns; scalability limits |
| STANCE | random | 1 distance kernel | two-stage | robust for minority cell types, but **calls marker genes ctSVGs** (no constant cross-cell-type effect); FDP 98%; runtime grows fast |
| C-SIDE | fixed | splines | Bonferroni | FDP 81%; QC filter removes true signals |
| CTSV | fixed | 10 kernels (x and y) | Cauchy | FDP 73%; scalability limits |
| spVC | fixed | bivariate spline basis | — | FDP 72% |

### 2.3 Simulation design (the idea we reuse)
- **Idealized**: classic synthetic patterns. Almost every method looks good here.
- **Realistic**: start from real **Xenium** single-cell data, fit **scDesign3** to produce ctSVGs and null genes with realistic marginals and spatial patterns, then aggregate into Visium-like pseudo-spots.
- **Decomposed**: move from idealized to realistic one factor at a time. The six factors are cell-type diversity, cell locations, composition patterns, null-gene distributions, capture efficiency, and the shape of the intra-cell-type spatial pattern. This shows **which** factor breaks each method.

### 2.4 Metrics
Early precision (EP), AUPRC, false-discovery proportion at nominal FDR 5%, type I error, KS distance of null p-values from Uniform(0,1), runtime and peak memory, and a qualitative usability score.

### 2.5 Main conclusions
1. Average EP falls from above 80% (idealized) to below 20% (realistic).
2. All methods inflate false discoveries in realistic settings.
3. Power falls with cell-type abundance. Random-effect models are more robust for rare types.
4. QC filters in CELINA and C-SIDE remove most true ctSVGs under realistic patterns.
5. No method dominates. The authors give a decision guide (their Fig. 6).
6. Open problems: power for rare cell types, better QC, scalability.

### 2.6 Companion benchmark: ctSVGbench
[FangZY-Lab/ctSVGbench](https://github.com/FangZY-Lab/ctSVGbench) is an independent public benchmark of ctSVG methods. Its design adds:
- **Rotation tests** (`fig4_rotate_sc.R`, `fig4_rotate_spot.R`): fixed-effect x/y bases are not rotation-invariant, so results change if the slide is rotated. STANCE was explicitly designed to fix this.
- **Deconvolution sensitivity** (`1-runsim_noRCTD.R`): how much ctSVG results depend on the RCTD cell-type proportions.
- **FPR vs. cell-type proportion** and **AUC vs. dropout**.
- **Permutation nulls** and **scalability** scripts.
- Real-data runs at both single-cell and spot resolution, including a lung cancer dataset.

Reading both benchmarks together, our checklist for any method is: *calibrated under permutation, rotation-invariant, robust to cell-type proportion error, powered for rare types, and scalable to 10^5–10^6 cells.*

## 3. Transferring ctSVG ideas to spatial proteomics

| Aspect | Spatial transcriptomics (benchmark) | Spatial proteomics (our target) | Implication |
|---|---|---|---|
| Resolution | spots (Visium) or cells (Xenium) | **single cells** (IMC/MIBI/CODEX); spots for CytAssist protein | cell-type proportions are one-hot and known from phenotyping; no deconvolution needed |
| Features | 10^2–10^4 genes | 30–60 markers | multiple-testing burden is small; per-feature runtime matters less, per-cell runtime matters more |
| Data type | UMI counts | continuous intensities (arcsinh/scaled); **counts** only for CytAssist / CITE-type ADTs | count models (C-SIDE, CTSV-ZINB, spVC-Poisson) need adaptation: Gaussian versions, or rounding/scaling with caution. Kernel/GP tests on transformed data (CELINA, STANCE, SPARK-X, nnSVG) carry over more naturally |
| Lineage markers | marker genes | marker proteins used **to define** cell types | circularity: a marker used for phenotyping should be excluded or flagged when testing within that type |
| Segmentation spillover | n/a (spots) / moderate (Xenium) | **lateral spillover** between neighbouring cells is a major artefact | a neighbouring cell's marker can leak in and create a false ctSVP near that cell type. Control with spillover compensation (e.g., REDSEA, CATALYST) or neighbour-composition covariates |
| Cohort structure | 1–few slides | many patients, often TMA cores | test per image, then meta-analyse across images (e.g., Cauchy/Fisher combination or mixed models) |

**Proposed protein-specific adaptation (project idea).** For each (marker *m*, cell type *k*) pair:
1. subset cells of type *k*;
2. regress out image/batch effects and the marker's lineage-level mean;
3. apply a spatial test (kernel score test à la SPARK-X/CELINA, NNGP à la nnSVG, or Moran's I) **within type *k***, optionally with neighbour-composition covariates to absorb spillover;
4. combine p-values across images;
5. validate with the framework in `algorithms/05_evaluation_framework`.

Statial's **SpatioMark** already does a version of steps 1–3. It regresses marker intensity in cell type *k* on distance to (or abundance of) cell type *j*. That makes it an interpretable, proteomics-native ctSV baseline.

## 4. Algorithms reviewed

### 4.1 Cell-type-specific spatially variable gene methods
- **CELINA**: Shang L, Wu P, Zhou X. *Statistical identification of cell type-specific spatially variable genes in spatial transcriptomics.* Nat Commun 16, 2025. doi:10.1038/s41467-025-56280-4. Spatially varying coefficient model with random effects. Handles both single-cell and spot resolution. Code: github.com/pekjoonwu/CELINA.
- **STANCE**: Su H, Wu Y, Chen B, Cui Y. *STANCE: a unified statistical model to detect cell-type-specific spatially variable genes in spatial transcriptomics.* Nat Commun 16, 2025. doi:10.1038/s41467-025-57117-w. Linear mixed model with a two-stage test (unified SVG test, then per-cell-type ctSVG test). Rotation-invariant. Code: github.com/Cui-STT-Lab/STANCE.
- **C-SIDE**: Cable DM et al. *Cell type-specific inference of differential expression in spatial transcriptomics.* Nat Methods 19, 2022. Built on RCTD deconvolution in `spacexr`. Supports parametric and non-parametric (spline) spatial covariates.
- **CTSV**: Yu S, Luo X. *Detecting cell-type-specific spatially variable genes with CTSV.* Bioinformatics, 2022. ZINB regression with kernel spatial terms. Bioconductor.
- **spVC**: Yu S, Li WV. *spVC for the detection and interpretation of spatial gene expression variation.* Genome Biol 25, 2024. doi:10.1186/s13059-024-03245-3. Generalized Poisson with constant and spatially varying coefficients.
- **ctSVG**: Zhuang H, Shang X, Hou W, Ji Z. *Identifying cell-type-specific spatially variable genes with ctSVG.* Genome Biol, 2025. doi:10.1186/s13059-025-03870-6. Built for Visium HD and Xenium (single-cell resolution): H&E segmentation → clustering → B-spline (PreTSA) per cluster → permutation null. Conceptually closest to our proteomics setting.

### 4.2 Global SVG methods (baselines)
- **SPARK-X**: Zhu J, Sun S, Zhou X. Genome Biol 2021. Non-parametric covariance test that scales linearly. It was the top performer in Li et al.'s 2025 Genome Biology benchmark of 14 SVG methods across 96 datasets ([Genome Biol 2025](https://link.springer.com/article/10.1186/s13059-025-03731-2)), which also found most methods poorly calibrated.
- **nnSVG**: Weber LM et al. Nat Commun 2023. NNGP with gene-specific length scales.
- **Moran's I**: a classic autocorrelation statistic. Competitive in benchmarks and the cheapest baseline (Squidpy, Giotto, Seurat).
- **SpatialDE / SpatialDE2**: Svensson et al. 2018; Kats et al. 2021. GP-based.
- **Hotspot**: DeTomaso & Yosef 2021. kNN-graph autocorrelation. Its `normal` model fits continuous protein intensities.
- Other benchmarks: [categorisation of 33 SVG methods (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11160866/); [OpenProblems SVG benchmark](https://openproblems.bio/benchmarks/spatially_variable_genes/); [pinellolab/SVG_Benchmarking](https://github.com/pinellolab/SVG_Benchmarking) (scDesign3-based simulations).

### 4.3 Spatial proteomics biomarker methods
- **Statial** (Ameen et al., Sydney Precision Data Science). **Kontextual** measures cell–cell localisation relative to a parent population, which makes it robust to tissue structure. **SpatioMark** finds marker changes associated with proximity to other cell types. Example data: `kerenSCE` (Keren 2018 MIBI).
- **spicyR** (Canete et al. 2022, *Bioinformatics*): tests whether co-localisation (L-function) differs between patient groups or relates to survival, with mixed effects for multiple images per patient.
- **lisaClust** (Patrick et al. 2023): local indicators of spatial association → spatial regions/niches.
- **mistyR** (Tanevski et al. 2022, *Genome Biol*): multi-view model separating intracellular, juxtacrine, and paracrine marker dependencies. Applied to IMC breast cancer.
- **imcRtools** (Windhager et al. 2023, *Nat Protoc*, "steinbock + imcRtools" workflow): spatial graphs, neighbourhood aggregation, histoCAT-style permutation interaction tests.
- **scimap** (Nirmal et al. 2024, *JOSS*): Python toolkit for CyCIF/CODEX/IMC: phenotyping, spatial distances, neighbourhoods, spatial LDA, interaction tests.
- **SPIAT** (Feng et al. 2023, *Nat Commun*): spatial metrics for multiplex IF, e.g., immune infiltration and tumour margins.
- **CytoCommunity** (Hu et al. 2024, *Nat Methods*): GNN that finds condition-specific tissue cellular neighbourhoods, linking niches to disease labels directly.
- **NaroNet** (Jiménez-Sánchez et al. 2022, *Med Image Anal*): weakly supervised patient-level biomarker discovery from multiplex images.
- **Cellular neighbourhoods** (Schürch et al. 2020, *Cell*): k-means on windowed cell-type composition, associated with CRC survival.

### 4.4 Pretrained / foundation models
- **KRONOS** (Mahmood lab, arXiv 2506.03373, 2025): ViT pretrained on more than 47 million patches, 175 markers, 16 tissues, 8 platforms. Tutorials cover cell/patch phenotyping, tissue search, and patient stratification. Weights are gated on Hugging Face (`MahmoodLab/kronos`). A successor, KRONOS2 (268 markers, CORAL toolkit), has been announced; check the repo.
- **VirTues** (Bunne lab, Nature 2026; arXiv 2501.06039): multi-scale transformer with ESM-2 protein-language-model marker embeddings, so it generalises to unseen marker panels. Supports virtual staining (inpainting a missing marker) and "virtual biomarker discovery". Weights `bunnelab/virtues` (sp31 MIT, sp32 CC BY-NC). Comes with **spora** (harmonized datasets) and **spora-bench**.
- **Mesmer** (Greenwald et al. 2022, *Nat Biotech*, DeepCell) and **Cellpose** (Stringer et al. 2021): pretrained whole-cell segmentation, the upstream step that defines single-cell protein quantification.
- **CellSighter** (Amitay et al. 2023) and **STELLAR** (Brbić et al. 2022): supervised cell-type annotation for multiplexed images / CODEX.

## 5. Gaps and opportunities for our project
1. **No published benchmark of ctSV methods on spatial proteomics.** Searches in Oct 2026 found SVG/ctSVG benchmarks only for transcriptomics. Porting the Li et al. framework to IMC/CODEX would be a genuine contribution.
2. **Count assumptions.** Several ctSVG methods assume counts. Measuring how they behave on transformed intensities vs. a Gaussian alternative is an easy, informative experiment.
3. **Spillover as a confounder.** It is specific to imaging proteomics and missing from transcriptomics benchmarks. It can be added as a seventh "decomposed" factor.
4. **Patient-level aggregation.** Most ctSVG tools assume one slide. Multi-image meta-analysis and its calibration are open questions.
5. **Foundation-model embeddings as covariates.** KRONOS/VirTues cell embeddings could define finer cell states, which could reduce confounding by cell type.

## 6. References (key)
- Li W, Ge X, Jiang Y. bioRxiv 2025. doi:10.1101/2025.11.26.690782. Code: https://github.com/weiqi-grace-li/ctsvg-benchmark
- ctSVGbench: https://github.com/FangZY-Lab/ctSVGbench
- Shang L, Wu P, Zhou X. Nat Commun 2025. doi:10.1038/s41467-025-56280-4 (CELINA)
- Su H et al. Nat Commun 2025. doi:10.1038/s41467-025-57117-w (STANCE)
- Cable DM et al. Nat Methods 2022 (C-SIDE)
- Yu S, Luo X. Bioinformatics 2022 (CTSV)
- Yu S, Li WV. Genome Biol 2024. doi:10.1186/s13059-024-03245-3 (spVC)
- Zhuang H et al. Genome Biol 2025. doi:10.1186/s13059-025-03870-6 (ctSVG)
- Zhu J, Sun S, Zhou X. Genome Biol 2021 (SPARK-X)
- Weber LM et al. Nat Commun 2023. doi:10.1038/s41467-023-39748-z (nnSVG)
- Li et al. Genome Biol 2025. doi:10.1186/s13059-025-03731-2 (SVG benchmark, 14 methods)
- Song D et al. Nat Biotechnol 2024 (scDesign3)
- Keren L et al. Cell 2018; Jackson HW et al. Nature 2020; Danenberg E et al. Nat Genet 2022; Ali HR et al. Nat Cancer 2020; Damond N et al. Cell Metab 2019; Hoch T et al. Sci Immunol 2022; Schürch CM et al. Cell 2020; Hickey JW et al. Nature 2023
- KRONOS: arXiv 2506.03373; VirTues: arXiv 2501.06039 / Nature 2026
