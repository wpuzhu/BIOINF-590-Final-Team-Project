# Proposed Project Plan

## Research question
Can ctSVG detection methods, built for spatial transcriptomics, find **cell-type-specific spatially variable proteins (ctSVPs)** in multiplexed imaging data reliably? Which method gives calibrated, biologically meaningful biomarkers?

## Aim 1: Protein-adapted simulation benchmark (after Li, Ge & Jiang 2025)
- **Base data**: Jackson 2020 IMC (D2) and Schürch 2020 CODEX (D8), single cells with phenotypes.
- **Generator**: scDesign3 (`family_use = "gaussian"` on arcsinh(x/5) intensities, or `"nb"` on rounded raw IMC counts). Spatial terms are fitted per cell type. Planted ctSVPs get Gaussian-bump, gradient, or niche-dependent ("near CD8 T cells") patterns.
- **Decomposed factors**: the paper's six (cell-type diversity, locations, composition patterns, null distribution, capture efficiency, pattern shape) **plus** two protein-specific factors:
  - (7) **lateral spillover**: mix x% of each neighbour's signal into the cell;
  - (8) **phenotyping error**: flip y% of cell-type labels.
- **Methods**: CELINA, STANCE, C-SIDE, CTSV, spVC, ctSVG, plus baselines (SPARK-X, nnSVG, Moran's I within cell type, Statial-SpatioMark).
- **Metrics**: EP, AUPRC, FDP@5%, type I error, KS, runtime and memory (`algorithms/05_evaluation_framework/metrics.py`).

## Aim 2: Real-data calibration and robustness (after ctSVGbench)
- Permutation null: shuffle coordinates within each cell type, then check type I error and KS.
- Rotation test: rotate coordinates by 30/60/90° and measure Jaccard overlap of the top-k calls.
- Phenotype-label perturbation, and including vs. excluding lineage markers.

## Aim 3: Biological validation and replication
- Positive controls: Ki-67 / pS6 gradients in tumour cells; PD-L1 / IDO1 in tumour cells near CD8⁺ T cells; HLA-DR in macrophages at the tumour–stroma interface.
- Replicate across cohorts: Jackson Basel ↔ Zurich; Danenberg ↔ Ali (METABRIC); Keren MIBI (TNBC subset).
- Outcome association: patient-level ctSVP scores → Cox model (or spicyR/Statial-style tests).
- Optional cross-modality check: Visium CytAssist tonsil (D13). Compare protein ctSV calls with RNA ctSV calls for the matched gene.

## Suggested team split (edit as needed)
| Member | Responsibility |
|---|---|
| A | Data wrangling: download D1/D2/D8, build a common `AnnData`/`SpatialExperiment` format |
| B | ctSVG methods (R): CELINA, STANCE, C-SIDE, CTSV, spVC, ctSVG |
| C | Baselines + proteomics methods: SPARK-X, nnSVG, Squidpy, Statial, spicyR |
| D | Simulation + evaluation framework; foundation-model embeddings (KRONOS/VirTues) |

## Timeline (approximate)
1. Week 1–2: data download and harmonization; run one method end-to-end on `kerenSCE`.
2. Week 3–4: simulation framework; run all methods.
3. Week 5: real-data calibration and replication.
4. Week 6: write-up and figures.
