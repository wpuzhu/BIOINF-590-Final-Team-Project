# Datasets

We do **not** store raw data in git. Files are large and each has its own license. This folder holds:
- `datasets_catalog.csv`: every dataset with link, technology, **pairing** (same cell / same section / adjacent / single modality) and its **v2 role**;
- `download_scripts/`: helpers that fetch data into `datasets/raw/` (git-ignored). Paired RNA + protein data: [`download_scripts/download_paired_multiomics.md`](download_scripts/download_paired_multiomics.md).

The v2 question needs RNA and protein from the same tissue, so datasets are ranked by **how tightly RNA and protein are paired**. Tighter pairing gives cleaner ablations; looser pairing is more realistic.

| Tier | Pairing | Datasets | Used for |
|---|---|---|---|
| **1a** | same cell (imaging) | **P1** Xenium Protein RCC, **P2** Xenium Protein ccRCC | main ablation dataset (E1–E6); H2 stage 1 |
| **1b** | same section, spot/pixel RNA | **P3** DBiTplus (whole transcriptome + ~50-plex CODEX), **D17a** spatial-CITE-seq (273 ADTs), **D13** Visium CytAssist Gene + Protein, **D17b/P5** SPOTS / Stereo-CITE-seq | spot-level E1–E4; E7 cross-dataset transfer; N1 |
| **2** | adjacent sections | **P4** SPATCH (Xenium 5K, CosMx 6K, Visium HD, Stereo-seq + 16-marker CODEX + scRNA-seq), D6 Hoch melanoma (IMC + RNAscope) | H1 main data; H2 stage 2; E8 registration error; N4 validation |
| **3** | protein only | D1–D5, D7–D12 (IMC, MIBI, CODEX, CyCIF; e.g. **D8** Schürch CRC CODEX, 56 markers) | H1 protein targets; H2 stage 3 (does the student transfer?) |
| **4** | RNA only | D14 Xenium, SPATCH SPT without CODEX, D15/D16 | H1 SPT side; N4 virtual protein |

> Check the license on each landing page before redistributing. 10x datasets need a short download form.

---

## Tier 1a. Same-cell RNA + protein (imaging)

### P1. 10x Xenium Protein: FFPE human renal cell carcinoma
Xenium In Situ Gene and Protein Expression run (Xenium Onboard Analysis v4.0). RNA: Human Multi-Tissue and Cancer panel (377 genes), plus 28 added targets for genes encoding the assayed proteins and for cell typing. Protein: 27-plex Cell Typing & Functional State marker set (six subpanels, e.g. CD4/CD20/CD8A/CD3E, CD138/HLA-DR/CD11c/CD68, β-catenin/CD31/PTEN/PanCK), plus E-cadherin/Vimentin/αSMA/CD45 segmentation markers. About 465k cells.
- https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma
- **Use**: the cleanest test bed. RNA and protein come from the *same segmented cell*, and cognate genes were added on purpose, so E2 (cognate vs random gene ablation) is well posed.

### P2. 10x Xenium Protein: FFPE human clear cell RCC
377-gene panel + 100 custom genes (protein-coding partners, kidney cell types, RCC genes) + all six protein subpanels. 10x integrated it with WNN. CC BY 4.0.
- https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc
- **Use**: a second same-cell slide for replication and for slide-level hold-out (train on P1, test on P2).

## Tier 1b. Same-section, spot/pixel-level RNA + protein (summary)

### P3. DBiTplus (Enninful et al., *Nat Methods* 2026)
DBiT-seq spatial whole transcriptome, then CODEX (~50 markers) and H&E on **the same section**. Human lymph node, FFPE lymphoma, mouse embryo. The authors integrate the modalities with a MaxFuse-derived method.
- GEO GSE308167 · code/data: Zenodo 10.5281/zenodo.16827046, https://github.com/Janezjz/DBiT-plus
- **Use**: E7 transfer target (different platform, tissue and resolution from Xenium Protein); H2 fallback.

### D17a. Spatial-CITE-seq · D13 Visium CytAssist Gene + Protein · D17b SPOTS · P5 SpatialGlue collection
Details in the sequencing-based section below. These are the standard benchmark datasets for spatial RNA + protein integration methods (SpatialGlue, MultiGATE, PRESENT).

## Tier 2. Adjacent sections

### P4. SPATCH (Ren et al., *Nat Commun* 2025)
Serial sections from colon adenocarcinoma, hepatocellular carcinoma and ovarian cancer, profiled with Stereo-seq v1.3, Visium HD FFPE, CosMx 6K and Xenium 5K. **CODEX (16 proteins) on adjacent sections** and scRNA-seq from the same samples. Manual cell-type annotations and segmentations.
- Web server (visualise / download): http://spatch.pku-genomics.org/ · raw: GSA-Human HRA011129
- **Use**: H1's main data (same tumours in both modalities); H2 stage 2 (register Xenium 5K → CODEX); E8 registration-error sweep; N4 virtual-protein validation against real adjacent CODEX.

---


## Tier 3. Protein-only cohorts (H1 protein targets, H2 stage-3 deployment, E7 transfer targets)

*Carried over from v1. IDs unchanged.*

### D1. Keren et al. 2018: TNBC, MIBI-TOF
*Cell* 174:1373. 41 triple-negative breast cancer patients, 36 markers, cell segmentation and phenotypes, survival. The paper classifies tumours as "compartmentalised" or "mixed", a ready-made spatial phenotype.
- Data: https://www.angelolab.com/mibi-data
- Quick access: `data("kerenSCE", package = "Statial")` (10-image subset, 57,811 cells, 36 proteins).
- **Use**: development and testing; survival association; SpatioMark/Kontextual baselines.

### D2. Jackson, Fischer et al. 2020: breast cancer, IMC
*Nature* 578:615. ~720 images from ~350 patients across two cohorts (Basel, Zurich), 35 markers, single-cell data, clinical data.
- Data: https://zenodo.org/records/3518284
- Bioconductor: `imcdatasets::JacksonFischer_2020_BreastCancer()`
- **Use**: **main benchmark dataset**. Basel and Zurich give a built-in replication cohort.

### D3. Danenberg et al. 2022: METABRIC breast cancer, IMC
*Nat Genet* 54:660. ~700 patients, 37 markers. Linked to METABRIC genomics and long-term survival. Defines 10 tumour microenvironment "structures".
- Data: https://zenodo.org/records/5850952 (full-stack TIFFs, masks, single-cell tables; check for newer versions)
- **Use**: large cohort for outcome association.

### D4. Ali et al. 2020: METABRIC breast cancer, IMC
*Nat Cancer* 1:163. 483 tumours, 37 markers.
- Data: https://idr.openmicroscopy.org/study/idr0076 (CC BY 4.0)
- **Use**: independent replication of D3 findings (image IDs are not shared with D3).

### D5. Damond et al. 2019: pancreas, type 1 diabetes, IMC
*Cell Metab* 29:755. 12 donors (non-diabetic, recent onset, long-duration), 38 markers.
- `imcdatasets::Damond_2019_Pancreas()`; Mendeley: https://data.mendeley.com/datasets/cydmwsfztj
- **Use**: non-cancer case; islet-centric spatial gradients as positive controls.

### D6. Hoch, Schulz et al. 2022: metastatic melanoma, IMC
*Sci Immunol* 7:eabk1692. Protein panel plus chemokine RNA panel on consecutive sections.
- `imcdatasets::HochSchulz_2022_Melanoma()`
- **Use**: chemokine-driven niches (CXCL9/10/13); protein–RNA consistency.

### D7. IMMUcan example data
Small multi-cancer IMC example used in the Bodenmiller-lab IMC workflow book.
- `imcdatasets::IMMUcan_2022_CancerExample()`; workflow book: https://bodenmillergroup.github.io/IMCDataAnalysis/
- **Use**: quick debugging.

### D8. Schürch et al. 2020: colorectal cancer, CODEX
*Cell* 182:1341. 35 advanced CRC patients (Crohn's-like reaction vs. diffuse inflammatory infiltration), 140 TMA spots, 56 markers, ~250k cells. Defined "cellular neighbourhoods".
- Data: https://data.mendeley.com/datasets/mpjzbtfgfr/1 (single-cell CSV with CN labels and clinical data; images on TCIA)
- **Use**: CODEX platform; neighbourhood-dependent ctSVPs (e.g., PD-1⁺ CD4 T cells in CN-9).

### D9. Hickey et al. 2023 / HuBMAP: healthy human intestine, CODEX
*Nature* 619:572. 64 sections from 8 donors, 57 markers, ~2.6M cells.
- HuBMAP portal: https://portal.hubmapconsortium.org/ (search "CODEX intestine")
- Processed tables are also linked from the STELLAR repo: https://github.com/snap-stanford/stellar
- **Use**: scalability tests at 10^6 cells; healthy reference.

### D10. HTAN multiplexed imaging (CyCIF / Orion / CODEX / MIBI)
Human Tumor Atlas Network. Includes Lin et al. 2023 (*Cell*) whole-slide CyCIF of CRC.
- https://humantumoratlas.org/explore (Level 3/4 single-cell tables are open-access for many atlases)
- **Use**: whole-slide scale; extra platforms.

### D11. spora [data]: harmonized spatial proteomics collection
More than 31 harmonized public datasets (IMC, MIBI, CODEX, …) in one format, with a Python loader. Released with VirTues.
- Datasets: https://spora.epfl.ch/datasets.html · loader: https://github.com/bunnelab/spora-io · benchmarks: https://github.com/bunnelab/spora-bench
- **Use**: the easiest way to get many cohorts in one format; also VirTues input.

### D12. Squidpy built-in examples
- `squidpy.datasets.imc()`: IMC breast cancer (Jackson et al. 2020 subset)
- `squidpy.datasets.mibitof()`: MIBI-TOF colorectal (Hartmann et al. 2021)
- **Use**: one-line loading in Python notebooks.

## Tier 1b. Sequencing-based spatial RNA + protein on the same section (details)

### D13. 10x Visium CytAssist Gene + Protein Expression: human tonsil (FFPE)
Whole transcriptome plus a 35-antibody immune panel (31 markers + 4 isotype controls), measured as **antibody UMI counts per spot**. CC BY 4.0.
- Standard: https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard
- With add-on antibodies (ICOS, STING, CD49B, MICA): https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard
- **Use**: the only dataset type where **spot-level ctSVG tools (CELINA, STANCE, C-SIDE, CTSV, spVC) run on protein as designed**. Deconvolve spots with the RNA modality (RCTD + a tonsil scRNA reference), then test the protein features.

### D17. Spatial-CITE-seq and SPOTS
- Spatial-CITE-seq (Liu et al. 2023, *Nat Biotechnol*): 189 ADTs (mouse) / 273 ADTs (human tonsil, skin) + whole transcriptome on the same section. GEO **GSE213264**.
- SPOTS (Ben-Chetrit et al. 2023, *Nat Biotechnol* 41:788): mouse spleen, breast tumour. GEO GSE198353 (verify).
- **Use**: higher-plex spatial protein counts.

## Legacy v1. Transcriptomics reference data for the ctSVG benchmarks

*Still used by `experiments/ctsvp_baseline` and for E6/N1 positive controls; not needed for H1/H2.*

### D14. 10x Xenium (source of the benchmark's realistic simulations)
- Human breast (Janesick et al. 2023): https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
- Ovarian cancer, healthy lymph node: browse https://www.10xgenomics.com/datasets (filter: Xenium)
- **Use**: reproduce the Li et al. realistic simulation before porting it to protein.

### D15. Datasets from the STANCE paper
- HER2+ breast cancer ST: https://zenodo.org/records/4751624
- Human kidney cancer: https://data.mendeley.com/datasets/g67bkbnhhg/1
- Mouse olfactory bulb ST: https://www.spatialresearch.org/resources-published-datasets/doi-10-1126science-aaf2403/

### D16. Datasets from the ctSVG paper (Zhuang et al. 2025)
Visium HD: [human CRC](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc), [human lung cancer](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-lung-cancer-if), [human pancreas](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-pancreas), [mouse brain](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he), [mouse embryo](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-embryo), [mouse intestine](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine), [mouse kidney](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-kidney).
Xenium: [human lung cancer](https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard), [human pancreas](https://www.10xgenomics.com/datasets/ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard), [mouse colon](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-colon-with-xenium-multimodal-cell-segmentation-1-standard).

---

### D18. Anchor benchmark data (Li, Ge & Jiang 2025)
Packaged realistic simulations, preprocessed Xenium-derived inputs, and all detection results: https://zenodo.org/records/21405877. Used by `run/realistic.R` in [weiqi-grace-li/ctsvg-benchmark](https://github.com/weiqi-grace-li/ctsvg-benchmark). The repo also bundles a breast-tumour simulation and the HER2+ ST real dataset.

---
