# Datasets

We do **not** store raw data in git. Files are large and each has its own license. This folder holds:
- `datasets_catalog.csv`: every dataset with link, technology, **pairing** (same cell / same section / adjacent / single modality) and its **v2 role**;
- `download_scripts/`: helpers that fetch data into `02_datasets/raw/` (git-ignored). Paired RNA + protein data: [`download_scripts/download_paired_multiomics.md`](download_scripts/download_paired_multiomics.md).
- [`dataset_playbook.md`](dataset_playbook.md): **what to do with each dataset**. For each one: input preparation, which experiments it supports, what can go wrong, and why its information may or may not transfer.

The v2 question needs RNA and protein from the same tissue, so datasets are ranked by **how tightly RNA and protein are paired**. Tighter pairing gives cleaner ablations; looser pairing is more realistic.

| Tier | Pairing | Datasets | Used for |
|---|---|---|---|
| **1a** | same cell (imaging) | **P1** Xenium Protein RCC, **P2** Xenium Protein ccRCC, **P6** CosMx Breast Multiomic (WTX + 64 proteins) | main ablation datasets (E1–E6); H2 stage 1 |
| **1b** | same section, spot/pixel RNA | **P3** DBiTplus (whole transcriptome + ~50-plex CODEX), **D17a** spatial-CITE-seq (273 ADTs), **D13** Visium CytAssist Gene + Protein, **D17b/P5** SPOTS / Stereo-CITE-seq | spot-level E1–E4; E7 cross-dataset transfer |
| **2** | adjacent sections | **P4** SPATCH (Xenium 5K, CosMx 6K, Visium HD, Stereo-seq + 16-marker CODEX + scRNA-seq), D6 Hoch melanoma (IMC + RNAscope) | H1 main data; H2 stage 2; E8 registration error; N4 validation |
| **3** | protein only | D1–D3, D8–D11 (IMC, MIBI, CODEX, CyCIF; e.g. **D8** Schürch CRC CODEX, 56 markers) | H1 protein targets; H2 stage 3 (does the student transfer?) |
| **4** | RNA only | public Xenium / CosMx / Visium HD sets; SPATCH SPT without CODEX | H1 SPT side; N4 virtual protein |

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

### P6. Bruker CosMx Human Breast Multiomic (same slide)
FFPE human breast profiled on the **same slide** with the CosMx Human Whole Transcriptome panel (~18,000 RNA targets) and the 64-plex Immuno-Oncology protein panel. Bruker also hosts a 64-plex protein-only tonsil set and RNA-only sets (6K human brain, pancreas, colon).
- https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/
- **Use**: the only same-cell dataset with a **whole transcriptome**. E2 can then separate "the cognate gene" from "a gene programme" without the panel-design bias of a 377-gene Xenium panel (where cognate genes were *added on purpose*). The larger protein panel helps N3 (minimal gene panel).

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


## Tier 1b. Sequencing-based spatial RNA + protein on the same section (details)

### D13. 10x Visium CytAssist Gene + Protein Expression: human tonsil (FFPE)
Whole transcriptome plus a 35-antibody immune panel (31 markers + 4 isotype controls), measured as **antibody UMI counts per spot**. CC BY 4.0.
- Standard: https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard
- With add-on antibodies (ICOS, STING, CD49B, MICA): https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard
- **Use**: quick spot-level E1–E4. The 4 isotype controls are built-in negative-control "proteins": RNA should explain none of their variance.

### D17. Spatial-CITE-seq and SPOTS
- Spatial-CITE-seq (Liu et al. 2023, *Nat Biotechnol*): 189 ADTs (mouse) / 273 ADTs (human tonsil, skin) + whole transcriptome on the same section. GEO **GSE213264**.
- SPOTS (Ben-Chetrit et al. 2023, *Nat Biotechnol* 41:788): mouse spleen, breast tumour. GEO **GSE198353** (also the SpatialGlue / MultiGATE benchmark data).
- **Use**: high-plex E2 (many proteins without a strong cognate signal); E7 transfer targets.

## Tier 3. Protein-only cohorts (deployment targets for H1/H2 and E7)

These cohorts never had RNA. They answer "does what we learned with RNA help where only protein exists?" Restrict to markers shared with the training panel, and report how many.

| ID | Dataset | Platform | Size / markers | Role | Access |
|---|---|---|---|---|---|
| D8 | Schürch et al. 2020, colorectal cancer | CODEX | 35 patients, 140 TMA spots, ~250k cells, 56 markers; cellular-neighbourhood labels | **H1 protein target** (tissue-matched to SPATCH colon) | [Mendeley mpjzbtfgfr](https://data.mendeley.com/datasets/mpjzbtfgfr/1) |
| D9 | Hickey et al. 2023, healthy intestine (HuBMAP) | CODEX | 8 donors, ~2.6M cells, 57 markers | H1 stretch target; STELLAR / Spatium example data | [HuBMAP](https://portal.hubmapconsortium.org/), STELLAR repo |
| D2 | Jackson, Fischer et al. 2020, breast cancer | IMC | ~350 patients (Basel + Zurich), 35 markers | H2 stage 3 (two cohorts: deployment + replication) | [Zenodo 3518284](https://zenodo.org/records/3518284), `imcdatasets` |
| D3 | Danenberg et al. 2022, METABRIC breast | IMC | ~700 patients, 37 markers | H2 stage 3 at scale | [Zenodo 5850952](https://zenodo.org/records/5850952) |
| D1 | Keren et al. 2018, TNBC | MIBI-TOF | 41 patients, 36 markers | small MIBI target (platform shift) | [Angelo lab](https://www.angelolab.com/mibi-data); `Statial::kerenSCE` |
| D6 | Hoch, Schulz et al. 2022, metastatic melanoma | IMC + RNAscope (consecutive sections) | protein panel + chemokine RNA panel | small paired-ish test of E3/E4 (chemokine RNA neighbourhoods → immune protein states) | `imcdatasets::HochSchulz_2022_Melanoma()` |
| D10 | HTAN atlases (e.g. Lin et al. 2023 CRC) | CyCIF / Orion / CODEX / MIBI | whole slides, 20–60 markers | scale tests | [HTAN portal](https://humantumoratlas.org/explore) |
| D11 | spora [data] (Bunne lab) | mixed | 31+ harmonised datasets | many cohorts in one format; VirTues input | [spora](https://spora.epfl.ch/datasets.html) |

## Tier 4. RNA only

Public Xenium / CosMx / Visium HD sets (10x and Bruker dataset pages) and the SPT sections of SPATCH. **Use:** the H1 SPT side, and N4 (virtual protein), only for proteins the information map marks as RNA-explainable and transferable.
