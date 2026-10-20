# Datasets

We do **not** store raw data in git. Files are large and each has its own license. This folder holds:
- `datasets_catalog.csv`: every dataset with its link, technology, size, and suggested role;
- `download_scripts/`: helpers that fetch the data into `datasets/raw/` (git-ignored).

| Script | What it fetches |
|---|---|
| `download_imc_bioconductor.R` | Jackson 2020, Damond 2019, Hoch/Schulz 2022, IMMUcan via **imcdatasets**; Keren 2018 via **Statial** (`kerenSCE`) |
| `download_squidpy_examples.py` | Squidpy built-in IMC and MIBI-TOF AnnData objects |
| `download_zenodo_mendeley.sh` | Zenodo/Mendeley deposits (Jackson raw, Danenberg METABRIC, STANCE data) |
| `download_10x.md` | step-by-step for 10x Visium CytAssist protein and Xenium datasets (10x requires a browser form) |

> Always check the license on each landing page before redistributing. Some spora datasets and HTAN controlled-access tiers have restrictions.

---

## A. Multiplexed protein imaging (primary validation data)

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

## B. Sequencing-based spatial protein (counts, which suits count-based ctSVG models)

### D13. 10x Visium CytAssist Gene + Protein Expression: human tonsil (FFPE)
Whole transcriptome plus a 35-antibody immune panel (31 markers + 4 isotype controls), measured as **antibody UMI counts per spot**. CC BY 4.0.
- Standard: https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard
- With add-on antibodies (ICOS, STING, CD49B, MICA): https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard
- **Use**: the only dataset type where **spot-level ctSVG tools (CELINA, STANCE, C-SIDE, CTSV, spVC) run on protein as designed**. Deconvolve spots with the RNA modality (RCTD + a tonsil scRNA reference), then test the protein features.

### D17. Spatial-CITE-seq and SPOTS
- Spatial-CITE-seq (Liu et al. 2023, *Nat Biotechnol*): ~200–300 ADTs + transcriptome, human tonsil/skin. GEO GSE213264 (verify the accession on the paper page).
- SPOTS (Ben-Chetrit et al. 2023, *Nat Methods*): mouse spleen, breast tumour. GEO GSE198353 (verify).
- **Use**: higher-plex spatial protein counts.

## C. Transcriptomics reference data used by the ctSVG benchmarks (method positive controls)

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

## Recommended minimal set for the project
1. **Develop** on D1 (`kerenSCE`) and D7: small, quick.
2. **Benchmark** on D2 (two cohorts) and D8 (CODEX).
3. **Replicate / outcome** on D3 + D4.
4. **Cross-modality sanity check** on D13.
5. **Method positive control** (does our pipeline reproduce the known ctSVGs?) on D15 / D14.
