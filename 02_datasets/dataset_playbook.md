# Dataset playbook: what to do with each dataset

For each dataset this page covers: **what it is → how to shape it into our input (`PairedSpatial` / AnnData) → which experiments it supports → what can go wrong → what we expect to transfer and what not.** Experiment IDs (E1–E8) are defined in [`../05_experiment_design/`](../05_experiment_design/); assumptions (A1–A8) in [`../03_assumptions/`](../03_assumptions/).

## Common target format

```
AnnData (one per slide or one concatenated with obs["sample"])
  X                      log1p-normalised RNA (cells/spots x genes)       raw counts kept in layers["counts"]
  obsm["protein"]        DataFrame cells x proteins: arcsinh(x / cofactor), z-scored per sample
  obsm["spatial"]        x, y in micrometres
  obs["cell_type"]       labels (say how they were derived: protein, RNA, both)
  obs["sample"]          patient / slide = the unit of splitting and resampling
  obs["region"]          optional: niche / domain / pathology annotation (for E5)
  uns["cognate"]         dict protein -> [genes]   (start from MaxFuse's protein_gene_conversion.csv)
```
Then `from_anndata(...)` in `07_code/ablation/xmodal_ablation.py`.

---

## Tier 1a: same cell

### P1 / P2: 10x Xenium Protein (RCC, ccRCC)
- **Shape:** cell_feature_matrix (gene + protein features) → split by `feature_types`; centroids from `cells.parquet`. Protein: arcsinh, cofactor per marker (start at the 20th percentile of positive cells). One slide each, so **splitting units are spatial blocks** (250–500 µm) or the two slides against each other.
- **Experiments:** E1–E6 (main), E8 (jitter sweep on exact pairs), H2 stage 1, N2.
- **Pitfalls:** the panel was *designed* so cognate genes are present (25 genes boosted in RCC, 100 custom in ccRCC), which inflates the cognate share in E2 relative to an unbiased panel. Segmentation spillover (A7) can create artificial RNA–protein correlation between neighbours, so check with a nucleus-only protein readout if available.
- **Expect:** lineage proteins (CD3E, CD8A, CD20, CD68, PanCK) transfer from the cognate gene. Activation/state proteins (HLA-DR, PD-L1-like, Ki-67) depend on programmes and context. Phospho/stability-regulated proteins (PTEN, β-catenin localisation) transfer poorly.

### P6: CosMx Human Breast Multiomic (WTX + 64 proteins, same slide)
- **Shape:** AtoMx export (expression matrices per modality + cell metadata with centroids). Whole transcriptome, so reduce with PCA (50) for "all RNA" features but keep cognate genes raw.
- **Experiments:** E2 with *unbiased* gene availability; N3 (minimal gene panel for 64 proteins); E1–E4 replication on a second platform and tissue.
- **Pitfalls:** WTX imaging has low counts per gene, so cognate genes may be sparse. Use neighbourhood-smoothed cognate expression as an extra input and report it separately.
- **Expect:** a larger share of protein variance from *programmes* than on Xenium, because the cognate gene is noisier.

## Tier 1b: same section, spot or pixel RNA

### P3: DBiTplus (lymph node, lymphoma)
- **Shape:** two resolutions. CODEX gives single cells; DBiT gives spots (~25–50 µm pixels). Either aggregate protein to spots (mean over cells in the spot) for a spot-level analysis, or use the authors' deconvolved single-cell transcriptomes (MaxFuse-derived pipeline in `integration/`).
- **Experiments:** E7 target (train on Xenium Protein, test here on shared genes × proteins); E1–E4 at spot level.
- **Pitfalls:** spot RNA is a mixture of cells, so "own RNA" already contains neighbourhood information, which blurs E3/E4. State the resolution in every result.
- **Expect:** high spot-level R² for lineage proteins (mixtures track composition), lower for per-cell states.

### D17a: Spatial-CITE-seq (GSE213264; tonsil, skin; 273 ADTs)
- **Shape:** ADT counts and mRNA counts per pixel. ADT: CLR or totalVI background-aware normalisation; RNA: log-norm.
- **Experiments:** E2 with a very large protein panel (many ADTs have no strong cognate signal); E7 target.
- **Pitfalls:** ADT background (non-specific binding) creates a shared floor across proteins. Include **isotype controls** as negative-control "proteins": RNA should explain ~0 of them.
- **Expect:** surface lineage markers transfer; intracellular or low-abundance ADTs barely do.

### D13: Visium CytAssist Gene + Protein (tonsil; 31 markers + 4 isotype controls)
- **Shape:** one `filtered_feature_bc_matrix.h5` with both "Gene Expression" and "Antibody Capture" features.
- **Experiments:** quick spot-level E1–E4; the isotype controls are built-in negative controls for the UQ calibration checks.

### D17b / P5: SPOTS (mouse spleen, GSE198353), Stereo-CITE-seq (mouse thymus), SpatialGlue collection
- **Experiments:** integration-method comparisons (SpatialGlue, MultiGATE, PRESENT were all shown on these), so we can reproduce published numbers before running our ablations. Mouse vs human also tests whether the information map transfers across species.

## Tier 2: adjacent sections

### P4: SPATCH (COAD, HCC, OV; Xenium 5K, CosMx 6K, Visium HD, Stereo-seq + 16-marker CODEX + scRNA-seq)
- **Shape:** per tumour, an SPT object and a CODEX object on *different* sections. Pair via (a) registration (STalign/PASTE2) + radius matching, and (b) MaxFuse matching. Keep both pairings.
- **Experiments:** H1 main data; H2 stage 2; **E8** (information versus measured registration error); N4 validation (virtual CODEX protein from Xenium vs real adjacent CODEX).
- **Pitfalls:** only 16 markers. Adjacent sections differ in cell content, so cell-level pairs are approximate; analyse at neighbourhood level (20–50 µm) and say so.
- **Expect:** neighbourhood-level and composition information transfers; cell-state information mostly does not survive the pairing error.

### D6: Hoch & Schulz melanoma (IMC protein + RNAscope chemokines on consecutive sections)
- **Use:** a small, biologically sharp test of E3/E4: do chemokine-RNA-rich neighbourhoods predict immune protein states in nearby cells?

## Tier 3: protein only (deployment targets)

Schürch CRC CODEX (D8, 56 markers), Jackson/Fischer and Danenberg IMC breast (D2/D3), Keren MIBI TNBC (D1), Hickey intestine CODEX (D9), HTAN CyCIF (D10), spora harmonised collection (D11).
- **Use:** H1 protein targets; H2 stage 3 ("does the RNA-taught student work on a cohort that never had RNA?"); E7 when only the *protein* side matters (e.g. does a niche definition learned with RNA help here?).
- **Pitfalls:** panel mismatch, so restrict to shared markers and report how many. Antibody clones differ across cohorts, so the same marker name is not the same measurement (A2).

## Tier 4: RNA only

Public Xenium/CosMx/Visium HD sets, and SPATCH SPT without CODEX.
- **Use:** H1 SPT side; N4 deployment of virtual protein (only for proteins the information map marks as RNA-explainable *and* transferable).

---

## Minimal path for a 4-week project

1. **Week 1:** P1 (Xenium Protein RCC) → E1–E4 + UQ. This is the information map.
2. **Week 2:** P2 or P6 → replicate; E5 region hold-out; E6 residuals.
3. **Week 3:** E7 transfer P1 → P3/D17a on the shared subset; E8 jitter sweep; SPATCH pairing for H2 stage 2.
4. **Week 4:** H1/H2 arms where the map predicts gains; write-up.
