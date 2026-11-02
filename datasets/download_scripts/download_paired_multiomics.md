# Downloading paired spatial RNA + protein data (v2)

Put everything under `datasets/raw/<dataset>/`. That folder is git-ignored.

## P1 / P2: 10x Xenium Protein (same cell)

1. Open the dataset page, fill in the short form, then copy the `curl` commands for the **Xenium Output Bundle**:
   - RCC: https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma
   - ccRCC: https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc
2. Unzip into `datasets/raw/xenium_protein_rcc/`.
3. Load (Python). The bundle's `cell_feature_matrix` holds both gene and protein features; split them by `feature_types`:

```python
import scanpy as sc, numpy as np, pandas as pd
ad = sc.read_10x_h5("datasets/raw/xenium_protein_rcc/cell_feature_matrix.h5", gex_only=False)
cells = pd.read_parquet("datasets/raw/xenium_protein_rcc/cells.parquet").set_index("cell_id")
ad = ad[cells.index.astype(str).intersection(ad.obs_names)].copy()
ad.obsm["spatial"] = cells.loc[ad.obs_names, ["x_centroid", "y_centroid"]].values
print(ad.var["feature_types"].value_counts())           # Gene Expression vs Protein Expression (check the exact label)
is_prot = ad.var["feature_types"].str.contains("Protein")
rna  = ad[:, ad.var["feature_types"] == "Gene Expression"].copy()
prot = pd.DataFrame(ad[:, is_prot].X.toarray(), index=ad.obs_names, columns=ad.var_names[is_prot])
sc.pp.normalize_total(rna); sc.pp.log1p(rna)
rna.obsm["protein"] = np.arcsinh(prot / 5)               # tune the cofactor per marker if needed
```

The file names and the protein `feature_types` label above follow Xenium Onboard Analysis v4 outputs. Check `analysis_summary.html` and the bundle listing if they differ. `spatialdata-io`'s `xenium()` reader is an alternative.

Then `experiments/ablation/xmodal_ablation.from_anndata(rna, protein_key="protein", celltype_key=..., sample_key=...)`.

## P3: DBiTplus

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308167 (spot × gene matrices)
- Zenodo 10.5281/zenodo.16827046 and https://github.com/Janezjz/DBiT-plus (processed CODEX tables + analysis code)

## P4: SPATCH

- Browse and download per platform and tumour at http://spatch.pku-genomics.org/ (processed objects, CODEX tables, annotations).
- Raw reads: GSA-Human HRA011129 (FTP/Aspera recommended by the host).
- Start with **colon (COAD)**: Xenium 5K + adjacent CODEX. That is the smallest complete pair, and it matches Schürch CRC CODEX (D8) for H1.

## D17a: Spatial-CITE-seq

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213264 (mouse tissues with 189 ADTs; human tonsil and skin with 273 ADTs)

## D13: Visium CytAssist Gene + Protein (tonsil)

See [`download_10x.md`](download_10x.md).

## P5: SPOTS / Stereo-CITE-seq / Visium-protein (SpatialGlue collection)

- Zenodo 10.5281/zenodo.10362607 (cited in the H2 plan; check its contents before relying on it). SpatialGlue and MultiGATE tutorials also link their preprocessed copies.
