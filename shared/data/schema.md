# Harmonised data schema (both stages)

Every dataset is converted once, by `shared/data/harmonise.py`, into one file per dataset with the layout below. All pipeline code reads only this layout, so adding a dataset means writing one loader.

Formats: `.h5ad` (AnnData; used on the cluster) or `.npz` (same keys; used for the simulated dry run where anndata is not installed). `shared/data/spt_io.py` reads both.

## Cells (`obs`)

| Column | Type | Meaning |
|---|---|---|
| `cell_id` | str | unique within dataset |
| `patient` | str | **the unit of resampling and of splitting**; never split a patient across train/test |
| `sample` / `slide` | str | tissue block / physical section |
| `x_um`, `y_um` | float | centroid in micrometres, in the slide's own frame |
| `cell_type_coarse`, `cell_type_fine` | str | labels; `NA` allowed |
| `label_source` | str | `protein_gating`, `rna_clustering`, `author`, `pathologist`: tells you if labels were themselves derived from one modality (a circularity risk in Stage 2) |
| `fold` | int | patient-level fold id (written by `splits.py`) |
| `spatial_block` | int | 250 µm block id for spatial-block resampling |

## Features

| Key | Shape | Meaning |
|---|---|---|
| `X_protein` | cells × markers | arcsinh(x/5) of mean intensity (imaging) or CLR of counts (antibody–oligo); `NaN` if not measured |
| `X_rna` | cells × genes | raw counts (integers); `NaN` rows if not measured on that cell |
| `protein_names`, `gene_names` | lists | HGNC symbols for genes; vendor marker names for proteins |
| `protein_cognate_gene` | list | HGNC symbol of the gene encoding each marker (`""` for PTMs, e.g. pS6, or multi-gene complexes) |

## Pairing (`uns["pairing"]`)

| Field | Values |
|---|---|
| `type` | `same_cell` (Xenium Protein, CosMx multiomic, DBiT-plus, same-section Xenium→PhenoCycler/COMET) · `adjacent_section` (SPATCH) · `unpaired` |
| `registration` | for `adjacent_section`: path to a 3×3 affine per (RNA slide, protein slide) plus its estimated error in µm |
| `rna_dataset` | for `adjacent_section`: the file holding the RNA cells |

Rule used by Stage 2: same-cell → cell-level teacher input; adjacent-section → the teacher sees the **neighbourhood RNA** around the registered location (radius ≥ 2× registration error), never a single "matched" cell; unpaired → no cell-level distillation, only Stage 1-style transfer.
