# Ablation + UQ code

Design: [`../../05_experiment_design/`](../../05_experiment_design/) · UQ: [`../../06_uncertainty_quantification/`](../../06_uncertainty_quantification/) · Results: [`../../08_results/simulation/`](../../08_results/simulation/)

| File | Contents |
|---|---|
| `xmodal_ablation.py` | `PairedSpatial`, `from_anndata`; `Features` (cell type, cognate, RNA PCs, neighbourhood RNA; frozen PCA; cell-type interactions); `evaluate`; experiments `variance_ladder` (E1), `gene_ablation` (E2), `radius_sweep` (E3), `context_controls` (E4), `region_holdout` (E5), `residual_map` (E6), `pairing_degradation` (E8) |
| `uq.py` | `spatial_blocks`, `unit_folds`, `fold_scores`, `cluster_bootstrap_ci`, `signflip_test`, `group_permutation_importance`, `toroidal_shift_test`, `mondrian_conformal`, `weighted_conformal_coverage`, `shift_decomposition`, `holm`, `bh` |
| `run_simulation_demo.py` | simulated paired tissue with known information sources; runs E1–E6, E8 |
| `run_uq_demo.py` | U1–U5 on simulation (8 samples) |
| `make_uq_figure.py` | 4-panel UQ figure |

## Use on real data

```python
import anndata as ad
from xmodal_ablation import from_anndata, Features, variance_ladder, context_controls
from uq import spatial_blocks, unit_folds, fold_scores, cluster_bootstrap_ci, signflip_test
a = ad.read_h5ad("xenium_protein_rcc.h5ad")   # X = log1p RNA; obsm["protein"] DataFrame; obsm["spatial"] in um
a.obs["block"] = spatial_blocks(a.obsm["spatial"], a.obs["sample"].values, block_um=300)
d = from_anndata(a, protein_key="protein", celltype_key="cell_type", sample_key="block",
                 cognate={"CD3E": ["CD3E"], "CD8A": ["CD8A"], "HLA-DR": ["HLA-DRA"], "PanCK": ["KRT8", "KRT18"]})
F = Features(d, n_pcs=50)
variance_ladder(d, F).to_csv("E1.csv")
```
- **Large slides (~500k cells):** subsample per cell type (e.g. ≤ 50k cells) for E2/E4. The radius graph is built per sample.
- **Cognate map:** start from MaxFuse's `docs/protein_gene_conversion.csv` (340 rows) and check by hand; record the mapping type (1:1, 1:many, isoform, PTM, none; see A2).
- **Design choices:** the predictor (ridge per protein) is fixed across ablations, so differences come from the inputs. PCA/scalers are fitted once on intact RNA per gene mask and reused for permuted inputs. Re-run key contrasts with a non-linear model before concluding.
