# Input-level ablation: which RNA information informs and transfers to spatial protein

The design is described in [`docs/02_ablation_experiments.md`](../../docs/02_ablation_experiments.md).

| File | What it is |
|---|---|
| `xmodal_ablation.py` | the framework: `PairedSpatial` container, `from_anndata`, feature builder with gene / radius / shuffle ablations, sample- or block-wise CV, and experiments E1–E6 (`variance_ladder`, `gene_ablation`, `radius_sweep`, `context_controls`, `region_holdout`, `residual_map`) |
| `run_simulation_demo.py` | self-check on simulated paired data where each protein's information source is known |
| `results_sim/` | CSVs (one per experiment) and `ablation_overview.png` from the self-check |

## Run

```bash
pip install numpy pandas scipy scikit-learn matplotlib
python run_simulation_demo.py --out results_sim      # ~35 s on a laptop
```

## Use on real data

```python
import anndata as ad
from xmodal_ablation import from_anndata, Features, variance_ladder, gene_ablation, context_controls
rna = ad.read_h5ad("xenium_protein_rcc.h5ad")          # X = log1p RNA, obsm["protein"] = DataFrame, obsm["spatial"] in um
d = from_anndata(rna, protein_key="protein", celltype_key="cell_type", sample_key="tile",
                 cognate={"CD3E": ["CD3E"], "CD8A": ["CD8A"], "HLA-DR": ["HLA-DRA"], "PanCK": ["KRT8", "KRT18"]})
F = Features(d, n_pcs=50)
variance_ladder(d, F).to_csv("E1.csv")
```

- **One slide?** Set `sample_key` to spatial tiles (e.g. 4–8 strips wider than the largest radius), or leave a single sample and the code falls back to spatial blocks.
- **Large slides (~500k cells):** subsample per cell type for E1–E4 (e.g. 50k cells). Build the radius graph on all cells, then index into it.
- **Cognate map:** start from MaxFuse's antibody-to-gene table and check each entry by hand (as planned in H1).

## Design choices

- The ridge predictor (one per protein, with own-cell-type interactions on the cognate and neighbourhood blocks) stays fixed across all ablations, so differences come from the inputs. Re-run the key contrasts with a nonlinear model (H2's GNN, or DGAT) before drawing conclusions.
- `context_controls` condition (c) swaps each cell's neighbourhood RNA for that of random cells *of the same type in the same sample*. It keeps local composition and breaks cell-state arrangement.
- `morans_i_within_type` is a lightweight stand-in for ctSVP testing. For final results, run the residuals through `experiments/ctsvp_baseline`.
