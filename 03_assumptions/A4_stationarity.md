# A4. Stationarity: the RNA → protein mapping is the same in training and test units

**Statement.** P(protein | RNA inputs) is the same across the regions, patients, tissues or platforms we transfer between. This is the assumption that *transferability* tests.

Two different things can change between source and target:

| Shift | What changes | Example | Fixable by reweighting? |
|---|---|---|---|
| **Covariate shift** | P(RNA inputs): which cells and contexts are present | tumour core has few T cells; target tissue has other cell types | yes: importance weighting, or per-cell-type models |
| **Concept shift** | P(protein \| RNA inputs): the mapping itself | hypoxia lowers translation in the core; a different antibody clone; a different platform's protein scale | no: needs target data or a model of the mechanism |
| **Label/prior shift** | P(cell type) with P(x \| type) fixed | different composition across patients | yes |
| **Support shift** | target contains inputs never seen in training | a cell type absent from the training slide | no: extrapolation; flag it |

**How we test it.** [`../05_experiment_design/transferability_diagnostics.md`](../05_experiment_design/transferability_diagnostics.md): region hold-out (E5), cross-dataset (E7), weighted conformal coverage and the covariate-vs-concept shift decomposition (`07_code/ablation/uq.py::shift_decomposition`).
