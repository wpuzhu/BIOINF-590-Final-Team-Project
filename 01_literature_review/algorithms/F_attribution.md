# F. Attribution and interpretability for trained models

**Use when** a model (H1 foundation model, H2 student/teacher, DGAT) works, and you need to know *which input* it used. These are model-level counterparts of the input-level ablations, and they come with their own failure modes.

| Tool | What it attributes | How to use here | Caveat |
|---|---|---|---|
| **Permutation importance** (grouped, conditional) | drop in held-out score when a feature group is permuted | implemented in `07_code/ablation/uq.py::group_permutation_importance`: permute a gene group jointly **within sample × cell type** at test time | marginal permutation of correlated features gives misleadingly large importances; permute *within strata* (conditional) |
| **Drop-and-refit (LOCO)** | drop in score when a group is removed and the model **retrained** | our E2 gene ablation is exactly this | measures *unique* information only: redundant groups each look unimportant |
| **Shapley values over groups** | fair split of the full-model gain among groups | for 4–8 gene groups, exact Shapley = refit on all 2^k subsets with the linear probe | cost grows 2^k; use groups, not genes |
| **Integrated gradients / SHAP** (Captum, `shap`) | per-cell, per-input-feature attributions for neural nets | which genes/markers drive the H2 teacher's predictions | baseline choice matters; noisy per cell, so aggregate per cell type |
| **GNNExplainer / PGExplainer** (PyG) | important edges and node features for a GNN prediction | do H2 students rely on neighbours? Compare with E3/E4 | unstable; report over seeds |
| **Attention weights** (SpatialGlue, MultiGATE, transformers) | where the model attends | modality attention as a descriptive map | attention ≠ importance; confirm with ablation |
| **Linear probes** on frozen embeddings | whether information is linearly decodable | N7: can Spatium embeddings decode cognate RNA? | a probe that is too strong learns the task itself; use ridge with a fixed small capacity and a random-embedding baseline |
| **Model-X knockoffs** (Candès et al. 2018) | feature selection with FDR control | which genes carry *conditional* information about protein, given all others | needs a model of the gene distribution; Gaussian knockoffs on RNA PCs are a practical approximation |

## Rule of thumb

Use **input ablation** (retrain without X) for claims of the form "X is needed". Use **permutation or attribution** (same model, X scrambled) for "the model uses X". Use **Shapley over groups** when groups are redundant (our simulation's IFNG vs T-cell markers). Use **knockoffs** when you need a list of genes with an FDR guarantee.
