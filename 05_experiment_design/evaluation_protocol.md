# Evaluation protocol (applies to every experiment)

## 1. Units and splits
- **Split unit = patient or slide.** With one slide, use **spatial blocks** whose side is ≥ 2–3 × the largest neighbourhood radius tested (E3 goes to 100 µm, so blocks ≥ 250–300 µm). Use `uq.spatial_blocks(coords, sample, block_um)`.
- Hyper-parameters (ridge α is chosen internally by `RidgeCV` on training folds; GNN settings) are chosen on **training units only**. For H1/H2, a separate validation unit chooses settings.
- Never let test-unit **protein** enter features. Test-unit **RNA** may enter neighbourhood features of other test cells, because RNA is an input.

## 2. Normalisation
| Modality | Transform | Fit on |
|---|---|---|
| RNA (imaging) | size-factor → log1p; PCA (30–50) for "all RNA" | training units (final runs) |
| RNA (sequencing spots) | as above; consider SCTransform | training units |
| Protein (imaging) | arcsinh(x / c), per-marker cofactor c; z-score per sample | per sample (no target leakage) |
| Protein (ADT) | CLR per cell, or totalVI background-corrected | per sample |

## 3. Metrics
| Metric | When | Why |
|---|---|---|
| R² (held-out) | primary for E1–E4 | variance explained; matches "information" |
| Pearson / Spearman r | always alongside R² | scale-free; robust when held-out variance is small or offsets differ (E5, E7) |
| MSE | shift decomposition | additive, decomposable |
| Conformal coverage and interval width | E5, E7, N4 | distribution-free check of whether uncertainty transfers |
| Macro-F1, balanced accuracy, ECE | H1/H2 classification | per H1/H2 plans |

**Always per cell type.** Pooled metrics mostly measure cell-type separation.

## 4. Minimum sample sizes
- ≥ 30 test cells per (unit, protein, cell type) for a fold score; ≥ 20 calibration cells per cell type for Mondrian conformal.
- ≥ 5 units for any CI. With n units, the exact sign-flip test's smallest p is 2⁻ⁿ (n = 8 gives 0.0039), which cannot survive Holm over many tests. Pre-register a **small set of primary contrasts** (Holm) and treat the rest as exploratory (BH).

## 5. Reporting template (one row per protein × cell type × contrast)
```
protein | cell_type | contrast | mean_dR2 | ci_lo | ci_hi | p_signflip | q_BH | n_units | n_cells | notes
```
plus one figure per experiment and the console log. See `08_results/simulation/` for examples.

## 6. Checklist before believing a result
- [ ] negative control (E4e, RNA shuffled) collapses to the cell-type baseline
- [ ] random-gene-set null for any gene-group claim
- [ ] effect present in ≥ 2 datasets or ≥ 2 platforms (E7)
- [ ] spillover check for neighbourhood effects on lineage markers (A7)
- [ ] non-linear predictor (DGAT / GNN) agrees in direction
