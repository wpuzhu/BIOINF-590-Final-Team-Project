# U4 · Model (epistemic) uncertainty

| Source | How to quantify | Where |
|---|---|---|
| Random seed / initialisation | ≥ 5 seeds; report mean ± SD; average within unit before unit-level tests | H1, H2, DGAT, SpatialGlue attention |
| Shuffle randomness in ablations | ≥ 5 shuffles (E4) or random gene sets (E2); average per unit | E2, E4, E8 |
| Model class | repeat key contrasts with ridge, a GNN (H2 encoder) and DGAT; conclusions must agree in sign | E1–E5 |
| Parameter uncertainty (linear probe) | Bayesian ridge posterior or bootstrap of coefficients; CIs on cell-type interaction terms | E1, E4 interpretation |
| Deep ensembles / MC dropout | spread of predictions across 5 models → per-cell epistemic SD; compare with conformal widths | N4 virtual protein |
| Hyper-parameter sensitivity | report results over a small grid (radius, n_pcs, k); flag conclusions that flip | E1–E4 |

**Epistemic vs aleatoric.** Ensemble spread (epistemic) should grow on out-of-distribution regions; residual variance (aleatoric) reflects measurement noise and RNA-independent protein. A protein with high aleatoric but low epistemic uncertainty is simply **not informed by RNA**: more data won't help. High epistemic uncertainty on a region points to a support or covariate shift.
