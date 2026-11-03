# 07 · Code

| Folder | Status | What |
|---|---|---|
| [`ablation/`](ablation/) | **runnable**, self-checked on simulation | information map (E1–E6, E8) + uncertainty quantification and transferability diagnostics |
| [`ctsvp_baseline/`](ctsvp_baseline/) | **runnable** (v1), results included | cell-type-specific spatially variable protein detection; reused for E6/N1/N2 |
| [`h2_privileged/`](h2_privileged/) | planned | H2 teacher → student (add Justin's `h2_scaffold.py` here) |
| [`h1_transfer/`](h1_transfer/) | planned | H1 Nicheformer ↔ Spatium |

```bash
cd 07_code/ablation
pip install numpy pandas scipy scikit-learn matplotlib
python run_simulation_demo.py      # E1-E6, E8   -> 08_results/simulation/ablation/   (~1 min)
python run_uq_demo.py              # U1-U5       -> 08_results/simulation/uq/         (~2 min)
python make_uq_figure.py           # uq_overview.png
```
