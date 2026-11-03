# 07 · Code

| Folder | Status | What |
|---|---|---|
| [`ablation/`](ablation/) | **runnable**, self-checked on simulation | information map (E1–E6, E8) + uncertainty quantification and transferability diagnostics |
| [`h2_privileged/`](h2_privileged/) | planned | H2 teacher → student (add Justin's `h2_scaffold.py` here) |
| [`h1_transfer/`](h1_transfer/) | planned; includes `protein_foundation_models_template.ipynb` (KRONOS, VirTues, Mesmer/Cellpose) | H1 Nicheformer ↔ Spatium; protein-only baselines |

```bash
cd 07_code/ablation
pip install numpy pandas scipy scikit-learn matplotlib
python run_simulation_demo.py      # E1-E6, E8   -> 08_results/simulation/ablation/   (~1 min)
python run_uq_demo.py              # U1-U5       -> 08_results/simulation/uq/         (~2 min)
python make_uq_figure.py           # uq_overview.png
```
