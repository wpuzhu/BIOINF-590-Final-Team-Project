# Experiments

| Folder | Status | What |
|---|---|---|
| [`ablation/`](ablation/) | **runnable**, self-checked on simulation | Input-level ablation framework (E1–E6): which RNA information predicts protein, at what spatial scale, and whether it holds across regions |
| [`h2_privileged/`](h2_privileged/) | planned (Justin's scaffold to be added) | RNA as privileged information: teacher (protein + RNA) → protein-only student |
| [`h1_transfer/`](h1_transfer/) | planned | Nicheformer ↔ Spatium cross-modal transfer benchmark |
| [`ctsvp_baseline/`](ctsvp_baseline/) | **runnable**, results included (v1) | Zero-shot cell-type-specific spatially variable protein detection on 10 public images; reused for E6 residual maps and idea N1 |

Suggested order: `ablation/` on Xenium Protein first (cheap, gives the information map) → `h2_privileged/` stage 0–1 → `h1_transfer/` Week-1 go/no-go.
