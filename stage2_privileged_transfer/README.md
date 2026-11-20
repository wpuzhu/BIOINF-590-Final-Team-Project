# Stage 2: RNA as privileged information

**Question.** RNA is available only for training. Can a **protein + spatial-graph student** become better on new patients that have no RNA? Through which mechanism, and which part of RNA's information makes it across?

**Design** (details and the reasoning on whether distillation makes sense: [`../docs/concept_and_design.md`](../docs/concept_and_design.md#stage-2-rna-as-privileged-information)):
- **Mechanisms compared on one student and identical splits**: none · KD soft labels · KD + features · contrastive alignment (InfoNCE; CCA in simulation) · auxiliary RNA reconstruction · RNA-specific reconstruction.
- **Controls**: self-training with a protein-only teacher (unlabelled-pool effect) · label smoothing (regularisation) · permuted dark knowledge.
- **Input-level ladder**: teacher-RNA graded nulls (within niche / type, type average, all) · neighbourhood-only RNA · pairing jitter 15–240 µm · graph off · panel size · dropped subtype markers · label fraction.
- **UQ**: patient → seed hierarchical bootstrap, sign-flip tests, TOST, variance split, conformal set size / coverage, ECE.

| Folder | Contents |
|---|---|
| `design/` | Justin's H2 plan summary, the original H2 idea note |
| `simulation/` | `run_stage2_sim.py` (4 scenarios × 10 arms; 18 ablation settings), `analyze_stage2.py` |
| `results/simulation/` | tables, figure and interpretation ([README](results/simulation/README.md)) |
| `pipeline/` | `train_stage2.py` (torch + PyG GraphSAGE student/teacher, all mechanisms; `--backend numpy` CPU smoke test; `--dry_run`), `models_torch.py`, `prepare_adjacent.py` (SPATCH neighbourhood RNA after registration), `configs/stage2.yaml` |

**Simulation in one line.** The designs separate the four truths. RNA that is informative but RNA-only transfers ~3%; transferable RNA keeps 24–36% of its informative gain through KD. Self-training and label smoothing give nothing. A useless teacher *hurts*, so distil only when Δ_info > 0. Graded nulls, jitter and graph removal identify cell-intrinsic vs niche-level routes, and niche-level transfer survives ~30–60 µm of pairing error.
