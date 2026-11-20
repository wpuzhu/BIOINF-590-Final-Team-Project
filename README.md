# BIOINF 590 Final Team Project

**Informative vs transferable information between spatial transcriptomics and spatial proteomics**

Team: Puzhu Wang · Meng-Jia (Justin) Lian · Nam Nguyen-Hoang

Spatial RNA and multiplexed protein imaging measure two views of the same tissue. We separate two questions that are often merged:
- **Informative**: does adding RNA improve a protein-based prediction on the same cells?
- **Transferable**: does a protein-only model trained with RNA's help do better on new patients that have no RNA?

The project runs in two stages, with patient-level uncertainty quantification throughout.

| Stage | Question | Models / mechanisms | Folder |
|---|---|---|---|
| **1 Cross-modal pretraining** | Does RNA pretraining help protein tasks (and vice versa)? If so, what carries over? | Nicheformer ↔ Spatium in a 2 × 2 with random-init, PCA, identity (A1), mapping (A2), overlap (A3) arms | [`stage1_cross_modal_pretraining/`](stage1_cross_modal_pretraining/) |
| **2 RNA as privileged information** | Can RNA seen only during training improve a protein + spatial-graph student? Which part of RNA's information makes it across? | none vs distillation vs contrastive alignment vs RNA reconstruction, with controls and an input-level ablation ladder | [`stage2_privileged_transfer/`](stage2_privileged_transfer/) |

**Start here:** [`docs/concept_and_design.md`](docs/concept_and_design.md) (definitions, designs, ablation ladder, UQ per experiment) → [`docs/related_work.md`](docs/related_work.md) → [`datasets/README.md`](datasets/README.md).

## Layout

```
docs/                          concept & design, related work, glossary, future ideas (N2–N7)
  background/                  v2 groundwork: literature review (with input formats), assumptions A1–A8,
                               information-map experiments E1–E8, UQ notes U1–U5
datasets/                      public datasets by stage (pairing, scale, leakage), catalog, registry.yaml, download scripts
stage1_cross_modal_pretraining/
  design/ simulation/ results/ pipeline/
stage2_privileged_transfer/
  design/ simulation/ results/ pipeline/
shared/                        data schema + harmonisers + registration · paired spatial simulator · numpy MLP ·
                               UQ library (hierarchical bootstrap, sign-flip, TOST, conformal) · information map
cluster/                       Great Lakes: conda env, env.sh, make_jobs.py, SLURM array templates, CPU smoke test
```

## Simulation results so far

- **Stage 1** (testbed, 5 labels/class): in-modality RNA pretraining +0.076 macro-F1 [0.066, 0.088]. Cross-modal ≈ 0 in both directions (RNA → protein +0.006 [−0.006, 0.018]; protein → RNA −0.010 to −0.018). Protein pretraining on 16 markers ≈ 0 even in-modality. → [results](stage1_cross_modal_pretraining/results/simulation_ablation_testbed/README.md)
- **Stage 2** (4 known-truth scenarios): RNA that is informative but RNA-only (+0.32 with RNA) transfers ~0. Transferable RNA keeps 24–36% of its gain via distillation (+0.04 to +0.06). Self-training and label smoothing give 0 or less, so the gain is RNA's. A useless teacher hurts (−0.036). Graded nulls, pairing jitter and graph removal tell cell-intrinsic from niche-level transfer. → [results](stage2_privileged_transfer/results/simulation/README.md)

## Reproduce

```bash
pip install numpy pandas scipy scikit-learn matplotlib pyyaml joblib
bash cluster/smoke/run_smoke.sh                                        # both pipelines end-to-end on simulated data (CPU, ~3 min)
cd stage2_privileged_transfer/simulation && python run_stage2_sim.py --experiment main && \
   python run_stage2_sim.py --experiment ablations --seeds 2 && python analyze_stage2.py     # ~40 min on 2 cores
cd stage1_cross_modal_pretraining/simulation/ablation_testbed && python summarize_2x2.py   # after run.py / run_reverse.py
```

Real data on UMich Great Lakes: [`cluster/README.md`](cluster/README.md).
