# H2: RNA as privileged information for spatial proteomics (summary)

*Summary of Justin's "H2 Experiment Plan" (Oct 23, 2026, 9 pp.), shared in the BIOINF 590 Project chat. The original PDF is the reference.*

**Hypothesis.** A protein-only model trained with RNA available **only during training** classifies cells and niches better on new protein-only slides than the same model trained without RNA. This is the generalized-distillation / learning-using-privileged-information setting (Lopez-Paz et al. 2016).

**Setup.** Unit = segmented cell with its protein vector, (x, y) and a spatial kNN graph (k = 10–15). Teacher = graph encoder on [protein, RNA]. Student = the same encoder on protein only, trained on labels + distillation toward the frozen teacher. Test = student only, on unseen slides, no RNA. Encoder: 2-layer GraphSAGE or GAT (PyG), hidden 128. Protein arcsinh (cofactor 5) + per-slide z-score; RNA log1p → 50 PCs or scVI latent.

**Student loss.** λ·CE(y, p_s) + (1−λ)·T²·KL(p_t^T ‖ p_s^T) + β·‖g(h_s) − h_t‖², with T = 2–4, tuned on validation patients only.

**Arms.**

| Arm | Teacher input | Student input | Question |
|---|---|---|---|
| A | none | protein | baseline |
| **E** | protein + RNA | protein | main test: does RNA help? |
| **F** | protein only (self-distillation) | protein | does distillation alone help? E must beat F |
| **G** | protein + RNA shuffled across cells | protein | is the gain from RNA content? G should match F |
| H | protein + RNA with registration jitter (0/20/50/100 µm) | protein | how much pairing error can H2 tolerate? |
| I | protein + RNA | protein, markers dropped to 100/50/25% | does RNA help more with small panels? |
| U | protein + RNA | protein + RNA | upper bound (needs RNA at test) |

**Data stages.** 0 simulation (scaffold `h2_scaffold.py`, not yet in this repo) → 1 Xenium Protein demos (same cell; e.g. FFPE RCC) → 2 SPATCH (adjacent sections; register Xenium 5K to CODEX with an affine fit on DAPI + STalign/PASTE2) → 3 an unpaired public CODEX cohort (deployment). Fallbacks: DBiTplus, SPOTS/Visium-protein (Zenodo 10.5281/zenodo.10362607). Split by patient, never by cell.

**Evaluation.** Macro-F1 (main), balanced accuracy, per-class F1 for rare types; niches ARI/NMI; ECE; label-efficiency curves at 1/5/10/25/100% labels; 5 seeds; one-sided paired Wilcoxon E vs F; bootstrap CIs over patients; Holm correction.

**Pre-registered success.** E − F ≥ 0.02 macro-F1 at 5% labels with the CI excluding 0 on Xenium Protein; the G − F CI includes 0; same sign on SPATCH; the student beats arm A on the unpaired CODEX cohort.

**Simulation result (stage 0, from the plan).** With informative RNA, E beats F by +0.05 to +0.07 macro-F1 at 1–5% labels, while G stays at F's level. With null RNA there is no gain, and slight negative transfer at 1% labels.

**How v2 extends this.** Arm G = our E4e; arm H = our E8; arm I relates to E2. The information map says which fine-grained types and proteins RNA can inform, which is where the E–F gap should be largest. Labels from RNA would make the teacher win trivially, so the main claim uses protein + morphology labels.
