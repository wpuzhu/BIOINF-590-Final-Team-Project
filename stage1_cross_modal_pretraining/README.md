# Stage 1: cross-modal pretraining transfer (Nicheformer ↔ Spatium)

**Question.** Does a model pretrained on one modality help on the other with few labels? Nicheformer (RNA-pretrained) is tested on protein, Spatium (protein-pretrained) on RNA, against the in-modality positive controls. If yes, *what* carries over: token identity, overall structure, or just "a pretrained network"?

**Framing.** A pre-registered **open test**: "no gain" (90% CI within ±0.02 macro-F1) is a result, not a failure. Full design: [`../docs/concept_and_design.md`](../docs/concept_and_design.md#stage-1-cross-modal-pretraining-transfer-nicheformer--spatium).

| Folder | Contents |
|---|---|
| `design/` | Justin's H1 plan summary, the original H1 idea note, testbed insights |
| `simulation/ablation_testbed/` | laptop-scale testbed: both directions (`run.py` RNA-pretrained, `run_reverse.py` protein-pretrained), A1–A3 ablations, `summarize_2x2.py` |
| `results/simulation_ablation_testbed/` | testbed results ([README](results/simulation_ablation_testbed/README.md)) |
| `pipeline/` | real-data pipeline: `adapters.py` (PCA / Nicheformer / Spatium behind one interface; cross-modal token mapping; arms), `run_stage1.py` (one job), `configs/stage1.yaml`, protein foundation-model notebook template |
| `results/real/` | written by `shared/uq/aggregate.py --stage 1` after the cluster run |

**Simulation so far.** In-modality RNA pretraining helps (+0.076 at 5 labels/class). Cross-modal transfer is ≈ 0 in both directions, and slightly negative for protein → RNA. Protein pretraining on a 16-plex panel gives ≈ 0 even in-modality. On real data, the Spatium in-modality control decides whether its cross-modal result can be interpreted at all.

**Run.** `python cluster/make_jobs.py --stage 1 --out jobs_stage1.tsv`, then the SLURM array in [`../cluster/README.md`](../cluster/README.md). Lines marked `VERIFY` in `pipeline/adapters.py` must be checked against the installed Nicheformer / Spatium versions.
