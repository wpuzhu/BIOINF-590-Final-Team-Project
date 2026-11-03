# H1: cross-modal transfer benchmark (planned)

Plan: [`04_ideas/H1_transfer_learning/justin_plan_summary.md`](../../04_ideas/H1_transfer_learning/justin_plan_summary.md).

To do (Week 1 go/no-go):
- [ ] Install Nicheformer and Spatium; embed 1,000 SPATCH colon cells (Xenium and CODEX) with each.
- [ ] Marker → HGNC gene map (MaxFuse antibody table, checked by hand); keep unmapped markers as new tokens.
- [ ] Patient split files shared with `07_code/ablation/` and `h2_privileged/`.
- [ ] Pretraining-leakage check against Nicheformer's corpus list and Spatium's study list.

`protein_foundation_models_template.ipynb`: how to load KRONOS and VirTues (and Mesmer/Cellpose for segmentation). Frozen embeddings + a linear probe are the protein-only baselines for H1 and the H2 student.
