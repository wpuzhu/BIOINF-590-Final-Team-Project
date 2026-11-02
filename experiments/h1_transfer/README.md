# H1: cross-modal transfer benchmark (planned)

Plan: [`docs/plans/H1_transfer_benchmark.md`](../../docs/plans/H1_transfer_benchmark.md).

To do (Week 1 go/no-go):
- [ ] Install Nicheformer and Spatium; embed 1,000 SPATCH colon cells (Xenium and CODEX) with each.
- [ ] Marker → HGNC gene map (MaxFuse antibody table, checked by hand); keep unmapped markers as new tokens.
- [ ] Patient split files shared with `ablation/` and `h2_privileged/`.
- [ ] Pretraining-leakage check against Nicheformer's corpus list and Spatium's study list.
