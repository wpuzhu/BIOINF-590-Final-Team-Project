# 05 · Experiment design

**Principle.** Keep **one fixed predictor** and change **only its inputs**. If held-out accuracy drops when a piece of input is removed, shuffled or moved, that piece carried information. If the accuracy holds on a region or dataset the model never saw, the information **transfers**.

```
paired cells ─┬─ RNA ──► [ablate: genes | radius | shuffles | jitter | regions] ─► features ─► fixed predictor ─► protein_j
              └─ protein ─────────────────────────────────────────────────────────────────────────────────────► target
evaluation: held-out patients / slides / spatial blocks · per protein × cell type · with UQ (06_uncertainty_quantification)
```

## The experiments

| ID | Name | Measures | Status |
|---|---|---|---|
| [E1](E1_variance_ladder.md) | variance ladder | how much each information source adds | code + sim ✔ |
| [E2](E2_gene_ablation.md) | gene-group ablation / conditional permutation | cognate vs programme vs diffuse information | code + sim ✔ |
| [E3](E3_radius_sweep.md) | radius sweep | spatial scale of informative context | code + sim ✔ |
| [E4](E4_context_controls.md) | context controls | arrangement vs composition vs leak | code + sim ✔ |
| [E5](E5_region_holdout.md) | region hold-out | spatial out-of-distribution transfer; covariate vs concept shift | code + sim ✔ |
| [E6](E6_residual_maps.md) | residual maps | protein spatial structure not explained by RNA | code + sim ✔ |
| [E7](E7_cross_dataset_transfer.md) | cross-dataset transfer | tissue/platform generality | design (needs real data) |
| [E8](E8_pairing_degradation.md) | pairing degradation | tolerance to registration error | code + sim ✔ |

Cross-cutting pages:
- [`evaluation_protocol.md`](evaluation_protocol.md): splits, normalisation, metrics, reporting template. **Read before running anything on real data.**
- [`transferability_diagnostics.md`](transferability_diagnostics.md): how to tell *informative* from *transferable*, and diagnose *why* something does not transfer.

## How the experiments connect to H1/H2

| This layer | H1 counterpart | H2 counterpart |
|---|---|---|
| E2 cognate vs random genes | A1 re-init mapped tokens; A2 random token mapping | arm I (marker drop) |
| E4e RNA shuffled | random-init control | arm G (shuffled-RNA teacher) |
| E8 jitter | – | arm H (registration jitter) |
| E5/E7 hold-outs | C1–C4 cross-modal transfer | stage 3 (unpaired CODEX cohort) |

**Pre-registered prediction:** H1/H2 gains concentrate in proteins and cell types where E1 shows high RNA-explainable variance *and* E5/E7 show it holds out of distribution.
