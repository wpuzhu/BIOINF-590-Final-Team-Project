# H2 · RNA as privileged information for protein-only models

| | |
|---|---|
| **Status** | planned (Justin's plan: [`justin_plan_summary.md`](justin_plan_summary.md)); stage-0 simulation exists (`h2_scaffold.py`, to be added) |
| **Question** | Does a protein-only student trained with an RNA+protein teacher beat the same student trained without RNA, on new protein-only slides? |
| **Assumptions** | A1 (pairing), A3 (label circularity), A4, A5 |
| **Data** | stage 1 Xenium Protein; stage 2 SPATCH; stage 3 unpaired CODEX cohort |
| **Methods** | GraphSAGE/GAT teacher–student distillation (Hinton KD + feature distillation) |
| **Experiments** | arms A, E, F, G, H, I, U; label-efficiency curves |
| **UQ** | paired one-sided Wilcoxon E vs F over patients × seeds; bootstrap CIs over patients; Holm across label fractions; ECE for calibration |
| **Cost** | small GNNs; hours |

## Motivation
RNA is available for a few training slides but not at deployment. Learning-using-privileged-information can transfer it into a protein-only model.

## Design
See the plan summary. v2 additions: (1) arm I and E8 share the jitter code; (2) the I1 map predicts which fine-grained types should gain; (3) Mondrian conformal on student predictions per type, to check that gains are not bought with miscalibration.

## Expected outcomes and how to read them
E > F ≈ G > A supports H2. E ≈ F means distillation alone; E ≈ G means no RNA content was used.

## Risks
Registration error on SPATCH (A1); few patients; teacher overfits RNA noise.
