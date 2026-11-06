# H1 · Transfer learning between SPT and spatial proteomics foundation models

| | |
|---|---|
| **Status** | frozen-mode benchmark run by Justin (5 Nov); random-init controls pending · ablation testbed with insights: [`ablation_insights.md`](ablation_insights.md) · plan: [`justin_plan_summary.md`](justin_plan_summary.md) |
| **Question** | Does pretraining on one modality (Nicheformer: RNA; Spatium: protein) beat from-scratch and best-native methods on the other modality? |
| **Assumptions** | A2 (token identity), A3, A5, A8 (pretraining leakage) |
| **Data** | SPATCH (CODEX + Xenium 5K/CosMx 6K), Schürch CRC CODEX |
| **Methods** | Nicheformer, Spatium, random-init twins; MAPS, scANVI, CellCharter, PCA baselines; Novae fallback ([D_foundation_models](../../01_literature_review/algorithms/D_foundation_models.md)) |
| **Experiments** | contrasts C1–C4; token ablations A1–A3; tasks T1 cell typing, T2 niches |
| **UQ** | 5 seeds × label levels; patient-level bootstrap CIs on macro-F1 differences; Holm correction across contrasts; leakage-stratified reporting |
| **Cost** | 20–40 GPU-hours (Justin's estimate) |

## Motivation
Both modalities share the data shape (cells, coordinates, graphs) and most markers have a gene partner. HEIST showed SPT pretraining can *run* on protein; nobody has shown it *wins* against proper controls.

## Design
See the plan summary. v2 addition: use the I1 information map to pre-register which cell types should gain most (high RNA-explainable, transferable variance) and test that prediction.

## Expected outcomes and how to read them
C1 > 0 and C2 ≥ 0 means transfer helps. If the gain disappears under A2 (random token mapping), transfer runs through feature identity. If gains concentrate where I1 predicts, the mechanism is the shared information.

## Risks
Models not running on the other modality (Week-1 go/no-go); SPATCH has only 16 markers; pretraining leakage.
