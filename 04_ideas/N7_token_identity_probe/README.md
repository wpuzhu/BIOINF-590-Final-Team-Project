# N7 · Token-identity probe for foundation models

| | |
|---|---|
| **Status** | design only; cheap add-on to H1 |
| **Question** | Do protein-pretrained embeddings (Spatium) linearly encode cognate RNA, and do RNA-pretrained embeddings (Nicheformer) encode protein? |
| **Assumptions** | A2, A8 |
| **Data** | same-cell data (P1/P6) |
| **Methods** | frozen embeddings + ridge probes |
| **Experiments** | probe accuracy vs PCA-of-input and random-init embeddings |
| **UQ** | patient/block bootstrap CIs on probe R²; permutation of targets for a probe null |
| **Cost** | minutes after embedding |

## Motivation
It shows whether pretraining stored cross-modal knowledge, without any fine-tuning.

## Design
Embed cells from one modality; probe the other modality's matched features.

## Expected outcomes and how to read them
Probe R² above the PCA baseline means cross-modal knowledge in the pretrained weights.

## Risks
Over-powerful probes; leakage via the pretraining corpus.
