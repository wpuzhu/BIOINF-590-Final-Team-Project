# N3 · Minimal gene panel that predicts the protein panel

| | |
|---|---|
| **Status** | design only |
| **Question** | How many and which transcripts are needed to reach X% of the full-RNA R² for each protein? |
| **Assumptions** | A2, A5 |
| **Data** | P6 CosMx WTX (unbiased gene availability), P1 |
| **Methods** | group LASSO / forward selection on the ridge probe; compare against cognate-only and random panels |
| **Experiments** | selection inside training folds; evaluate on held-out units |
| **UQ** | **stability selection** (bootstrap units, report each gene's selection frequency); CI of R² at each panel size |
| **Cost** | minutes |

## Motivation
It answers a design question (what to measure) and cross-checks E2: non-cognate selections mean programme-level information.

## Design
Nested CV: selection on inner folds only.

## Expected outcomes and how to read them
A per-protein curve of R² vs panel size, plus a stable gene list.

## Risks
Selection on full data leaks (A8).
