# N6 · Counterfactual neighbourhoods

| | |
|---|---|
| **Status** | design only |
| **Question** | How would a cell's protein change if its neighbourhood RNA changed? |
| **Assumptions** | A4, A6, no unmeasured confounding of neighbourhood and state |
| **Data** | same-cell data |
| **Methods** | Celcomen, Cellina; our probe with edited neighbourhood inputs ([E_neighbourhood_context](../../01_literature_review/algorithms/E_neighbourhood_context.md)) |
| **Experiments** | node/edge perturbations; validation against natural contrasts (matched cells) |
| **UQ** | matching estimator with patient-bootstrap CI; agreement between model families |
| **Cost** | moderate |

## Motivation
It moves from 'context predicts protein' to 'what happens if context changes'.

## Design
Remove T cells / IFNG from a tumour cell's neighbourhood in silico; compare predicted protein change with matched real cells that have fewer T cells nearby.

## Expected outcomes and how to read them
Direction and size of neighbourhood effects with CIs.

## Risks
Counterfactuals are unvalidated extrapolations; confounding by region.
