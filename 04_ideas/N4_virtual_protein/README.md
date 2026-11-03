# N4 · Virtual protein for RNA-only cohorts

| | |
|---|---|
| **Status** | design only |
| **Question** | Can predicted protein from RNA-only slides be trusted, and where? |
| **Assumptions** | A4 (stationarity) critically |
| **Data** | train on P1/P2/P6; deploy on SPATCH Xenium 5K; validate on SPATCH adjacent CODEX |
| **Methods** | DGAT, ridge probe, totalVI floor |
| **Experiments** | neighbourhood-level comparison of predicted vs real adjacent-section protein |
| **UQ** | **conformal prediction intervals** per cell type; flag cells outside the training support (domain classifier); weighted conformal under covariate shift |
| **Cost** | moderate |

## Motivation
It turns the information map into a usable tool, with honest uncertainty.

## Design
Only report proteins the I1 map marks as transferable; attach intervals to every prediction.

## Expected outcomes and how to read them
Coverage on the adjacent CODEX near nominal for lineage proteins; under-coverage flags non-transferable proteins.

## Risks
Platform shift Xenium → CosMx; only 16 CODEX markers to validate.
