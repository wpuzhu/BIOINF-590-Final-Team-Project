# Future ideas (outside the two-stage plan)

These ideas came from the v2 brainstorm. They are kept for the discussion section or follow-up work. The full original write-ups are in git history (commit `bcfa562`, folder `docs/future_ideas.md/`).

| ID | Question | Data | Methods | Status / link to the stages |
|---|---|---|---|---|
| N2 | After removing what RNA (own + neighbourhood) predicts, is protein still spatially structured within a cell type? | same-cell (P1/P2/P6/P7) | cross-fitted RNA → protein predictor; within-type Moran's I on residuals; toroidal shift tests | code in `shared/info_map` (residual map) and `shared/uq/uq_core.py`; explains *informative-only* signal in Stage 2 |
| N3 | How many and which transcripts are needed to reach X% of the full-RNA R² for each protein? | P6 CosMx WTX, P1 | group LASSO / forward selection; cognate-only and random panels as baselines | design only; a gene-level version of the Stage 2 teacher-input ablation |
| N4 | Can protein predicted from RNA-only slides be trusted, and where? | train P1/P2/P6; deploy SPATCH Xenium 5K; validate on adjacent CODEX | DGAT, ridge probe, totalVI; conformal under shift | design only |
| N5 | Which variance components are shared by RNA and protein, which are specific, which are spatially smooth? | same-cell | MOFA+ / MEFISTO, APOLLO-style partially shared latents | design only; gives `align` vs `recon_specific` a factor-level reading |
| N6 | How would a cell's protein change if its neighbourhood RNA changed? | same-cell | Celcomen, Cellina; edited neighbourhood inputs | design only; needs no-unmeasured-confounding assumptions |
| N7 | Do Spatium embeddings linearly encode cognate RNA, and Nicheformer embeddings protein? | same-cell (P1/P6/P7) | frozen embeddings + ridge probes | cheap add-on to Stage 1 (run on the same embeddings) |
