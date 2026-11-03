# I1 · The information map: which RNA information informs protein, and does it transfer?

| | |
|---|---|
| **Status** | core idea of v2 · code runnable (`07_code/ablation/`) · validated on simulation |
| **Question** | For each protein × cell type: how much held-out variance comes from cell type, the cognate gene, other genes and neighbourhood RNA? Does that hold on held-out regions and datasets? |
| **Assumptions** | A1, A2, A3, A5, A6, A7 ([03_assumptions](../../03_assumptions/)) |
| **Data** | P1/P2 Xenium Protein (main), P6 CosMx multiomic, P3 DBiTplus / D17a spatial-CITE-seq (transfer targets) |
| **Methods** | fixed ridge probe (B0); re-check with DGAT/MISTy ([01_literature_review/algorithms](../../01_literature_review/algorithms/)) |
| **Experiments** | E1–E8 ([05_experiment_design](../../05_experiment_design/)) |
| **UQ** | cluster bootstrap CIs, sign-flip tests, conditional group permutation, Mondrian/weighted conformal, shift decomposition ([06_uncertainty_quantification](../../06_uncertainty_quantification/)) |
| **Cost** | laptop; minutes per dataset (subsample to ≤50k cells for E2/E4) |

## Motivation
H1 and H2 ask *whether* RNA helps a protein model. They do not say *which part* of RNA helps or *where*. The information map is the measurement layer underneath both: keep one simple predictor fixed and change only its inputs.

## Design
1. Variance ladder (E1). 2. Gene-group ablation (E2). 3. Radius sweep (E3). 4. Context controls (E4). 5. Region hold-out (E5). 6. Residual maps (E6). 7. Cross-dataset transfer (E7). 8. Pairing degradation (E8). All with the UQ protocol: unit-level splits, bootstrap CIs, permutation p-values, Holm/BH correction.

## Expected outcomes and how to read them
A figure with proteins × {cell type, cognate, programme, neighbourhood, residual-SV}, coloured by ΔR² with CIs, plus a second panel for transfer (region and dataset). It predicts where H1/H2 should help: proteins with high, transferable RNA-explainable variance.

## Risks
Linear probe underfits, so confirm with a non-linear model. Panel design inflates cognate share (Xenium boosted cognate genes), so replicate on CosMx WTX. Spillover fakes neighbourhood effects (A7).
