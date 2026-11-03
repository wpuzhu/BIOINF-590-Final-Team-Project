# E1 · Variance ladder

**Question.** How much held-out protein variance does each information source add, in a fixed order: cell type → cognate gene → all RNA → neighbourhood RNA?

**Code.** `variance_ladder(d, F, radius=30)` in `07_code/ablation/xmodal_ablation.py`; per-unit scores with `uq.fold_scores`

## Manipulation (input level only)
Nested input sets:
1. one-hot cell type
2. + cognate gene(s) of the protein (with cell-type interactions)
3. + PCs of all RNA (own cell)
4. + neighbourhood-averaged RNA PCs at radius r (with cell-type interactions)

Variant orders (run both): put *all RNA* before *cognate* to measure the cognate gene's **unique** contribution beyond the transcriptome.

## Controls
RNA shuffled across cells within sample (all steps above cell type should collapse to step 1); RNA-derived vs protein-derived cell-type labels (A3).

## Read-out
ΔR² per step, per protein × cell type, mean over held-out units. Also Pearson r (scale-free).

## Uncertainty quantification
Paired per-unit differences → cluster bootstrap 95% CI; one-sided sign-flip test (H0: no gain); BH across proteins × steps (Holm for pre-registered primary contrasts).

## Interpretation
Large cognate step = cell-intrinsic, transcript-determined protein (transferable if the cognate is measured in the target). Large neighbourhood step = context-dependent protein (transfer depends on context being similar). Small total = RNA-independent protein (see B1–B4 in `03_assumptions/why_information_does_not_transfer.md`).

## Failure modes
Order dependence: gains are path-dependent when inputs are correlated (use Shapley over steps if needed). Cell-type labels derived from protein make step 1 contain protein information.

## Simulation self-check
CD3 0.55 → 0.84 at the cognate step (+0.31, CI 0.29–0.32, sign-flip p = 0.004 with 8 units); PD-L1 gains only at the neighbourhood step (+0.11, CI 0.08–0.13). See `08_results/simulation/`.
