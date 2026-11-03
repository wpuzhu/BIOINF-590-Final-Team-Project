# E2 · Gene-group ablation (drop-and-refit) and conditional permutation importance

**Question.** Is a protein's RNA information in its own gene, in a gene programme, or spread across the transcriptome?

**Code.** `gene_ablation(...)` (drop-and-refit); `uq.group_permutation_importance(...)` (same model, permuted group)

## Manipulation (input level only)
Remove from **all** inputs (own and neighbourhood): own cognate gene(s); all cognate genes; named programmes (MHC-II, IFN response, hypoxia, proliferation, ligand–receptor sets); random gene sets of equal size (≥5 draws).
Permutation variant: at test time, permute the group's genes jointly **within sample × cell type** (conditional permutation keeps cell-type information intact).

## Controls
Random gene sets of matched size and matched mean expression (null distribution); a group of genes known to be irrelevant (housekeeping).

## Read-out
ΔR² vs full model (drop-and-refit) and importance = R²(intact) − R²(permuted).

## Uncertainty quantification
Bootstrap CIs over units; empirical p-value of the named group's ΔR² against the random-set distribution; BH across groups.

## Interpretation
Drop-and-refit measures **unique** information: two redundant groups each look unimportant. Permutation importance measures what **this model** uses. If they disagree, there is redundancy → use joint ablation or Shapley over groups.

## Failure modes
Panel design bias (cognate genes added on purpose); 1:many markers (A2); correlated programmes.

## Simulation self-check
Dropping the MHC-II programme costs HLA-DR 0.51 R², its cognate alone 0.16; random sets cost 0. *IFNG* alone costs PD-L1 ~0 under drop-and-refit (redundant with T-cell markers), but conditional permutation importance is small and significant (0.004, CI 0.003–0.004).
