# U5 · Multiplicity, power and pre-registration

## How many tests?
E1 alone: proteins (≈ 27) × cell types (≈ 8) × steps (3) ≈ 650 tests per dataset. Plan for it:
1. **Primary contrasts (pre-registered, Holm):** a handful, e.g. "cognate gain > 0 for lineage proteins", "neighbourhood gain > 0 for PD-L1-like proteins in tumour cells", "arrangement > 0 for any protein".
2. **Everything else (exploratory, BH-FDR 10%)**, reported as a map with q-values.
3. **Hierarchical testing:** test at protein level first (pooled over types); test types only within significant proteins. This keeps FDR control with fewer tests.

`uq.holm(p)`, `uq.bh(p)`.

## Power with few units
With n patients, the exact sign-flip test cannot give p < 2⁻ⁿ. With n = 8 that is 0.0039, so Holm over 30 tests cannot reach 0.05 (our simulation: smallest Holm-adjusted p = 0.117 despite large effects). Options:
- more units (spatial blocks within patients, used with a hierarchical bootstrap);
- fewer primary tests;
- report **effect sizes with CIs** as the main result and treat p-values as secondary.

## Knockoffs for gene-level claims
Model-X knockoffs (Candès et al., *JRSS-B* 2018) give FDR-controlled selection of genes that carry **conditional** information about a protein. With strongly correlated genes, construct knockoffs on gene **groups** or on PCs.

## Pre-registration (lightweight)
Before running real data, commit in `04_ideas/I1_information_map/` a list of: primary contrasts, the split unit, block size, radius grid, α levels, and the expected direction for each protein class (table in `03_assumptions/why_information_does_not_transfer.md`).
