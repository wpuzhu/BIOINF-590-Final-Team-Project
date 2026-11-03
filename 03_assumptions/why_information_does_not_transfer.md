# Why some information (and some regions) do not transfer

When RNA predicts a protein well in one place but not another, there are two families of explanation. **Biological** reasons mean the RNA → protein relationship really differs. **Technical/statistical** reasons mean it only looks different. A good experiment has to separate them, so each mechanism below comes with a **signature** in our experiments and a **diagnostic** to confirm it.

## Biological mechanisms

| # | Mechanism | Which proteins | Signature in our experiments | Diagnostic |
|---|---|---|---|---|
| B1 | **Translational control** (mTOR, hypoxia, stress granules) | proliferation and growth proteins (pS6, Ki-67 partly, MYC), stress proteins | E5: transfers between similar regions but fails into hypoxic/necrotic regions (tumour core) at fixed mRNA | shift decomposition shows a **concept-shift** component; residual spatial structure (E6) correlates with a hypoxia gene score of neighbours (toroidal test) |
| B2 | **Protein stability and degradation** (half-lives of hours to days vs mRNA minutes to hours) | PTEN, β-catenin, p53, cyclins | weak cognate gain in E1 despite good cognate detection; residual not explained by any RNA | cognate R² ≪ cross-platform protein reproducibility; residual uncorrelated with neighbourhood |
| B3 | **Time lag** between transcription and protein accumulation | induced proteins (PD-L1 after IFN-γ, HLA-DR after IFN-γ) | protein correlates better with *neighbour* inducer RNA (IFNG) than with its own mRNA; E3 peak at paracrine radius | E1: neighbourhood gain > cognate gain; E4: arrangement effect if induction is local |
| B4 | **Post-translational modification / cleavage** | phospho-epitopes, cleaved caspase-3, cleaved collagen | no RNA input explains them (as for pS6 in our simulation) | expected ~0 R² across E1; strong residual spatial structure (E6) |
| B5 | **Secretion and extracellular deposition** | collagens, cytokines, ECM, granzyme (secreted) | protein is located away from the producing cell; neighbourhood RNA at larger radius informs it | E3 peak at larger radius; protein of non-producer types predicted by neighbours' RNA |
| B6 | **Membrane trafficking and shedding** | receptors (e.g. surface vs internal pools) | cognate mRNA tracks total, antibody detects surface only | cognate gain varies by cell type → **cell-type interaction** term matters (E1 with vs without interactions) |
| B7 | **Cell-type–specific coupling** | proteins expressed in several types with different regulation | pooled R² looks fine, but per-cell-type R² differs; transfer fails where composition differs | Mondrian conformal coverage per type; per-type E1 |
| B8 | **Region-specific programmes** (tertiary lymphoid structures, invasive margin, necrosis) | context-induced states | model learns region identity (large-radius features), and E5 fails for regions absent from training | E3 + E5 together; support check (is the region's context seen in training?) |

## Technical and statistical mechanisms

| # | Mechanism | Signature | Diagnostic / fix |
|---|---|---|---|
| T1 | **Covariate shift** (different cell/context mix) | error rises in target, but per-cell-type error is unchanged | importance-weighted source error ≈ target error; weighted conformal restores coverage (`uq.py`) |
| T2 | **Support shift** (unseen cell types or contexts) | very large errors for a subset of target cells | domain classifier AUC ≈ 1 for those cells; flag as extrapolation, don't report as "non-transferable biology" |
| T3 | **Pairing error** (adjacent sections, matching) | state information lost, composition kept | E8 decay curves; compare MaxFuse vs registration pairing |
| T4 | **Platform and antibody differences** | systematic offsets or scale differences per protein | per-protein calibration on a few anchor cells; rank-based metrics (Spearman) vs R² |
| T5 | **Segmentation and spillover** | neighbourhood "helps" lineage proteins of neighbour types | A7 check; nucleus-only intensities |
| T6 | **Panel mismatch** | cognate gene missing in target | restrict to shared genes; report coverage of the cognate map |
| T7 | **Small held-out variance** | R² very negative even when predictions are reasonable (as for PD-L1 in the core in our simulation: R² −1.5, Pearson 0.18) | report Pearson / MSE alongside R²; conformal coverage is scale-free |
| T8 | **Model misspecification** | a linear model fails out of distribution where a non-linear one would extrapolate | oracle refit in target (shift decomposition): a big "target-specific gain" with good in-target fit means the *model class* needs region-specific terms |
| T9 | **Batch effects between slides** | everything is worse on new slides, uniformly | per-sample normalisation; multi-sample integration (PRESENT) before E7 |

## Decision tree for a failed transfer

```
transfer error on target  >>  in-source held-out error ?
├─ no  → transferable (report CI on the gap)
└─ yes → is the target supported? (domain classifier, per-cell-type counts)
         ├─ no  → T2 support shift: extrapolation, not biology
         └─ yes → does importance weighting reproduce the target error?  (shift_decomposition)
                  ├─ yes → T1 covariate shift: transferable *conditional on cell type/context*;
                  │        fix with per-type models / weighting
                  └─ no  → concept shift: P(protein | RNA) differs
                           ├─ consistent with pairing/platform? (E8 curve, T4 offsets) → technical
                           └─ residual spatially structured within type, linked to a biological
                              covariate (toroidal test)? → biological (B1–B8): a finding
```

## What we expect, protein by protein (hypotheses to pre-register)

| Protein class | Expected main information source | Expected transferability |
|---|---|---|
| Lineage surface markers (CD3, CD8, CD20, CD68, PanCK) | cell type + cognate | high across regions and datasets |
| Activation/antigen presentation (HLA-DR, PD-L1) | programme + neighbourhood (IFN-γ) | medium; region-dependent (covariate shift) |
| Proliferation (Ki-67) | cognate + cell cycle programme | medium-high |
| Signalling/phospho (pS6, pERK), stability-regulated (PTEN, β-catenin) | little RNA information | low; residual is spatial (B1/B2/B4) |
| ECM/secreted (collagen, granzyme B) | neighbourhood RNA at larger radius | low-medium; depends on A6 |
