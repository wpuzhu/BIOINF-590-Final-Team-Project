# E. Neighbourhood-context and counterfactual models

**Use when** the question is about *spatial context*: does a cell's protein depend on who its neighbours are, or on what state they are in, and over what distance? Each method below is a model-level counterpart to one of our input-level ablations (E3 radius, E4 shuffles).

---

## E1. MISTy / mistyR (Tanevski et al., *Genome Biol* 2022)

**Idea.** Multi-view learning. Each target marker is predicted from separate **views**: the *intraview* (other markers in the same cell), the *juxtaview* (direct neighbours within a distance threshold), the *paraview* (all cells weighted by a distance kernel with range `l`), and any **custom view** you supply. A meta-model combines the views and reports, per target:
- `multi.R2`: variance explained by all views;
- `gain.R2`: gain from adding the spatial views over the intraview alone;
- per-view contributions and per-predictor importances.

**Input (R).**
```r
library(mistyR)
views <- create_initial_view(protein_df) %>%                       # targets: protein markers per cell
         add_paraview(coords, l = 30, family = "gaussian")          # spatial context of protein
rna_para <- create_initial_view(rna_pc_df) %>% add_paraview(coords, l = 30)
views <- views %>% add_views(create_view("para.rna.30", rna_para[["paraview.30"]]$data, "para.rna.30"))
run_misty(views, "misty_results"); res <- collect_results("misty_results")
```
The vignette `FunctionalPipelinePathwayActivityLigands.Rmd` builds exactly this kind of **cross-feature-set paraview** (ligand expression → pathway activity).

**How we use it.** An independent replication of E1/E3 with a different model class (random forests, view stacking). A significant "para.rna" contribution to a protein means neighbouring RNA informs that protein. Running the **same** MISTy on RNA targets and on protein targets compares how much spatial context matters in each modality.

---

## E2. NCEM (Fischer et al., *Nat Biotechnol* 2023)

Node-centric expression models: a cell's expression is regressed on its own cell type plus the cell types in its neighbourhood (linear, or GNN variants). Coefficients read as "type B within r µm changes gene/protein g in type A by β". **Use:** the interaction coefficients are a parametric version of E4's "composition" effect. Fit on protein (target) with RNA-defined neighbour types to keep the analysis cross-modal.

## E3. BANKSY (Singhal et al., *Nat Genet* 2024)

Augments each cell's features with the mean neighbourhood expression and an azimuthal Gabor filter (local gradient), weighted by λ, then clusters. It works on RNA and protein. **Use:** define niches the same way in both modalities for E5 region hold-out. The λ sweep (0 = cell types → 0.8 = domains) is a cheap "how spatial is this structure" dial.

## E4. CellCharter (Varrone et al., *Nat Genet* 2024)

Aggregates latent features over 1…L neighbourhood hops, then fits a Gaussian mixture; it selects stable cluster numbers automatically. The latent comes from scVI (RNA) or a trVAE with **MSE loss** (protein). **Use:** H1's niche baseline (T2), and region labels for E5.

## E5. STELLAR (Brbić et al., *Nat Methods* 2022)

GNN that transfers cell-type annotations from an annotated spatial dataset to an unannotated one, and can discover novel types. Shown on CODEX (healthy intestine, Hickey et al.) and MERFISH. **Input:** cell × marker matrix, coordinates, distance-threshold graph, labels for the reference. **Use:** a spatial label-transfer baseline in H2 and for the cell-typing task in H1.

---

## Counterfactual models (idea N6)

| Model | Intervention it supports | Notes |
|---|---|---|
| **Celcomen** (Megas et al., ICLR 2025) | knock out / alter a gene **at a location** and generate the counterfactual tissue | causally identifiable disentanglement of intra- vs inter-cellular gene–gene interactions; a global interaction matrix shared across the tissue |
| **Cellina** (arXiv 2606.08493) | **edge perturbation** (rewire who neighbours whom) and **node perturbation** (change neighbours' expression) | formalises "tissue graph counterfactuals"; compares with MintFlow, Concert, Celcomen, SpatialProp; >2.5M cells in benchmarks |
| **MintFlow, Concert, SpatialProp** | label-conditioned or graph-based in-silico perturbations | additional baselines listed by Cellina |

**Link to our ablations.** E4's "shuffle neighbours within cell type" is a crude **node perturbation**. Changing the radius graph is an **edge perturbation**. Counterfactual models do this inside a generative model rather than on the input of a fixed predictor. Agreement between the two is evidence that a neighbourhood effect is real.

**Caution.** Counterfactual predictions are not validated by held-out accuracy. Validate the *direction* against naturally occurring contrasts: cells of the same type and similar own-RNA whose neighbourhoods differ. That is a matching estimator; use bootstrap CIs over patients.
