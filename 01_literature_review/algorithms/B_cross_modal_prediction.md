# B. Cross-modal prediction: spatial RNA → spatial protein

**Use when** you have paired data to train on and want protein estimates (or prediction accuracy as a *measurement of information*) for cells where only RNA is used.
**This is the core family for our project.** "How well can protein be predicted from RNA input X, on held-out units?" is our operational definition of *informative* (E1–E4), and "…on held-out regions/datasets?" is our definition of *transferable* (E5, E7).

---

## B0. Our fixed linear probe (in this repo)

`07_code/ablation/xmodal_ablation.py`: ridge regression per protein on [cell-type one-hot, cognate gene(s), RNA PCs, neighbourhood-averaged RNA PCs], with cell-type interactions.

**Input.** `PairedSpatial(rna, protein, coords, cell_type, sample, genes, proteins, cognate, region)`, or `from_anndata(adata, protein_key="protein")`.
**Why linear?** It is fast enough to run hundreds of ablations; its gains cannot come from capacity; and if a deep model beats it, the extra must come from non-linear structure. Every conclusion should be re-checked with a non-linear predictor (below).

---

## B1. DGAT (Wang, Cody & Osmanbeyoglu, bioRxiv 2025)

**Idea.** Dual-graph attention network. An mRNA encoder uses gene expression plus spatial context; a protein encoder models spatial co-expression of proteins. A shared latent feeds an mRNA reconstruction decoder and a protein prediction decoder. At inference only the mRNA encoder and protein decoder are used, so RNA-only slides get predicted protein.

**Input.** Paired gene and protein matrices with spatial coordinates (training); RNA-only spatial data (prediction). Pretrained weights, demo notebooks and data are on Hugging Face (`Carlwanghy1128/DGAT`): Demo1_Train, Demo2_Predict, Demo3_Predict_ST.

**Evaluation in paper.** Leave-one-sample-out, which matches our split rule.

**How we use it.** The non-linear replacement for B0 in E1–E5 and the engine for N4 (virtual protein). Its attention over neighbours can be compared with our radius sweep (E3).

---

## B2. STProtein (Jiang et al., arXiv 2602.05811; AAAI 2026 SPARTA workshop, oral)

**Idea.** GNN with multi-task learning to predict spatial protein from spatial multi-omics (mainly transcriptomics); attention over neighbouring spots. It is benchmarked against totalVI and scArches, and evaluated by downstream clustering of predicted protein (e.g. mouse spleen).

**How we use it.** A second non-linear predictor, so conclusions don't hinge on one architecture.

---

## B3. Non-spatial predictors (floors)

| Method | Input | Note |
|---|---|---|
| **totalVI** | CITE-seq-style AnnData: RNA counts in `X`, protein counts in `obsm["protein_expression"]` | joint VAE; predicts protein for RNA-only cells. One of the two best protein predictors in Hu et al. 2024 |
| **scArches** (totalVI reference mapping) | reference CITE-seq model + query RNA | the other top protein predictor in Hu et al. 2024; transfer via architectural surgery |
| **sciPENN, Seurat bridge/label transfer** | CITE-seq reference | additional baselines if time permits |

**Why include them.** If a spatial model barely beats totalVI, spatial context adds little for that protein. That is a result about *what information matters*, not a failure.

---

## B4. Spatial "imputation" methods from the RNA side

Tangram, gimVI and SpaGE (best at predicting unmeasured spatial genes in Li et al., *Nat Methods* 2022) map scRNA-seq onto spatial RNA. They are not RNA → protein, but are useful for N4: impute a missing *cognate gene* on a small-panel SPT slide before predicting protein.

---

## Design notes for using predictors as information meters

1. **The metric is held-out R² / Pearson per protein × cell type,** never pooled across cell types. Pooled R² mostly measures cell-type separation (the v1 confounding lesson).
2. **Split by patient/slide or by spatial block** wider than the largest neighbourhood radius ([`../../06_uncertainty_quantification/`](../../06_uncertainty_quantification/)).
3. **Count predictors are a mismatch for protein intensities.** Use Gaussian/MSE losses on arcsinh-transformed protein, or totalVI's background-aware model for ADT counts. (Justin's slides note that Spatial-DC found count-likelihood methods did poorly on proteomics.)
4. **Predicted protein inherits RNA's biases.** A virtual protein cannot show post-transcriptional regulation; report it only where the information map says RNA explains the protein *and* the mapping transfers.
