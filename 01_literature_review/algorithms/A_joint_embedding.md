# A. Joint embedding of paired spatial RNA + protein

**Use when** RNA and protein were measured on the *same spot or cell*: Tier 1 data in [`../../02_datasets/`](../../02_datasets/).
**Answers:** is a tissue domain defined more by RNA or by protein? Which signal is shared and which is modality-specific?

---

## A1. SpatialGlue (Long et al., *Nat Methods* 2024)

**Idea.** A graph neural network with *dual attention*. For each modality, one encoder runs on a **spatial graph** (k nearest spots) and one on a **feature graph** (k nearest spots in expression space). Within-modality attention mixes the two. A second, between-modality attention mixes RNA and protein into one embedding. The between-modality attention weights say, per spot, how much each modality contributed.

**Input.**
```
adata_omics1  AnnData  RNA: spots x genes (counts), obsm["spatial"]
adata_omics2  AnnData  protein (ADT): spots x proteins (counts), obsm["spatial"]   # same spots, same order
```
Pre-processing (repo `SpatialGlue/preprocess.py`):
- RNA: filter → HVG (~3,000) → normalise → log → scale → PCA (often n_comps = number of proteins − 1) into `obsm["feat"]`
- protein: CLR per cell (`clr_normalize_each_cell`) → scale → PCA into `obsm["feat"]`
- graphs: `construct_neighbor_graph(adata_omics1, adata_omics2, datatype=..., n_neighbors=3 or 6)` for spatial; feature kNN with k = 20 (correlation metric)

**Output.** `obsm["SpatialGlue"]` (joint latent), `alpha` attention weights (within-modality: spatial vs feature graph; between-modality: RNA vs protein), clusters.

**Data in paper.** SPOTS mouse spleen (GEO GSE198353), Stereo-CITE-seq mouse thymus, spatial ATAC–RNA mouse brain, in-house Visium CytAssist human lymph node, plus simulations. Benchmarked against Seurat WNN, totalVI, MultiVI, MOFA+, MEFISTO, scMM and StabMap.

**How we use it.**
- **Modality-attention map** = a learned, unsupervised counterpart of our variance ladder (E1). Regions where protein attention ≫ RNA attention are candidates for *protein-specific* (non-transferable) organisation.
- **Ablation test:** retrain with the RNA input shuffled within cell type, then with the protein input shuffled. If domains survive the RNA shuffle, they are protein-driven.
- **UQ:** attention weights vary by seed. Report the mean ± SD over ≥5 seeds, and compare against attention obtained with the shuffled-RNA control.

**Pitfalls.** It needs spots with both modalities. Attention is not causal importance (see [F](F_attribution.md)). Default hyper-parameters are tuned for spot data, so increase `n_neighbors` for single-cell data.

---

## A2. MultiGATE (Miao et al., *Nat Commun* 2025)

**Idea.** A two-level graph-attention autoencoder. Level 1 attends over the spatial graph. Level 2 attends **across modalities between features**, e.g. a protein and its coding gene, or an ATAC peak and a gene, which also yields an inferred cross-modal regulatory link.

**Input.** Two AnnData objects (RNA + protein, or RNA + ATAC) with coordinates. The SPOTS mouse-spleen tutorial is the RNA + protein template ([docs](https://multigate.readthedocs.io/en/latest)).

**Output.** Joint embedding; feature-level cross-modal attention scores.

**How we use it.** Its protein ↔ gene attention is a *learned cognate map*. Compare it with the hand-curated protein → gene map used in E2. A protein whose attention concentrates on *non-cognate* genes is informed by a programme rather than its own transcript (cf. HLA-DR in our simulation).

---

## A3. PRESENT (Li et al., *Brief Bioinform*; bioRxiv 2024)

**Idea.** Contrastive multi-view autoencoder: per-omics graph-attention encoders with count-likelihood decoders, plus an inter-omics alignment loss. It extends to **multiple samples** with batch alignment and a cross-sample graph.

**Input.** AnnData per modality (RNA, ADT, ATAC), coordinates. **Output.** Joint embedding; batch-integrated embeddings across slides.

**How we use it.** Multi-slide integration before cross-dataset transfer (E7), so that "transfer failure" isn't just a batch effect.

---

## A4. Other integrators worth knowing

| Method | One-liner | Why it matters for us |
|---|---|---|
| **COSMOS** (*Nat Commun* 2025) | cooperative GNN integration of spatial multi-omics | alternative integrator for robustness checks |
| **ARISE** (*Bioinformatics*) | RNA-anchored shared-edge graph: RNA similarity ∩ spatial proximity; auxiliary modalities encoded on that scaffold | its "RNA as anchor" design matches our RNA → protein framing ([code](https://github.com/XiangxiangWang-code/ARISE)) |
| **totalVI** (Gayoso et al., *Nat Methods* 2021) | VAE for CITE-seq; protein = foreground/background NB mixture | non-spatial floor; models protein background, which matters for low-signal ADTs |
| **Seurat WNN / muon** | per-cell modality weights from kNN prediction accuracy | simplest baseline; 10x used it on the ccRCC Xenium Protein demo |
| **MOFA+ / MEFISTO** | linear factor models; MEFISTO adds smooth (spatial) factors | interpretable **shared vs modality-specific** factors → idea N5 |

## Experiments these feed

E1 (does protein-specific organisation exist?), E5 (does a domain learned on one region separate the others?), N5 (shared vs specific factors). See [`../../05_experiment_design/`](../../05_experiment_design/).
