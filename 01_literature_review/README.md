# 01 · Literature review: algorithms for spatial transcriptomics ↔ spatial proteomics

**Goal of this section.** Know every method we might use well enough to (i) prepare its **input**, (ii) know what its **output** means, and (iii) know which **experiment** in [`../05_experiment_design/`](../05_experiment_design/) it serves.

## How the field is organised (our taxonomy)

The methods sort by what they assume about pairing and what they output:

| Family | Pairing it needs | Output | Role in our project | File |
|---|---|---|---|---|
| **A. Joint embedding** of paired RNA + protein | same spot/cell | one latent per spot/cell, plus modality weights | which modality drives each domain; shared vs specific signal | [`algorithms/A_joint_embedding.md`](algorithms/A_joint_embedding.md) |
| **B. Cross-modal prediction** (RNA → protein) | same spot/cell for training | predicted protein for new RNA-only data | predictor backbone for ablations; "virtual protein" | [`algorithms/B_cross_modal_prediction.md`](algorithms/B_cross_modal_prediction.md) |
| **C. Matching and registration** | none (unpaired) or adjacent sections | cell–cell matches, aligned coordinates | building pairs from SPATCH-like data; registration error for E8 | [`algorithms/C_matching_registration.md`](algorithms/C_matching_registration.md) |
| **D. Foundation models** | none (pretrained on one modality) | frozen embeddings, fine-tuned classifiers | H1 transfer; linear probes (N7) | [`algorithms/D_foundation_models.md`](algorithms/D_foundation_models.md) |
| **E. Neighbourhood-context and counterfactual models** | within one modality | per-view importance, neighbour effects, simulated perturbations | model-level counterparts of our spatial ablations (E3/E4); N6 | [`algorithms/E_neighbourhood_context.md`](algorithms/E_neighbourhood_context.md) |
| **F. Attribution and interpretability** | any trained model | feature/edge importances | explaining H1/H2 models; checking E2 at model level | [`algorithms/F_attribution.md`](algorithms/F_attribution.md) |
| **G. Spatial-variability detectors** (v1) | one modality | per-feature (and per-cell-type) spatial variability p-values | residual maps (E6), N1, N2 | [`algorithms/G_spatial_variability.md`](algorithms/G_spatial_variability.md) |

A one-page version of the same catalogue is [`methods_catalog_short.md`](methods_catalog_short.md).

## Input formats at a glance

Almost everything consumes an **AnnData** object (Python) or a **SpatialExperiment / Seurat** object (R). The differences are in *what goes in `X`* and *what has to be pre-computed*.

| Method | Language | Container | What goes in | Pre-processing it expects | Key hyper-parameters |
|---|---|---|---|---|---|
| SpatialGlue | Py (torch) | two AnnData, same spots, same order | RNA counts; protein ADT counts; `obsm["spatial"]` | RNA: HVG → log-norm → scale → PCA; protein: CLR → scale → PCA (stored in `obsm["feat"]`) | spatial kNN `n_neighbors` (3 for SPOTS, 6 for others), feature kNN `k=20` |
| MultiGATE | Py | two AnnData (RNA + protein/ATAC) | normalised matrices + coordinates | spatial graph; cross-modality feature attention (for RNA + protein, seed with protein ↔ coding-gene links; check the SPOTS tutorial) | spatial radius, attention heads |
| PRESENT | Py | AnnData per modality | raw counts per modality (count-likelihood decoders per omics) | spatial graph | latent dim, contrastive weight |
| totalVI | Py (scvi-tools) | one AnnData | RNA **raw counts** in `X`, protein **raw counts** in `obsm["protein_expression"]` | none (models counts) | latent dim, batch key |
| DGAT | Py | paired AnnData (RNA + protein, coords) | gene + protein matrices, coordinates | spatial graph + feature graph | GAT layers |
| MaxFuse | Py | four NumPy arrays | `shared_arr1/2` (linked features: protein ↔ cognate gene), `active_arr1/2` (all features per modality) | log-norm RNA, scaled protein; protein→gene name table (340 rows shipped in the repo) | `metacell_size`, `matching_ratio`, smoothing method |
| Tangram | Py | AnnData sc + AnnData spatial | shared genes | `tg.pp_adatas` intersects genes | cell vs cluster mode |
| STalign / PASTE2 | Py | coordinates (+ images) / AnnData per slice | cell positions or rasterised density; expression for PASTE2 | – | LDDMM smoothness / overlap fraction `s` |
| Nicheformer | Py | AnnData, genes aligned to model vocabulary | **raw counts** | library-size normalise → divide by technology mean (Xenium/CosMx/MERFISH/dissociated) → rank → token IDs; context tokens for species/assay/modality | max sequence length, context tokens |
| Spatium | Py | AnnData (cell × protein) saved to disk | protein intensities + metadata (Technology, Tissue, Health status, Disease) | `DataSchema.set()`, `FineTuneTokenizer` ranks proteins → token IDs (234-protein vocabulary) | cell-type column; niche tokenizer for spatial tasks |
| HEIST | Py | spatial graph + per-cell-type gene co-expression graphs | normalised HVG expression, coordinates, cell types | Voronoi spatial graph; MAGIC denoising; mutual-information co-expression edges per cell type | MI threshold τ |
| Novae | Py | AnnData | normalised expression + `obsm["spatial"]` | `novae.spatial_neighbors(adata)` | zero-shot vs fine-tune |
| KRONOS / VirTues | Py | multiplex **images** (marker channels) + marker names | per-channel images; VirTues uses marker embeddings from a protein language model | patching / tiling | patch size, marker list |
| MISTy (mistyR) | R | data frames per view | target matrix (e.g. protein) + any predictor view (e.g. RNA of neighbours) + coordinates | `create_initial_view` → `add_juxtaview` / `add_paraview(l, family)` → custom views via `create_view` | paraview radius `l`, kernel family |
| BANKSY | Py/R | AnnData / SpatialExperiment | normalised expression + coordinates | neighbourhood-augmented matrix | λ (neighbour weight), k |
| CellCharter | Py | AnnData | latent from scVI (RNA) or trVAE-MSE (protein) | neighbourhood aggregation over hops | number of hops, cluster K |
| STELLAR | Py (PyG) | annotated reference + unannotated target | cell × marker matrix + coordinates + labels | distance-threshold graph | distance threshold |
| ctSVG/ctSVP detectors (v1) | R/Py | counts or intensities + coordinates + cell types | see `notebooks_v1/` | per method | – |

## Folder contents

```
01_literature_review/
├── README.md                ← this page (taxonomy + input formats)
├── methods_catalog_short.md ← one-page catalogue with links
├── algorithms/              ← one file per family, each method: idea · input · output · how we use it · pitfalls
├── notebooks_v1/            ← v1 template notebooks (ctSVG, protein-native methods, KRONOS/VirTues)
└── archive_v1/              ← v1 literature review and plan (ctSVP biomarker project)
```

## Key references that frame the whole project

- Hu et al., *Nat Methods* 2024: benchmark of 14 prediction and 18 integration algorithms for single-cell multi-omics. totalVI and scArches were best at predicting protein from RNA. This is the non-spatial baseline our spatial predictors must beat.
- Chen et al., *Nat Biotechnol* 2024 (MaxFuse): formalises **weak linkage**. Protein panels share few, weakly correlated features with RNA, which is exactly why RNA → protein transfer is hard.
- Ren et al., *Nat Commun* 2025 (SPATCH): measured the same tumours on four SPT platforms plus adjacent-section CODEX, and reported transcript–protein alignment. Our main source of realistic cross-modal pairs.
- Madhu et al., ICLR 2026 (HEIST): an SPT-pretrained graph foundation model, evaluated on CODEX/MIBI. The closest existing work to H1.
- Tejada-Lapuerta et al., *Nat Methods* 2025 (Nicheformer) and Spatium (bioRxiv 2026): the two foundation models in H1.
