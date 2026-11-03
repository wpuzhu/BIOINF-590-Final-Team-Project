# Methods catalog (one page): spatial transcriptomics ↔ proteomics

*Short version. Detailed write-ups with input formats are in [`algorithms/`](algorithms/); the taxonomy and input table are in [`README.md`](README.md).*

Grouped by the role each method can play in this project. **Role** says where it plugs in: *measure* (information map, I1), *H1*, *H2*, *baseline*, or *explain*. Check each repo's README for current install instructions before running anything.

## A. Joint embedding of paired spatial RNA + protein (same section)

Use these to build a shared representation of paired data. For ablations, use their per-modality attention or weights to see which modality drives each domain.

| Method | Idea | Data shown on | Role | Link |
|---|---|---|---|---|
| **SpatialGlue** (Long et al., *Nat Methods* 2024) | GNN with dual attention: within-modality (spatial + features), then between-modality | SPOTS spleen, Stereo-CITE-seq thymus, spatial ATAC–RNA | baseline integrator; its modality attention weights are a built-in "which modality matters here" readout | [doi:10.1038/s41592-024-02316-4](https://doi.org/10.1038/s41592-024-02316-4) |
| **MultiGATE** (Miao et al., *Nat Commun* 2025) | two-level graph-attention autoencoder; also infers cross-modality feature links | SPOTS, spatial ATAC–RNA | cross-modal feature attention as a learned cognate map (compare with our E2) | [docs](https://multigate.readthedocs.io/en/latest) |
| **PRESENT** (Li et al., *Brief Bioinform*) | contrastive multi-view autoencoder; multi-sample integration | spatial RNA + ADT + ATAC | multi-slide integration before E7 transfer | [docs](https://bio-present.readthedocs.io) |
| **COSMOS** (*Nat Commun* 16:27, 2025) | cooperative GNN integration of spatial multi-omics | spatial multi-omics | alternative integrator | search "Cooperative integration of spatially resolved multi-omics data with COSMOS" |
| **ARISE** (*Bioinformatics*) | RNA-anchored shared-edge graph (RNA similarity ∩ spatial proximity), hierarchical fusion | bi-/tri-modal benchmarks | "RNA as anchor" design fits our RNA→protein framing | [github](https://github.com/XiangxiangWang-code/ARISE) |
| **totalVI** (Gayoso et al., *Nat Methods* 2021) | VAE for CITE-seq RNA + protein; protein background/foreground model | dissociated CITE-seq | non-spatial baseline; also predicts protein from RNA | [scvi-tools](https://scvi-tools.org) |
| **WNN** (Seurat / muon) | per-cell modality weights | any paired data | the simplest baseline; 10x used it on the ccRCC Xenium Protein demo | Seurat v4+ |

## B. Predicting spatial protein from spatial RNA (cross-modal prediction)

These are the closest relatives of our ablation predictor. Swap one in for ridge to check that conclusions hold with a nonlinear model.

| Method | Idea | Role | Link |
|---|---|---|---|
| **DGAT** (Wang, Cody, Osmanbeyoglu, bioRxiv 2025) | dual-graph attention; RNA encoder + protein decoder trained on paired spatial data, then applied to RNA-only slides; evaluated leave-one-sample-out | strong spatial predictor for E1–E5 and the "virtual protein" idea | [bioRxiv](https://doi.org/10.1101/2025.07.05.662121) · [HF](https://huggingface.co/Carlwanghy1128/DGAT) |
| **STProtein** (Jiang et al., arXiv 2602.05811, AAAI-SPARTA 2026) | GNN + multi-task learning to predict spatial protein from spatial multi-omics | second spatial predictor; compares against totalVI and scArches | [arXiv](https://arxiv.org/abs/2602.05811) |
| **totalVI / scArches** | best protein predictors in the Hu et al. 2024 *Nat Methods* benchmark of multi-omics prediction | non-spatial floor: does adding space help at all? | [benchmark](https://doi.org/10.1038/s41592-024-02429-w) |

## C. Linking unpaired or adjacent-section data

Needed for H2 stage 2 (SPATCH adjacent sections) and E8.

| Method | Idea | Role | Link |
|---|---|---|---|
| **MaxFuse** (Chen et al., *Nat Biotechnol* 2024) | iterative co-embedding, smoothing and matching that works with **weak linkage** (few, weakly correlated shared features), e.g. CODEX ↔ scRNA-seq | match CODEX cells to RNA cells; its antibody-to-gene table gives the protein → gene mapping (H1 harmonisation) | [doi:10.1038/s41587-023-01935-0](https://doi.org/10.1038/s41587-023-01935-0) |
| **STalign**, **PASTE2** | image/landmark and optimal-transport alignment of serial sections | registration for SPATCH; gives the µm error that sets E8 and H2 arm H | GitHub: JEFworks-Lab/STalign, raphael-group/paste2 |
| **Tangram** | maps single cells onto space by optimal transport | alternative RNA→space mapping when only scRNA-seq is paired | GitHub: broadinstitute/Tangram |

## D. Foundation models (H1)

| Model | Modality | Notes from Justin's H1 plan and sources | Link |
|---|---|---|---|
| **Nicheformer** (Tejada-Lapuerta et al., *Nat Methods* 2025) | SPT + dissociated | pretrained on SpatialCorpus-110M (~57M dissociated + 53M spatial cells); a cell is a sequence of rank-ordered gene tokens over a 20,310-gene vocabulary | [doi:10.1038/s41592-025-02814-z](https://doi.org/10.1038/s41592-025-02814-z) |
| **Spatium** | protein | protein foundation model for spatial proteomics (H1 plan: 51.8M cells, 70 studies, ~239-marker vocabulary, rank tokens) | [docs](https://spatium.readthedocs.io/en/latest/) · [github](https://github.com/ploughhh/Spatium) |
| **HEIST** (Madhu et al., arXiv 2506.11152; ICLR 2026 per H1 plan) | SPT → protein | hierarchical graph foundation model; pretrained on SPT, evaluated on CODEX and MIBI. Justin's 29 Oct link | [arXiv](https://arxiv.org/abs/2506.11152) |
| **KRONOS** (Mahmood lab) | protein images | self-supervised on 47M patches, 175 markers, 8 platforms; segmentation-free | [github](https://github.com/mahmoodlab/KRONOS) |
| **VirTues** (Bunne lab) | protein images | marker-aware multi-scale transformer; protein-language-model marker tokens; arbitrary panels | [arXiv](https://arxiv.org/abs/2501.06039) |
| **Novae** | SPT | graph-based foundation model for spatial domains; H1 fallback if Nicheformer will not run | see H1 plan sources |

## E. Neighbourhood-context and perturbation models (explain / ablate)

These make the "spatial context" part of our question explicit and offer model-level counterparts to our input-level shuffles.

| Method | Idea | Role | Link |
|---|---|---|---|
| **MISTy / mistyR** | multi-view model: intrinsic vs juxta vs para views; reports each view's contribution | the closest analogue of our E1/E3 at the model level; runs on protein *or* RNA, so compare the view importances across modalities | Bioconductor `mistyR` |
| **NCEM** (Fischer et al. 2023) | node-centric expression models; neighbour cell-type effects on expression | interpretable neighbourhood effects | GitHub: theislab/ncem |
| **Celcomen** (Megas et al., ICLR 2025) | causally identifiable generative GNN; disentangles intra- vs inter-cellular gene programmes; spatial counterfactuals ("knock out a gene here") | in-silico perturbation of the RNA neighbourhood, then read predicted protein | [arXiv](https://arxiv.org/abs/2409.05804) |
| **Cellina** (arXiv 2606.08493) | tissue-graph counterfactuals by edge perturbation (rewire neighbours) and node perturbation (change neighbour expression) | formalises our E4 shuffles as interventions; compares to MintFlow, Concert, SpatialProp | [arXiv](https://arxiv.org/html/2606.08493v2) |
| **BANKSY**, **CellCharter**, **STELLAR** | neighbour-augmented features / clustering / GNN annotation; each runs on RNA and protein data | shared niche definitions for E5 region hold-out on both modalities | GitHub: prabhakarlab/Banksy, CSOgroup/cellcharter, snap-stanford/stellar |

## F. Attribution tools (to explain trained H1/H2 models)

- **Integrated gradients / SHAP** (Captum, `shap`) on input genes or markers: model-level counterparts of E2.
- **GNNExplainer / PGExplainer** (PyTorch Geometric) for which neighbours and edges a GNN uses: model-level counterpart of E3/E4.
- **Linear probes** on frozen embeddings (as in H1) for cell type, niche, and *cognate-gene expression*. If a protein-pretrained embedding linearly predicts gene expression, the knowledge is shared.

## G. Baselines carried over from v1 (still useful)

- ctSVG / ctSVP detectors (CELINA, STANCE, SPARK-X, within-type Moran's I) in [`notebooks_v1/`](notebooks_v1/) and [`../07_code/ctsvp_baseline/`](../07_code/ctsvp_baseline/): used for **E6 residual maps** and for comparing within-cell-type spatial variability across modalities (idea N1 in [`../04_ideas/`](../04_ideas/)).
- Protein-native methods from v1 (mistyR, Statial, imcRtools, scimap, CytoCommunity): baselines for protein-side tasks.
