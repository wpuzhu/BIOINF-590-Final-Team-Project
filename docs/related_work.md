# Related work: informative and transferable information between spatial transcriptomics and proteomics

No paper we found asks our exact question: which RNA information is *informative* for protein-based tasks, and which part *transfers* to a protein-only model. The question sits where five lines of work meet. For each line, this page says what it gives us and what it leaves open.

## 1. RNA–protein concordance in the same cell or section (descriptive)

- **Same-section / same-cell platforms**: CosMx multiomic (980 RNA + 108 proteins, NSCLC and breast; He et al. 2022, *Nat Biotechnol*), Xenium In Situ Gene and Protein (2025), Xenium followed by IMC, PhenoCycler or COMET on the same slide (vendor workflows; Duchini et al. 2026 open resource), DBiT-plus (*Nat Methods* 2026).
- **Typical findings**: lineage markers agree well (CD3E/CD3, MS4A1/CD20). State, post-translational and secreted markers agree poorly. Examples: Bruker's 2026 CellScape + CosMx study found FOXP3 *protein* levels splitting Treg states by neighbourhood in a way RNA alone confounds. A Frontiers 2025 pipeline co-registers ST, SP and H&E to compute per-cell transcript–protein correlation.
- **Leaves open**: correlation is not usefulness. These studies do not ask whether RNA improves a *protein-based* prediction, or whether that improvement survives without RNA.

## 2. Cross-modal prediction (RNA → protein): the "informative" side

- **Single-cell**: totalVI, sciPENN, scTranslator, Seurat WNN / bridge integration; cross-modality translators such as scButterfly, JAMIE and scPair.
- **Spatial**: scProSpatial (2025; spatial protein inferred from scRNA-seq), SpaTranslator (bioRxiv 2025; spatial multi-omics translation within and across slices), STProtein (arXiv 2026; GNN multi-task prediction of spatial protein from ST), the STP challenge (ST + H&E → protein in glioma), SpatialGlue / MultiGATE / PRESENT / COSMOS (joint embeddings on spatial-CITE-seq and SPOTS).
- **Gives us**: how much protein variance RNA (own + neighbourhood) explains. That is our information map (`shared/info_map`) and Δ_info.
- **Leaves open**: these methods need RNA at test time. Their accuracy measures *informativeness*, not *transferability*. They are also rarely evaluated with patient-level uncertainty.

## 3. Shared vs modality-specific factors

- APOLLO (partially overlapping latent spaces; CITE-seq, SHARE-seq, imaging; 2024–2025), IndiSeek (information-guided disentanglement; 2025), UniVI (*Genome Res*), MOFA+ / MEFISTO, MultiVI, scMM.
- **Gives us**: the formal split that motivates our arms. The *shared* part is what alignment (`align_*`) can transfer. The *RNA-specific* part is informative but should not transfer (`recon_specific` negative control). IndiSeek's framing (shared features carry the cross-modal dependence; specific features should be independent of them) is close to ours.
- **Leaves open**: no spatial component, and no test of what a single-modality model gains downstream.

## 4. Foundation models and cross-modal pretraining transfer (Stage 1)

- **RNA**: Nicheformer (*Nat Methods* 2025; SpatialCorpus-110M), scGPT, Geneformer, CellPLM, SToFM (SToCorpus-88M).
- **Protein**: Spatium (> 51M cells, 70 studies; CODEX / IMC / CyCIF), KRONOS (2025), VirTues (*Nature* 2026), CA-MAE.
- **Both**: HEIST (graph foundation model pretrained mainly on MERFISH/Xenium, evaluated on proteomics cohorts) is the closest published claim that SPT pretraining helps on proteomics.
- **Leaves open**: no random-init controls at matched architecture, no token-identity ablations, no equivalence testing. Justin's H1 benchmark and our testbed suggest the cross-modal gain is ≈ 0, which makes this a real open question.

## 5. Privileged information and cross-modal distillation (Stage 2)

- **Theory**: learning using privileged information (Vapnik & Vashist 2009); generalized distillation (Lopez-Paz et al. 2016); modality hallucination and missing-modality distillation (e.g. depth → RGB, multimodal → mono-modal segmentation); multi-teacher privileged KD; the modality focusing hypothesis (Xue et al. 2023): what transfers is the decisive, *shared* part.
- **Closest biology analogues**: these distil **omics → histology**, not RNA → protein.
  - *Cross-Modal Knowledge Distillation from Spatial Transcriptomics to Histology* (Hizmi, Bakulin, Bagon & Yosef, arXiv 2026): a frozen ST niche model (NOLAN on scVI features) teaches an H&E-only student on matched neighbourhoods. This is our Stage 2 with H&E in place of protein.
  - MINT (ST-supervised fine-tuning of pathology foundation models, 2026); PEKA (MICCAI 2025); multi-teacher distillation for gene-expression prediction from H&E; MoMKD (genomics → histology, 2026); subspace distillation for transcriptome-agnostic WSI inference.
  - MKDR (multi-omics modality completion + KD; *Front Oncol* 2025).
- **Gives us**: the training recipe (soft labels, temperature, feature matching, neighbourhood-level matching when pairing is imperfect).
- **Leaves open**: (a) RNA → *protein* students; (b) controls that separate "RNA taught it" from "unlabelled cells / regularisation taught it" (our `self_train`, `label_smooth`, `kd_permuted`); (c) graded nulls and a pairing-distance sweep that locate *where* and *at what spatial scale* the transferable information lives; (d) patient-level CIs and equivalence tests.

## Where this project is new

1. Explicit **informative vs transferable** decomposition (Δ_info, Δ_trans, η, RNA-specific transfer) with patient-level CIs.
2. **Stage 1** as a pre-registered 2 × 2 (Nicheformer/Spatium × RNA/protein target) with random-init, identity, mapping and overlap ablations, and equivalence testing.
3. **Stage 2** comparing distillation, alignment and reconstruction on one protein + graph student. The controls and the input-level ladder (nulls, neighbourhood-only RNA, jitter, graph off, panel, labels) show *why* information does or does not transfer.

Search date: Nov 2026. Preprints are marked as such. Re-check for newer RNA → protein distillation work before writing the final report.
