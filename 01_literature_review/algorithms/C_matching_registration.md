# C. Matching unpaired cells and registering adjacent sections

**Use when** RNA and protein come from **different cells**: adjacent sections (SPATCH), or dissociated scRNA-seq versus CODEX.
**Why it matters:** every cross-modal analysis on Tier-2 data is only as good as the pairing. Pairing error is itself a reason information "does not transfer" (assumption A1 in [`../../03_assumptions/`](../../03_assumptions/)), so we quantify it explicitly (E8).

---

## C1. MaxFuse (Chen et al., *Nat Biotechnol* 2024)

**Problem it solves.** *Weak linkage*: the linked features (protein ↔ its gene) are few and weakly correlated. A prototypical example is spatial proteomics versus single-cell sequencing.

**Idea (three stages).**
1. **Fuzzy smoothing.** Each modality's features are smoothed over its own kNN graph built on *all* features, which raises the signal-to-noise of the few linked features. Initial matching uses the smoothed linked features (linear assignment, in batches).
2. **Iterative refinement.** Matched pairs ("pivots") give a CCA co-embedding of *all* features. The cells are re-matched in that space; repeat.
3. **Propagation and filtering.** Drop low-quality matches, then propagate matches to unmatched cells.

**Input (exact API, `maxfuse.model.Fusor`).**
```python
fusor = Fusor(shared_arr1, shared_arr2,   # linked features only (protein  <-> cognate gene), same column order
              active_arr1, active_arr2,   # all informative features of each modality
              method="centroid_shrinkage" | "graph_smoothing", labels1=None, labels2=None)
fusor.split_into_batches(...); fusor.construct_graphs(...); fusor.find_initial_pivots(...)
fusor.refine_pivots(...); fusor.filter_bad_matches(...); fusor.propagate(...)
matching = fusor.get_matching(); emb1, emb2 = fusor.get_embedding(...)
```
The repo ships `docs/protein_gene_conversion.csv` (340 antibody → gene rows, e.g. CD273 → *PDCD1LG2*). Use it as the starting **cognate map** for E2 and H1, then check entries by hand.

**How we use it.**
- Build cell-level RNA ↔ protein pairs for SPATCH (Xenium 5K ↔ CODEX) as an *alternative* to coordinate registration. The two pairing methods disagree where pairing is uncertain, which gives a natural UQ signal.
- **Matching confidence as a filter:** run E1–E4 on high-confidence matches only, then on all matches. If the information map changes, conclusions depend on pairing quality.
- MaxFuse's own benchmark shows that label-transfer accuracy under weak linkage depends on the linked features. Our E2 asks the mirror question: how much of protein is in those linked features?

---

## C2. STalign (Clifton et al., *Nat Commun* 2023)

Diffeomorphic (LDDMM) alignment of spatial data: section ↔ section, section ↔ histology, across technologies, and to 3D atlases. **Input:** cell coordinates rasterised into density images, or images (e.g. DAPI). **Output:** a smooth deformation that maps one section's coordinates onto another's. **Use:** register Xenium/CosMx sections onto the adjacent CODEX section in SPATCH. Landmarks (manually picked) improve robustness.

## C3. PASTE2 (Liu, Zeira & Raphael, *Genome Res* 2023; now maintained in PASTE3)

Optimal-transport **partial** alignment of adjacent slices that only partly overlap. **Input:** AnnData per slice (expression + coordinates), overlap fraction `s` (or select it with `select_overlap_fraction`). **Output:** a partial alignment matrix and projected coordinates. **Use:** a second registration method. Disagreement between STalign and PASTE2 gives a registration uncertainty in µm.

## C4. Tangram (Biancalani et al., *Nat Methods* 2021)

Maps dissociated cells to spatial voxels through shared genes (`tg.pp_adatas` → `tg.map_cells_to_space` → `tg.project_genes`). It cannot map to protein-only data directly. **Use:** add whole-transcriptome RNA to SPT spots/cells before pairing with protein.

---

## Registration error budget (what E8 needs)

| Source | Typical size | How to estimate |
|---|---|---|
| Section-to-section distance | 5–10 µm (one section) | known from protocol |
| Deformation after alignment | tens of µm | landmark residuals; STalign vs PASTE2 disagreement |
| Segmentation mismatch between platforms | ~cell diameter | compare nucleus centroids on DAPI |

E8 sweeps artificial jitter (0–100 µm) on **same-cell** data (Xenium Protein) to learn how fast information decays with pairing error. Comparing that curve with SPATCH's measured error tells us whether H2 stage 2 can work at all.
