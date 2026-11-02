# Additional ideas (for team discussion)

Ranked by how cheaply they fit next to H1/H2. N1–N3 reuse the ablation code and the v1 ctSVP pipeline, so they are the cheapest.

## N1. Spatial-variability transfer: do RNA ctSVGs become protein ctSVPs?

Our v1 work detects **cell-type-specific spatially variable** features. On paired data, run the same detectors (CELINA, SPARK-X, within-type Moran's I) on each protein *and* on its cognate gene, within the same cell type. Then sort every (protein, cell type) pair into one of four boxes:

| | protein ctSV | protein not ctSV |
|---|---|---|
| **gene ctSV** | concordant spatial programme (transfers) | RNA pattern buffered at protein level |
| **gene not ctSV** | protein-only spatial regulation (does not transfer) | neither |

This is a direct, interpretable answer to "which *spatial* information transfers". It also gives H1/H2 a list of proteins where RNA should and should not help.

## N2. Residual ctSVPs: protein spatial structure RNA cannot explain

This is E6 extended. Predict protein from RNA (own + neighbourhood), then run the full ctSVP battery on the **residual**. Residual ctSVPs are candidates for post-transcriptional or microenvironment-driven control: hypoxia-driven translation, phospho-signalling such as pS6, or protein secreted by neighbours. Phospho-epitopes and secreted proteins are the expected positives. Lineage markers are the expected negatives (their residual should be flat).

## N3. Minimal gene panel for protein: what would you need to measure?

Greedy or LASSO-group selection of genes that best predicts the protein panel (held-out R²), compared with (a) the cognate genes only and (b) random panels of equal size. It answers a practical design question: how many and which transcripts make a "virtual protein panel"? It also cross-checks E2: if selected genes are mostly non-cognate, the information is programme-level.

## N4. Virtual protein for RNA-only cohorts

Train the best predictor (from B in the methods catalog, e.g. DGAT or our ridge+context) on paired data. Apply it to RNA-only slides (Xenium 5K / CosMx 6K in SPATCH) and score it against SPATCH's **adjacent-section CODEX**. That is a real, held-out protein readout without same-cell pairing. Restrict to proteins the information map flags as RNA-explainable and region-stable.

## N5. Shared vs modality-specific latent factors

Fit a two-view factor model (MOFA+, or MEFISTO with spatial smoothness) on paired data, then decompose variance into **shared** factors (RNA and protein both load) and **modality-specific** factors. Then ask which factors are spatially smooth, and which transfer across slides (by projecting held-out slides). This gives a complementary, model-based version of E1.

## N6. Counterfactual neighbourhoods

Use Celcomen-style or Cellina-style interventions (methods catalog E). Rewire a tumour cell's neighbourhood (e.g. remove nearby T cells, or knock out *IFNG* in neighbours) and see how predicted protein changes. Validate the direction against real cells with naturally different neighbourhoods (matched by cell type and own RNA). This turns E3/E4 from "does context help prediction" into "what happens to protein when context changes".

## N7. Token-identity probe for H1 (cheap add-on)

Freeze Nicheformer and Spatium embeddings and fit linear probes from each to (a) the other modality's matched feature and (b) cell type. If Spatium embeddings predict cognate gene expression on paired data better than PCA on protein does, pretraining has encoded cross-modal knowledge. This sits between H1's A1 and A2 ablations and needs no fine-tuning.

## Pitfalls to plan for

- **Composition confounding:** report everything per cell type. Pooled R² mostly measures cell-type differences (the v1 benchmark's main lesson).
- **Leakage:** split by patient or slide. With a single slide, use spatial blocks wider than the largest radius tested.
- **Label circularity (from H2):** if cell types are defined from RNA, a protein model that "recovers" them is doing cross-modal prediction, not typing. Use protein + morphology labels for the main claim.
- **Panel mismatch across datasets:** for E7 and H1, restrict to the shared gene × protein subset and say so.
- **Segmentation differences:** same-cell data (Xenium Protein) avoids most of this; for spot data (DBiTplus, spatial-CITE-seq, Visium CytAssist), run at spot level, or deconvolve first and say which.
