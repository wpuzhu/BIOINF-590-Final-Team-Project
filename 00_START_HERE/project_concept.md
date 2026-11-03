# Project concept (v2): what RNA knows about protein, in space

*BIOINF 590 Final Team Project. Builds on Justin's H1/H2 plans (Oct 23 and Oct 28, 2026; summaries in [`04_ideas/H1_transfer_learning/`](../04_ideas/H1_transfer_learning/) and [`04_ideas/H2_privileged_information/`](../04_ideas/H2_privileged_information/)).*

## The one-sentence version

Spatial transcriptomics (SPT) and spatial proteomics measure the same tissues in the same data shape. We want to find out **which information in spatial RNA is also informative about spatial protein, and whether that information transfers** to new slides, tissues and platforms where only protein is measured.

## Why the two modalities can share models

Once cells are segmented, both modalities produce the same object (from Justin's "Why SPT and Spatial Proteomics Can Share Models"):

| Part | SPT | Spatial proteomics | Same? |
|---|---|---|---|
| row | one segmented cell | one segmented cell | yes |
| columns | genes (hundreds to ~6,000) | markers (~16–60) | **differs** |
| values | sparse transcript counts | dense mean intensities | **differs** |
| location, neighbours, labels, container | (x, y), kNN graph, cell type/niche, AnnData | same | yes |

Most protein markers have a gene partner (CD3 → *CD3E*, CD8 → *CD8A*, Ki-67 → *MKI67*, E-cadherin → *CDH1*), and tissue structure (who sits next to whom) does not depend on the assay. That shared structure is why transfer is plausible.

They also differ in ways that limit transfer. mRNA and protein only partly agree, because translation, degradation, secretion and post-translational modification decouple them. Protein panels are small, and phospho-epitopes, pan-cytokeratin and DNA stains have no clean gene partner. Noise models differ too: count models built for RNA do poorly on protein intensities.

## Two ways RNA could help a protein model (Justin's H1 and H2)

In both, the finished model sees **only protein** when it is used on new slides.

| | **H1: transfer learning** ("textbook") | **H2: privileged information** ("tutor") |
|---|---|---|
| Source of RNA | public corpora, other tissues and patients | our own tissue: same cells (Xenium Protein) or the adjacent section (SPATCH) |
| When RNA is used | before training, through pretrained weights | during training only, through a teacher |
| Model | foundation models: Nicheformer (SPT) ↔ Spatium (protein) | small graph network (GraphSAGE/GAT); the teacher sees protein + RNA, the student sees protein only |
| Decisive control | the same architecture trained from scratch, and the best protein-native method | a teacher fed RNA shuffled across cells; self-distillation (protein teacher) |
| Main risk | getting each model to run on the other modality | registration error between adjacent sections |
| Data | SPATCH CODEX + Xenium 5K/CosMx 6K; Schürch CRC CODEX | Xenium Protein demos, then SPATCH, then an unpaired CODEX cohort |

The two can be combined: start H2's student from H1's pretrained weights, then distil from the RNA teacher. That tests whether general (corpus) RNA knowledge and sample-specific RNA knowledge add up. Justin's suggested options for a 4-week project are A (H1 as the main study), B (H2 as the main study) and C (H2 plus one H1 arm).

## What this version adds: measure the information first, at the input level

H1 and H2 both answer *"does RNA help?"* with a yes/no on a downstream task. Neither, by itself, says **which part of the RNA** carried the help. The team asked for that, so we add a measurement layer that sits underneath both hypotheses:

1. **Measure.** On paired data, fit one fixed, simple predictor from RNA-derived inputs to each protein. Then ablate the inputs: genes (the cognate gene, gene programmes, random genes), spatial context (neighbourhood radius, neighbours shuffled within cell type), and tissue regions (train on some regions, test on others). The result is an **information map**: for each protein and each cell type, how much comes from cell identity, its own transcript, other genes, and the neighbourhood, and how much of that holds up on held-out regions and patients. See [`04_ideas/I1_information_map/`](../04_ideas/I1_information_map/) and [`05_experiment_design/`](../05_experiment_design/); code is in [`07_code/ablation/`](../07_code/ablation/).
2. **Transfer.** H1 and H2 then test whether a protein-only model can actually use that information. The information map predicts where they should win: proteins and cell types with high RNA-explainable, region-stable variance.
3. **Explain what does not transfer.** The protein signal RNA cannot explain (the residual) is biology in its own right: post-transcriptional or microenvironment-driven regulation. We test whether the residual is still spatially structured *within* a cell type. See idea [N2](../04_ideas/N2_residual_spatial_structure/) and experiment [E6](../05_experiment_design/E6_residual_maps.md).

Layers 1 and 3 need only paired data and a laptop, so they make a safe core result. H1 and H2 are the higher-risk, higher-reward modelling arms.

## Research questions

- **Q1 (informative).** For each protein × cell type, what fraction of held-out variance is explained by (a) cell type, (b) the cognate transcript, (c) the rest of the transcriptome, (d) the RNA neighbourhood? How does this change with neighbourhood radius?
- **Q2 (spatial or compositional).** When the neighbourhood helps, is it because of *which cell types* are nearby, or *what state* those cells are in? This is tested by shuffling neighbours within cell type.
- **Q3 (transferable).** Does the RNA→protein mapping hold across regions (tumour core → margin → stroma), patients, tissues and platforms (Xenium Protein → DBiTplus/spatial-CITE-seq → SPATCH)?
- **Q4 (useful).** Do protein-only models gain from RNA through pretraining (H1) or distillation (H2), and is the gain largest where Q1–Q3 say the information is?
- **Q5 (residual biology).** Which proteins keep within-cell-type spatial structure after RNA is accounted for?

## How results will be judged

- Splits are always by patient or slide, never by cell. With one slide, use spatial blocks.
- Every positive result needs its negative control: RNA shuffled across cells (should lose everything), neighbours shuffled within cell type (separates composition from state), and random gene sets of equal size (separates "those genes" from "any genes").
- Effects are reported as differences with bootstrap CIs over patients, with Holm correction across conditions (as in the H2 plan).
