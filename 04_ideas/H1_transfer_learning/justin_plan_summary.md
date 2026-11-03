# H1: Cross-modal transfer between spatial transcriptomics and proteomics (summary)

*Summary of Justin's "H1 Benchmark Plan" (Oct 28, 2026, 13 pp.), shared in the BIOINF 590 Project chat. The original PDF is the reference; this page is a short index of it.*

**Hypothesis.** A model pretrained on one spatial modality, then fine-tuned on the other, beats (i) the same model trained from scratch and (ii) the best method built for the target modality. The test runs in both directions: SPT → protein and protein → SPT.

**Primary contrasts (pre-registered, held-out patients).**

| ID | Comparison | Question |
|---|---|---|
| C1 | SPT-pretrained vs same architecture from scratch, on protein | Does SPT pretraining help protein tasks at all? |
| C2 | SPT-pretrained vs best protein method (native or protein-pretrained), on protein | Is the help large enough to matter? |
| C3 | Protein-pretrained vs from scratch and vs best SPT method, on SPT | Does transfer work in reverse? |
| C4 | Transfer efficiency SPT → protein vs protein → SPT | Is transfer asymmetric? |

**Models.** Nicheformer (SPT; rank tokens, 110M-cell corpus) and Spatium (protein; ~239-marker vocabulary, 51.8M cells). Both are 12-layer transformers pretrained with masked token prediction, with public weights. Random-init twins of each are the scratch controls. Novae is the fallback if Nicheformer will not run by the end of Week 1.

**Tasks.** T1 cell-type annotation (macro-F1, balanced accuracy); T2 niche/domain identification (ARI, NMI). Modes: frozen embeddings + linear probe / Leiden, and fine-tuned at 5 and 20 labelled cells per class.

**Baselines.** MAPS (protein typing), scANVI (SPT typing), CellCharter (niches), PCA + logistic regression / Leiden floor.

**Data.** SPATCH (same tumours measured with Xenium 5K, CosMx 6K and 16-marker CODEX on adjacent sections; colon, liver, ovarian) and Schürch CRC CODEX (56 markers). Start with colon only if time is short. Protein markers are mapped to HGNC genes with MaxFuse's antibody-to-gene table. Pretraining-leakage check against Nicheformer's corpus list and Spatium's 70-study list.

**Ablations of *what* transfers.** A1 re-initialise mapped token embeddings (feature identity); A2 map markers to random gene tokens (negative control); A3 drop mapped markers to 100/75/50/25% (panel overlap).

**Timeline.** Week 1 both models embed SPATCH colon (go/no-go) → Week 2 first matrix heatmap (frozen) → Week 3 fine-tuning + A2/A3 → Week 4 A1 (if a cross-modal gain appeared), Holm correction, report. Roughly 100 short fine-tuning runs, 20–40 GPU-hours.

**Stretch goals.** Masked-feature imputation (T3); Sorin LUAD IMC with CosMx NSCLC; more models (Novae, HEIST, KRONOS2); target-only pretraining; controlled pretraining of one small architecture.

**How v2 extends this.** The input-level information map ([`../I1_information_map/`](../I1_information_map/) and [`05_experiment_design`](../../05_experiment_design/)) predicts *which* proteins and cell types should show C1/C2 gains. A1–A3 are the model-level counterparts of our E2 gene ablations.
