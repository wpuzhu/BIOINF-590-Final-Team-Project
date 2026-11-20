# Concept and design: two stages, informative vs transferable

## The question

Spatial transcriptomics (SPT) and spatial proteomics measure two views of the same tissue. We want to know **which information in one view helps with the other**, and whether that help **survives when only one view is available**. These are two different properties:

| Concept | Plain meaning | Needs RNA at test time? | Measured as |
|---|---|---|---|
| **Informative** | Adding RNA to protein improves the prediction on the same cells | **yes** | Δ_info = perf(protein + RNA model) − perf(protein-only model) |
| **Transferable** | A protein-only model trained *with RNA's help* is better on new slides/patients that have no RNA | **no** | Δ_trans = perf(protein-only student taught by RNA) − perf(protein-only baseline) |
| Transfer efficiency | Share of the informative gain that a protein-only model can keep | — | η = Δ_trans / Δ_info (defined only when Δ_info's CI is above 0) |
| RNA-specific transfer | Part of Δ_trans that needs RNA at all | no | perf(RNA-taught student) − perf(self-taught student, same unlabelled cells, protein-only teacher) |

Information can be informative but not transferable when it lives only in RNA: no function of protein and space can reproduce it, so a student cannot learn it, however it is taught. It can be transferable but not *needed* when protein already contains it and RNA only makes it easier to find with few labels. The simulation below builds all four combinations to check that the designs tell them apart.

Stage 1 uses the same split for **pretraining**: in-modality pretraining gain (RNA→RNA, protein→protein) is how *informative* pretraining is; the cross-modal gain (RNA→protein, protein→RNA) is how much of it *transfers*.

> **About the premise "ST → proteomics will obviously give a large gain".** The evidence so far points the other way. Justin's H1 benchmark had PCA ≥ every foundation model on protein, with Spatium reading 135/5,001 Xenium genes. Our testbed shows ~0 cross-modal gain in both directions while in-modality RNA pretraining gains +0.076 macro-F1 at 5 labels/class. Stage 1 is therefore run as a **pre-registered open test** in which "no gain" (equivalence within ±0.02 macro-F1) is a reportable result. The positive in-modality control is required so that a null is not just a broken pipeline.

---

## Stage 1: cross-modal pretraining transfer (Nicheformer ↔ Spatium)

**Design: 2 × 2 × arms × labels.**

| | target = RNA | target = protein |
|---|---|---|
| **Nicheformer** (RNA-pretrained) | in-modality (positive control) | **cross-modal** |
| **Spatium** (protein-pretrained) | **cross-modal** | in-modality (positive control) |

- Arms: pretrained · scratch (same architecture, random weights) · PCA probe · A1 re-init mapped token embeddings (removes *identity*) · A2 random token map (removes *correct* identity) · A3 overlap q ∈ {0.25, 0.5, 0.75} (vocabulary overlap at constant information).
- Modes: frozen embedding + logistic probe; full fine-tuning.
- Labels: 5, 20, 100 per class, and all. Gains from pretraining are expected only at few labels; at "all" they should vanish (sanity check).
- Datasets: matched-tissue pairs with low pretraining leakage first (SPATCH Xenium 5K ↔ adjacent CODEX; Xenium Protein same-cell), then high-plex (CosMx 108-protein), then patient-powered protein cohorts (Schürch CRC). See [`../datasets/README.md`](../datasets/README.md).
- Decision rules (pre-registered): cross-modal transfer "exists" if C1 (pretrained − scratch) has a 95% CI > 0 **and** beats PCA, at ≥ 1 label level, on ≥ 2 datasets. "Absent" if the 90% CI lies within ±0.02 (TOST). Anything else is "inconclusive" and is reported with its CI.
- What carries over (if anything): identity = pretrained − A1; correct mapping = pretrained − A2; overlap curve from A3.

## Stage 2: RNA as privileged information

### Does a distillation framework make sense?

**Yes, as one of several ways to use training-only RNA, and only under conditions we can test.** A teacher that sees protein + RNA and a student that sees protein + space is *learning using privileged information* (Vapnik & Vashist 2009). Lopez-Paz et al. (2016) unify that with distillation as *generalized distillation*. It can help the student for three reasons only:

1. **More labels in effect**: the teacher labels the unlabelled training cells, which matters only when labels are scarce.
2. **Better labels**: soft labels carry class similarity ("dark knowledge"), and RNA can separate subtypes that protein mixes up.
3. **Representation shaping**: matching the teacher's features or aligning with RNA moves the student's features toward directions that protein *does* contain but that few labels would not find.

It **cannot** give the student information that protein and space do not contain. Xue et al. (ICLR 2023) call this the modality focusing hypothesis: what transfers is the part of the teacher's knowledge that is decisive *and* shared with the student's modality. So distillation is right when the question is "is RNA's extra information recoverable from protein + neighbourhood?". It is not a way to impute RNA-only information.

Spatial context matters twice. **The student's graph** lets it read niche information from neighbours' protein: RNA-defined subtypes that sit in particular niches can then be learned through space. **The pairing** sets what the teacher can see:

| Pairing | Teacher's privileged input | What can transfer |
|---|---|---|
| same cell (CosMx 980/108, Xenium Protein) | the cell's own RNA | cell-intrinsic + niche information |
| same section, re-stained (Xenium → COMET / PhenoCycler) | own RNA after registration (µm error) | as above, minus registration noise |
| adjacent section (SPATCH) | **neighbourhood RNA** around the registered position | niche-scale information only (space is the bridge) |
| unpaired | none at cell level | Stage 1-style transfer only |

### Mechanisms compared (the same student and splits throughout)

| Arm | What it does | Identifies |
|---|---|---|
| `none` | protein + graph, labelled cells only | baseline |
| `teacher` | protein + graph + RNA, **evaluated with RNA** | Δ_info (information ceiling) |
| `kd` | student matches teacher soft labels on all training cells (λ CE + (1−λ) T² KL) | Δ_trans via outputs |
| `kd_feat` | + match teacher hidden features (β MSE) | does representation transfer add anything |
| `align_infonce` / `align_cca` | student embedding aligned with RNA embedding (InfoNCE; CCA-shared target in the simulation) | transfers only the RNA part *shared* with protein |
| `recon_rna` | auxiliary head reconstructs RNA PCs | multi-task, no teacher |
| `recon_specific` | auxiliary head reconstructs the RNA part *not* predictable from protein | negative control: should not help |
| **controls** `self_train` | distil from a **protein-only** teacher | removes RNA: is the gain just the unlabelled pool? |
| `label_smooth` | label smoothing 0.1 | is the gain just regularisation? |
| `kd_permuted` | teacher probabilities permuted among cells with the same teacher label | keeps pseudo-labels, removes dark knowledge |

Pre-registered: pick the winner on `default` by the RNA-specific contrast (`mechanism − self_train`), then run all ablations with the winner (KD in the simulation).

### Step-by-step ablation ladder (input level)

Each step answers one question; later steps are interpreted only if earlier ones pass.

| Step | Manipulation | Question | Reading |
|---|---|---|---|
| 0 | positive controls (in-modality pretraining, `teacher` on a known RNA-defined label) | is the pipeline able to see a gain? | must pass |
| 1 | `teacher − none` | is RNA **informative** for these labels? | CI > 0 needed to continue |
| 2 | `kd − none` | does any of it **transfer** to protein-only? | Δ_trans, η |
| 3 | `kd − self_train`, `label_smooth`, `kd_permuted` | is the transfer from **RNA**, or from unlabelled cells / regularisation / pseudo-labels? | RNA-specific part |
| 4 | teacher RNA nulls: shuffle **within niche** / **within coarse type** / type average / all | **where** is the transferable information: cell-intrinsic (killed by within-type shuffle) vs niche-level (survives within-niche shuffle)? | graded |
| 5 | teacher sees **neighbourhood RNA only**; **pairing jitter** 0–240 µm | the **spatial scale** at which RNA information still transfers (what adjacent sections can deliver) | jitter curve |
| 6 | student **graph off** | does transfer go **through space**? | spatial vs protein route |
| 7 | student **panel** 5/10/20 markers; **drop subtype markers** | how much does the student's own view limit transfer? | ceiling |
| 8 | **label fraction** 1–50% | regime where RNA help matters | should shrink with labels |
| 9 | (real data) different label source (protein-gated vs RNA-clustered) | circularity check | report both |

### Uncertainty quantification per experiment

The unit is **the patient**. Cells within a patient are not independent, so cell-level CIs are far too narrow. Seeds are a second, inner source of variance.

| Quantity | Method | Notes / literature |
|---|---|---|
| Paired gain (any contrast) | **two-stage hierarchical bootstrap**: resample patients, then seeds within patient; 95% percentile CI | cluster bootstrap (Davison & Hinkley 1997; Field & Welsh 2007); hierarchical bootstrap (Saravanan et al. 2020) |
| Test of gain ≠ 0 | **sign-flip permutation** on patient-mean paired differences (exact for ≤ 15 patients) | smallest p = 2^(1−n): with 3 patients p ≥ 0.25 → report as descriptive |
| "No gain" claims | **TOST equivalence**, margin 0.02 macro-F1, 90% CI inside the margin | Schuirmann 1987; Lakens 2017 |
| Many contrasts | Holm within a family (per dataset/scenario), BH across ablation settings | |
| Variance sources | random-effects ANOVA: patient vs seed vs residual share | shows whether more seeds or more patients help; Bouthillier et al. 2021 (variance in ML benchmarks) |
| Transfer efficiency η | patient bootstrap of the ratio of means; undefined if Δ_info's CI includes 0 | ratio CIs are unstable near 0 |
| Jitter / label curves | pointwise CI bands from the same bootstrap; the "transfer distance" is where the lower band crosses 0 | |
| Per-cell confidence | **split-conformal APS** sets, calibrated on held-out *patients*; coverage and set size; Mondrian (per cell type) for class-conditional coverage | Romano, Sesia & Candès 2020; Vovk et al. 2005 |
| Calibration | expected calibration error (ECE) per arm | KD often changes calibration |
| Cohort shift (deployment) | weighted conformal under covariate shift; covariate vs concept shift decomposition (shared/uq/uq_core.py) | Tibshirani et al. 2019 |
| CV caveat | k-fold test sets overlap training sets across folds, so naive t-tests over folds are anti-conservative; that is why we resample patients, not folds | Nadeau & Bengio 2003; Bates, Hastie & Tibshirani 2024 |

Simulation results and their reading: [`../stage2_privileged_transfer/results/simulation/README.md`](../stage2_privileged_transfer/results/simulation/README.md) and [`../stage1_cross_modal_pretraining/results/simulation_ablation_testbed/README.md`](../stage1_cross_modal_pretraining/results/simulation_ablation_testbed/README.md).

## References (methods)

Vapnik & Vashist 2009, *Neural Networks* (LUPI) · Lopez-Paz, Bottou, Schölkopf & Vapnik 2016, ICLR (generalized distillation) · Hinton, Vinyals & Dean 2015 (KD) · Xue, Gao, Shao & Zhao 2023, ICLR (modality focusing hypothesis) · Romano, Sesia & Candès 2020, NeurIPS (APS) · Tibshirani, Barber, Candès & Ramdas 2019, NeurIPS (weighted conformal) · Lakens 2017, *Soc Psychol Personal Sci* (TOST) · Saravanan, Berman & Sober 2020, *Neurons Behav Data Anal Theory* (hierarchical bootstrap) · Bouthillier et al. 2021, MLSys · Nadeau & Bengio 2003, *Mach Learn* · Bates, Hastie & Tibshirani 2024, *JASA* (cross-validation CIs).
