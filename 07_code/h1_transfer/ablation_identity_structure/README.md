# H1 ablation testbed: identity vs structure vs panel overlap

A controlled, laptop-scale version of the H1 ablations (A1–A3 in Justin's plan). It answers the question: *when a model pretrained on RNA is used on protein, which kind of pretrained knowledge, if any, carries over?*

Insights: [`../../../04_ideas/H1_transfer_learning/ablation_insights.md`](../../../04_ideas/H1_transfer_learning/ablation_insights.md) · Results: [`../../../08_results/simulation/h1_ablation/`](../../../08_results/simulation/h1_ablation/)

## Why a testbed (and not Nicheformer/Spatium directly)

The real models need a GPU and the SPATCH data, and their pretraining corpora are fixed, so the truth about what *should* transfer is unknown. Here a single simulated latent cell state generates both modalities, and we **set** how strongly each protein follows its gene. That lets us check whether the ablations recover a known answer before trusting them on real data.

## Components

| File | Contents |
|---|---|
| `sim.py` | simulated biology: 12 fine cell types (4 lineages × 3 subtypes), 24 latent programmes, 300 genes (sparse Poisson counts), a protein panel of 16 (or 60) markers. Marker *j* = ρⱼ · (its cognate gene's noise-free expression) + (1−ρⱼ) · (protein-specific readout of the same latent state) + noise + patient batch scale. |
| `model.py` | a small "foundation model": **layer 1 = token-embedding table** (what each gene token means → *identity*), **layer 2 = shared encoder** (how evidence combines → *structure*). Rank-based tokenisation like Nicheformer/Spatium. Self-supervised pretraining, the arms, fine-tuning, frozen probing. |
| `run.py` | runs one scenario × 5 seeds × {5, 20, 100} labels/class × all arms, in fine-tuned and frozen modes; optional in-modality RNA positive control |
| `analyze.py` | paired contrasts → decomposition, bootstrap CIs, sign-flip tests, Holm; A3 overlap curves; figure |

## Arms

| Arm | What the protein model starts from | Removes |
|---|---|---|
| `scratch` | random weights (fine-tuned) / random features (frozen) | everything (H1 control C1) |
| `noise_pre` | pretrained on RNA with each gene **shuffled across cells**: same marginals, no co-expression | biology; keeps "a pretrained network" |
| `pretrained` | RNA-pretrained; markers mapped to their cognate gene tokens | nothing |
| `A1_reinit` | pretrained encoder kept; embedding rows of the mapped tokens re-initialised | gene **identity** |
| `A2_random_map` | pretrained; each marker mapped to a **random wrong** gene token | **correct** identity (misleading prior) |
| `A3_overlap_q` | pretrained; only a fraction q of markers mapped to cognate tokens, the rest get new tokens | vocabulary **overlap**, with information content held constant (no markers dropped) |
| `PCA_probe` | PCA of the protein + logistic regression | reference bar, as in Justin's benchmark |

**Decomposition** (exact): total = pretrained − scratch = identity (pretrained − A1) + structure (A1 − noise_pre) + generic (noise_pre − scratch).

## Scenarios

| Scenario | Coupling ρ | Panel | Truth |
|---|---|---|---|
| `coupled` | 0.9 all markers | 16 | protein ≈ its gene: identity *could* transfer |
| `realistic` | lineage 0.85, state 0.3; 4 markers without a gene partner | 16 | CODEX-like |
| `decoupled` | 0.15 all markers | 16 | the gene is a misleading guide; only latent state is shared |
| `coupled_P60` | 0.9 | 60 | high-plex panel (like Schürch 56, spatial-CITE-seq 273) |

Positive controls (in `coupled`): **RNA target** (in-modality, 300 genes, where pretraining should help) and **RNA restricted to the 16 cognate genes** (same modality, protein-sized panel).

## Run

```bash
pip install numpy pandas scipy scikit-learn matplotlib
python run.py --scenario coupled --seeds 5 --rna_target
python run.py --scenario realistic --seeds 5
python run.py --scenario decoupled --seeds 5
python run.py --scenario coupled --P 60 --tag coupled_P60 --seeds 5
python analyze.py                       # -> 08_results/simulation/h1_ablation/
```
About 2 minutes per scenario on 2 CPU cores (the `coupled` run with the RNA controls takes ~3 minutes).

## Design decisions and what was tried first

- **scikit-learn, not PyTorch:** PyTorch couldn't be installed in the build environment. An embedding-bag MLP keeps exactly the two ingredients the ablations need: a token table (identity) and a shared encoder (structure). It is not a transformer; see limitations.
- **Pretraining objective:** a raw masked-token reconstruction target was tried first. It is dominated by Poisson noise, stopped improving after a few epochs, and gave **no** in-modality benefit. The final objective predicts each cell's *denoised* (top-30-PC) expression from 70% of its tokens. It is still label-free, and it gives a clear in-modality few-label gain (+0.076 macro-F1 at 5 labels), which is the precondition for the ablations to mean anything.
- **Simulation strength** was tuned so RNA tokens classify 12 fine types at ~0.70 accuracy (fine types are hard, as in real data).
- **A3 keeps all markers** and only changes which ones use pretrained tokens. Dropping markers would confound "less overlap" with "less information".
- **The pretraining-loss line in the logs** compares losses on different targets (the shuffled corpus has a lower-variance target), so the two numbers aren't comparable.

## Limitations

- Toy model and toy biology. Effect *sizes* will not match Nicheformer/Spatium; use the *pattern* of contrasts.
- **The model is single-cell (no spatial graph),** so "structure" here means cell-profile structure, not tissue organisation. Spatial-structure transfer needs a graph model (HEIST, Novae).
- **Units are seed × held-out patient (15)** but share patients, so the CIs are somewhat optimistic. The sign-flip test over 5 seeds has a minimum two-sided p of 0.0625, so use the CIs for inference.
- Simulated protein is close to linear in the latent state, which favours the PCA probe.
