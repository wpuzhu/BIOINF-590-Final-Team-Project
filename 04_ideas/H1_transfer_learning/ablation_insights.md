# Which part of RNA-pretraining knowledge transfers to protein? Insights from the A1–A3 testbed

*Sub-branch `h1-ablation-identity-structure`. Code: [`07_code/h1_transfer/ablation_identity_structure/`](../../07_code/h1_transfer/ablation_identity_structure/). Numbers: [`08_results/simulation/h1_ablation/`](../../08_results/simulation/h1_ablation/). Everything below except the "real data" section comes from **simulation**: it shows what the design can detect and what to expect, not what is true of Nicheformer and Spatium.*

## The four hypotheses

When an RNA-pretrained model helps on protein, the help could come from:

| | Hypothesis | Test |
|---|---|---|
| **H-identity** | the model knows what each gene *means*, and the protein inherits it through its cognate token | A1 (re-init mapped tokens), A2 (wrong map) |
| **H-overlap** | the help scales with how much of the panel the model's vocabulary covers | A3 (fraction mapped, information held constant) |
| **H-structure** | the encoder learned how cell profiles are organised, independent of names | A1 − shuffled-corpus pretraining |
| **H-generic** | any pretrained network is a better starting point or feature map than a random one, with no biology involved | shuffled-corpus pretraining − scratch |
| **H-null** | nothing transfers | total gain ≈ 0 |

## 1. The design works: it attributes a real gain correctly

On the **in-modality positive control** (pretrained on RNA, fine-tuned on new RNA, 5 labels per class), pretraining helps by +0.076 macro-F1 [0.067, 0.086]. The decomposition assigns essentially all of it to **identity** (+0.077 [0.068, 0.086]). Structure contributes +0.008 (CI includes 0) and generic −0.009. A wrong gene map costs 0.107.

So when transfer exists, A1/A2 find where it comes from: here, the token table. The gain fades with more labels (+0.077 → +0.044 → +0.015 at 5 → 20 → 100 labels), as expected for prior knowledge.

## 2. On protein targets, the dominant outcome is H-null, with a small H-generic component

Fine-tuned, 16-marker protein panel, 5 labels per class:

| Scenario | total | identity | generic |
|---|---|---|---|
| coupled (protein ≈ gene, ρ = 0.9) | +0.006 (CI includes 0) | −0.006 | +0.013 |
| realistic (CODEX-like) | +0.001 | −0.009 | +0.014 |
| decoupled (ρ = 0.15) | −0.001 | **−0.021** | +0.021 |

- **Even when protein follows its gene almost perfectly, the pretrained gene meanings do not help.** The identity component is ≈ 0 or slightly negative in every protein scenario.
- **When coupling is weak, the correct gene map is actively harmful** (decoupled: identity −0.021, and the true map is worse than a *random* map by 0.017). The cognate prior points the model the wrong way and fine-tuning has to undo it.
- **The only positive component is generic (+0.01 to +0.02).** A pretrained network is a slightly better-conditioned starting point than a random one, and the identity loss cancels it out.
- **A3 shows no dose-response.** Mapping 25%, 50%, 75% or 100% of markers to pretrained tokens makes no difference (gain −0.003 to +0.011). There is no overlap effect because there is no identity effect to scale.

## 3. Why identity fails on protein: panel size, not modality

Two controls isolate the cause:

- **Same modality, protein-sized panel.** Restrict the RNA target to the 16 cognate genes only, and the in-modality pretraining gain (+0.076 with 300 genes) drops to −0.002. Pretrained knowledge helps when the input is **high-dimensional and noisy relative to the labels**. A 16-feature panel gives it nothing to denoise or fill in.
- **Bigger protein panel.** With 60 coupled markers, identity knowledge becomes measurable, but **only in frozen mode**: pretrained embeddings beat re-initialised ones by +0.104 [0.086, 0.123], and the true map beats a wrong map by +0.095. Fine-tuning removes the advantage (identity −0.006), and even the best frozen embedding stays **0.12 below a PCA probe** on the raw markers.

**Mechanism:** a cognate token's pretrained embedding encodes how that gene behaves among ~40 expressed genes out of 300. A protein cell presents 16 dense, rank-weighted markers, a different input distribution. Small, clean protein panels are already near-optimally summarised by PCA, so prior knowledge has little to add.

## 4. How much of each hypothesis is likely true for the real H1 benchmark?

Combining the testbed with Justin's frozen-mode results ("H1 Benchmark Results" deck, 5 Nov 2026; plan in [`justin_plan_summary.md`](justin_plan_summary.md)):

| Real-data observation (Justin) | Testbed counterpart | Agreement |
|---|---|---|
| PCA probe 0.82 on CODEX and Xenium; no foundation model beats it with all labels | frozen pretrained 0.10–0.12 below PCA in every scenario | ✔ same direction |
| Nicheformer on CODEX: −0.30 (5 labels), −0.12 (all labels) vs best baseline | frozen 16-marker cross-modal: −0.10 vs PCA, identity ≈ 0 | ✔ (real gap larger) |
| Nicheformer on Xenium (in-modality): +0.039 [+0.005, +0.076] at 5 labels, gone with all labels | in-modality gain +0.076 at 5 labels → +0.015 at 100, carried by identity | ✔ same shape |
| Spatium reads only 135 of 5,001 Xenium genes and scores 0.56 vs 0.80 | protein-sized panel kills pretraining benefit even in-modality | ✔ consistent explanation |

**Best current estimate for 16-marker CODEX targets** (to be confirmed with the random-init controls Justin listed as next steps):

| Hypothesis | Likelihood | Expected size |
|---|---|---|
| **H-null** (no meaningful cross-modal gain after fine-tuning) | **most likely** | total gain ≈ 0 ± 0.02 macro-F1 |
| **H-generic** (pretrained = better initialisation) | likely, small | +0.01–0.02 at 5 labels, ~0 with all labels |
| **H-identity** | unlikely on 16 markers; possible on high-plex panels in frozen mode | ≈ 0 (could be *negative* for weakly coupled state/phospho markers); +0.05–0.10 frozen on ≥50-marker panels |
| **H-overlap** | no independent effect expected | flat A3 curve on 16 markers |
| **H-structure** (cell-profile) | unlikely for Nicheformer | ≈ 0 |
| **H-structure** (spatial organisation) | **untested**; the one route that could still work | Justin's niche results (Nicheformer on CODEX niches 0.29 vs PCA 0.31) suggest small |

The asymmetry Justin saw (Nicheformer → CODEX 0.70 vs Spatium → Xenium 0.56) is best explained by **vocabulary coverage**, not a deeper asymmetry in transfer. Spatium sees 2.7% of the Xenium genes (135 of 5,001), whereas most CODEX markers have a gene token in Nicheformer's 20,310-gene vocabulary.

## 5. Recommendations for the real experiments

1. **Run the random-init controls first** (Justin's next step 1). If total ≈ 0, report H-null with CIs, and A1–A3 become moot for 16-marker CODEX.
2. **Move the identity test to a high-plex panel.** Run A1/A2 on **Schürch CRC CODEX (56 markers)** or **spatial-CITE-seq (273 ADTs)**, in **frozen mode**, where the testbed says identity effects are detectable (+0.10). SPATCH CODEX (16 markers) is under-powered for this question.
3. **Stratify A2 by marker coupling.** Use our information map (I1, ablation E2) to split markers into strongly vs weakly gene-coupled. The testbed predicts the true map helps only for strongly coupled markers and *hurts* for weakly coupled ones. That would be a direct cross-check between I1 and H1.
4. **Test spatial structure explicitly** with a graph model (HEIST or Novae): random-init vs pretrained, with and without shuffled cell coordinates. It is the only hypothesis this testbed cannot address.
5. **Statistics:** with 3–5 seeds, sign-flip p-values cannot reach significance after Holm. Pre-register 2–3 primary contrasts (total gain at 5 labels; identity on the high-plex panel) and report effects with patient-level bootstrap CIs.

## Caveats

The testbed uses an embedding-bag MLP, not a transformer, and simulated, near-linear biology. Effect sizes will differ from Nicheformer/Spatium; the *pattern* of contrasts is the transferable lesson. The CIs come from 15 seed × patient units that share 3 test patients, so they are somewhat optimistic.
