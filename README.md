# BIOINF 590 Final Team Project

**What does spatial RNA know about spatial protein? Finding the informative and transferable information between spatial transcriptomics and spatial proteomics**

Spatial transcriptomics (SPT) and multiplexed protein imaging produce the same kind of object: segmented cells with coordinates, neighbours and labels. They differ in the features measured (thousands of sparse transcript counts vs ~16–60 dense protein intensities). This project asks:

1. **Which information in spatial RNA predicts spatial protein?** The candidates are cell identity, the protein's own transcript, gene programmes, and the RNA of neighbouring cells.
2. **Does that information transfer** to new regions, patients, tissues and platforms?
3. **Can protein-only models benefit from it?** Either through pretraining on RNA (H1) or through RNA used only during training (H2).
4. **What protein spatial biology does RNA *not* explain?**

Start with [`docs/01_project_concept.md`](docs/01_project_concept.md).

---

## Three layers

| Layer | Question | How | Where |
|---|---|---|---|
| **Measure** | Which RNA inputs carry information about each protein, at what spatial scale, and does it hold on held-out regions/slides? | One fixed predictor; ablate only the **inputs**: genes, neighbourhood radius, within-type neighbour shuffles, region hold-out, all-RNA shuffle control | [`docs/02_ablation_experiments.md`](docs/02_ablation_experiments.md) · [`experiments/ablation/`](experiments/ablation/) (runnable) |
| **Transfer** | Can protein-only models use RNA knowledge? | **H1** pretrained Nicheformer (SPT) ↔ Spatium (protein), both directions · **H2** teacher (protein + RNA) → protein-only student, with shuffled-RNA and self-distillation controls | [`docs/plans/`](docs/plans/) · [`experiments/h1_transfer/`](experiments/h1_transfer/), [`experiments/h2_privileged/`](experiments/h2_privileged/) |
| **Explain** | What is left? | Protein minus RNA-predicted protein → within-cell-type spatial variability (reusing our v1 ctSVP pipeline) | [`experiments/ctsvp_baseline/`](experiments/ctsvp_baseline/) · ideas N1–N2 in [`docs/04_new_ideas.md`](docs/04_new_ideas.md) |

The measure layer predicts where H1/H2 should win: proteins and cell types whose variance RNA explains and whose RNA→protein mapping holds on held-out regions.

---

## Repository layout

```
BIOINF-590-Final-Team-Project/
├── README.md
├── docs/
│   ├── 01_project_concept.md          ← start here: question, H1 vs H2, the measurement layer
│   ├── 02_ablation_experiments.md     ← input-level ablation design (E1–E8) + simulation results
│   ├── 03_methods_catalog.md          ← SPT ↔ protein methods: integration, prediction, matching,
│   │                                     foundation models, neighbourhood/counterfactual models
│   ├── 04_new_ideas.md                ← N1–N7: spatial-variability transfer, residual ctSVPs, ...
│   ├── plans/                         ← summaries of Justin's H1 and H2 plans
│   └── archive/                       ← v1 literature review and plan (ctSVP biomarkers)
├── datasets/
│   ├── README.md                      ← datasets ranked by RNA–protein pairing (same cell → none)
│   ├── datasets_catalog.csv           ← every dataset: link, pairing, v2 role
│   └── download_scripts/              ← incl. download_paired_multiomics.md (Xenium Protein, DBiTplus, SPATCH)
├── experiments/
│   ├── ablation/                      ← RUNNABLE: xmodal_ablation.py + simulation self-check + results
│   ├── h2_privileged/                 ← planned (H2 scaffold goes here)
│   ├── h1_transfer/                   ← planned
│   └── ctsvp_baseline/                ← RUNNABLE v1 pipeline (ctSVP detection), reused for residual maps
└── algorithms/                        ← v1 method notebooks (ctSVG, protein-native methods, KRONOS/VirTues)
```

---

## Key datasets

| Pairing | Dataset | RNA | Protein | Role |
|---|---|---|---|---|
| same cell | **Xenium Protein RCC / ccRCC** (10x, 2025) | 377-gene panel + added cognate genes | 27-plex | main ablation data; H2 stage 1 |
| same section | **DBiTplus** (Nat Methods 2026) | whole transcriptome (spots) | ~50-plex CODEX | cross-dataset transfer (E7) |
| same section | **Spatial-CITE-seq** (GSE213264) | whole transcriptome | 273 ADTs | high-plex transfer target |
| same section | Visium CytAssist Gene + Protein tonsil; SPOTS; Stereo-CITE-seq | transcriptome (spots) | ~20–35 ADTs | spot-level checks |
| adjacent | **SPATCH** (Nat Commun 2025): colon, liver, ovarian | Xenium 5K, CosMx 6K, Visium HD, Stereo-seq | 16-marker CODEX | H1 main data; H2 stage 2; registration sweep |
| protein only | Schürch CRC CODEX (56 markers), Jackson/Danenberg IMC, Keren MIBI, … | — | 30–60 | H1 targets; H2 deployment |

Full list: [`datasets/README.md`](datasets/README.md).

## Key methods

SpatialGlue, MultiGATE and PRESENT for joint RNA + protein embedding. DGAT and STProtein for spatial RNA → protein prediction; totalVI and scArches as non-spatial floors. MaxFuse for weak-linkage matching, with STalign/PASTE2 for registration. Nicheformer, Spatium, HEIST, KRONOS and VirTues as foundation models. MISTy, NCEM, Celcomen and Cellina for neighbourhood effects and counterfactuals. Details: [`docs/03_methods_catalog.md`](docs/03_methods_catalog.md).

## Quick start

```bash
cd experiments/ablation
pip install numpy pandas scipy scikit-learn matplotlib
python run_simulation_demo.py          # ~35 s; writes results_sim/ (CSVs + ablation_overview.png)
```

## History

- **v1 (Oct 2026):** cell-type-specific spatially variable protein (ctSVP) biomarker detection, anchored on Li, Ge & Jiang's ctSVG benchmark. Kept: `experiments/ctsvp_baseline/`, `algorithms/`, protein cohorts in `datasets/`. The literature review and plan are in `docs/archive/`.
- **v2 (Nov 2026):** RNA ↔ protein information and transfer, from the team's H1/H2 discussion.

## Team

BIOINF 590, University of Michigan. Puzhu Wang, Meng-Jia "Justin" Lian, Nam Nguyen-Hoang. Add roles here.

Each dataset and algorithm keeps its original license; please cite the original papers.
