# BIOINF 590 Final Team Project

**What does spatial RNA know about spatial protein? Finding informative and transferable information between spatial transcriptomics and spatial proteomics**

Spatial transcriptomics (SPT) and multiplexed protein imaging produce the same kind of object: segmented cells with coordinates, neighbours and labels. They differ in what is measured: thousands of sparse transcript counts vs ~16–60 dense protein intensities. We ask:

1. **Informative:** which RNA information predicts each protein: cell identity, its own transcript, gene programmes, or neighbouring cells' RNA? At what spatial scale?
2. **Transferable:** does that information hold in new regions, patients, tissues and platforms, and under imperfect RNA–protein pairing?
3. **Why not:** when it fails to transfer, is it biology (translation, stability, PTMs, local signalling) or technique (shift, pairing error, platform, spillover)?
4. **Useful:** can protein-only models gain from RNA, through pretraining (**H1**) or RNA used only during training (**H2**)?

All with **uncertainty quantification** at the right unit (patients/blocks): cluster bootstrap, sign-flip and toroidal permutation tests, conformal prediction, shift decomposition.

**→ Start at [`00_START_HERE/`](00_START_HERE/).**

---

## Repository layout (read in order)

```
00_START_HERE/              concept, roadmap, 4-week plan, glossary
01_literature_review/       algorithms by family, each with INPUT FORMAT, output, role, pitfalls
│  ├── algorithms/          A joint embedding · B RNA→protein prediction · C matching/registration ·
│  │                        D foundation models · E neighbourhood/counterfactual · F attribution · G spatial variability
│  ├── notebooks_v1/        v1 template notebooks (ctSVG, protein-native methods, KRONOS/VirTues)
│  └── archive_v1/          v1 literature review + plan
02_datasets/                datasets by RNA–protein pairing tier + dataset_playbook.md + catalog + download scripts
03_assumptions/             A1–A8, one file each + why_information_does_not_transfer.md (mechanisms, decision tree)
04_ideas/                   one folder per idea: I1 information map · H1 transfer · H2 privileged info · N1–N7
05_experiment_design/       E1–E8 one file each + evaluation protocol + transferability diagnostics
06_uncertainty_quantification/  U1 bootstrap · U2 permutation · U3 conformal · U4 model uncertainty · U5 multiplicity/power
07_code/                    ablation/ (information map + UQ, runnable) · ctsvp_baseline/ (v1, runnable) · h1/h2 stubs
08_results/                 simulation self-checks (figures, CSVs, logs); real-data results go here
```

## Quick start

```bash
cd 07_code/ablation
pip install numpy pandas scipy scikit-learn matplotlib
python run_simulation_demo.py   # experiments E1-E6, E8 on simulated data with known answers  (~1 min)
python run_uq_demo.py           # bootstrap, permutation, conformal, shift decomposition       (~2 min)
python make_uq_figure.py
```

The simulation self-check recovers every planted information source, for example:
- cognate gain for CD3: +0.31 R² (95% CI 0.29–0.32);
- PD-L1's neighbourhood effect, peaking at ~50 µm and compositional;
- pS6's RNA-independent spatial residual.

It also shows two pitfalls. Naive cell-level permutation nulls are up to 13× too narrow. Same-section pairing error above ~10–20 µm erases cell-state information. Details: [`08_results/`](08_results/).

## Key datasets (full list in [`02_datasets/`](02_datasets/))

| Pairing | Dataset | RNA | Protein |
|---|---|---|---|
| same cell | 10x **Xenium Protein** RCC / ccRCC | 377-gene panel + added cognate genes | 27-plex |
| same cell | Bruker **CosMx Breast Multiomic** | whole transcriptome (~18k) | 64-plex |
| same section | **DBiTplus** (GSE308167), **spatial-CITE-seq** (GSE213264), Visium CytAssist Gene + Protein, SPOTS (GSE198353) | transcriptome (spots/pixels) | 30–273 |
| adjacent | **SPATCH** (COAD/HCC/OV; Xenium 5K, CosMx 6K, Visium HD, Stereo-seq) | 5–6k genes | 16-marker CODEX |
| protein only | Schürch CRC CODEX, Jackson/Danenberg IMC, Keren MIBI, Hickey CODEX, HTAN, spora | – | 30–60 |

## Key methods (full write-ups in [`01_literature_review/`](01_literature_review/))

- **Integration:** SpatialGlue, MultiGATE, PRESENT, totalVI, WNN.
- **Prediction:** DGAT, STProtein, totalVI/scArches.
- **Matching/registration:** MaxFuse, STalign, PASTE2, Tangram.
- **Foundation models:** Nicheformer, Spatium, HEIST, Novae, KRONOS, VirTues.
- **Context/counterfactual:** MISTy, NCEM, BANKSY, CellCharter, STELLAR, Celcomen, Cellina.
- **Spatial variability:** CELINA, SPARK-X, STANCE, and others (v1).

## History

- **v1 (Oct 2026):** cell-type-specific spatially variable protein (ctSVP) detection, anchored on Li, Ge & Jiang's ctSVG benchmark. Kept as `07_code/ctsvp_baseline/`, `01_literature_review/notebooks_v1/` and `archive_v1/`.
- **v2 (Nov 2026):** RNA ↔ protein information and transfer (team H1/H2 discussion), the information map, UQ, and this progressive structure.

## Team

BIOINF 590, University of Michigan: Puzhu Wang, Meng-Jia "Justin" Lian, Nam Nguyen-Hoang. Each dataset and algorithm keeps its own license; please cite the original papers.
