# Datasets

Raw data is **not** stored in git. This folder holds:
- `datasets_catalog.csv`: every dataset with link, technology, pairing, size and role;
- `registry.yaml`: how each dataset is converted into the shared schema ([`../shared/data/schema.md`](../shared/data/schema.md)) by [`../shared/data/harmonise.py`](../shared/data/harmonise.py);
- `download_scripts/`: fetch helpers (write into `$SCRATCH/raw`, never into the repo);
- [`dataset_playbook.md`](dataset_playbook.md): per-dataset input preparation, pitfalls, and why its information may or may not transfer.

Selection rules (checked Nov 2026):
1. **Public**: downloadable without a data-access committee (10x / Bruker forms are fine).
2. **High quality**: peer-reviewed or vendor reference data, single-cell segmentation, ≥ 20 protein markers, or a whole-transcriptome RNA side.
3. **Patients ≥ 4** when possible. Patients are the unit for confidence intervals. Datasets with 1–3 patients are labelled *descriptive*: no patient-level CI.
4. **Pretraining leakage**: Nicheformer was pretrained on SpatialCorpus-110M (57M dissociated + 53M spatial cells, mostly public MERFISH / Xenium / CosMx). Spatium was pretrained on > 51M cells from 70 public proteomics studies (CODEX, IMC, CyCIF). A test set that one model saw during pretraining favours that model. Each row below has a leakage column. **Check the released corpus manifests before the final run**; until then, treat older public sets as *possibly seen*.

---

## Stage 1: cross-modal pretraining transfer (Nicheformer ↔ Spatium)

Stage 1 needs **labelled targets in each modality from the same tissues**. Pairing is not needed, but matched tissues make the comparison fair.

| ID | Dataset | Modality / panel | Tissue · scale | Role in Stage 1 | Leakage risk |
|---|---|---|---|---|---|
| **P4** | SPATCH (Ren et al., *Nat Commun* 2025) | Xenium 5K, CosMx 6K, Visium HD, Stereo-seq **+ CODEX (16 markers) on adjacent sections** + scRNA-seq | colon, HCC, ovarian; manual annotations | **Main**: same tumours in both modalities. RNA target for Spatium→RNA; protein target for Nicheformer→protein | low (released 2025, after both corpora were built) |
| **P1/P2** | 10x Xenium Protein RCC / ccRCC (2025) | 377 (+28/+100) genes + 27 proteins, same cell | kidney cancer, ~465k cells; 1 patient each | RNA and protein targets on identical cells (panel effect without tissue confound) | low (released July 2025); **descriptive** (1 patient each) |
| **P7** | CosMx multiomic NSCLC + breast (He et al. 2022) | 980 RNA + **108 proteins**, same cell | > 800k cells released | **High-plex protein target**: the best test of whether Nicheformer's gene vocabulary helps on a large panel | medium (CosMx 2022 release could be in SpatialCorpus) |
| **D8** | Schürch 2020 CRC CODEX | 56 markers | 35 patients, 140 TMA spots, ~250k cells | Protein target matched to SPATCH colon; the main **patient-powered** protein target | **high for Spatium** (classic CODEX set) |
| **D17a** | Spatial-CITE-seq (Liu 2023) | 273 ADTs + whole transcriptome | human tonsil, skin | High-plex protein target; spot-level | medium |
| **D2/D3** | Jackson 2020 / Danenberg 2022 IMC | 35–37 markers | ~350 / ~700 patients | Patient-powered protein targets | **high for Spatium** |
| RNA | public Xenium 5K (e.g. ovarian cancer, lymph node; CC BY 4.0) | 5,001 genes | 1 sample each | RNA target for Spatium→RNA | **high for Nicheformer** if released before 2024 |

**Priority**: SPATCH (P4) + Schürch (D8) as the main pair, P7 for high plex, P1/P2 as same-cell anchors. Report Nicheformer and Spatium separately on a *low-leakage* subset (P4, P1/P2) and on the full set.

## Stage 2: RNA as privileged information (teacher sees RNA, student sees protein + space)

Stage 2 needs RNA and protein **from the same cell or the same section**. Adjacent sections work only at neighbourhood level.

| Tier | ID | Dataset | RNA / protein | Pairing | Scale | Use |
|---|---|---|---|---|---|---|
| **1a same cell** | **P7** | CosMx multiomic NSCLC + breast | 980 / 108 | same cell | > 800k cells, several patients | **Main Stage 2 dataset** (most proteins, most cells) |
| 1a | **P6** | CosMx Human Breast Multiomic | WTX (> 18,000) / 64 IO | same slide | 1 section | whole-transcriptome teacher (cognate vs programme ablations) |
| 1a | **P8** | CosMx Same-Cell Multiomics | > 19,000 / 72 | same section | see Bruker page | replication of P6 |
| 1a | **P1/P2** | Xenium Protein RCC / ccRCC | 377+ / 27 | same cell | ~465k cells | same-cell, cognate genes added on purpose (A1/A2-style teacher ablations) |
| **1b same section** | **P9** | Xenium + COMET + H&E (Zenodo 21713660, 2026) | Xenium panel / cyclic IF | same section (re-stained) | 4 cores (HCC, tonsil, adenoma), 3.9 GB, CC BY 4.0 | replication on a different protein chemistry; *descriptive* |
| 1b | **P10** | oral cavity atlas: Xenium 300-plex → PhenoCycler-Fusion 40-plex | 300 / 40 | same section | atlas scale | replication with a 40-plex CODEX-like panel |
| 1b | **P11** | TACIT: Xenium 280-plex + PhenoCycler 36-plex | 280 / 36 | same slide | see paper | replication |
| 1b | **P3** | DBiT-plus (*Nat Methods* 2026) | whole transcriptome (spots) / ~26–50-plex CODEX | same section | lymph node, lymphoma | spot-level teacher |
| **2 adjacent** | **P4** | SPATCH Xenium 5K / CosMx 6K → adjacent CODEX | 5,000+ / 16 | adjacent (registration) | 3 tumours | **"space is the bridge"**: teacher gets neighbourhood RNA after registration |
| **3 deployment** | D8, D2, D3 | protein-only cohorts | — / 35–56 | none | 35–700 patients | Does the RNA-taught student also do better on cohorts that never had RNA? (marker intersection only) |

Notes:
- **Label circularity**: if cell types were called from RNA clustering, an RNA teacher "agrees with the labels" by construction. Record `label_source` and run Stage 2 with protein-gated labels where available (P7 author annotations, SPATCH manual labels).
- Xenium Protein's RCC release was re-analysed with Xenium Onboard Analysis v4.0 and replaced on Aug 27 2025: download the current version.

---

## Download and harmonise (Great Lakes)

```bash
export RAW=/scratch/$SLURM_ACCOUNT/$USER/raw HARM=/scratch/$SLURM_ACCOUNT/$USER/harmonised
bash download_scripts/download_zenodo_mendeley.sh $RAW      # Zenodo / Mendeley sets (P9, D2, D3, D8)
# 10x (P1, P2, Xenium 5K) and Bruker (P6, P7, P8): fill the vendor form, then see download_scripts/*.md
sbatch ../cluster/slurm/harmonise.sbatch                    # all datasets in registry.yaml -> $HARM/*.h5ad
python ../shared/data/register_adjacent.py --rna $HARM/spatch_colon_xenium5k.h5ad \
       --protein $HARM/spatch_colon_codex.h5ad --out $HARM/spatch_colon_registration.json
```

Earlier per-dataset notes (input formats, what can go wrong) are in [`dataset_playbook.md`](dataset_playbook.md); the full v2 literature context is in [`../docs/background/literature_review/`](../docs/background/literature_review/).
