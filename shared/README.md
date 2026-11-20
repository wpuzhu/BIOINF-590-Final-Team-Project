# shared/: code used by both stages

| Folder | Contents |
|---|---|
| `data/` | harmonised schema (`schema.md`), readers/writers (`spt_io.py`: .h5ad or .npz), raw → schema loaders (`harmonise.py`: Xenium incl. Xenium Protein, CosMx flat files, generic cell tables), adjacent-section registration (`register_adjacent.py`), simulated export for smoke tests |
| `sim/` | `paired_spatial_sim.py`: paired spatial RNA + protein simulator with four scenarios (null, informative-only, transferable via protein, transferable via space), plus generic array operations used on real data too: neighbour means, pairing jitter, graded nulls |
| `models/` | `np_mlp.py`: numpy MLP with hard labels + distillation + auxiliary target in one loss (simulation and CPU smoke tests) |
| `uq/` | `uq_stats.py`: paired-improvement UQ (hierarchical bootstrap, sign-flip, TOST, variance components, conformal APS, ECE, Holm/BH, `improvement_report`) · `uq_core.py`: information-map UQ (cluster bootstrap, toroidal shift, Mondrian/weighted conformal, shift decomposition) · `aggregate.py`: cluster results → tables |
| `info_map/` | the v2 information map (which RNA information predicts each protein; E1–E8) with its simulation checks: background for both stages |

Every comparison is paired (same patients, labelled cells and seeds) and resampled at the **patient** level. See [`../docs/concept_and_design.md`](../docs/concept_and_design.md#uncertainty-quantification-per-experiment).
