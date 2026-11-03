# A1. Cell correspondence: the RNA and protein we pair really belong together

**Statement.** For each analysis unit (cell, spot or neighbourhood), the RNA vector and the protein vector describe the same biological material.

| Data tier | How well A1 holds | Main violation |
|---|---|---|
| Same cell (Xenium Protein, CosMx multiomic) | strong | segmentation errors and lateral spillover mix neighbours (see A7) |
| Same section, spot RNA (DBiTplus, spatial-CITE-seq, Visium protein) | at spot level | a spot mixes several cells, so "own RNA" partly contains neighbour RNA |
| Adjacent sections (SPATCH) | only at neighbourhood scale | 5–10 µm section offset plus registration error of tens of µm |
| Unpaired (scRNA-seq ↔ CODEX via MaxFuse) | statistical only | matches are probabilistic; errors concentrate in ambiguous cell states |

**If violated.** Information about **cell states** is destroyed first (it varies at the scale of one cell). **Composition** information survives longest (it varies at the scale of niches). A violation therefore looks like "state information does not transfer", even when it would with perfect pairing.

**How we test it.**
- **E8:** jitter same-cell pairs by 0–100 µm and record how fast each protein's R² decays. This gives a decay curve for each protein.
- Compare with the measured registration error on SPATCH to predict which proteins can still be informed there.
- Pairing-method sensitivity: registration-based vs MaxFuse-based pairs. Disagreement = uncertainty.
