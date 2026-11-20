# E3 · Neighbourhood radius sweep

**Question.** At what spatial scale does neighbouring RNA inform a protein?

**Code.** `radius_sweep(d, F, radii=(0,10,20,30,50,100))`

## Manipulation (input level only)
Neighbourhood RNA features computed at radius r ∈ {0, 10, 20, 30, 50, 100, 200} µm (or k-NN with k ∈ {5, 10, 20, 50}).

## Controls
Same sweep with neighbours shuffled within cell type (E4c) to separate state from composition at each scale.

## Read-out
R²(r) curve per protein × cell type; r* = argmax.

## Uncertainty quantification
Bootstrap the **argmax** over units (CI on r*); report the curve with pointwise bootstrap bands. Blocks for splitting must be ≥ 2–3 × max r (A5).

## Interpretation
r* ≈ 10–20 µm: juxtacrine/contact or spillover (check A7). 30–100 µm: paracrine. Monotone increase up to the largest r: the model is learning **region identity**, so expect poor E5 transfer.

## Failure modes
Large radii leak region identity; spot data already integrate ~25–55 µm.

## Simulation self-check
PD-L1 R² rises from 0.17 (0 µm) to 0.49 (50 µm), matching the planted induction kernel (σ = 30 µm).
