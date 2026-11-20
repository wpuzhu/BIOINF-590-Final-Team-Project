# E8 · Pairing degradation (registration jitter)

**Question.** How fast does RNA → protein information decay with pairing error, and can adjacent-section data (SPATCH) still carry it?

**Code.** `pairing_degradation(d, F, deltas, same_type)` and `jitter_index(...)` in `xmodal_ablation.py`; mirrors H2 arm H

## Manipulation (input level only)
On same-cell data, replace each cell's RNA with the RNA of a random cell within distance δ ∈ {0, 10, 20, 50, 100} µm (optionally restricted to the same type), and rebuild all features.

## Controls
δ = ∞ (random cell in sample) equals E4e.

## Read-out
R²(δ) per protein; the δ at which half the RNA-explained variance is lost (δ½).

## Uncertainty quantification
Bootstrap CI on δ½ over units; repeat the jitter draws.

## Interpretation
Compare δ½ with the measured SPATCH registration error: proteins with δ½ > error can be informed by adjacent sections; others cannot (H2 stage-2 feasibility).

## Failure modes
Same-type-restricted jitter preserves composition, so state proteins decay faster than lineage proteins (expected, and informative).

## Simulation self-check
With any-cell jitter, Ki-67 falls from 0.79 (0 µm) to 0.14 (20 µm) and HLA-DR from 0.50 to 0.12, so cell-state information is gone at one cell diameter. Composition-driven PD-L1 keeps 0.40 → 0.38 under same-type jitter up to 100 µm. Lineage CD3 drops to its cell-type floor (0.54). So adjacent-section data can only carry composition information unless registration is better than ~10 µm.
