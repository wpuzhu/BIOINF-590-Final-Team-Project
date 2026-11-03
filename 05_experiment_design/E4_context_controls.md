# E4 · Context controls: arrangement vs composition vs nothing

**Question.** When the neighbourhood helps, is it because of *which cell types* are around (composition), *what state* they are in (arrangement), or is it a leak?

**Code.** `context_controls(d, F, radius)`

## Manipulation (input level only)
(a) own RNA only; (b) + real neighbourhood; (c) neighbourhood RNA taken from random cells **of the same type in the same sample** (keeps composition, destroys state arrangement); (d) neighbourhood from random cells of any type (destroys composition); (e) **all** RNA shuffled across cells (negative control).

## Controls
(e) is the control: everything above the cell-type baseline must vanish.

## Read-out
R² for (a)–(e); contrasts: arrangement = b − c; composition = c − d; neighbourhood = b − a; RNA content = b − e.

## Uncertainty quantification
Per-unit contrasts → bootstrap CI + sign-flip p; repeat shuffles ≥5 times and average within unit before testing.

## Interpretation
b > c: the **state** of neighbours matters (a mechanistic, cell-state signal). b ≈ c > d: only **composition** matters, which transfers wherever composition is similar but is not a state signal. e < baseline: something is wrong (leak or bug).

## Failure modes
Spillover (A7) mimics arrangement for lineage markers of neighbour types.

## Simulation self-check
PD-L1: composition effect 0.11 (CI 0.10–0.13, p = 0.004), arrangement ~0, as planted (induction depends on T-cell presence).
