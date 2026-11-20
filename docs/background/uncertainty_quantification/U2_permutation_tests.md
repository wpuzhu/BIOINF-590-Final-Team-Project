# U2 · Permutation tests: choose the null that matches H0

| H0 | Permute what | Keeps | Code |
|---|---|---|---|
| "adding input S does not improve held-out R²" (paired, per unit) | **sign** of each unit's difference | unit structure | `signflip_test(diffs)` |
| "S carries no information beyond cell type" | S among cells **of the same type within sample** | cell-type signal, composition | `shuffle_within(sample, cell_type)`; `group_permutation_importance` |
| "neighbour *states* don't matter, only composition" | neighbourhood RNA among same-type cells | composition | E4c |
| "two spatial fields are unrelated" (residual vs hypoxia, vs immune density) | **toroidal shift** of one field (or block rotation) | each field's own autocorrelation | `toroidal_shift_test(x, y, coords, sample)` |
| "this gene group is no more important than a random group of the same size" | gene-set labels (random draws) | group size, expression level | E2 random sets |
| "the residual has no spatial structure within type" | values among same-type cells | marginal distribution | within-type Moran's I |
| "cell-level predictions are uncorrelated with truth" (probe null, N7) | target vector across units | – | permute targets at the unit level |

## Why toroidal shifts?
Two smooth spatial fields are often correlated **by chance**: both vary slowly, so there are few effective independent observations. A naive cell-level shuffle destroys that smoothness, so the null is far too narrow and p-values come out too small. Shifting one field on a torus (wrapping at the slide's bounding box) keeps its autocorrelation and breaks only the alignment.

In our simulation the null SD of the correlation was 0.017 (naive) vs 0.22 (toroidal) for the smooth pS6 residual vs a smooth covariate, a 13× difference. For PD-L1 residual vs the y-coordinate, the naive p was 0.055 vs 0.18 toroidal.

**Caveats.** Edges: the wrap joins opposite slide edges; for irregular tissue, prefer **random rotations + translations within the tissue mask**, or Moran spectral randomisation. With several slides, shift each slide independently (as implemented).

## Exactness and resolution
With n units, the sign-flip test has 2ⁿ possible sign patterns. Our code enumerates all of them when 2ⁿ ≤ 10,000, so the test is exact. With n = 8 the minimum p is 1/256 ≈ 0.0039.
