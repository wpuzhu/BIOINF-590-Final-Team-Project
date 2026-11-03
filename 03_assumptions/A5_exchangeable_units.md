# A5. Exchangeable, independent units for validation and inference

**Statement.** The units we split on and resample (patients, slides, spatial blocks) are approximately independent and exchangeable. Cells are not: neighbouring cells share environment, segmentation errors and imaging artefacts.

**If violated.**
- Cell-level train/test splits leak neighbourhood information, which inflates R², especially for spatial features.
- Cell-level bootstrap or permutation gives confidence intervals and p-values that are far too narrow. Our simulation shows the null SD of a correlation is up to **13× too small** under naive cell permutation (`08_results/simulation/uq/`).

**Rules.**
1. Split by **patient/slide**; with one slide, split by **spatial blocks** at least 2–3× the largest neighbourhood radius (E3 goes to 100 µm, so blocks ≥ 250–300 µm).
2. Bootstrap **units**, not cells (cluster bootstrap).
3. Permutation nulls must keep what H0 keeps: within-cell-type shuffles for "no within-type information"; toroidal shifts or block permutations for "no spatial association".
4. With few units (e.g. 8 patients), an exact sign-flip test has a minimum p of 1/256. Plan the number of tests (and Holm vs BH correction) around that.
