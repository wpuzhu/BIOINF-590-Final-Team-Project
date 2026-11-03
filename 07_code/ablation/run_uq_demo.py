"""
Uncertainty quantification + transferability diagnostics on simulated paired data.

  U1  cluster bootstrap CIs and sign-flip permutation p-values for ablation contrasts (unit = sample)
  U2  conditional group-permutation importance (within sample x cell type), with bootstrap CIs
  U3  toroidal-shift test vs naive cell permutation for spatial associations of residuals
  U4  Mondrian (per-cell-type) split conformal: coverage in-distribution vs on a held-out region
  U5  weighted conformal + shift decomposition: covariate shift vs concept shift

Usage: python run_uq_demo.py [--out DIR]
"""
import argparse
import os

import numpy as np
import pandas as pd

from run_simulation_demo import simulate
from uq import (bh, cluster_bootstrap_ci, fold_scores, group_permutation_importance, holm,
                mondrian_conformal, shift_decomposition, signflip_test, toroidal_shift_test,
                unit_folds, weighted_conformal_coverage)
from xmodal_ablation import Features, radius_graph, shuffle_within


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="../../08_results/simulation/uq")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)

    d = simulate(n_samples=8, n_cells=900, seed=7)
    F = Features(d, n_pcs=30)
    folds = unit_folds(d.sample)
    rng = np.random.default_rng(0)
    perm_ct = shuffle_within(d.sample, d.cell_type, rng)
    perm_any = shuffle_within(d.sample, np.zeros(len(d.sample)), rng)
    conds = {
        "celltype": lambda p: F.build(),
        "cognate": lambda p: F.build(cognate_for=p),
        "own_rna": lambda p: F.build(cognate_for=p, all_rna=True),
        "nbr_real": lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=30),
        "nbr_shuf_within_type": lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=30, nbr_perm=perm_ct),
        "nbr_shuf_any": lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=30, nbr_perm=perm_any),
        "rna_shuffled": lambda p: F.build(cognate_for=p, all_rna=True, nbr_radius=30, rna_perm=perm_any),
    }
    fs = pd.concat([fold_scores(d, f, folds, lab) for lab, f in conds.items()])
    fs.to_csv(f"{a.out}/U0_fold_scores.csv", index=False)
    wide = fs.pivot_table(index=["protein", "unit"], columns="condition", values="r2")

    # ---------------- U1 contrasts
    contrasts = {
        "cognate_gain (cognate - celltype)": ("cognate", "celltype", "greater"),
        "neighbourhood_gain (nbr_real - own_rna)": ("nbr_real", "own_rna", "greater"),
        "arrangement (nbr_real - nbr_shuf_within_type)": ("nbr_real", "nbr_shuf_within_type", "greater"),
        "composition (nbr_shuf_within_type - nbr_shuf_any)": ("nbr_shuf_within_type", "nbr_shuf_any", "greater"),
        "rna_content (nbr_real - rna_shuffled)": ("nbr_real", "rna_shuffled", "greater"),
    }
    rows = []
    for prot in d.proteins:
        w = wide.loc[prot]
        for name, (c1, c0, alt) in contrasts.items():
            diffs = (w[c1] - w[c0]).values
            est, lo, hi = cluster_bootstrap_ci(diffs)
            rows.append({"protein": prot, "contrast": name, "mean_dR2": est, "ci_lo": lo, "ci_hi": hi,
                         "p_signflip": signflip_test(diffs, alternative=alt), "n_units": len(diffs)})
    u1 = pd.DataFrame(rows)
    u1["p_holm"] = holm(u1.p_signflip.values)
    u1["q_bh"] = bh(u1.p_signflip.values)
    u1.to_csv(f"{a.out}/U1_contrasts_bootstrap_signflip.csv", index=False)
    print("\nU1 contrasts (unit = sample, 8 units; exact sign-flip has min p = 1/256)")
    print(u1.round(3).to_string(index=False))

    # ---------------- U2 conditional group permutation importance
    groups = {"MHCII_programme": ["CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"], "IFNG": ["IFNG"],
              "CD3E": ["CD3E"], "MKI67": ["MKI67"],
              "random_4": list(rng.choice([g for g in d.genes if g.startswith("G")], 4, replace=False))}
    imp = pd.concat([group_permutation_importance(d, F, p, groups, folds, n_rep=5) for p in ["HLADR", "PDL1", "Ki67"]])
    per_unit = imp.groupby(["protein", "group", "unit"]).importance.mean().reset_index()
    u2 = []
    for (p, g), sub in per_unit.groupby(["protein", "group"]):
        est, lo, hi = cluster_bootstrap_ci(sub.importance.values)
        u2.append({"protein": p, "group": g, "importance_dR2": est, "ci_lo": lo, "ci_hi": hi,
                   "p_signflip": signflip_test(sub.importance.values)})
    u2 = pd.DataFrame(u2)
    u2.to_csv(f"{a.out}/U2_group_permutation_importance.csv", index=False)
    print("\nU2 conditional group permutation importance")
    print(u2.round(3).to_string(index=False))

    # ---------------- U3 toroidal shift vs naive permutation
    X = conds["nbr_real"]
    resid = {}
    for p in ["PDL1", "pS6", "Ki67"]:
        j = d.proteins.index(p)
        pred = np.zeros(len(d.sample))
        for name, tr, te in folds:
            from uq import _fit_predict
            pred[te], _, _ = _fit_predict(X(p), d.protein[:, j], tr, te)
        resid[p] = d.protein[:, j] - pred
    W = radius_graph(d.coords, d.sample, 60)
    t_density = W @ (d.cell_type == "Tcell").astype(float)
    y_coord = d.coords[:, 1]
    tumour = d.cell_type == "Tumor"
    u3 = []
    for p in resid:
        for cov_name, cov in [("local_Tcell_density_60um", t_density), ("y_coordinate", np.sin(y_coord / 160))]:
            r = toroidal_shift_test(resid[p][tumour], cov[tumour], d.coords[tumour], d.sample[tumour], n_perm=299)
            u3.append({"residual_of": p, "covariate": cov_name, **r})
    u3 = pd.DataFrame(u3)
    u3.to_csv(f"{a.out}/U3_toroidal_vs_naive.csv", index=False)
    print("\nU3 residual associations in tumour cells: toroidal-shift vs naive cell permutation")
    print(u3.round(4).to_string(index=False))

    # ---------------- U4 Mondrian conformal: in-distribution vs held-out region
    samples = np.unique(d.sample)
    tr_s, cal_s, te_s = samples[:4], samples[4:6], samples[6:]
    in_src = np.isin(d.region, ["margin", "stroma"])
    train = np.isin(d.sample, tr_s) & in_src
    calib = np.isin(d.sample, cal_s) & in_src
    test_in = np.isin(d.sample, te_s) & in_src
    test_core = np.isin(d.sample, te_s) & (d.region == "core")
    u4 = []
    for p in d.proteins:
        Xp = X(p)
        for lab, te in [("test_same_regions", test_in), ("test_core_region", test_core)]:
            r = mondrian_conformal(d, Xp, p, train, calib, te, alpha=0.1)
            r["test_set"] = lab
            u4.append(r)
    u4 = pd.concat(u4)
    u4.to_csv(f"{a.out}/U4_mondrian_conformal.csv", index=False)
    print("\nU4 per-cell-type 90% conformal coverage (target 0.90)")
    print(u4.pivot_table(index=["protein", "cell_type"], columns="test_set", values="coverage_test").round(3).to_string())

    # ---------------- U5 weighted conformal + shift decomposition (source: margin+stroma, target: core)
    # shift is measured on low-dimensional CONTEXT features: own cell type + neighbourhood composition.
    # Density ratios on 100+ dimensional RNA features are unstable (see docs).
    Zctx = np.hstack([F.ct, radius_graph(d.coords, d.sample, 30) @ F.ct])
    u5a, u5b = [], []
    for p in d.proteins:
        Xp = X(p)
        u5a.append(weighted_conformal_coverage(d, Xp, p, train, calib, test_core, alpha=0.1, shift_features=Zctx))
        u5b.append(shift_decomposition(d, Xp, p, source=in_src, target=d.region == "core", shift_features=Zctx))
    u5a, u5b = pd.DataFrame(u5a), pd.DataFrame(u5b)
    u5a.to_csv(f"{a.out}/U5a_weighted_conformal.csv", index=False)
    u5b.to_csv(f"{a.out}/U5b_shift_decomposition.csv", index=False)
    print("\nU5a weighted vs unweighted conformal coverage on the held-out core region")
    print(u5a.round(3).to_string(index=False))
    print("\nU5b shift decomposition (MSE units; protein z-scored overall)")
    print(u5b.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
