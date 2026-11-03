"""Figure for the UQ self-check (reads the CSVs written by run_uq_demo.py)."""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d = sys.argv[1] if len(sys.argv) > 1 else "../../08_results/simulation/uq"
u1 = pd.read_csv(f"{d}/U1_contrasts_bootstrap_signflip.csv")
u3 = pd.read_csv(f"{d}/U3_toroidal_vs_naive.csv")
u4 = pd.read_csv(f"{d}/U4_mondrian_conformal.csv")
u5 = pd.read_csv(f"{d}/U5a_weighted_conformal.csv")

fig, ax = plt.subplots(2, 2, figsize=(14, 10))
# A forest plot
u1 = u1.copy()
u1["label"] = u1.protein + " | " + u1.contrast.str.split(" ").str[0]
u1 = u1[u1.contrast.str.startswith(("cognate", "neighbourhood", "arrangement", "composition"))]
y = np.arange(len(u1))
sig = u1.q_bh < 0.05
ax[0, 0].errorbar(u1.mean_dR2, y, xerr=[u1.mean_dR2 - u1.ci_lo, u1.ci_hi - u1.mean_dR2], fmt="none", ecolor="grey")
ax[0, 0].scatter(u1.mean_dR2, y, c=np.where(sig, "C3", "C0"), zorder=3)
ax[0, 0].axvline(0, color="k", lw=0.8)
ax[0, 0].set_yticks(y); ax[0, 0].set_yticklabels(u1.label, fontsize=6)
ax[0, 0].set_xlabel("mean ΔR² over held-out samples (95% cluster-bootstrap CI)")
ax[0, 0].set_title("A. Ablation contrasts (red: BH q < 0.05, sign-flip test)", fontsize=9)
# B toroidal vs naive
lab = u3.residual_of + " ~ " + u3.covariate.str.replace("_", " ")
x = np.arange(len(u3))
ax[0, 1].bar(x - 0.2, u3.null_sd_naive, 0.4, label="naive cell permutation")
ax[0, 1].bar(x + 0.2, u3.null_sd_toroidal, 0.4, label="toroidal shift")
for i, r in u3.iterrows():
    ax[0, 1].text(i, max(r.null_sd_naive, r.null_sd_toroidal) + 0.005,
                  f"p={r.p_naive_cell_perm:.3f}\np={r.p_toroidal:.3f}", ha="center", fontsize=6)
ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(lab, rotation=30, ha="right", fontsize=7)
ax[0, 1].set_ylabel("SD of null correlation"); ax[0, 1].legend(fontsize=7)
ax[0, 1].set_title("B. Spatial nulls: naive permutation is too narrow (p: naive / toroidal)", fontsize=9)
# C conformal coverage
pv = u4.pivot_table(index=["protein", "cell_type"], columns="test_set", values="coverage_test")
ax[1, 0].scatter(pv["test_same_regions"], pv["test_core_region"])
for (p, c), r in pv.iterrows():
    if r["test_core_region"] < 0.85:
        ax[1, 0].annotate(f"{p}/{c}", (r["test_same_regions"], r["test_core_region"]), fontsize=7)
ax[1, 0].axhline(0.9, ls="--", c="grey"); ax[1, 0].axvline(0.9, ls="--", c="grey")
ax[1, 0].plot([0.7, 1], [0.7, 1], c="k", lw=0.5)
ax[1, 0].set_xlabel("coverage, held-out samples, same regions"); ax[1, 0].set_ylabel("coverage, held-out core region")
ax[1, 0].set_title("C. Per-cell-type 90% conformal coverage: in-distribution vs transferred", fontsize=9)
# D weighted conformal
x = np.arange(len(u5))
ax[1, 1].bar(x - 0.2, u5.coverage_unweighted, 0.4, label="split conformal")
ax[1, 1].bar(x + 0.2, u5.coverage_weighted, 0.4, label="covariate-shift weighted")
ax[1, 1].axhline(0.9, ls="--", c="grey")
ax[1, 1].set_xticks(x); ax[1, 1].set_xticklabels(u5.protein)
ax[1, 1].set_ylim(0.5, 1.02); ax[1, 1].legend(fontsize=7)
ax[1, 1].set_title("D. Core-region coverage: weighting fixes covariate shift (PDL1), not all", fontsize=9)
fig.tight_layout()
fig.savefig(f"{d}/uq_overview.png", dpi=150)
print("saved", f"{d}/uq_overview.png")
