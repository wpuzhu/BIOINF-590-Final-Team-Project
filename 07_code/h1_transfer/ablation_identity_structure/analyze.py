"""
Decompose the transfer gain into identity / structure / generic parts, with uncertainty.

Paired contrasts (same seed, same test patient, same labelled cells):
  total      = pretrained    - scratch           (frozen mode: - random_features)
  identity   = pretrained    - A1_reinit         gene-token meaning carried to protein
  mapping    = pretrained    - A2_random_map     correct vs WRONG gene->token map
  structure  = A1_reinit     - noise_pre         encoder knowledge from real co-expression
  generic    = noise_pre     - scratch           pretraining with no biology (optimisation/regularisation)
  vs_PCA     = pretrained    - PCA_probe         (frozen only) the bar foundation models must clear
  total = identity + structure + generic  (exactly, by construction)

Uncertainty: units are (seed x held-out patient) pairs (5 x 3 = 15). Percentile bootstrap
(10,000 resamples) of the mean paired difference; sign-flip test over seeds (patient-averaged,
exact with 5 seeds -> min one-sided p = 1/32). Holm across the contrasts of each scenario x mode.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "ablation"))
from uq import holm, signflip_test  # noqa: E402

D = sys.argv[1] if len(sys.argv) > 1 else "../../../08_results/simulation/h1_ablation"
raw = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{D}/raw_*.csv"))])
rng = np.random.default_rng(0)

CONTRASTS = {
    "finetune": {"total": ("pretrained", "scratch"), "identity": ("pretrained", "A1_reinit"),
                 "mapping": ("pretrained", "A2_random_map"), "structure": ("A1_reinit", "noise_pre"),
                 "generic": ("noise_pre", "scratch")},
    "frozen": {"total": ("pretrained", "random_features"), "identity": ("pretrained", "A1_reinit"),
               "mapping": ("pretrained", "A2_random_map"), "structure": ("A1_reinit", "noise_pre"),
               "generic": ("noise_pre", "random_features"), "vs_PCA": ("pretrained", "PCA_probe")},
    "finetune_RNA_target": {"total": ("pretrained", "scratch"), "identity": ("pretrained", "A1_reinit"),
                            "mapping": ("pretrained", "A2_random_map"), "structure": ("A1_reinit", "noise_pre"),
                            "generic": ("noise_pre", "scratch")},
    "finetune_RNA_panel": {"total": ("pretrained", "scratch")},
}


def boot_ci(v, B=10000):
    idx = rng.integers(0, len(v), (B, len(v)))
    m = v[idx].mean(1)
    return float(v.mean()), float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


rows = []
for (sc, mode, k), g in raw.groupby(["scenario", "mode", "labels"]):
    wide = g.pivot_table(index=["seed", "test_patient"], columns="arm", values="macro_f1")
    for name, (a1, a0) in CONTRASTS.get(mode, {}).items():
        if a1 not in wide or a0 not in wide:
            continue
        d = (wide[a1] - wide[a0]).dropna()
        est, lo, hi = boot_ci(d.values)
        per_seed = d.groupby(level="seed").mean().values
        rows.append(dict(scenario=sc, mode=mode, labels=k, contrast=name, mean=est, ci_lo=lo, ci_hi=hi,
                         p_two_sided=signflip_test(per_seed, alternative="two-sided"), n_units=len(d)))
res = pd.DataFrame(rows)
res["p_holm"] = np.nan
for _, idx in res.groupby(["scenario", "mode", "labels"]).groups.items():
    res.loc[idx, "p_holm"] = holm(res.loc[idx, "p_two_sided"].values)
res.to_csv(f"{D}/decomposition.csv", index=False)

# mean F1 table
f1 = raw.groupby(["scenario", "mode", "arm", "labels"]).macro_f1.mean().unstack().round(3)
f1.to_csv(f"{D}/mean_macro_f1.csv")

# A3 overlap curves (gain over scratch as mapped fraction shrinks, information held constant)
a3 = []
for (sc, k), g in raw[raw["mode"] == "finetune"].groupby(["scenario", "labels"]):
    wide = g.pivot_table(index=["seed", "test_patient"], columns="arm", values="macro_f1")
    for q, arm in [(1.0, "pretrained"), (0.75, "A3_overlap_0.75"), (0.5, "A3_overlap_0.50"), (0.25, "A3_overlap_0.25")]:
        est, lo, hi = boot_ci((wide[arm] - wide["scratch"]).values)
        a3.append(dict(scenario=sc, labels=k, overlap=q, gain=est, ci_lo=lo, ci_hi=hi))
a3 = pd.DataFrame(a3)
a3.to_csv(f"{D}/A3_overlap.csv", index=False)

# ---------------------------------------------------------------- figure
scen = ["coupled", "realistic", "decoupled", "coupled_P60"]
parts = ["identity", "structure", "generic"]
colors = {"identity": "#3b6ea8", "structure": "#d08c3c", "generic": "#9aa0a6", "total": "black"}
fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))


def deco_panel(ax, mode, k, title, scen_list):
    sub = res[(res["mode"] == mode) & (res["labels"] == k)]
    x = np.arange(len(scen_list))
    for i, part in enumerate(parts + ["total"]):
        vals = [sub[(sub.scenario == s) & (sub.contrast == part)] for s in scen_list]
        m = [v["mean"].iloc[0] if len(v) else np.nan for v in vals]
        lo = [v["ci_lo"].iloc[0] if len(v) else np.nan for v in vals]
        hi = [v["ci_hi"].iloc[0] if len(v) else np.nan for v in vals]
        off = (i - 1.5) * 0.19
        ax.bar(x + off, m, 0.18, color=colors[part], label=part, alpha=0.9 if part != "total" else 0.35,
               edgecolor="black" if part == "total" else None)
        ax.errorbar(x + off, m, yerr=[np.subtract(m, lo), np.subtract(hi, m)], fmt="none", ecolor="k", lw=0.8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(scen_list, fontsize=8)
    ax.set_ylabel("Δ macro-F1 (95% bootstrap CI)"); ax.set_title(title, fontsize=9)


deco_panel(axes[0, 0], "finetune", 5, "A. Protein target, FINE-TUNED, 5 labels/class", scen)
deco_panel(axes[0, 1], "frozen", 5, "B. Protein target, FROZEN + probe, 5 labels/class", scen)
deco_panel(axes[0, 2], "finetune_RNA_target", 5, "C. Positive control: RNA target (in-modality), 5 labels", ["coupled"])
axes[0, 0].legend(fontsize=7, loc="lower left")

for sc, c in zip(scen, ["C0", "C1", "C2", "C3"]):
    s = a3[(a3.scenario == sc) & (a3.labels == 5)].sort_values("overlap")
    axes[1, 0].errorbar(s.overlap, s.gain, yerr=[s.gain - s.ci_lo, s.ci_hi - s.gain], marker="o", color=c, label=sc, capsize=2)
axes[1, 0].axhline(0, color="k", lw=0.8)
axes[1, 0].set_xlabel("fraction of markers mapped to their gene token (rest: new tokens)")
axes[1, 0].set_ylabel("fine-tuned gain over scratch"); axes[1, 0].legend(fontsize=7)
axes[1, 0].set_title("D. A3 panel overlap (5 labels): no dose-response", fontsize=9)

fz = raw[(raw["mode"] == "frozen") & (raw.labels == 5)].groupby(["scenario", "arm"]).macro_f1.mean().unstack()
order = ["PCA_probe", "pretrained", "A2_random_map", "A1_reinit", "noise_pre", "random_features"]
fz[order].loc[scen].plot.bar(ax=axes[1, 1], width=0.8)
axes[1, 1].set_ylabel("macro-F1"); axes[1, 1].set_title("E. Frozen embeddings vs PCA probe (5 labels)", fontsize=9)
axes[1, 1].tick_params(axis="x", rotation=0, labelsize=8); axes[1, 1].legend(fontsize=6)

ft = raw[(raw["mode"].isin(["finetune_RNA_target", "finetune_RNA_panel"]))]
for (mode, arm), g in ft.groupby(["mode", "arm"]):
    if arm not in ("scratch", "pretrained"):
        continue
    m = g.groupby("labels").macro_f1.mean()
    axes[1, 2].plot(m.index, m.values, marker="o", ls="-" if arm == "pretrained" else "--",
                    color="C0" if mode == "finetune_RNA_target" else "C3",
                    label=f"{'300-gene RNA' if mode == 'finetune_RNA_target' else '16 cognate genes only'} · {arm}")
axes[1, 2].set_xscale("log"); axes[1, 2].set_xlabel("labels per class"); axes[1, 2].set_ylabel("macro-F1")
axes[1, 2].legend(fontsize=7); axes[1, 2].set_title("F. Same modality, panel cut to protein size: gain vanishes", fontsize=9)
fig.suptitle("H1 ablation testbed: where does RNA-pretraining knowledge go when the target is protein?", fontsize=12)
fig.tight_layout()
fig.savefig(f"{D}/h1_ablation_overview.png", dpi=150)

pd.set_option("display.width", 220)
print(res[res.labels == 5].round(3).to_string(index=False))
print()
print(a3[a3.labels == 5].round(3).to_string(index=False))
