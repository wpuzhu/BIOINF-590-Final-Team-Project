"""
Summarise the Stage 2 simulation with uncertainty.

Outputs (in results/simulation/):
  main_improvements.csv     every pre-registered contrast x scenario: estimate, 95% hierarchical
                            bootstrap CI (patients -> seeds), sign-flip p, Holm p, TOST, variance shares
  transfer_efficiency.csv   Delta_trans / Delta_info with a patient-bootstrap CI
  information_ceiling.csv   NLL gap protein-only vs teacher, and the fraction a student realises
  conformal.csv             90% APS coverage and set size per arm (cross-conformal by patient)
  ablation_deltas.csv       Delta_info (teacher - none) and Delta_trans (kd - none) per ablation setting
  stage2_overview.png       four-panel figure
"""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "shared", "uq"))
from uq_stats import hierarchical_bootstrap, improvement_report, paired_table  # noqa: E402

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "stage2_privileged_transfer", "results", "simulation")
SC = ["null", "informative_only", "transferable_protein", "transferable_spatial"]
SC_LABEL = {"null": "Null", "informative_only": "Informative only",
            "transferable_protein": "Transferable (protein)", "transferable_spatial": "Transferable (spatial)"}
COMPARISONS = {
    "informative: teacher - none": ("teacher", "none"),
    "transferable: kd - none": ("kd", "none"),
    "transferable: kd_feat - none": ("kd_feat", "none"),
    "transferable: align_shared - none": ("align_shared", "none"),
    "transferable: recon_rna - none": ("recon_rna", "none"),
    "transferable: recon_specific - none": ("recon_specific", "none"),
    "control: self_train - none (unlabelled pool)": ("self_train", "none"),
    "control: label_smooth - none (regularisation)": ("label_smooth", "none"),
    "RNA-specific: kd - self_train": ("kd", "self_train"),
    "dark knowledge: kd - kd_permuted": ("kd", "kd_permuted"),
}
MARGIN = 0.02
COL = {"null": "#8a8a85", "informative_only": "#eb6834", "transferable_protein": "#2a78d6",
       "transferable_spatial": "#1baf7a"}


def load(kind):
    fs = sorted(glob.glob(os.path.join(D, f"raw_{kind}_*.csv")))
    return pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs]) if fs else pd.DataFrame()


def ratio_ci(d_num, d_den, B=5000, seed=0):
    a = d_num.groupby("patient")["diff"].mean()
    b = d_den.groupby("patient")["diff"].mean().reindex(a.index)
    rng = np.random.default_rng(seed)
    est = a.mean() / b.mean()
    bs = []
    for _ in range(B):
        i = rng.integers(0, len(a), len(a))
        den = b.values[i].mean()
        bs.append(a.values[i].mean() / den if den > 1e-6 else np.nan)
    bs = np.array(bs)
    lo, hi = np.nanquantile(bs, [0.025, 0.975])
    return est, lo, hi


def main():
    m = load("main")
    out = []
    rep = improvement_report(m, "arm", "macro_f1", COMPARISONS, margin=MARGIN, alternative="two-sided",
                             strata=("scenario",))
    rep.to_csv(os.path.join(D, "main_improvements.csv"), index=False)
    print(rep[["scenario", "comparison", "estimate", "ci_lo", "ci_hi", "p_signflip", "p_holm",
               "tost_equivalent", "share_patient", "share_seed"]].round(3).to_string(index=False))

    eff = []
    for sc in SC:
        g = m[m.scenario == sc]
        if g.empty:
            continue
        di = paired_table(g, "arm", "macro_f1", "teacher", "none")
        for arm in ("kd", "kd_feat", "align_shared", "recon_rna", "recon_specific"):
            dt = paired_table(g, "arm", "macro_f1", arm, "none")
            e, lo, hi = ratio_ci(dt, di)
            info = hierarchical_bootstrap(di)
            eff.append(dict(scenario=sc, arm=arm, efficiency=e, ci_lo=lo, ci_hi=hi,
                            delta_info=info["estimate"], info_ci_lo=info["ci_lo"],
                            defined=bool(info["ci_lo"] > 0)))
    eff = pd.DataFrame(eff)
    eff.to_csv(os.path.join(D, "transfer_efficiency.csv"), index=False)

    ce = (m.groupby(["scenario", "arm"]).nll.mean().unstack())
    ceil = pd.DataFrame({"nll_none": ce["none"], "nll_self_train": ce["self_train"], "nll_teacher": ce["teacher"],
                         "nll_kd": ce["kd"]})
    ceil["ceiling_nats"] = ceil.nll_none - ceil.nll_teacher
    ceil["realised_fraction_kd"] = (ceil.nll_none - ceil.nll_kd) / ceil.ceiling_nats
    # distillation also fixes over-confidence; measure against the self-trained student, which
    # gets the same soft-label calibration benefit without RNA
    ceil["realised_fraction_kd_vs_self_train"] = ((ceil.nll_self_train - ceil.nll_kd)
                                                  / (ceil.nll_self_train - ceil.nll_teacher))
    ceil.to_csv(os.path.join(D, "information_ceiling.csv"))
    conf = m.groupby(["scenario", "arm"])[["coverage", "set_size"]].mean().round(3)
    ss = improvement_report(m, "arm", "set_size", {"set size: kd - none": ("kd", "none"),
                                                   "set size: kd - self_train": ("kd", "self_train")},
                            margin=0.1, alternative="two-sided", strata=("scenario",))
    ss.to_csv(os.path.join(D, "conformal_set_size_contrasts.csv"), index=False)
    conf.to_csv(os.path.join(D, "conformal.csv"))
    print(ceil.round(3).to_string())
    print(conf.unstack(0).round(3).to_string())

    # ---------------- ablations
    a = load("ablations")
    base = m[m.arm.isin(["none", "teacher", "kd"]) & (m.seed < (a.seed.max() + 1 if len(a) else 3))]
    rows = []
    if len(a):
        a = pd.concat([a, base.assign(experiment="default", setting="default")])
        for (sc, ex, st), g in a.groupby(["scenario", "experiment", "setting"]):
            for name, (x1, x0) in {"delta_info": ("teacher", "none"), "delta_trans": ("kd", "none")}.items():
                r = hierarchical_bootstrap(paired_table(g, "arm", "macro_f1", x1, x0), B=4000)
                rows.append(dict(scenario=sc, experiment=ex, setting=st, quantity=name, **r,
                                 f1_none=g[g.arm == "none"].macro_f1.mean()))
    ab = pd.DataFrame(rows)
    if len(ab):
        ab.to_csv(os.path.join(D, "ablation_deltas.csv"), index=False)
    figure(rep, ab)


def _series(ab, ex, sc, q, default_x, cast=float):
    g = ab[(ab.scenario == sc) & (ab.quantity == q) & (ab.experiment.isin([ex, "default"]))].copy()
    g["x"] = [default_x if s == "default" else cast(s) for s in g.setting]
    return g.sort_values("x")


def figure(rep, ab):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#3a3a36", "xtick.color": "#5a5a55",
                         "ytick.color": "#5a5a55"})
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.6))
    # A: forest of main contrasts
    A = ax[0, 0]
    show = ["informative: teacher - none", "transferable: kd - none", "transferable: kd_feat - none",
            "transferable: align_shared - none", "transferable: recon_rna - none",
            "transferable: recon_specific - none", "control: self_train - none (unlabelled pool)",
            "RNA-specific: kd - self_train", "dark knowledge: kd - kd_permuted"]
    labels = ["Teacher (sees RNA)", "KD soft labels", "KD + features", "Shared-RNA alignment",
              "RNA reconstruction", "RNA-specific reconstruction", "Self-training (no RNA)",
              "KD minus self-training", "KD minus permuted KD"]
    for j, sc in enumerate(SC):
        g = rep[rep.scenario == sc].set_index("comparison")
        ys = np.arange(len(show)) + (j - 1.5) * 0.18
        est = [g.loc[c, "estimate"] if c in g.index else np.nan for c in show]
        lo = [g.loc[c, "ci_lo"] if c in g.index else np.nan for c in show]
        hi = [g.loc[c, "ci_hi"] if c in g.index else np.nan for c in show]
        A.errorbar(est, ys, xerr=[np.subtract(est, lo), np.subtract(hi, est)], fmt="o", ms=4,
                   color=COL[sc], ecolor=COL[sc], elinewidth=1.5, capsize=0, label=SC_LABEL[sc])
    A.axvline(0, color="#8a8a85", lw=1)
    A.set_yticks(range(len(show)), labels)
    A.invert_yaxis()
    A.set_xlabel("Change in macro-F1 vs comparator (95% CI, patient bootstrap)")
    A.set_title("A  Gains by mechanism and scenario", loc="left", fontweight="bold")
    A.legend(frameon=False, fontsize=8, loc="lower right")
    if len(ab):
        panels = [(ax[0, 1], "jitter_um", 0, "Pairing offset (µm), teacher RNA from a nearby cell",
                   "B  Adjacent-section pairing: how far can space carry RNA?"),
                  (ax[1, 0], "label_frac", 0.05, "Fraction of training cells labelled (log scale)",
                   "C  Label budget")]
        for axx, ex, dx, xl, tt in panels:
            for sc in SC[1:]:
                for q, ls in (("delta_trans", "-"), ("delta_info", ":")):
                    g = _series(ab, ex, sc, q, dx)
                    axx.plot(g.x, g.estimate, ls, marker="o", ms=4, lw=2, color=COL[sc])
                    if q == "delta_trans":
                        axx.fill_between(g.x, g.ci_lo, g.ci_hi, color=COL[sc], alpha=0.15, lw=0)
            axx.axhline(0, color="#8a8a85", lw=1)
            axx.set_xlabel(xl)
            axx.set_ylabel("Gain in macro-F1")
            axx.set_title(tt, loc="left", fontweight="bold")
            if ex == "label_frac":
                axx.set_xscale("log")
        ax[0, 1].plot([], [], "-", color="#3a3a36", label="Transferable (KD - none), band = 95% CI")
        ax[0, 1].plot([], [], ":", color="#3a3a36", label="Informative (teacher - none)")
        for sc in SC[1:]:
            ax[0, 1].plot([], [], "s", color=COL[sc], label=SC_LABEL[sc])
        ax[0, 1].legend(frameon=False, fontsize=8)
        # D: graded nulls
        Dx = ax[1, 1]
        order = ["default", "within_niche", "within_type", "type_average", "all"]
        names = ["Real RNA", "Shuffled\nwithin niche", "Shuffled\nwithin type", "Type\naverage", "Shuffled\nall"]
        for j, sc in enumerate(SC[1:]):
            g = ab[(ab.scenario == sc) & (ab.quantity == "delta_trans") & ab.experiment.isin(["teacher_null", "default"])]
            g = g.set_index("setting").reindex(order)
            xs = np.arange(len(order)) + (j - 1) * 0.25
            Dx.bar(xs, g.estimate, width=0.23, color=COL[sc], label=SC_LABEL[sc])
            Dx.errorbar(xs, g.estimate, yerr=[g.estimate - g.ci_lo, g.ci_hi - g.estimate], fmt="none",
                        ecolor="#3a3a36", elinewidth=1)
        Dx.axhline(0, color="#8a8a85", lw=1)
        Dx.set_xticks(range(len(order)), names)
        Dx.set_ylabel("Transferable gain (KD - none), macro-F1")
        Dx.set_title("D  Graded nulls on the teacher's RNA", loc="left", fontweight="bold")
        Dx.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fn = os.path.join(D, "stage2_overview.png")
    fig.savefig(fn, dpi=150)
    print("saved", fn)


if __name__ == "__main__":
    main()
