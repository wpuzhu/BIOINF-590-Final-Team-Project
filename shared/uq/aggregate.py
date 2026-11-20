"""
Collect per-job result CSVs from the cluster and produce the pre-registered UQ tables.

  python aggregate.py --stage 2 --results $SCRATCH/results/stage2 --out ../../stage2_privileged_transfer/results/real
  python aggregate.py --stage 1 --results $SCRATCH/results/stage1 --out ../../stage1_cross_modal_pretraining/results/real

Stage 2 tables (per dataset): main contrasts (same as the simulation), ablation deltas,
transfer efficiency, conformal coverage/set size from the saved probabilities
(calibration = held-out calibration patients, never trained on).
Stage 1 tables (per dataset x modality x model x mode x labels): C1 gain, vs-PCA, identity,
correct-map, A3 curve; transfer efficiency against the in-modality model on the same target.
Unit = patient throughout; seeds are the inner level of the hierarchical bootstrap.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uq_stats import conformal_summary, ece, hierarchical_bootstrap, improvement_report, paired_table  # noqa: E402

S2_MAIN = {
    "informative: teacher - none": ("teacher", "none"),
    "transferable: kd - none": ("kd", "none"),
    "transferable: kd_feat - none": ("kd_feat", "none"),
    "transferable: align_infonce - none": ("align_infonce", "none"),
    "transferable: align_cca - none": ("align_cca", "none"),
    "transferable: recon_rna - none": ("recon_rna", "none"),
    "transferable: recon_specific - none": ("recon_specific", "none"),
    "control: self_train - none": ("self_train", "none"),
    "control: label_smooth - none": ("label_smooth", "none"),
    "RNA-specific: kd - self_train": ("kd", "self_train"),
    "dark knowledge: kd - kd_permuted": ("kd", "kd_permuted"),
}


def stage2(res, out, margin):
    df = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(res, "*", "*", "fold*_seed*.csv"))])
    main = df[df.setting == "default"]
    rep = improvement_report(main, "arm", "macro_f1", S2_MAIN, margin=margin, alternative="two-sided",
                             strata=("dataset",))
    rep.to_csv(os.path.join(out, "stage2_main_improvements.csv"), index=False)
    rows = []
    for (ds, st), g in df.groupby(["dataset", "setting"]):
        for q, (a1, a0) in {"delta_info": ("teacher", "none"), "delta_trans": ("kd", "none")}.items():
            d = paired_table(g, "arm", "macro_f1", a1, a0)
            if len(d):
                rows.append(dict(dataset=ds, setting=st, quantity=q, **hierarchical_bootstrap(d)))
    pd.DataFrame(rows).to_csv(os.path.join(out, "stage2_ablation_deltas.csv"), index=False)
    conf = []
    for f in glob.glob(os.path.join(res, "*", "default", "fold*_seed*_probs.npz")):
        z = np.load(f, allow_pickle=True)
        ds = f.split(os.sep)[-3]
        cal, te = z["split"] == "calib", z["split"] == "test"
        ok = z["y"] >= 0
        for arm in [k for k in z.files if k not in ("y", "patient", "split")]:
            P = z[arm]
            cs = conformal_summary(P[cal & ok], z["y"][cal & ok], P[te & ok], z["y"][te & ok], alpha=0.1)
            conf.append(dict(dataset=ds, file=os.path.basename(f), arm=arm, coverage=cs.coverage[0],
                             set_size=cs.set_size[0], ece=ece(P[te & ok], z["y"][te & ok])))
    if conf:
        c = pd.DataFrame(conf)
        c.to_csv(os.path.join(out, "stage2_conformal_by_job.csv"), index=False)
        c.groupby(["dataset", "arm"])[["coverage", "set_size", "ece"]].agg(["mean", "std"]).to_csv(
            os.path.join(out, "stage2_conformal_summary.csv"))
    print(rep.round(3).to_string(index=False))


def stage1(res, out, margin):
    df = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(res, "*", "*", "*.csv"))])
    df["arm_q"] = np.where(df.arm == "A3_overlap", "A3_" + df.q.round(2).astype(str), df.arm)
    pca = df[df.model == "pca"].assign(arm_q="pca")
    rows = []
    for (ds, mod, model, mode, k), g in df[df.model != "pca"].groupby(["dataset", "modality", "model", "mode", "labels"]):
        g = pd.concat([g, pca[(pca.dataset == ds) & (pca.modality == mod) & (pca.labels == k)]])
        comps = {"C1 gain: pretrained - scratch": ("pretrained", "scratch"),
                 "C2: pretrained - pca": ("pretrained", "pca"),
                 "identity: pretrained - A1_reinit": ("pretrained", "A1_reinit"),
                 "correct map: pretrained - A2_random_map": ("pretrained", "A2_random_map")}
        for q in sorted(set(g.arm_q) - {"pretrained", "scratch", "pca", "A1_reinit", "A2_random_map"}):
            comps[f"overlap: {q} - scratch"] = (q, "scratch")
        r = improvement_report(g, "arm_q", "macro_f1", comps, margin=margin, alternative="two-sided")
        if len(r):
            r.insert(0, "labels", k); r.insert(0, "mode", mode); r.insert(0, "model", model)
            r.insert(0, "modality", mod); r.insert(0, "dataset", ds)
            rows.append(r)
    if not rows:
        print("no foundation-model results yet (PCA only):")
        print(pca.groupby(["dataset", "modality", "labels"]).macro_f1.mean().round(3).to_string())
        return
    rep = pd.concat(rows)
    native = {"nicheformer": "rna", "spatium": "protein"}
    rep["direction"] = [("in-modality" if native.get(m) == mod else "cross-modal") for m, mod in zip(rep.model, rep.modality)]
    rep.to_csv(os.path.join(out, "stage1_contrasts.csv"), index=False)
    g = rep[rep.comparison.str.startswith("C1")].set_index(["dataset", "modality", "mode", "labels", "direction"]).estimate
    eff = [dict(dataset=ds, modality=mod, mode=mo, labels=k, cross_gain=v, in_gain=g.get((ds, mod, mo, k, "in-modality")),
                efficiency=v / g[(ds, mod, mo, k, "in-modality")] if (ds, mod, mo, k, "in-modality") in g.index
                and g[(ds, mod, mo, k, "in-modality")] > 0.01 else np.nan)
           for (ds, mod, mo, k, d), v in g.items() if d == "cross-modal"]
    pd.DataFrame(eff).to_csv(os.path.join(out, "stage1_transfer_efficiency.csv"), index=False)
    print(rep.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, choices=[1, 2], required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--margin", type=float, default=0.02, help="pre-registered equivalence margin (macro-F1)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    (stage1 if a.stage == 1 else stage2)(a.results, a.out, a.margin)


if __name__ == "__main__":
    main()
