"""
Stage 2 simulation: can RNA seen only at training time improve a protein + spatial-graph student?

Compares, on identical splits (patients x seeds x labelled cells):
  none            student on labelled cells only (protein + 1-hop/2-hop neighbour protein)
  label_smooth    none + label smoothing 0.1            (control: is KD just regularisation?)
  self_train      student distilled from a PROTEIN-ONLY teacher on all training cells
                                                          (control: is the gain the unlabelled pool?)
  teacher         protein + graph + RNA model, evaluated WITH RNA    -> informative gain
  kd              soft-label distillation from the RNA teacher         -> transferable gain
  kd_feat         kd + feature distillation (match teacher hidden layer)
  kd_permuted     kd with teacher probabilities permuted among cells of the same teacher argmax
                  (keeps pseudo-labels, destroys dark knowledge)
  align_shared    hard labels + auxiliary head regressing the RNA part shared with protein
                  (top CCA variates; closed-form stand-in for contrastive alignment)
  recon_rna       hard labels + auxiliary head reconstructing all RNA PCs
  recon_specific  hard labels + auxiliary head reconstructing the RNA-specific residual
                  (RNA minus its linear prediction from protein)

Ablations (none / teacher / kd per setting): label fraction, graded nulls on the teacher's RNA,
teacher sees only neighbourhood RNA, adjacent-section pairing jitter, graph removal, panel size.

Unit of analysis = patient (5 folds x 2 test patients). Usage:
  python run_stage2_sim.py --experiment main
  python run_stage2_sim.py --experiment ablations --seeds 2
  python run_stage2_sim.py --quick            # 1 scenario, 1 seed, smoke test
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.metrics import f1_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path[:0] = [os.path.join(ROOT, "shared", "sim"), os.path.join(ROOT, "shared", "models"),
                os.path.join(ROOT, "shared", "uq")]
from np_mlp import NumpyMLP  # noqa: E402
from paired_spatial_sim import (SCENARIOS, SimConfig, graded_null, jitter_pairing,  # noqa: E402
                                neighbour_mean, simulate)
from uq_stats import conformal_summary  # noqa: E402

K = 8
MAIN_ARMS = ["none", "label_smooth", "self_train", "teacher", "kd", "kd_feat", "kd_permuted",
             "align_shared", "recon_rna", "recon_specific"]
ABLATIONS = ([("label_frac", v) for v in (0.01, 0.02, 0.2, 0.5)]
             + [("teacher_null", v) for v in ("all", "within_type", "within_niche", "type_average")]
             + [("teacher_input", "rna_neighbourhood_only")]
             + [("jitter_um", v) for v in (15, 30, 60, 120, 240)]
             + [("graph", "off")]
             + [("panel", v) for v in (5, 10)])


# --------------------------------------------------------------------------------------
def cca_shared(X, R, k=8, reg=1e-3):
    """RNA canonical weights for the top-k CCA directions shared with student features X."""
    Xc, Rc = X - X.mean(0), R - R.mean(0)
    n = len(X)
    def inv_sqrt(C):
        w, V = np.linalg.eigh(C + reg * np.eye(len(C)))
        return V @ np.diag(w ** -0.5) @ V.T
    Wx, Wr = inv_sqrt(Xc.T @ Xc / n), inv_sqrt(Rc.T @ Rc / n)
    U, s, Vt = np.linalg.svd(Wx @ (Xc.T @ Rc / n) @ Wr)
    return Wr @ Vt[:k].T, s[:k]


def specific_residual(X, R, reg=1.0):
    Xc = np.c_[X, np.ones(len(X))]
    B = np.linalg.solve(Xc.T @ Xc + reg * np.eye(Xc.shape[1]), Xc.T @ R)
    return R - Xc @ B


def permute_within_argmax(P, seed):
    rng = np.random.default_rng(seed)
    Q, a = P.copy(), P.argmax(1)
    for c in np.unique(a):
        idx = np.where(a == c)[0]
        Q[idx] = P[rng.permutation(idx)]
    return Q


# --------------------------------------------------------------------------------------
def build(d, setting=None):
    """Student and teacher feature matrices for one ablation setting."""
    kind, val = setting if setting else (None, None)
    prot = d.protein[:, :val] if kind == "panel" else d.protein
    S = [prot]
    if kind != "graph":
        S += [neighbour_mean(prot, d.coords, d.patient, 40), neighbour_mean(prot, d.coords, d.patient, 100)]
    S = np.hstack(S)
    rna = d.rna
    if kind == "teacher_null":
        rna = graded_null(rna, val, d.coarse, d.niche, d.patient, seed=1)
    elif kind == "teacher_input":
        rna = neighbour_mean(rna, d.coords, d.patient, 40)
    elif kind == "jitter_um":
        rna = rna[jitter_pairing(d.coords, d.patient, val, seed=1)]
    return S, rna


def metrics(y, P, patient):
    out = []
    for p in np.unique(patient):
        m = patient == p
        out.append((p, f1_score(y[m], P[m].argmax(1), average="macro"),
                    float(-np.log(P[m, y[m]] + 1e-9).mean())))
    return out


def run_job(scenario, d, S, rna, fold, seed, arms, label_frac, tag, setting, conformal):
    pats = np.unique(d.patient)
    test_p = pats[fold * 2: fold * 2 + 2]
    te = np.isin(d.patient, test_p)
    tr = ~te
    rng = np.random.default_rng(1000 * seed + fold)
    lab = rng.random(tr.sum()) < label_frac
    y_tr = d.fine[tr]
    pca = PCA(30, random_state=0).fit(rna[tr])
    R = pca.transform(rna)
    T_feat = np.c_[S, R]
    E = 200
    mk = lambda: NumpyMLP(K, epochs=E, seed=seed)
    probs, models = {}, {}

    def need(a):
        return a in arms

    base = mk().fit(S[tr], y_tr, lab)
    models["none"] = base
    if any(need(a) for a in ("teacher", "kd", "kd_feat", "kd_permuted")):
        models["teacher"] = mk().fit(T_feat[tr], y_tr, lab)
        soft = models["teacher"].predict_proba(T_feat[tr], T=2.0)
    if need("label_smooth"):
        models["label_smooth"] = mk().fit(S[tr], y_tr, lab, smooth=0.1)
    if need("self_train"):
        models["self_train"] = mk().fit(S[tr], y_tr, lab, soft=base.predict_proba(S[tr], T=2.0), lam=0.3)
    if need("kd"):
        models["kd"] = mk().fit(S[tr], y_tr, lab, soft=soft, lam=0.3)
    if need("kd_feat"):
        models["kd_feat"] = mk().fit(S[tr], y_tr, lab, soft=soft, lam=0.3,
                                     target=models["teacher"].hidden(T_feat[tr]), beta=1.0)
    if need("kd_permuted"):
        models["kd_permuted"] = mk().fit(S[tr], y_tr, lab, soft=permute_within_argmax(soft, seed), lam=0.3)
    if need("align_shared"):
        W, _ = cca_shared(S[tr], R[tr])
        models["align_shared"] = mk().fit(S[tr], y_tr, lab, target=(R[tr] - R[tr].mean(0)) @ W, beta=1.0)
    if need("recon_rna"):
        models["recon_rna"] = mk().fit(S[tr], y_tr, lab, target=R[tr], beta=1.0)
    if need("recon_specific"):
        models["recon_specific"] = mk().fit(S[tr], y_tr, lab, target=specific_residual(S[tr], R[tr]), beta=1.0)

    rows = []
    for a in arms:
        X = T_feat if a == "teacher" else S
        P = models[a].predict_proba(X[te])
        probs[a] = P
        for p, f1, nll in metrics(d.fine[te], P, d.patient[te]):
            row = dict(scenario=scenario, experiment=tag, setting=str(setting), arm=a, patient=int(p),
                       seed=seed, fold=fold, macro_f1=f1, nll=nll, label_frac=label_frac)
            if conformal:
                # cross-conformal by patient: calibrate on the other test patient (exchangeable unit)
                yt, pt = d.fine[te], d.patient[te]
                cal, tst = pt != p, pt == p
                cs = conformal_summary(P[cal], yt[cal], P[tst], yt[tst], alpha=0.1)
                row.update(coverage=float(cs.coverage[0]), set_size=float(cs.set_size[0]))
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", choices=["main", "ablations"], default="main")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--label_frac", type=float, default=0.05)
    ap.add_argument("--cells", type=int, default=800)
    ap.add_argument("--n_jobs", type=int, default=2)
    ap.add_argument("--scenarios", nargs="*", default=list(SCENARIOS))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=os.path.join(ROOT, "stage2_privileged_transfer", "results", "simulation"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.quick:
        a.scenarios, a.seeds = ["transferable_protein"], 1
    t0 = time.time()
    for sc in a.scenarios:
        d = simulate(SimConfig(scenario=sc, cells_per_patient=a.cells, seed=0))
        if a.experiment == "main":
            S, rna = build(d)
            jobs = [delayed(run_job)(sc, d, S, rna, f, s, MAIN_ARMS, a.label_frac, "main", "default", True)
                    for f in range(5) for s in range(a.seeds)]
        else:
            jobs = []
            for st in ABLATIONS:
                S, rna = build(d, st)
                lf = st[1] if st[0] == "label_frac" else a.label_frac
                jobs += [delayed(run_job)(sc, d, S, rna, f, s, ["none", "teacher", "kd"], lf,
                                          st[0], st[1], False) for f in range(5) for s in range(a.seeds)]
        rows = sum(Parallel(n_jobs=a.n_jobs)(jobs), [])
        df = pd.DataFrame(rows)
        fn = os.path.join(a.out, f"raw_{a.experiment}_{sc}.csv")
        df.to_csv(fn, index=False)
        print(f"[{time.time() - t0:6.0f}s] {sc}: {len(df)} rows -> {fn}", flush=True)
        print(df.groupby(["experiment", "setting", "arm"]).macro_f1.mean().unstack().round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
