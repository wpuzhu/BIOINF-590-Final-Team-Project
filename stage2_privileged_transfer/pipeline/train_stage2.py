"""
Stage 2 on real data: one job = (dataset, ablation setting, test fold, seed); it trains every
mechanism on the same split so that all paired differences share labels and seeds.

  python train_stage2.py --config configs/stage2.yaml --dataset cosmx_multiomic_nsclc \
         --setting default --fold 0 --seed 0 --out $SCRATCH/results/stage2
  python train_stage2.py ... --dry_run        # data, split, ablation transforms; no torch

Split per fold: test patients = fold f; calibration patients = fold f+1 (conformal only, never
trained on); training = the rest. Mechanisms (same names as the simulation):
  none, label_smooth, self_train, teacher, kd, kd_feat, kd_permuted, align_infonce, align_cca,
  recon_rna, recon_specific
Ablation settings are defined in the YAML (label fraction, teacher-RNA nulls, neighbourhood-only
RNA, pairing jitter, graph off, panel size, dropped subtype markers).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path[:0] = [os.path.join(ROOT, "shared", "data"), os.path.join(ROOT, "shared", "sim"),
                os.path.join(ROOT, "stage2_privileged_transfer", "simulation")]
from paired_spatial_sim import graded_null, jitter_pairing  # noqa: E402  (generic array ops)
from spt_io import load, radius_edges  # noqa: E402


# ------------------------------------------------------------------ data preparation (numpy)
def lognorm(counts, target=1e4):
    c = np.nan_to_num(np.asarray(counts, float))
    return np.log1p(c / (c.sum(1, keepdims=True) + 1e-9) * target)


def neighbourhood_niche(coarse, edges, n, k=10, seed=0):
    from sklearn.cluster import KMeans
    cats, y = np.unique(coarse, return_inverse=True)
    comp = np.zeros((n, len(cats)))
    np.add.at(comp, (edges[1], y[edges[0]]), 1.0)
    comp /= comp.sum(1, keepdims=True) + 1e-9
    return KMeans(k, n_init=4, random_state=seed).fit_predict(comp)


def prepare(cfg, ds_name, setting, fold, seed, verbose=True):
    ds = cfg["datasets"][ds_name]
    st = cfg["settings"][setting]
    c = load(os.path.expandvars(ds["path"]))
    n = c.n
    rng = np.random.default_rng(seed)
    y_str = c["cell_type_fine"].astype(str)
    classes = np.array(sorted(set(y_str) - {"NA", "nan", ""}))
    y = np.array([np.searchsorted(classes, s) if s in set(classes) else -1 for s in y_str])
    group = c["slide"] if "slide" in c else c["patient"]
    edges = radius_edges(c.coords, group, radius=st.get("radius_um", cfg.get("radius_um", 40)))

    # student input: protein panel (with panel / marker-drop ablations)
    names = list(map(str, c["protein_names"]))
    keep = np.ones(len(names), bool)
    if st.get("drop_markers"):
        keep &= ~np.isin(names, st["drop_markers"])
    if st.get("panel"):
        idx = np.where(keep)[0]
        keep[:] = False
        keep[rng.choice(idx, size=min(st["panel"], len(idx)), replace=False)] = True
    Xp = np.nan_to_num(c["X_protein"][:, keep].astype(np.float32))

    # teacher privileged input: RNA (cell, neighbourhood, or registered neighbourhood)
    pairing = c.get("pairing", {}).get("type", "unpaired")
    src = st.get("teacher_rna", "cell" if pairing == "same_cell" else "registered_neighbourhood")
    if src == "registered_neighbourhood":
        R = np.load(os.path.expandvars(ds["registered_rna"]))       # from prepare_adjacent.py
    else:
        R = lognorm(c["X_rna"])
        if st.get("jitter_um", 0) > 0:
            R = R[jitter_pairing(c.coords, group, st["jitter_um"], seed=seed + 1)]
        if src == "neighbourhood":
            M = np.zeros_like(R); deg = np.zeros(n)
            np.add.at(M, edges[1], R[edges[0]]); np.add.at(deg, edges[1], 1)
            R = np.where(deg[:, None] > 0, M / np.maximum(deg, 1)[:, None], R)
    if st.get("teacher_null", "none") != "none":
        coarse = c["cell_type_coarse"].astype(str)
        niche = neighbourhood_niche(coarse, edges, n, seed=seed)
        R = graded_null(R, st["teacher_null"], coarse, niche, c["patient"], seed=seed + 2)

    # split
    folds = c["fold"].astype(int)
    k = folds.max() + 1
    test = folds == fold
    calib = folds == (fold + 1) % k
    train = ~test & ~calib
    lab_frac = st.get("label_frac", cfg.get("label_frac", 0.05))
    labelled = train & (y >= 0) & (rng.random(n) < lab_frac)
    from sklearn.decomposition import PCA
    pca = PCA(min(50, R.shape[1] - 1), random_state=0).fit(R[train])
    Rp = pca.transform(R).astype(np.float32)
    info = dict(dataset=ds_name, setting=setting, fold=fold, seed=seed, n_cells=int(n),
                n_classes=len(classes), n_markers=int(keep.sum()), teacher_rna=src, pairing=pairing,
                n_train=int(train.sum()), n_labelled=int(labelled.sum()), n_test=int(test.sum()),
                n_calib=int(calib.sum()), n_edges=int(edges.shape[1]),
                test_patients=sorted(map(str, np.unique(c["patient"][test]))))
    if verbose:
        print(json.dumps(info, indent=1))
    return dict(Xp=Xp, R=Rp, y=y, edges=edges, train=train, test=test, calib=calib, labelled=labelled,
                patient=c["patient"].astype(str), group=np.asarray(group).astype(str), classes=classes,
                use_graph=st.get("graph", True), info=info)


# ------------------------------------------------------------------ training (torch)
def train_all(D, cfg, mechanisms, seed, device):
    import torch
    import torch.nn.functional as F
    from models_torch import RNAProjector, Student, Teacher, info_nce, kd_loss
    torch.manual_seed(seed)
    hp = cfg["train"]
    T, lam, beta = hp["T"], hp["lambda"], hp["beta"]
    K = len(D["classes"])
    x = torch.tensor(D["Xp"], device=device)
    r = torch.tensor(D["R"], device=device)
    y = torch.tensor(np.maximum(D["y"], 0), device=device)
    # per-slide subgraphs = mini-batches (cells only connect within a slide)
    slides = []
    for g in np.unique(D["group"]):
        idx = np.where(D["group"] == g)[0]
        loc = -np.ones(len(D["y"]), int); loc[idx] = np.arange(len(idx))
        e = D["edges"][:, np.isin(D["edges"][1], idx)]
        slides.append((torch.tensor(idx, device=device), torch.tensor(loc[e], device=device),
                       bool(D["train"][idx].any())))
    lab = torch.tensor(D["labelled"], device=device)
    trn = torch.tensor(D["train"], device=device)

    def run(model, extra=None, rna_in=False):
        out = {}
        for idx, e, _ in slides:
            o = model(x[idx], e, r[idx]) if rna_in else model(x[idx], e)
            for k_, v in o.items():
                if v is not None:
                    out.setdefault(k_, []).append((idx, v))
        res = {}
        for k_, parts in out.items():
            full = torch.zeros(len(D["y"]), parts[0][1].shape[1], device=device)
            for idx, v in parts:
                full[idx] = v
            res[k_] = full
        return res

    def fit(model, loss_fn, params, rna_in=False):
        opt = torch.optim.AdamW(params, lr=hp["lr"], weight_decay=hp["wd"])
        for ep in range(hp["epochs"]):
            model.train()
            for idx, e, has_train in slides:
                if not has_train:
                    continue
                o = model(x[idx], e, r[idx]) if rna_in else model(x[idx], e)
                loss = loss_fn(o, idx)
                if loss is None:
                    continue
                opt.zero_grad(); loss.backward(); opt.step()
        model.eval()
        return model

    def ce_part(o, idx, smooth=0.0):
        m = lab[idx]
        if m.sum() == 0:
            return None
        return F.cross_entropy(o["logits"][m], y[idx][m], label_smoothing=smooth)

    def add(*terms):
        terms = [t for t in terms if t is not None]
        return sum(terms) if terms else None

    d_p, d_r = x.shape[1], r.shape[1]
    G = D["use_graph"]
    models, probs = {}, {}
    s = Student(d_p, K, use_graph=G).to(device)
    models["none"] = fit(s, lambda o, i: ce_part(o, i), s.parameters())

    if "label_smooth" in mechanisms:
        s = Student(d_p, K, use_graph=G).to(device)
        models["label_smooth"] = fit(s, lambda o, i: ce_part(o, i, 0.1), s.parameters())

    need_teacher = any(m in mechanisms for m in ("teacher", "kd", "kd_feat", "kd_permuted"))
    if need_teacher:
        t = Teacher(d_p, d_r, K, use_graph=G).to(device)
        models["teacher"] = fit(t, lambda o, i: ce_part(o, i), t.parameters(), rna_in=True)
        with torch.no_grad():
            to = run(models["teacher"], rna_in=True)
            pT, hT = F.softmax(to["logits"] / T, 1), to["h"]

    def distil(soft, feat_target=None):
        s = Student(d_p, K, use_graph=G, aux_dim=feat_target.shape[1] if feat_target is not None else 0).to(device)
        def lf(o, i):
            m = trn[i]
            kd = kd_loss(o["logits"][m], soft[i][m], T) if m.any() else None
            ft = beta * F.mse_loss(o["aux"][m], feat_target[i][m]) if feat_target is not None and m.any() else None
            ce = ce_part(o, i)
            return add(lam * ce if ce is not None else None, (1 - lam) * kd if kd is not None else None, ft)
        return fit(s, lf, s.parameters())

    if "self_train" in mechanisms:
        with torch.no_grad():
            p0 = F.softmax(run(models["none"])["logits"] / T, 1)
        models["self_train"] = distil(p0)
    if "kd" in mechanisms:
        models["kd"] = distil(pT)
    if "kd_feat" in mechanisms:
        models["kd_feat"] = distil(pT, (hT - hT.mean(0)) / (hT.std(0) + 1e-6))
    if "kd_permuted" in mechanisms:
        a = pT.argmax(1).cpu().numpy()
        perm = np.arange(len(a))
        g = np.random.default_rng(seed)
        for cl in np.unique(a):
            ii = np.where(a == cl)[0]; perm[ii] = g.permutation(ii)
        models["kd_permuted"] = distil(pT[torch.tensor(perm, device=device)])

    def aux_model(target):
        tg = (target - target[trn].mean(0)) / (target[trn].std(0) + 1e-6)
        s = Student(d_p, K, use_graph=G, aux_dim=tg.shape[1]).to(device)
        def lf(o, i):
            m = trn[i]
            return add(ce_part(o, i), beta * F.mse_loss(o["aux"][m], tg[i][m]) if m.any() else None)
        return fit(s, lf, s.parameters())

    if "recon_rna" in mechanisms:
        models["recon_rna"] = aux_model(r)
    if "recon_specific" in mechanisms or "align_cca" in mechanisms:
        from run_stage2_sim import cca_shared, specific_residual
        tr_np = D["train"]
        Xs = D["Xp"]
        if "recon_specific" in mechanisms:
            res = specific_residual(Xs[tr_np], D["R"][tr_np])
            full = np.zeros_like(D["R"]); full[tr_np] = res
            models["recon_specific"] = aux_model(torch.tensor(full, device=device))
        if "align_cca" in mechanisms:
            W, _ = cca_shared(Xs[tr_np], D["R"][tr_np], k=16)
            models["align_cca"] = aux_model(torch.tensor((D["R"] - D["R"][tr_np].mean(0)) @ W,
                                                         dtype=torch.float32, device=device))
    if "align_infonce" in mechanisms:
        s = Student(d_p, K, use_graph=G).to(device)
        rp = RNAProjector(d_r).to(device)
        def lf(o, i):
            m = trn[i]
            nce = beta * info_nce(o["z"][m], rp(r[i][m]), tau=hp.get("tau", 0.1)) if m.sum() > 8 else None
            return add(ce_part(o, i), nce)
        models["align_infonce"] = fit(s, lf, list(s.parameters()) + list(rp.parameters()))

    with torch.no_grad():
        for name in ["none"] + [m for m in mechanisms if m != "none"]:
            o = run(models[name], rna_in=(name == "teacher"))
            probs[name] = F.softmax(o["logits"], 1).cpu().numpy()
    return probs


def train_all_numpy(D, cfg, mechanisms, seed):
    """CPU smoke-test backend (no torch): numpy MLP on [protein, neighbour-mean protein]; runs
    none / label_smooth / self_train / teacher / kd / kd_permuted. Same outputs as train_all."""
    sys.path.insert(0, os.path.join(ROOT, "shared", "models"))
    from np_mlp import NumpyMLP
    from run_stage2_sim import permute_within_argmax
    n = len(D["y"])
    M = np.zeros_like(D["Xp"]); deg = np.zeros(n)
    np.add.at(M, D["edges"][1], D["Xp"][D["edges"][0]]); np.add.at(deg, D["edges"][1], 1)
    S = np.c_[D["Xp"], M / np.maximum(deg, 1)[:, None]] if D["use_graph"] else D["Xp"]
    Tf = np.c_[S, D["R"]]
    tr, y, lab = D["train"], np.maximum(D["y"], 0), D["labelled"][D["train"]]
    K, hp = len(D["classes"]), cfg["train"]
    mk = lambda: NumpyMLP(K, epochs=150, seed=seed)
    models = {"none": mk().fit(S[tr], y[tr], lab)}
    if "label_smooth" in mechanisms:
        models["label_smooth"] = mk().fit(S[tr], y[tr], lab, smooth=0.1)
    if "self_train" in mechanisms:
        models["self_train"] = mk().fit(S[tr], y[tr], lab, soft=models["none"].predict_proba(S[tr], T=hp["T"]), lam=hp["lambda"])
    if any(m in mechanisms for m in ("teacher", "kd", "kd_permuted")):
        models["teacher"] = mk().fit(Tf[tr], y[tr], lab)
        soft = models["teacher"].predict_proba(Tf[tr], T=hp["T"])
        if "kd" in mechanisms:
            models["kd"] = mk().fit(S[tr], y[tr], lab, soft=soft, lam=hp["lambda"])
        if "kd_permuted" in mechanisms:
            models["kd_permuted"] = mk().fit(S[tr], y[tr], lab, soft=permute_within_argmax(soft, seed), lam=hp["lambda"])
    return {k: v.predict_proba(Tf if k == "teacher" else S) for k, v in models.items() if k in mechanisms or k == "none"}


def evaluate(D, probs):
    from sklearn.metrics import f1_score
    rows = []
    y = D["y"]
    for arm, P in probs.items():
        for part in ("test",):
            m = D[part] & (y >= 0)
            for p in np.unique(D["patient"][m]):
                mm = m & (D["patient"] == p)
                rows.append(dict(**{k: D["info"][k] for k in ("dataset", "setting", "fold", "seed")}, arm=arm,
                                 patient=p, n_cells=int(mm.sum()),
                                 macro_f1=f1_score(y[mm], P[mm].argmax(1), average="macro",
                                                   labels=np.arange(len(D["classes"])), zero_division=0),
                                 nll=float(-np.log(P[mm, y[mm]] + 1e-9).mean())))
    return pd.DataFrame(rows)


def main():
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(HERE, "configs", "stage2.yaml"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--setting", default="default")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.expandvars("$SCRATCH/results/stage2"))
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--backend", choices=["torch", "numpy"], default="torch",
                    help="numpy = CPU smoke test of the whole pipeline (subset of mechanisms)")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    t0 = time.time()
    D = prepare(cfg, a.dataset, a.setting, a.fold, a.seed)
    if a.dry_run:
        print(f"[dry run] data and split OK ({time.time() - t0:.1f}s); mechanisms: "
              f"{cfg['settings'][a.setting].get('mechanisms', cfg['mechanisms'])}")
        return
    mechs = cfg["settings"][a.setting].get("mechanisms", cfg["mechanisms"])
    if a.backend == "numpy":
        dev = "cpu (numpy backend)"
        probs = train_all_numpy(D, cfg, mechs, a.seed)
    else:
        import torch
        dev = "cuda" if torch.cuda.is_available() else "cpu"
        probs = train_all(D, cfg, mechs, a.seed, dev)
    od = os.path.join(a.out, a.dataset, a.setting)
    os.makedirs(od, exist_ok=True)
    tag = f"fold{a.fold}_seed{a.seed}"
    evaluate(D, probs).to_csv(os.path.join(od, f"{tag}.csv"), index=False)
    keep = D["test"] | D["calib"]
    np.savez_compressed(os.path.join(od, f"{tag}_probs.npz"), y=D["y"][keep], patient=D["patient"][keep],
                        split=np.where(D["test"][keep], "test", "calib"), **{k: v[keep] for k, v in probs.items()})
    json.dump(D["info"], open(os.path.join(od, f"{tag}_info.json"), "w"))
    print(f"done in {time.time() - t0:.0f}s on {dev}")


if __name__ == "__main__":
    main()
