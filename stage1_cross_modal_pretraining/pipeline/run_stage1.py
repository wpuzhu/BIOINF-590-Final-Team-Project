"""
Stage 1 on real data. One job = (dataset, target modality, model, arm, mode, fold, seed);
every label level runs inside the job on the same split.

  python run_stage1.py --config configs/stage1.yaml --dataset spatch_colon_codex --modality protein \
         --model nicheformer --arm pretrained --mode frozen --fold 0 --seed 0 --out $SCRATCH/results/stage1
  python run_stage1.py ... --model pca --dry_run

The 2 x 2 design: model in {nicheformer (RNA-pretrained), spatium (protein-pretrained)} x
target modality in {rna, protein}. In-modality cells (nicheformer-rna, spatium-protein) are the
positive controls; cross-modal cells are the question. Pre-registered contrasts (aggregate.py):
  C1 gain       pretrained - scratch
  C2 vs PCA     pretrained - pca
  identity      pretrained - A1_reinit
  correct map   pretrained - A2_random_map
  overlap curve A3_overlap(q) for q in config
  efficiency    gain(cross) / gain(in-modality, same target)
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
sys.path[:0] = [HERE, os.path.join(ROOT, "shared", "data")]
from adapters import get_adapter  # noqa: E402
from spt_io import load  # noqa: E402


def labels_per_class(y, pool, k, rng):
    idx = []
    for c in np.unique(y[pool]):
        cand = pool[y[pool] == c]
        idx.extend(rng.choice(cand, size=min(k, len(cand)), replace=False))
    return np.array(idx)


def main():
    import yaml
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import f1_score
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(HERE, "configs", "stage1.yaml"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--modality", choices=["rna", "protein"], required=True)
    ap.add_argument("--model", default="pca")
    ap.add_argument("--arm", default="pretrained")
    ap.add_argument("--q", type=float, default=1.0, help="A3 overlap fraction")
    ap.add_argument("--mode", choices=["frozen", "finetune"], default="frozen")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.expandvars("$SCRATCH/results/stage1"))
    ap.add_argument("--dry_run", action="store_true")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    t0 = time.time()
    c = load(os.path.expandvars(cfg["datasets"][a.dataset]["path"]))
    key = "X_rna" if a.modality == "rna" else "X_protein"
    names = list(map(str, c["gene_names" if a.modality == "rna" else "protein_names"]))
    ok = ~np.isnan(c[key]).any(1) & (c["cell_type_fine"].astype(str) != "NA")
    X = c[key][ok]
    y_str = c["cell_type_fine"].astype(str)[ok]
    classes, y = np.unique(y_str, return_inverse=True)
    patient, folds = c["patient"].astype(str)[ok], c["fold"].astype(int)[ok]
    test = np.where(folds == a.fold)[0]
    pool = np.where(folds != a.fold)[0]
    rng = np.random.default_rng(a.seed)
    info = dict(dataset=a.dataset, modality=a.modality, model=a.model, arm=a.arm, q=a.q, mode=a.mode,
                fold=a.fold, seed=a.seed, n_cells=int(ok.sum()), n_features=len(names), n_classes=len(classes),
                test_patients=sorted(set(patient[test])))
    print(json.dumps(info, indent=1))
    if a.dry_run and a.model != "pca":
        print("[dry run] data and split OK; model not loaded")
        return
    ad = get_adapter(a.model, cfg.get("models", {}).get(a.model, {}))
    kw = {}
    if a.model == "nicheformer" and a.modality == "protein":
        kw["cognate"] = list(map(str, c["protein_cognate_gene"]))
    if a.model == "spatium" and a.modality == "rna":
        kw["marker_gene"] = json.load(open(os.path.expandvars(cfg["marker_gene_map_json"])))
    ids, mask = ad.encode(X, names, a.modality, a.arm, a.q, rng, patient=patient, **kw)
    rows = []
    Z = ad.embed(ids, mask, a.arm, **({"train_idx": pool} if a.model == "pca" else {})) if a.mode == "frozen" else None
    for k in cfg["labels_per_class"]:
        tr = labels_per_class(y, pool, k, rng) if k != "all" else pool
        if a.mode == "frozen":
            clf = LogisticRegression(C=1.0, max_iter=3000).fit(Z[tr], y[tr])
            yhat = clf.predict(Z[test])
        else:
            predict = ad.finetune(ids, mask, y, tr, a.arm, len(classes), **cfg.get("finetune", {}))
            yhat = predict(test)
        for p in np.unique(patient[test]):
            m = patient[test] == p
            rows.append(dict(**{k_: info[k_] for k_ in ("dataset", "modality", "model", "arm", "q", "mode", "fold", "seed")},
                             labels=k, patient=p,
                             macro_f1=f1_score(y[test][m], yhat[m], average="macro", labels=np.arange(len(classes)),
                                               zero_division=0)))
    df = pd.DataFrame(rows)
    od = os.path.join(a.out, a.dataset, a.modality)
    os.makedirs(od, exist_ok=True)
    fn = os.path.join(od, f"{a.model}_{a.arm}{'' if a.arm != 'A3_overlap' else f'_{a.q:.2f}'}_{a.mode}_fold{a.fold}_seed{a.seed}.csv")
    df.to_csv(fn, index=False)
    print(df.groupby("labels").macro_f1.mean().round(3).to_string())
    print(f"wrote {fn} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
