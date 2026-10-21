"""Step 3: aggregate detection results across methods, images and datasets.

Outputs (results/aggregate/):
    pairs_<mode>.csv           one row per (dataset, image, cell_type, feature): p / q / sig per method,
                               votes, ACAT-combined p, consensus call
    method_agreement.csv       pairwise Jaccard of significant sets and Spearman of -log10 p
    replication.csv            per (dataset, cell_type, feature): images tested / images with consensus
    cross_dataset_common.csv   consensus themes (coarse cell type x harmonised marker) in >= 2 datasets
Consensus rule: ACAT q < FDR AND a strict majority of the ct methods that tested the pair call it (BH q < FDR).
"""
import re
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from config import COARSE_TYPE, FDR, RES
from methods_py import acat

AGG = RES / "aggregate"
KEY = ["dataset", "image", "cell_type", "feature"]


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n:
        o = np.argsort(p[ok])
        r = p[ok][o] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
        tmp = np.empty(n); tmp[o] = r; out[ok] = tmp
    return out


def load(mode):
    files = [RES / "raw" / f"py_{mode}.csv"] + sorted((RES / "raw").glob(f"r_{mode}*.csv"))
    parts = [pd.read_csv(f) for f in files if f.exists() and f.stat().st_size > 5]
    d = pd.concat(parts, ignore_index=True)
    d["q_value"] = d.groupby(["dataset", "method"])["p_value"].transform(bh)   # BH per method x dataset
    return d


def pair_table(d):
    ct = d[d.method.str.startswith("ct_")]
    methods = sorted(ct.method.unique())
    P = ct.pivot_table(index=KEY, columns="method", values="p_value", aggfunc="first")
    Q = ct.pivot_table(index=KEY, columns="method", values="q_value", aggfunc="first")
    S = (Q < FDR).where(Q.notna())
    out = pd.concat([P.add_prefix("p_"), Q.add_prefix("q_"), S.add_prefix("sig_")], axis=1).reset_index()
    out["n_methods"] = P.notna().sum(1).to_numpy()
    out["votes"] = S.sum(1).to_numpy()
    out["acat_p"] = [acat(row[~np.isnan(row)]) if (~np.isnan(row)).any() else np.nan for row in P.to_numpy()]
    out["acat_q"] = out.groupby("dataset")["acat_p"].transform(bh)
    out["consensus"] = (out.acat_q < FDR) & (out.votes > out.n_methods / 2)
    out["unanimous"] = out.votes == out.n_methods
    return out, methods


def agreement(pairs, methods):
    rows = []
    for a, b in combinations(methods, 2):
        both = pairs[[f"p_{a}", f"p_{b}", f"sig_{a}", f"sig_{b}"]].dropna()
        if both.empty:
            continue
        sa, sb = both[f"sig_{a}"].astype(bool), both[f"sig_{b}"].astype(bool)
        rho = spearmanr(-np.log10(both[f"p_{a}"].clip(1e-300)), -np.log10(both[f"p_{b}"].clip(1e-300)))[0]
        rows.append(dict(method_a=a, method_b=b, n_pairs=len(both),
                         jaccard=(sa & sb).sum() / max((sa | sb).sum(), 1), spearman=rho))
    return pd.DataFrame(rows)


def coarse_type(t):
    t = t.lower()
    return next((lab for pat, lab in COARSE_TYPE if re.search(pat, t)), "other")


def harmonise_marker(m):
    m = re.sub(r"_[A-Z][a-z]?\d{2,3}$", "", m)                      # drop IMC metal tag
    m = m.upper().replace("-", "").replace(" ", "").replace("_", "")
    return {"ASMA": "SMA", "CD3EPSI": "CD3", "PANCYTO": "PANCK", "CYTOKERATIN": "PANCK", "CK": "PANCK",
            "ECADHE": "ECAD", "VIMENTI": "VIMENTIN", "S6P": "PS6", "S6": "PS6"}.get(m, m)


def main(mode="observed"):
    AGG.mkdir(parents=True, exist_ok=True)
    d = load(mode)
    pairs, methods = pair_table(d)
    pairs.to_csv(AGG / f"pairs_{mode}.csv", index=False)
    print(f"[{mode}] ct methods: {methods}; pairs: {len(pairs)}; consensus: {pairs.consensus.sum()}")
    if mode != "observed":
        return
    agreement(pairs, methods).to_csv(AGG / "method_agreement.csv", index=False)

    rep = (pairs.groupby(["dataset", "cell_type", "feature"])
           .agg(n_images=("image", "nunique"), n_consensus=("consensus", "sum"),
                best_acat_q=("acat_q", "min")).reset_index())
    rep.to_csv(AGG / "replication.csv", index=False)

    cons = pairs[pairs.consensus].assign(coarse_type=lambda x: x.cell_type.map(coarse_type),
                                          marker=lambda x: x.feature.map(harmonise_marker))
    common = (cons.groupby(["coarse_type", "marker"])
              .agg(n_datasets=("dataset", "nunique"), datasets=("dataset", lambda s: ";".join(sorted(set(s)))),
                   n_images=("image", "nunique"), cell_types=("cell_type", lambda s: ";".join(sorted(set(s)))))
              .reset_index().query("n_datasets >= 2").sort_values(["n_datasets", "n_images"], ascending=False))
    common.to_csv(AGG / "cross_dataset_common.csv", index=False)
    print(f"  cross-dataset themes (>=2 datasets): {len(common)}")


if __name__ == "__main__":
    for mode in (sys.argv[1:] or ["observed", "permuted"]):
        main(mode)
