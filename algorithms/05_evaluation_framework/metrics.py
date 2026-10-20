"""Evaluation metrics for ctSVG / ctSVP detection, following Li, Ge & Jiang (2025).

Input format (long tables):
    results : columns [method, cell_type, feature, p_value]   (p_value = -1 or NaN if the method did not test it)
    truth   : columns [cell_type, feature, is_true]            (one row per tested pair)

Metrics per method:
    EP        early precision: precision among the top-R ranked pairs, R = number of true positives
    AUPRC     area under the precision-recall curve (untested pairs ranked last)
    FDP@q     false-discovery proportion among BH-adjusted p < q
    power@q   recall among BH-adjusted p < q
    typeI     fraction of null pairs with raw p < alpha
    KS        Kolmogorov-Smirnov distance of null p-values from Uniform(0, 1)
    n_untested_true  true pairs dropped by the method's QC filters
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import kstest
from sklearn.metrics import average_precision_score


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    tmp = np.empty(n)
    tmp[order] = adj
    out[ok] = tmp
    return out


def _one_method(df: pd.DataFrame, q: float, alpha: float) -> dict:
    p = df["p_value"].to_numpy(float)
    p = np.where(p < 0, np.nan, p)                       # -1 = not tested
    y = df["is_true"].to_numpy(bool)
    tested = ~np.isnan(p)
    score = np.where(tested, -np.log10(np.clip(p, 1e-300, 1)), -1.0)   # untested ranked last

    R = int(y.sum())
    top = np.argsort(-score, kind="stable")[:R]
    ep = y[top].mean() if R else np.nan
    auprc = average_precision_score(y, score) if 0 < R < len(y) else np.nan

    padj = bh(p)
    called = np.nan_to_num(padj, nan=1.0) < q
    fdp = (called & ~y).sum() / max(called.sum(), 1)
    power = (called & y).sum() / max(R, 1)

    null_p = p[~y & tested]
    type1 = (null_p < alpha).mean() if null_p.size else np.nan
    ks = kstest(null_p, "uniform").statistic if null_p.size else np.nan
    return dict(EP=ep, AUPRC=auprc, **{f"FDP@{q}": fdp, f"power@{q}": power,
                f"typeI@{alpha}": type1}, KS=ks,
                n_tested=int(tested.sum()), n_untested_true=int((~tested & y).sum()))


def evaluate(results: pd.DataFrame, truth: pd.DataFrame, q: float = 0.05, alpha: float = 0.05,
             by_cell_type: bool = False) -> pd.DataFrame:
    """Score every method in `results` against `truth`. Pairs missing from a method count as untested."""
    rows = []
    for m, r in results.groupby("method"):
        df = truth.merge(r[["cell_type", "feature", "p_value"]], on=["cell_type", "feature"], how="left")
        df["p_value"] = df["p_value"].fillna(-1)
        groups = df.groupby("cell_type") if by_cell_type else [("all", df)]
        for ct, d in groups:
            rows.append(dict(method=m, cell_type=ct, **_one_method(d, q, alpha)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stress-test helpers
def rotate_coords(coords: np.ndarray, degrees: float) -> np.ndarray:
    """Rotate (n, 2) coordinates about their centroid (ctSVGbench rotation test)."""
    t = np.deg2rad(degrees)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    c = coords.mean(0)
    return (coords - c) @ R.T + c


def permute_within_type(coords: np.ndarray, ctype: np.ndarray, rng=None) -> np.ndarray:
    """Shuffle locations among cells of the same type: keeps composition, destroys within-type patterns."""
    rng = rng or np.random.default_rng()
    out = coords.copy()
    for t in np.unique(ctype):
        idx = np.flatnonzero(ctype == t)
        out[idx] = coords[rng.permutation(idx)]
    return out


def jaccard_topk(a: pd.DataFrame, b: pd.DataFrame, k: int = 20) -> float:
    """Overlap of the top-k (cell_type, feature) calls between two result tables."""
    top = lambda d: set(map(tuple, d.nsmallest(k, "p_value")[["cell_type", "feature"]].to_numpy()))
    A, B = top(a), top(b)
    return len(A & B) / max(len(A | B), 1)


if __name__ == "__main__":                              # tiny self-check
    rng = np.random.default_rng(0)
    truth = pd.DataFrame({"cell_type": ["T"] * 100, "feature": [f"m{i}" for i in range(100)],
                          "is_true": [True] * 10 + [False] * 90})
    good = truth.assign(method="good", p_value=np.r_[rng.uniform(0, 1e-4, 10), rng.uniform(size=90)])
    bad = truth.assign(method="random", p_value=rng.uniform(size=100))
    print(evaluate(pd.concat([good, bad]), truth).round(3).to_string(index=False))
