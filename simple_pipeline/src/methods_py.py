"""Step 2a: Python detection methods (zero-shot: fixed defaults, no tuning).

Global baselines (ignore cell type, i.e. the confounded analysis the benchmark warns about):
    global_moran, global_sparkx
Cell-type-specific (ct) methods, run on the cells of ONE type in ONE image, one per model class
from the benchmark:
    ct_moran         spatial autocorrelation (kNN Moran's I)
    ct_sparkx_py     fixed-effect projection kernels (Python port of the SPARK-X statistic:
                     1 linear + 5 Gaussian + 5 cosine coordinate transforms, Cauchy combination;
                     the same model class as CTSV/C-SIDE x-y bases)
    ct_kernel_score  random-effect Gaussian-kernel variance-component score test (the class used
                     by CELINA/STANCE), 3 bandwidths, Satterthwaite null, Cauchy combination
Output: results/raw/py_<mode>.csv with columns dataset,image,cell_type,feature,method,stat,p_value,n_cells
<mode> = "observed" or "permuted" (coordinates shuffled within each cell type: calibration null).
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors

from config import FDR, KNN, MAX_CELLS_PER_TYPE, MIN_CELLS_PER_TYPE, PROC, RES, SEED


# ----------------------------------------------------------------------------- helpers
def acat(p: np.ndarray, axis: int = 0) -> np.ndarray:
    """Cauchy combination of p-values (Liu & Xie 2020), as used by SPARK-X/CTSV/CELINA."""
    p = np.clip(np.asarray(p, float), 1e-15, 1 - 1e-15)
    t = np.tan((0.5 - p) * np.pi).mean(axis=axis)
    return np.clip(0.5 - np.arctan(t) / np.pi, 0, 1)


def zscore(Y):
    Y = Y - Y.mean(0)
    sd = Y.std(0)
    sd[sd == 0] = 1
    return Y / sd


def knn_W(coords, k=KNN):
    k = min(k, len(coords) - 1)
    _, ind = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords)
    rows = np.repeat(np.arange(len(coords)), k)
    return csr_matrix((np.full(rows.size, 1.0 / k), (rows, ind[:, 1:].ravel())),
                      shape=(len(coords),) * 2)


# ----------------------------------------------------------------------------- methods
def morans_i(Y, coords):
    """Moran's I per column, one-sided p (positive autocorrelation) under the normality assumption."""
    n = len(Y)
    W = knn_W(coords)
    Z = Y - Y.mean(0)
    den = (Z ** 2).sum(0)
    den[den == 0] = np.nan
    I = (Z * (W @ Z)).sum(0) / den            # S0 = n for row-standardised W
    S0 = n
    Wsym = W + W.T
    S1 = 0.5 * Wsym.multiply(Wsym).sum()
    S2 = ((np.asarray(W.sum(1)).ravel() + np.asarray(W.sum(0)).ravel()) ** 2).sum()
    EI = -1.0 / (n - 1)
    VI = (n * n * S1 - n * S2 + 3 * S0 ** 2) / ((n * n - 1) * S0 ** 2) - EI ** 2
    p = stats.norm.sf((I - EI) / np.sqrt(VI))
    return I, np.nan_to_num(p, nan=1.0)


def sparkx_py(Y, coords):
    """Python port of the SPARK-X test statistic with projection kernels.
    For each coordinate transform E (n x 2): T = n * ||Q'y||^2 / ||y||^2 ~ chi2(rank E)."""
    n = len(Y)
    C = zscore(coords.astype(float))
    lengths = np.quantile(np.abs(C), [0.2, 0.4, 0.6, 0.8, 1.0])
    transforms = [C]
    transforms += [np.exp(-C ** 2 / (2 * l ** 2)) for l in lengths]
    transforms += [np.cos(2 * np.pi * C / l) for l in lengths]
    Yc = Y - Y.mean(0)
    yy = (Yc ** 2).sum(0)
    yy[yy == 0] = np.nan
    P, T = [], []
    for E in transforms:
        E = E - E.mean(0)
        Q, R = np.linalg.qr(E)
        r = int((np.abs(np.diag(R)) > 1e-8 * np.abs(R).max()).sum())
        if r == 0:
            continue
        stat = n * ((Q[:, :r].T @ Yc) ** 2).sum(0) / yy
        T.append(stat)
        P.append(stats.chi2.sf(stat, r))
    P = np.nan_to_num(np.array(P), nan=1.0)
    return np.nanmax(np.array(T), 0), acat(P, axis=0)


def kernel_score(Y, coords, quantiles=(0.1, 0.25, 0.5)):
    """Variance-component score test y = mu + g + e, g ~ N(0, tau K), H0: tau = 0.
    Q = y'PKPy / sigma^2; null approximated by a scaled chi-square (Satterthwaite)."""
    n = len(Y)
    D = squareform(pdist(coords))
    Yc = Y - Y.mean(0)
    s2 = (Yc ** 2).sum(0) / (n - 1)
    s2[s2 == 0] = np.nan
    J = np.eye(n) - 1.0 / n
    P, T = [], []
    dvec = D[np.triu_indices(n, 1)]
    for q in quantiles:
        h = np.quantile(dvec, q)
        K = np.exp(-(D ** 2) / (2 * h ** 2))
        PKP = J @ K @ J
        lam = np.clip(np.linalg.eigvalsh(PKP), 0, None)
        mean, var = lam.sum(), 2 * (lam ** 2).sum()
        a, df = var / (2 * mean), 2 * mean ** 2 / var
        Qs = (Yc * (K @ Yc)).sum(0) / s2
        T.append(Qs)
        P.append(stats.chi2.sf(Qs / a, df))
    P = np.nan_to_num(np.array(P), nan=1.0)
    return np.nanmax(np.array(T), 0), acat(P, axis=0)


CT_METHODS = {"ct_moran": morans_i, "ct_sparkx_py": sparkx_py, "ct_kernel_score": kernel_score}
GLOBAL_METHODS = {"global_moran": morans_i, "global_sparkx_py": sparkx_py}


# ----------------------------------------------------------------------------- driver
def iter_images():
    for meta_f in sorted(PROC.glob("*/*__meta.csv")):
        dataset, image = meta_f.parent.name, meta_f.name.replace("__meta.csv", "")
        meta = pd.read_csv(meta_f)
        expr = pd.read_csv(meta_f.with_name(f"{image}__expr.csv"), index_col=0)
        yield dataset, image, meta, expr


def permute_within_type(coords, ctype, rng):
    out = coords.copy()
    for t in np.unique(ctype):
        idx = np.flatnonzero(ctype == t)
        out[idx] = coords[rng.permutation(idx)]
    return out


def run(mode="observed"):
    rng = np.random.default_rng(SEED)
    rows = []
    for dataset, image, meta, expr in iter_images():
        coords = meta[["x", "y"]].to_numpy(float)
        ctype = meta["cell_type"].to_numpy()
        if mode == "permuted":
            coords = permute_within_type(coords, ctype, rng)
        X = expr.to_numpy(float)
        feats = expr.columns
        if mode == "observed":                                   # global baselines
            for m, f in GLOBAL_METHODS.items():
                Xs, cs = X, coords
                if len(X) > 3000 and m != "global_moran":
                    pick = rng.choice(len(X), 3000, replace=False); Xs, cs = X[pick], coords[pick]
                st, p = f(Xs, cs)
                rows += [dict(dataset=dataset, image=image, cell_type="ALL", feature=g, method=m,
                              stat=s, p_value=pv, n_cells=len(Xs)) for g, s, pv in zip(feats, st, p)]
        for t in sorted(np.unique(ctype)):
            idx = np.flatnonzero(ctype == t)
            if len(idx) < MIN_CELLS_PER_TYPE:
                continue
            if len(idx) > MAX_CELLS_PER_TYPE:
                idx = np.sort(rng.choice(idx, MAX_CELLS_PER_TYPE, replace=False))
            Xt, ct = X[idx], coords[idx]
            keep = Xt.std(0) > 0
            for m, f in CT_METHODS.items():
                st, p = f(Xt[:, keep], ct)
                rows += [dict(dataset=dataset, image=image, cell_type=t, feature=g, method=m,
                              stat=s, p_value=pv, n_cells=len(idx))
                         for g, s, pv in zip(feats[keep], st, p)]
        print(f"  [{mode}] {dataset}/{image} done")
    out = RES / "raw"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / f"py_{mode}.csv", index=False)
    return df


if __name__ == "__main__":
    for mode in (sys.argv[1:] or ["observed", "permuted"]):
        run(mode)
