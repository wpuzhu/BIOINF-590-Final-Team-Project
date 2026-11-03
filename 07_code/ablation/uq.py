"""
Uncertainty quantification and transferability diagnostics for the input-level ablation framework.

Everything here respects one rule: the unit of resampling / permutation is the unit of
independence. In spatial omics that is the patient (or slide), and inside a slide it is a
spatial block much larger than the neighbourhood radius - never the individual cell.

Contents
  splitting       spatial_blocks, unit_folds
  per-fold scores fold_scores                        (R2 / Pearson per protein x held-out unit)
  bootstrap       cluster_bootstrap_ci               (CI of a mean / mean difference over units)
  permutation     signflip_test                      (paired contrast between two input conditions)
                  group_permutation_importance       (conditional, within sample x cell type)
                  toroidal_shift_test                (association of two spatial fields, keeps autocorrelation)
  conformal       mondrian_conformal                 (per-cell-type split conformal intervals + coverage)
                  weighted_conformal_coverage        (covariate-shift-weighted conformal, Tibshirani et al. 2019)
  transfer        shift_decomposition                (is a transfer failure covariate shift or concept shift?)
  multiplicity    holm, bh
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler

from xmodal_ablation import ALPHAS, PairedSpatial, r2

# --------------------------------------------------------------------------------------
# splitting
# --------------------------------------------------------------------------------------
def spatial_blocks(coords, sample, block_um=250.0):
    """Label each cell with sample:bx_by on a square grid of side `block_um` (per sample).
    Choose block_um >= 2-3x the largest neighbourhood radius so blocks are ~independent."""
    lab = np.empty(len(sample), dtype=object)
    for s in np.unique(sample):
        m = sample == s
        xy = coords[m] - coords[m].min(0)
        b = np.floor(xy / block_um).astype(int)
        lab[m] = [f"{s}:{i}_{j}" for i, j in b]
    return lab.astype(str)


def unit_folds(units):
    """Leave-one-unit-out folds: (name, train_mask, test_mask)."""
    return [(u, units != u, units == u) for u in np.unique(units)]


def _fit_predict(X, y, tr, te):
    sc = StandardScaler().fit(X[tr])
    m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[tr]), y[tr])
    return m.predict(sc.transform(X[te])), m, sc


def fold_scores(d: PairedSpatial, make_X, folds, label="cond", min_cells=30):
    """Held-out R2 and Pearson per (protein, fold). Folds with < min_cells test cells are skipped."""
    rows = []
    for j, p in enumerate(d.proteins):
        X, y = make_X(p), d.protein[:, j]
        for name, tr, te in folds:
            if te.sum() < min_cells:
                continue
            yh, _, _ = _fit_predict(X, y, tr, te)
            pr = np.corrcoef(y[te], yh)[0, 1] if yh.std() > 0 and y[te].std() > 0 else np.nan
            rows.append({"protein": p, "unit": name, "condition": label,
                         "r2": r2(y[te], yh), "pearson": pr, "n": int(te.sum())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------------------
def cluster_bootstrap_ci(values, B=2000, alpha=0.05, seed=0, stat=np.mean):
    """Percentile CI for stat(values) where each value is one independent unit (patient/block).
    For a contrast, pass per-unit DIFFERENCES (paired), not two separate vectors."""
    v = np.asarray(values, float)
    v = v[~np.isnan(v)]
    if len(v) < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.array([stat(v[rng.integers(0, len(v), len(v))]) for _ in range(B)])
    return float(stat(v)), float(np.quantile(boots, alpha / 2)), float(np.quantile(boots, 1 - alpha / 2))


# --------------------------------------------------------------------------------------
# permutation tests
# --------------------------------------------------------------------------------------
def signflip_test(diffs, n_perm=10000, seed=0, alternative="greater"):
    """Paired permutation test over units: under H0 each unit's difference is symmetric around 0,
    so its sign can be flipped. Exact when 2^n_units <= n_perm."""
    d = np.asarray(diffs, float)
    d = d[~np.isnan(d)]
    n = len(d)
    if n == 0:
        return np.nan
    obs = d.mean()
    if 2 ** n <= n_perm:
        signs = np.array(np.meshgrid(*[[-1, 1]] * n)).reshape(n, -1).T
    else:
        signs = np.random.default_rng(seed).choice([-1, 1], size=(n_perm, n))
    null = (signs * d).mean(1)
    if alternative == "greater":
        return float((np.sum(null >= obs - 1e-12)) / len(null))
    if alternative == "less":
        return float((np.sum(null <= obs + 1e-12)) / len(null))
    return float((np.sum(np.abs(null) >= abs(obs) - 1e-12)) / len(null))


def group_permutation_importance(d: PairedSpatial, F, protein, groups: dict, folds, radius=30,
                                 n_rep=20, seed=0):
    """Conditional permutation importance of RNA gene GROUPS for one protein.
    The model is trained once per fold on intact inputs; at test time the group's genes are
    permuted jointly among cells of the same (sample, cell type) - this keeps the cell-type
    signal (which is not the group's to claim) and destroys only within-type information.
    Neighbourhood features are rebuilt from the permuted RNA, so the group's spatial
    contribution is removed too. Returns per (fold, group, rep) drop in R2."""
    from xmodal_ablation import shuffle_within
    rng = np.random.default_rng(seed)
    j = d.proteins.index(protein)
    y = d.protein[:, j]
    X_full = F.build(cognate_for=protein, all_rna=True, nbr_radius=radius)
    rows = []
    rna0 = d.rna.copy()
    for name, tr, te in folds:
        yh, m, sc = _fit_predict(X_full, y, tr, te)
        base = r2(y[te], yh)
        for g, genes in groups.items():
            gi = [d.gidx[x] for x in genes if x in d.gidx]
            for r in range(n_rep):
                perm = shuffle_within(d.sample, d.cell_type, rng)
                d.rna = rna0.copy()
                d.rna[:, gi] = rna0[perm][:, gi]
                Xp = F.build(cognate_for=protein, all_rna=True, nbr_radius=radius)
                rows.append({"protein": protein, "unit": name, "group": g, "rep": r,
                             "importance": base - r2(y[te], m.predict(sc.transform(Xp[te])))})
            d.rna = rna0
    d.rna = rna0
    return pd.DataFrame(rows)


def toroidal_shift_test(x, y, coords, sample, n_perm=999, seed=0):
    """Association (Pearson) between two spatial fields x and y measured on the same cells,
    with a null that preserves each field's own spatial autocorrelation: field y is
    translated by a random (dx, dy) on a torus (wrapping at the sample's bounding box) and
    re-read at the original cell positions by nearest neighbour. A naive cell-level
    permutation would destroy autocorrelation and give far too small p-values."""
    rng = np.random.default_rng(seed)
    obs = np.corrcoef(x, y)[0, 1]
    null = np.empty(n_perm)
    groups = [np.where(sample == s)[0] for s in np.unique(sample)]
    trees = []
    for idx in groups:
        lo, hi = coords[idx].min(0), coords[idx].max(0)
        trees.append((idx, lo, hi - lo + 1e-9))
    for k in range(n_perm):
        y_shift = np.empty_like(y)
        for idx, lo, span in trees:
            shift = rng.uniform(0, 1, 2) * span
            moved = (coords[idx] - lo + shift) % span
            _, nn = cKDTree(moved).query(coords[idx] - lo)
            y_shift[idx] = y[idx][nn]
        null[k] = np.corrcoef(x, y_shift)[0, 1]
    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (n_perm + 1)
    # naive (wrong) permutation p-value, for teaching contrast
    naive = np.array([np.corrcoef(x, rng.permutation(y))[0, 1] for _ in range(199)])
    p_naive = (1 + np.sum(np.abs(naive) >= abs(obs))) / 200
    return {"r": float(obs), "p_toroidal": float(p), "p_naive_cell_perm": float(p_naive),
            "null_sd_toroidal": float(null.std()), "null_sd_naive": float(naive.std())}


# --------------------------------------------------------------------------------------
# conformal prediction
# --------------------------------------------------------------------------------------
def mondrian_conformal(d: PairedSpatial, X, protein, train, calib, test, alpha=0.1):
    """Split conformal with a separate calibration quantile per cell type (Mondrian), giving
    per-cell-type ~(1-alpha) coverage IF calibration and test cells are exchangeable.
    Put calibration and test in DIFFERENT units (samples/blocks) from training. Report the
    empirical coverage on test: a drop below 1-alpha on a held-out region/dataset is a
    distribution-free signal that the RNA->protein mapping did not transfer."""
    j = d.proteins.index(protein)
    y = d.protein[:, j]
    sc = StandardScaler().fit(X[train])
    m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[train]), y[train])
    pred = m.predict(sc.transform(X))
    res = np.abs(y - pred)
    rows = []
    for ct in np.unique(d.cell_type):
        c = calib & (d.cell_type == ct)
        t = test & (d.cell_type == ct)
        if c.sum() < 20 or t.sum() < 10:
            continue
        n = c.sum()
        q = np.quantile(res[c], min(1, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")
        rows.append({"protein": protein, "cell_type": ct, "q_halfwidth": float(q),
                     "coverage_test": float(np.mean(res[t] <= q)), "n_calib": int(n), "n_test": int(t.sum())})
    return pd.DataFrame(rows)


def _density_ratio(Xs, Xt, seed=0):
    """w(x) = p_target(x) / p_source(x) via a probabilistic classifier (source=0, target=1)."""
    Z = np.vstack([Xs, Xt])
    lab = np.r_[np.zeros(len(Xs)), np.ones(len(Xt))]
    sc = StandardScaler().fit(Z)
    clf = LogisticRegression(C=0.1, max_iter=2000, random_state=seed)
    p = cross_val_predict(clf, sc.transform(Z), lab, cv=5, method="predict_proba")[:, 1]
    p = np.clip(p, 1e-3, 1 - 1e-3)
    w = p / (1 - p) * (len(Xs) / len(Xt))
    try:
        from sklearn.metrics import roc_auc_score
        auc = roc_auc_score(lab, p)
    except Exception:
        auc = np.nan
    return w[: len(Xs)], w[len(Xs):], auc


def weighted_conformal_coverage(d, X, protein, train, calib, test, alpha=0.1, shift_features=None):
    """Weighted split conformal (Tibshirani et al., NeurIPS 2019): calibration residuals are
    re-weighted by the estimated density ratio between test and calibration features. If the
    only difference between regions is WHICH cells are present (covariate shift), weighted
    conformal restores coverage; if it does not, P(protein | RNA) itself changed (concept shift)."""
    j = d.proteins.index(protein)
    y = d.protein[:, j]
    sc = StandardScaler().fit(X[train])
    m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[train]), y[train])
    res = np.abs(y - m.predict(sc.transform(X)))
    Z = X if shift_features is None else shift_features   # low-dim context (cell type + neighbourhood composition) is more stable
    wc, wt, auc = _density_ratio(Z[calib], Z[test])
    rc = res[calib]
    order = np.argsort(rc)
    rc, wc = rc[order], wc[order]
    covered = []
    # per test point weighted quantile (vectorised approximation with the mean test weight)
    w_test = np.median(wt)
    p = wc / (wc.sum() + w_test)
    cum = np.cumsum(p)
    k = np.searchsorted(cum, 1 - alpha)
    q_w = rc[min(k, len(rc) - 1)] if k < len(rc) else np.inf
    n = len(rc)
    q_u = np.quantile(res[calib], min(1, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")
    return {"protein": protein, "coverage_unweighted": float(np.mean(res[test] <= q_u)),
            "coverage_weighted": float(np.mean(res[test] <= q_w)),
            "q_unweighted": float(q_u), "q_weighted": float(q_w), "domain_auc": float(auc)}


# --------------------------------------------------------------------------------------
# transferability: covariate shift vs concept shift
# --------------------------------------------------------------------------------------
def shift_decomposition(d, X, protein, source, target, n_inner=5, seed=0, shift_features=None):
    """Decompose the error of a source-trained model on a target domain (region/dataset).

    mse_target_transfer : source model evaluated on target                 (what we observe)
    mse_source_weighted : source model on held-out source cells, importance-weighted to look
                          like the target feature distribution             (what covariate shift alone predicts)
    mse_target_oracle   : model refit on target by inner CV                (best achievable with these inputs)

    covariate-shift part = mse_source_weighted - mse_source_unweighted
    concept-shift part   = mse_target_transfer - mse_source_weighted  (P(protein|RNA) changed)
    target-specific gain = mse_target_transfer - mse_target_oracle   (non-transferable information)
    """
    rng = np.random.default_rng(seed)
    j = d.proteins.index(protein)
    y = d.protein[:, j]
    src_idx = np.where(source)[0]
    rng.shuffle(src_idx)
    half = len(src_idx) // 2
    s_tr, s_va = src_idx[:half], src_idx[half:]
    sc = StandardScaler().fit(X[s_tr])
    m = RidgeCV(alphas=ALPHAS).fit(sc.transform(X[s_tr]), y[s_tr])
    e_va = (y[s_va] - m.predict(sc.transform(X[s_va]))) ** 2
    tgt = np.where(target)[0]
    e_tg = (y[tgt] - m.predict(sc.transform(X[tgt]))) ** 2
    Z = X if shift_features is None else shift_features
    w_va, _, auc = _density_ratio(Z[s_va], Z[tgt], seed)
    w_va = np.clip(w_va, 0, np.quantile(w_va, 0.99))
    mse_src = e_va.mean()
    mse_src_w = np.sum(w_va * e_va) / np.sum(w_va)
    # oracle: refit on target with inner K-fold over contiguous chunks
    folds = np.array_split(tgt[np.argsort(d.coords[tgt, 0])], n_inner)
    e_or = []
    for k in range(n_inner):
        te = folds[k]
        tr = np.concatenate([folds[i] for i in range(n_inner) if i != k])
        sc2 = StandardScaler().fit(X[tr])
        m2 = RidgeCV(alphas=ALPHAS).fit(sc2.transform(X[tr]), y[tr])
        e_or.append((y[te] - m2.predict(sc2.transform(X[te]))) ** 2)
    mse_or = np.concatenate(e_or).mean()
    var_t = y[tgt].var()
    return {"protein": protein, "domain_auc": float(auc), "var_target": float(var_t),
            "mse_source": float(mse_src), "mse_source_weighted": float(mse_src_w),
            "mse_target_transfer": float(e_tg.mean()), "mse_target_oracle": float(mse_or),
            "covariate_shift_part": float(mse_src_w - mse_src),
            "concept_shift_part": float(e_tg.mean() - mse_src_w),
            "target_specific_gain": float(e_tg.mean() - mse_or)}


# --------------------------------------------------------------------------------------
# multiplicity
# --------------------------------------------------------------------------------------
def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    m = len(p)
    run = 0
    for rank, i in enumerate(o):
        run = max(run, (m - rank) * p[i])
        adj[i] = min(1, run)
    return adj


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    adj = np.empty_like(p)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * m / (rank + 1))
        adj[i] = prev
    return adj
