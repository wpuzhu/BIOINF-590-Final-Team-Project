"""
Model-level uncertainty quantification for Stage 1 and Stage 2 comparisons.

The central quantity in both stages is a MARGINAL IMPROVEMENT: the paired difference in a
held-out metric between two arms run on identical splits (same patients, same labelled cells,
same seeds), e.g.  Stage 1: pretrained - random-init;  Stage 2: RNA-taught student - protein-only.

Functions (numpy only, no statsmodels):
  paired_table            long results -> wide table of paired differences per (patient, seeds)
  hierarchical_bootstrap  CI of the mean improvement: resample patients, then seeds within patients
  signflip_patients       exact/MC sign-flip test on patient-mean differences
  tost_equivalence        two one-sided tests: is |improvement| < margin? (claims "no gain")
  variance_components     ANOVA (method-of-moments) split of variance into patient / seed / residual
  aps_sets, conformal_summary   adaptive prediction sets (Romano et al. 2020), Mondrian by cell type
  ece                     expected calibration error
  holm, bh                multiplicity
Each function documents the unit of resampling. Never resample cells.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------------------
def paired_table(df, arm_col, metric, a1, a0, keys=("patient", "seed")):
    """Wide table of metric[a1] - metric[a0] for every combination of `keys` present in both."""
    keys = list(keys)
    w = df.pivot_table(index=keys, columns=arm_col, values=metric)
    d = (w[a1] - w[a0]).dropna().rename("diff").reset_index()
    return d


def hierarchical_bootstrap(d, value="diff", outer="patient", inner=("seed",), B=10000, alpha=0.05, seed=0):
    """Two-stage cluster bootstrap of the mean paired difference.

    Stage 1 resamples OUTER units (patients) with replacement; stage 2 resamples rows (seed
    replicates) within each drawn patient. This propagates both biological (between-patient)
    and training (between-seed) variability. With < 6 patients percentile CIs under-cover:
    report n and the sign-flip p alongside.
    """
    rng = np.random.default_rng(seed)
    groups = [g[value].to_numpy() for _, g in d.groupby(outer)]
    n = len(groups)
    est = float(np.mean([g.mean() for g in groups]))
    boots = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, n, n)
        boots[b] = np.mean([groups[i][rng.integers(0, len(groups[i]), len(groups[i]))].mean() for i in pick])
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return {"estimate": est, "ci_lo": float(lo), "ci_hi": float(hi), "n_patients": n,
            "n_rows": int(len(d)), "boot_sd": float(boots.std())}


def signflip_patients(d, value="diff", outer="patient", alternative="greater", n_perm=20000, seed=0):
    """Sign-flip permutation test on patient-mean differences (H0: symmetric around 0).
    Exact when 2^n_patients <= n_perm. Smallest attainable one-sided p = 2^-n_patients."""
    m = d.groupby(outer)[value].mean().to_numpy()
    n = len(m)
    obs = m.mean()
    if 2 ** n <= n_perm:
        signs = np.array(list(itertools.product([-1, 1], repeat=n)))
    else:
        signs = np.random.default_rng(seed).choice([-1, 1], size=(n_perm, n))
    null = (signs * m).mean(1)
    if alternative == "greater":
        p = np.mean(null >= obs - 1e-12)
    elif alternative == "less":
        p = np.mean(null <= obs + 1e-12)
    else:
        p = np.mean(np.abs(null) >= abs(obs) - 1e-12)
    return float(p)


def tost_equivalence(d, margin, value="diff", outer="patient", alpha=0.05, B=10000, seed=0):
    """Equivalence by the CI-inclusion form of TOST: the improvement is declared practically
    zero if the (1 - 2*alpha) hierarchical-bootstrap CI lies entirely inside (-margin, margin).
    Pre-register `margin` (e.g. 0.02 macro-F1)."""
    r = hierarchical_bootstrap(d, value, outer, B=B, alpha=2 * alpha, seed=seed)
    return {**r, "margin": margin, "equivalent": bool(r["ci_lo"] > -margin and r["ci_hi"] < margin)}


def variance_components(d, value="diff", factors=("patient", "seed")):
    """Additive random-effects ANOVA on a balanced table: y = mu + a_patient + b_seed + e.
    sigma2_f = max(0, (MS_f - MS_resid) / (N / n_levels_f)). Reports the share of variance
    due to biology (patients) vs training randomness (seeds) vs residual (interaction)."""
    y = d[value].to_numpy()
    N, mu = len(y), y.mean()
    ss_tot = ((y - mu) ** 2).sum()
    out, ss_used, df_used = {}, 0.0, 0
    ms = {}
    factors = [f for f in factors if d[f].nunique() > 1]       # a factor with one level has no variance
    for f in factors:
        lv = d.groupby(f)[value]
        k = lv.ngroups
        ss = (lv.count() * (lv.mean() - mu) ** 2).sum()
        ms[f] = (ss / (k - 1), k)
        ss_used += ss
        df_used += k - 1
    df_res = N - 1 - df_used
    ms_res = (ss_tot - ss_used) / df_res if df_res > 0 else np.nan
    for f, (m_f, k) in ms.items():
        out[f"var_{f}"] = max(0.0, (m_f - ms_res) / (N / k)) if np.isfinite(ms_res) else np.nan
    out["var_residual"] = ms_res
    tot = sum(v for v in out.values() if np.isfinite(v))
    for kname in list(out):
        out[kname.replace("var_", "share_")] = out[kname] / tot if tot > 0 else np.nan
    return out


# --------------------------------------------------------------------------------------
def aps_scores(probs, y):
    """Adaptive prediction set score (Romano, Sesia & Candes 2020), randomised version."""
    rng = np.random.default_rng(0)
    order = np.argsort(-probs, 1)
    sorted_p = np.take_along_axis(probs, order, 1)
    cum = np.cumsum(sorted_p, 1)
    rank = np.argmax(order == y[:, None], 1)
    u = rng.random(len(y))
    return cum[np.arange(len(y)), rank] - u * sorted_p[np.arange(len(y)), rank]


def aps_sets(probs, qhat):
    order = np.argsort(-probs, 1)
    sorted_p = np.take_along_axis(probs, order, 1)
    cum = np.cumsum(sorted_p, 1)
    keep = cum - sorted_p < qhat           # include classes until the cumulative mass passes qhat
    sets = np.zeros_like(probs, bool)
    np.put_along_axis(sets, order, keep, 1)
    return sets


def conformal_summary(probs_cal, y_cal, probs_test, y_test, alpha=0.1, groups_cal=None, groups_test=None):
    """Split-conformal APS sets. Calibration cells must come from patients disjoint from the test
    patients AND from training (exchangeability holds between patients, not between cells).
    With groups (cell types): Mondrian calibration per group (class-conditional coverage).
    Returns coverage and mean set size overall and per group: at equal coverage, a SMALLER set
    means the model is genuinely more certain."""
    def qhat(scores):
        n = len(scores)
        return np.quantile(scores, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")
    rows = []
    if groups_cal is None:
        q = qhat(aps_scores(probs_cal, y_cal))
        S = aps_sets(probs_test, q)
        cov = S[np.arange(len(y_test)), y_test]
        return pd.DataFrame([{"group": "all", "coverage": cov.mean(), "set_size": S.sum(1).mean(), "n": len(y_test)}])
    for g in np.unique(groups_test):
        c, t = groups_cal == g, groups_test == g
        if c.sum() < 20 or t.sum() < 5:
            continue
        q = qhat(aps_scores(probs_cal[c], y_cal[c]))
        S = aps_sets(probs_test[t], q)
        cov = S[np.arange(t.sum()), y_test[t]]
        rows.append({"group": g, "coverage": cov.mean(), "set_size": S.sum(1).mean(), "n": int(t.sum())})
    return pd.DataFrame(rows)


def ece(probs, y, n_bins=15):
    conf = probs.max(1)
    pred = probs.argmax(1)
    bins = np.linspace(0, 1, n_bins + 1)
    e = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (conf > lo) & (conf <= hi)
        if m.any():
            e += m.mean() * abs((pred[m] == y[m]).mean() - conf[m].mean())
    return float(e)


# --------------------------------------------------------------------------------------
def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj, run, m = np.empty_like(p), 0.0, len(p)
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i])
        adj[i] = min(1.0, run)
    return adj


def bh(p):
    p = np.asarray(p, float)
    m, o = len(p), np.argsort(p)
    adj, prev = np.empty_like(p), 1.0
    for r in range(m - 1, -1, -1):
        prev = min(prev, p[o[r]] * m / (r + 1))
        adj[o[r]] = prev
    return adj


def improvement_report(df, arm_col, metric, comparisons, margin=0.02, alternative="greater",
                       outer="patient", inner=("seed",), strata=()):
    """One row per (stratum, comparison): estimate, 95% hierarchical-bootstrap CI, sign-flip p,
    TOST verdict, variance shares; Holm within each stratum."""
    rows = []
    strata = list(strata)
    groups = df.groupby(strata) if strata else [((), df)]
    for key, g in groups:
        key = key if isinstance(key, tuple) else (key,)
        for name, (a1, a0) in comparisons.items():
            if a1 not in g[arm_col].values or a0 not in g[arm_col].values:
                continue
            d = paired_table(g, arm_col, metric, a1, a0, keys=(outer, *inner))
            if d.empty:
                continue
            r = hierarchical_bootstrap(d, outer=outer)
            r["p_signflip"] = signflip_patients(d, outer=outer, alternative=alternative)
            r["tost_equivalent"] = tost_equivalence(d, margin, outer=outer)["equivalent"]
            vc = variance_components(d, factors=(outer, *inner))
            r.update({k: v for k, v in vc.items() if k.startswith("share_")})
            rows.append({**dict(zip(strata, key)), "comparison": name, **r})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["p_holm"] = np.nan
        grp = out.groupby(strata).groups if strata else {(): out.index}
        for _, idx in grp.items():
            out.loc[idx, "p_holm"] = holm(out.loc[idx, "p_signflip"].values)
    return out
