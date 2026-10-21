"""Step 4: validate aggregated calls against the ORIGINAL dataset annotations
(cell-type labels, the authors' cellular-neighbourhood labels for Schurch CODEX) and against nulls.
These annotations were never used for detection beyond defining which cells form a "cell type".

V1  Lineage confounding (the benchmark's central warning): are top global SVPs / ct calls just the
    markers that DEFINE cell types? Lineage markers per type come from the annotated labels.
V1b Variance explained by annotated cell type (eta^2) vs spatial score, global vs within-type.
V2  Annotation support: within cell type k, is marker m associated with (a) the authors' cellular
    neighbourhood label (Schurch) or (b) the annotated cell-type composition of each cell's kNN neighbourhood?
V3  Spillover audit: a call (k, m) is flagged as a likely segmentation-spillover artefact (not a
    cell-intrinsic biomarker) if some annotated neighbour type j expresses m >= LINEAGE_DELTA SD more than
    type k does AND j's local abundance is significantly positively correlated (Spearman p < 0.01) with m
    inside type k, e.g. CD3 in tumour cells next to T cells.
V3b Expression plausibility: a cell-intrinsic biomarker should be "native" to the cell type, i.e. the
    type's mean is within NATIVE_GAP SD of the best-expressing annotated type. Spatial variation of a
    marker the type does not express (e.g. FOXP3 in tumour cells) points to background / imaging-tile
    gradients or spillover.
V3c Shared background: within type k, regress m on the cell's mean signal over all other markers (a
    per-cell brightness / regional-background score, the protein analogue of library-size normalisation)
    and re-test the residual with Moran's I. If spatial structure disappears (BH q >= FDR), the call is
    "background-driven" (e.g. imaging-tile or tissue-region intensity gradients shared by all markers).
V3d Canonical-lineage mismatch (Schurch only): the authors annotate every antibody channel with the
    population it marks ("CD5 - T cells", "LAG-3 - checkpoint"). A call whose marker is documented for a
    different lineage than the tested cell type (e.g. T-bet in tumour cells) is flagged. Functional
    channels (checkpoint, proliferation, metabolism, ...) are allowed in any cell type.
V4  Calibration: type-I error / KS on coordinates permuted within cell type (all calls there are false).
V5  Replication across images of the same dataset, vs. the same statistic on the permuted null.
V6  Python vs. official R SPARK-X concordance.
Outputs: results/validation/*.csv
"""
import re

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors

from aggregate import AGG, bh, load
from methods_py import morans_i
from config import FDR, KNN, LINEAGE_DELTA, NATIVE_GAP, MAX_LINEAGE_PER_TYPE, MIN_CELLS_PER_TYPE, PROC, RES

VAL = RES / "validation"


def load_images(dataset):
    for meta_f in sorted((PROC / dataset).glob("*__meta.csv")):
        image = meta_f.name.replace("__meta.csv", "")
        yield image, pd.read_csv(meta_f), pd.read_csv(meta_f.with_name(f"{image}__expr.csv"), index_col=0)


def lineage_markers(dataset):
    """Annotation-derived lineage markers. Each marker gets at most one owner: the cell type with the
    highest mean (z-scored over the dataset), and only if that mean beats the rest of the cells by
    >= LINEAGE_DELTA. Each type keeps its strongest MAX_LINEAGE_PER_TYPE markers."""
    metas, exprs = zip(*[(m, e) for _, m, e in load_images(dataset)])
    meta, expr = pd.concat(metas), pd.concat(exprs)
    Z = (expr - expr.mean()) / expr.std().replace(0, 1)
    counts = meta.cell_type.value_counts()
    types = counts[counts >= MIN_CELLS_PER_TYPE].index
    means = Z.groupby(meta.cell_type.to_numpy()).mean().loc[types]
    out = {t: [] for t in meta.cell_type.unique()}
    for m in Z.columns:
        owner = means[m].idxmax()
        delta = means.at[owner, m] - Z.loc[(meta.cell_type != owner).to_numpy(), m].mean()
        if delta >= LINEAGE_DELTA:
            out[owner].append((delta, m))
    return {t: [m for _, m in sorted(v, reverse=True)[:MAX_LINEAGE_PER_TYPE]] for t, v in out.items()}


def type_means(dataset):
    """Dataset-level mean z-scored intensity of every marker in every annotated cell type."""
    metas, exprs = zip(*[(m, e) for _, m, e in load_images(dataset)])
    meta, expr = pd.concat(metas), pd.concat(exprs)
    Z = (expr - expr.mean()) / expr.std().replace(0, 1)
    counts = meta.cell_type.value_counts()
    M = Z.groupby(meta.cell_type.to_numpy()).mean()
    M.attrs["top"] = M.loc[counts[counts >= MIN_CELLS_PER_TYPE].index].max()   # best-expressing common type
    return M


# Authors' channel descriptions (Schurch CSV headers) -> lineage they mark; function terms are lineage-free
LINEAGE_TERMS = {"t cell": "T", "th1": "T", "th2": "T", "regulatory t": "T", "cytotoxic t": "T", "naive t": "T",
                 "t helper": "T", "memory": "T", "b cell": "B", "plasma": "B", "macrophage": "myeloid",
                 "monocyt": "myeloid", "dcs": "myeloid", "granulocyte": "granulocyte", "nk": "NK",
                 "epithelia": "epithelial", "vascul": "endothelial", "lymphatic": "endothelial",
                 "smooth muscle": "stromal", "stroma": "stromal", "nerve": "nerve", "neuroendocrine": "epithelial",
                 "hematopoietic": "immune"}
TYPE_LINEAGE = [("tumor", "epithelial"), ("treg|t cell", "T"), ("b cell|plasma", "B"),
                ("macroph|monocyt|dc", "myeloid"), ("granulocyte", "granulocyte"), ("nk", "NK"),
                ("vascul|lymphatic", "endothelial"), ("stroma|smooth|adipocyte", "stromal"), ("nerve", "nerve")]


def canonical_mismatch(dataset):
    """(cell_type, marker) -> True if the authors' channel description names a different lineage."""
    f = PROC / dataset / "marker_annotations.csv"
    if not f.exists():
        return None
    ann = pd.read_csv(f)
    mk = {r.marker: next((lab for k, lab in LINEAGE_TERMS.items() if k in r.description.lower()), None)
          for r in ann.itertuples()}
    def check(t, m):
        lin = mk.get(m)
        tl = next((lab for pat, lab in TYPE_LINEAGE if re.search(pat, t.lower())), None)
        if lin is None or tl is None:
            return False
        if lin == "immune":
            return tl in ("epithelial", "endothelial", "stromal", "nerve")
        return lin != tl
    return check


def v1b_eta2(raw):
    """Fraction of each marker's variance explained by annotated cell type (eta^2), per image, vs the
    spatial score of global methods and of within-type methods. Composition confounding predicts a strong
    positive correlation for global methods and a much weaker one within cell types."""
    rows = []
    for ds in raw.dataset.unique():
        for image, meta, expr in load_images(ds):
            ct = meta.cell_type.to_numpy()
            eta2 = {m: 1 - expr[m].groupby(ct).transform("mean").sub(expr[m]).pow(2).sum()
                    / max(((expr[m] - expr[m].mean()) ** 2).sum(), 1e-12) for m in expr.columns}
            r = raw[(raw.dataset == ds) & (raw.image == image)]
            for method, g in r.groupby("method"):
                score = (-np.log10(g.p_value.clip(lower=1e-300))).groupby(g.feature).mean()
                e = pd.Series(eta2).reindex(score.index)
                rows.append(dict(dataset=ds, image=image, method=method,
                                 spearman_with_eta2=stats.spearmanr(score, e, nan_policy="omit")[0]))
    return pd.DataFrame(rows)


def v1_lineage(raw, pairs, L):
    rows = []
    for ds, Ld in L.items():
        all_lin = set().union(*Ld.values())
        g = raw[(raw.dataset == ds) & raw.method.str.startswith("global_")]
        base = g.feature.isin(all_lin).mean()
        for m, gm in g.groupby("method"):
            top = gm.sort_values("p_value").groupby("image").head(10)
            rows.append(dict(dataset=ds, analysis=m, set="top-10 per image", n=len(top),
                             frac_lineage_any=top.feature.isin(all_lin).mean(), base_rate=base))
        p = pairs[pairs.dataset == ds]
        cls = p.apply(lambda r: "own-lineage" if r.feature in Ld.get(r.cell_type, [])
                      else ("foreign-lineage" if r.feature in all_lin else "non-lineage"), axis=1)
        for name, mask in [("all tested pairs", np.ones(len(p), bool)), ("consensus ctSVPs", p.consensus.to_numpy())]:
            c = cls[mask].value_counts(normalize=True)
            rows.append(dict(dataset=ds, analysis="ct", set=name, n=int(mask.sum()),
                             frac_lineage_any=c.get("own-lineage", 0) + c.get("foreign-lineage", 0),
                             frac_own=c.get("own-lineage", 0), frac_foreign=c.get("foreign-lineage", 0),
                             base_rate=base))
    return pd.DataFrame(rows)


def v2_v3_annotation(pairs, L):
    """Per tested pair: niche-label test (Schurch) and kNN-composition regression; driver neighbour type."""
    rows = []
    for ds in pairs.dataset.unique():
        M = type_means(ds)
        mismatch = canonical_mismatch(ds)
        for image, meta, expr in load_images(ds):
            coords = meta[["x", "y"]].to_numpy(float)
            types = meta.cell_type.to_numpy()
            levels = np.unique(types)
            _, ind = NearestNeighbors(n_neighbors=KNN + 1).fit(coords).kneighbors(coords)
            comp = (types[ind[:, 1:], None] == levels[None, None, :]).mean(1)       # cells x types
            Zimg = (expr - expr.mean()) / expr.std().replace(0, 1)
            tested = pairs[(pairs.dataset == ds) & (pairs.image == image)]
            for t, tp in tested.groupby("cell_type"):
                idx = np.flatnonzero(types == t)
                C = comp[idx][:, levels != t]
                keepc = C.std(0) > 0
                C, cnames = C[:, keepc], levels[levels != t][keepc]
                Xd = np.column_stack([np.ones(len(idx)), C[:, :-1]]) if C.shape[1] > 1 else np.ones((len(idx), 1))
                niche = meta.niche.astype(str).to_numpy()[idx] if "niche" in meta else None
                for _, r in tp.iterrows():
                    y = expr[r.feature].to_numpy()[idx]
                    rec = dict(dataset=ds, image=image, cell_type=t, feature=r.feature,
                               consensus=r.consensus, acat_q=r.acat_q,
                               type_mean_z=M.at[t, r.feature],
                               gap_to_top_type=M.attrs["top"][r.feature] - M.at[t, r.feature],
                               expressed_in_type=M.attrs["top"][r.feature] - M.at[t, r.feature] < NATIVE_GAP,
                               canonical_mismatch=mismatch(t, r.feature) if mismatch else False)
                    # (b) composition regression F-test + driver = neighbour type most correlated with y
                    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
                    rss1, rss0 = ((y - Xd @ beta) ** 2).sum(), ((y - y.mean()) ** 2).sum()
                    k1, n = Xd.shape[1] - 1, len(y)
                    if k1 > 0 and rss1 > 0 and n - k1 - 1 > 0:
                        F = ((rss0 - rss1) / k1) / (rss1 / (n - k1 - 1))
                        rec["p_composition"] = stats.f.sf(F, k1, n - k1 - 1)
                        cor = [stats.spearmanr(C[:, j], y) if C[:, j].std() > 0 else (0.0, 1.0)
                               for j in range(C.shape[1])]
                        rho = np.array([c[0] for c in cor]); prho = np.array([c[1] for c in cor])
                        j = int(np.nanargmax(rho))
                        rec["driver_type"], rec["driver_rho"] = cnames[j], rho[j]
                        # spillover: some neighbour type expressing m >= LINEAGE_DELTA SD more than type t
                        # is significantly POSITIVELY associated with m inside type t
                        excess = np.array([M.at[c, r.feature] - M.at[t, r.feature] for c in cnames])
                        hit = (excess >= LINEAGE_DELTA) & (rho > 0) & (prho < 0.01)
                        rec["spillover_suspect"] = bool(hit.any())
                        if hit.any():
                            jj = int(np.argmax(np.where(hit, rho, -np.inf)))
                            rec["spillover_source"], rec["spillover_rho"] = cnames[jj], rho[jj]
                    # V3c shared background: remove the cell's mean signal over all OTHER markers, re-test
                    if r.consensus:
                        bg = Zimg.drop(columns=r.feature).to_numpy()[idx].mean(1)
                        Xb = np.column_stack([np.ones(len(idx)), bg])
                        res = y - Xb @ np.linalg.lstsq(Xb, y, rcond=None)[0]
                        rec["p_after_background"] = morans_i(res[:, None], coords[idx])[1][0]
                    # (a) authors' cellular-neighbourhood labels (Schurch only)
                    if niche is not None:
                        grp = [y[niche == c] for c in np.unique(niche) if c != "nan" and (niche == c).sum() >= 10]
                        if len(grp) >= 2:
                            rec["p_niche_label"] = stats.kruskal(*grp).pvalue
                    rows.append(rec)
    out = pd.DataFrame(rows)
    for col in ["p_composition", "p_niche_label", "p_after_background"]:
        if col in out:
            out[col.replace("p_", "q_")] = out.groupby("dataset")[col].transform(bh)
    out["spillover_suspect"] = out["spillover_suspect"].astype("boolean").fillna(False).astype(bool)
    return out


def enrichment(ann):
    rows = []
    for ds, a in ann.groupby("dataset"):
        for col in ["q_niche_label", "q_composition"]:
            if col not in a or a[col].isna().all():
                continue
            b = a.dropna(subset=[col])
            sup = b[col] < FDR
            tab = pd.crosstab(b.consensus, sup).reindex(index=[False, True], columns=[False, True], fill_value=0)
            orr, p = stats.fisher_exact(tab.to_numpy())
            rows.append(dict(dataset=ds, annotation=col.replace("q_", ""), n_pairs=len(b),
                             support_rate_consensus=sup[b.consensus].mean(),
                             support_rate_non_consensus=sup[~b.consensus].mean(), odds_ratio=orr, fisher_p=p))
    return pd.DataFrame(rows)


def v4_calibration(perm):
    ct = perm[perm.method.str.startswith("ct_")]
    return (ct.groupby(["dataset", "method"])
            .apply(lambda g: pd.Series(dict(n=len(g), typeI_005=(g.p_value < 0.05).mean(),
                                            ks=stats.kstest(g.p_value.dropna(), "uniform").statistic,
                                            n_bh_calls=(g.q_value < FDR).sum())), include_groups=False)
            .reset_index())


def v5_replication(pairs, pairs_perm):
    rows = []
    for name, p in [("observed", pairs), ("permuted null", pairs_perm)]:
        r = p.groupby(["dataset", "cell_type", "feature"]).agg(n_img=("image", "nunique"), n_cons=("consensus", "sum"))
        r = r[r.n_img >= 2]
        for ds, g in r.groupby(level=0):
            called = g[g.n_cons >= 1]
            rows.append(dict(dataset=ds, data=name, pairs_seen_in_2plus_images=len(g),
                             called_in_any=len(called), replicated_2plus=(called.n_cons >= 2).sum(),
                             replication_rate=(called.n_cons >= 2).mean() if len(called) else np.nan))
    return pd.DataFrame(rows)


def v6_sparkx(pairs):
    if "p_ct_sparkx_R" not in pairs:
        return pd.DataFrame()
    b = pairs[["dataset", "p_ct_sparkx_py", "p_ct_sparkx_R"]].dropna()
    return (b.groupby("dataset").apply(lambda g: pd.Series(dict(
        n=len(g), spearman=stats.spearmanr(g.p_ct_sparkx_py, g.p_ct_sparkx_R)[0],
        agree_at_005=((g.p_ct_sparkx_py < .05) == (g.p_ct_sparkx_R < .05)).mean())), include_groups=False)
        .reset_index())


def main():
    VAL.mkdir(parents=True, exist_ok=True)
    raw = load("observed")
    perm = load("permuted")
    pairs = pd.read_csv(AGG / "pairs_observed.csv")
    pairs_perm = pd.read_csv(AGG / "pairs_permuted.csv")
    L = {ds: lineage_markers(ds) for ds in pairs.dataset.unique()}
    pd.DataFrame([(ds, t, ";".join(m)) for ds, d in L.items() for t, m in d.items()],
                 columns=["dataset", "cell_type", "lineage_markers"]).to_csv(VAL / "lineage_markers.csv", index=False)

    v1 = v1_lineage(raw, pairs, L); v1.to_csv(VAL / "v1_lineage_confounding.csv", index=False)
    v1b = v1b_eta2(raw); v1b.to_csv(VAL / "v1b_celltype_variance_vs_spatial_score.csv", index=False)
    ann = v2_v3_annotation(pairs, L); ann.to_csv(VAL / "v2_annotation_per_pair.csv", index=False)
    enr = enrichment(ann); enr.to_csv(VAL / "v2_annotation_enrichment.csv", index=False)
    v4 = v4_calibration(perm); v4.to_csv(VAL / "v4_calibration.csv", index=False)
    v5 = v5_replication(pairs, pairs_perm); v5.to_csv(VAL / "v5_replication.csv", index=False)
    v6 = v6_sparkx(pairs); v6.to_csv(VAL / "v6_sparkx_py_vs_R.csv", index=False)

    # Final validated list: consensus + annotation-supported + not spillover-suspect + replicated (if >1 image)
    rep = pd.read_csv(AGG / "replication.csv")
    a = ann[ann.consensus].copy()
    a["annotation_supported"] = (a.get("q_niche_label", pd.Series(np.nan, index=a.index)) < FDR) | (a.q_composition < FDR)
    a = a.merge(rep, on=["dataset", "cell_type", "feature"], how="left")
    a["replicated"] = (a.n_images == 1) | (a.n_consensus >= 2)
    a["background_driven"] = ~(a.q_after_background < FDR)
    a["validated"] = (a.annotation_supported & ~a.spillover_suspect & a.expressed_in_type
                      & ~a.background_driven & ~a.canonical_mismatch & a.replicated)
    a.to_csv(VAL / "consensus_with_validation.csv", index=False)
    summary = (a.groupby("dataset").agg(consensus=("feature", "size"),
                                        annotation_supported=("annotation_supported", "sum"),
                                        spillover_suspect=("spillover_suspect", "sum"),
                                        not_expressed_in_type=("expressed_in_type", lambda x: (~x.astype(bool)).sum()),
                                        background_driven=("background_driven", "sum"),
                                        canonical_mismatch=("canonical_mismatch", "sum"),
                                        replicated=("replicated", "sum"), validated=("validated", "sum"))
               .reset_index())
    summary.to_csv(VAL / "validation_funnel.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
