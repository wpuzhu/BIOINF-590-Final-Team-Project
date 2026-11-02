"""Step 5: figures + results/SUMMARY.md from the aggregate and validation tables."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import FDR, RES

AGG, VAL, FIG = RES / "aggregate", RES / "validation", RES / "figures"
# Reference categorical palette (slots 1-3 validate for all pairs; fixed order = fixed dataset identity)
DS_COLOR = {"schurch_codex": "#2a78d6", "jackson_imc": "#eb6834", "hartmann_mibi": "#1baf7a"}
DS_LABEL = {"schurch_codex": "Schürch CODEX (CRC)", "jackson_imc": "Jackson IMC (breast)",
            "hartmann_mibi": "Hartmann MIBI-TOF (CRC)"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False})


def grouped_barh(ax, table, xlabel, ref=None):
    """table: index = method, columns = dataset."""
    methods, dss = list(table.index), [d for d in DS_COLOR if d in table.columns]
    h = 0.8 / len(dss)
    for i, ds in enumerate(dss):
        y = np.arange(len(methods)) + (i - (len(dss) - 1) / 2) * h
        ax.barh(y, table[ds], height=h - 0.04, color=DS_COLOR[ds], label=DS_LABEL[ds], edgecolor=SURF, lw=1)
    ax.set_yticks(np.arange(len(methods)), methods)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    if ref is not None:
        ax.axvline(ref, color=INK2, lw=1, ls="--")
    ax.legend(frameon=False, loc="lower right", fontsize=8)


def fig_confounding():
    d = pd.read_csv(VAL / "v1b_celltype_variance_vs_spatial_score.csv")
    t = d.groupby(["method", "dataset"]).spearman_with_eta2.mean().unstack()
    t = t.loc[sorted(t.index, key=lambda m: (not m.startswith("global"), m))]
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(t) + 1.2))
    grouped_barh(ax, t, "Spearman ρ: spatial score vs. variance explained by cell type (η²)", ref=0)
    ax.set_title("Global methods track cell-type identity; within-type methods much less",
                 loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "fig1_confounding_eta2.png", dpi=160); plt.close(fig)
    return t


def fig_calibration():
    d = pd.read_csv(VAL / "v4_calibration.csv")
    t = d.pivot(index="method", columns="dataset", values="typeI_005")
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(t) + 1.2))
    grouped_barh(ax, t, "Type-I error at p < 0.05 (coordinates permuted within cell type)", ref=0.05)
    ax.set_title("Calibration under within-type permutation (dashed = nominal 0.05)",
                 loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "fig2_calibration.png", dpi=160); plt.close(fig)
    return d


def fig_agreement():
    a = pd.read_csv(AGG / "method_agreement.csv")
    ms = sorted(set(a.method_a) | set(a.method_b))
    J = pd.DataFrame(np.eye(len(ms)), index=ms, columns=ms)
    for r in a.itertuples():
        J.loc[r.method_a, r.method_b] = J.loc[r.method_b, r.method_a] = r.jaccard
    fig, ax = plt.subplots(figsize=(1.0 * len(ms) + 2.2, 0.8 * len(ms) + 1.4))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ)
    im = ax.imshow(J.to_numpy(), cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(ms)), [m.replace("ct_", "") for m in ms], rotation=35, ha="right")
    ax.set_yticks(range(len(ms)), [m.replace("ct_", "") for m in ms])
    for i in range(len(ms)):
        for j in range(len(ms)):
            v = J.iat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8, color=INK if v < 0.55 else SURF)
    ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.046, label="Jaccard of BH-significant (q < 0.05) pairs")
    ax.set_title("Agreement between ct methods", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "fig3_method_agreement.png", dpi=160); plt.close(fig)
    return J


def fig_funnel():
    f = pd.read_csv(VAL / "validation_funnel.csv").set_index("dataset")
    a = pd.read_csv(VAL / "consensus_with_validation.csv")
    steps = [("Consensus calls", lambda x: np.ones(len(x), bool)),
             ("+ annotation-supported", lambda x: x.annotation_supported),
             ("+ no spillover", lambda x: x.annotation_supported & ~x.spillover_suspect),
             ("+ native to cell type", lambda x: x.annotation_supported & ~x.spillover_suspect & x.expressed_in_type),
             ("+ not background-driven", lambda x: x.annotation_supported & ~x.spillover_suspect
              & x.expressed_in_type & ~x.background_driven),
             ("+ lineage-consistent", lambda x: x.annotation_supported & ~x.spillover_suspect
              & x.expressed_in_type & ~x.background_driven & ~x.canonical_mismatch),
             ("+ replicated = validated", lambda x: x.validated)]
    t = pd.DataFrame({ds: [int(np.asarray(fn(g)).sum()) for _, fn in steps] for ds, g in a.groupby("dataset")},
                     index=[s for s, _ in steps])
    fig, ax = plt.subplots(figsize=(7, 4))
    grouped_barh(ax, t, "Number of (image, cell type, marker) calls remaining")
    ax.set_title("Consensus ctSVPs filtered by annotation-based checks", loc="left", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "fig4_validation_funnel.png", dpi=160); plt.close(fig)
    return t


def md(df, floatfmt="{:.3g}"):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(floatfmt.format(v) if isinstance(v, (float, np.floating)) else str(v)
                                     for v in r) + " |")
    return "\n".join(out)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    conf = fig_confounding(); cal = fig_calibration(); J = fig_agreement(); funnel = fig_funnel()
    pairs = pd.read_csv(AGG / "pairs_observed.csv")
    a = pd.read_csv(VAL / "consensus_with_validation.csv")
    enr = pd.read_csv(VAL / "v2_annotation_enrichment.csv")
    rep = pd.read_csv(VAL / "v5_replication.csv")
    sx = pd.read_csv(VAL / "v6_sparkx_py_vs_R.csv") if (VAL / "v6_sparkx_py_vs_R.csv").stat().st_size > 5 else None
    common = pd.read_csv(AGG / "cross_dataset_common.csv")
    methods = sorted(c[2:] for c in pairs.columns if c.startswith("p_ct_"))

    v = a[a.validated]
    top = (v.groupby(["dataset", "cell_type", "feature"])
           .agg(images=("image", "nunique"), best_acat_q=("acat_q", "min")).reset_index()
           .sort_values(["dataset", "images", "best_acat_q"], ascending=[True, False, True])
           .groupby("dataset").head(12))
    flagged = (a[~a.validated & a.spillover_suspect].groupby(["dataset", "cell_type", "feature", "spillover_source"])
               .size().rename("images").reset_index().sort_values("images", ascending=False).head(12))

    from aggregate import coarse_type, harmonise_marker
    vc = v.assign(coarse_type=v.cell_type.map(coarse_type), marker=v.feature.map(harmonise_marker))
    vcommon = (vc.groupby(["coarse_type", "marker"]).agg(datasets=("dataset", lambda s: ";".join(sorted(set(s)))),
                                                         n_datasets=("dataset", "nunique"),
                                                         images=("image", "nunique"))
               .reset_index().query("n_datasets >= 2").sort_values(["n_datasets", "images"], ascending=False))

    s = [f"# Simple pipeline: results summary\n",
         "Auto-generated by `src/report.py`. Zero-shot run: fixed defaults, no tuning; validation uses the "
         "original dataset annotations, which detection never sees beyond the cell-type labels.\n",
         f"**ct methods aggregated:** {', '.join(methods)}  \n"
         f"**Tested (image, cell type, marker) pairs:** {len(pairs)}  \n"
         f"**Consensus calls** (ACAT q < {FDR} and a strict majority of methods at BH q < {FDR}): "
         f"{int(pairs.consensus.sum())}  \n**Validated candidates:** {int(a.validated.sum())}\n",
         "## 1. Composition confounding (benchmark claim reproduced)\n",
         "Mean Spearman ρ between each marker's spatial score and the share of its variance explained by "
         "annotated cell type (η²). Global SVG methods mostly rediscover cell-type identity; within-type "
         "testing removes most of that.\n",
         md(conf.reset_index().rename(columns=DS_LABEL)), "\n![confounding](figures/fig1_confounding_eta2.png)\n",
         "## 2. Calibration (within-type permutation null)\n",
         md(cal[["dataset", "method", "typeI_005", "ks", "n_bh_calls"]]),
         "\n![calibration](figures/fig2_calibration.png)\n",
         "## 3. What the methods have in common\n",
         "Pairwise Jaccard of significant sets / Spearman of −log10 p:\n",
         md(pd.read_csv(AGG / "method_agreement.csv")), "\n![agreement](figures/fig3_method_agreement.png)\n"]
    if sx is not None and len(sx):
        s += ["Python SPARK-X port vs official R SPARK-X (same pairs; R subsamples cells independently):\n", md(sx), ""]
    s += ["Consensus themes shared by ≥ 2 datasets **before** validation (coarse cell type × marker), top 15:\n",
          md(common.head(15)[["coarse_type", "marker", "n_datasets", "n_images", "datasets"]]),
          "\nThe same, **after** validation:\n",
          md(vcommon.head(15)) if len(vcommon) else "_none_",
          "\n## 4. Validation against original annotations\n",
          "Enrichment of annotation support (authors' cellular-neighbourhood labels; annotated neighbour "
          "composition) among consensus vs. non-consensus pairs:\n", md(enr),
          "\nReplication across images (consensus in ≥ 2 images), observed vs. permuted null:\n", md(rep),
          "\nValidation funnel:\n", md(funnel.reset_index().rename(columns={"index": "step", **DS_LABEL})),
          "\n![funnel](figures/fig4_validation_funnel.png)\n",
          "Most frequent calls removed as likely **spillover** (marker leaking from the named neighbour type):\n",
          md(flagged) if len(flagged) else "_none_",
          "\n## 5. Validated candidate ctSVPs (top 12 per dataset)\n", md(top),
          "\nFull lists: `validation/consensus_with_validation.csv` (every flag per call), "
          "`aggregate/pairs_observed.csv` (every method's p/q per pair).\n"]
    (RES / "SUMMARY.md").write_text("\n".join(s))
    print("wrote", RES / "SUMMARY.md")


if __name__ == "__main__":
    main()
