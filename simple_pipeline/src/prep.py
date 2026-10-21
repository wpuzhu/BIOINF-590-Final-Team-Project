"""Step 1: load the three annotated datasets and write one standard table pair per image:
    data/processed/<dataset>/<image>__expr.csv   cells x markers (transformed intensities)
    data/processed/<dataset>/<image>__meta.csv   cell_id, x, y, cell_type, [niche], [group]
The detection methods only ever see expr + x/y + cell_type.
`niche` and `group` are original-paper annotations, kept for validation only.
"""
import re

import anndata as ad
import numpy as np
import pandas as pd

from config import DATASETS, PROC, RAW, SEED


def dense(X):
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def transform_raw(X: pd.DataFrame) -> pd.DataFrame:
    """Same rule for every raw-intensity dataset (no tuning): arcsinh with a per-marker cofactor
    (20% of the marker's 99th percentile), then clip at the 99.9th percentile."""
    q99 = X.quantile(0.99).replace(0, np.nan).fillna(X.max()).clip(lower=1e-6)
    Y = np.arcsinh(X / (0.2 * q99))
    return Y.clip(upper=Y.quantile(0.999), axis=1)


def write(dataset: str, image: str, expr: pd.DataFrame, meta: pd.DataFrame) -> None:
    out = PROC / dataset
    out.mkdir(parents=True, exist_ok=True)
    meta = meta.reset_index(drop=True)
    expr = expr.reset_index(drop=True)
    meta.insert(0, "cell_id", [f"{image}_c{i}" for i in range(len(meta))])
    expr.index = meta["cell_id"]
    expr.to_csv(out / f"{image}__expr.csv")
    meta.to_csv(out / f"{image}__meta.csv", index=False)
    print(f"  {dataset}/{image}: {len(meta)} cells, {expr.shape[1]} markers, "
          f"{meta.cell_type.nunique()} cell types")


def prep_schurch(cfg):
    d = pd.read_csv(RAW / cfg["file"])
    d = d[~d["ClusterName"].isin(cfg["drop_types"])]
    marker_cols = [c for c in d.columns if " - " in c]                 # 56 antibody channels
    names = {c: c.split(" - ")[0] for c in marker_cols}
    # Small batch: per clinical group (1 = CLR, 2 = DII), the spots closest to the median size
    sizes = d.groupby(["groups", "spots"]).size().rename("n").reset_index()
    picked = []
    for g, s in sizes.groupby("groups"):
        s = s.assign(dev=(s.n - s.n.median()).abs()).sort_values(["dev", "spots"])
        picked += s.spots.head(cfg["n_images_per_group"]).tolist()
    expr_all = transform_raw(d[marker_cols].rename(columns=names))
    (PROC / "schurch_codex").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"marker": [names[c] for c in marker_cols],                 # authors' channel annotation
                  "description": [c.split(" - ")[1].split(":")[0] for c in marker_cols]}
                 ).to_csv(PROC / "schurch_codex" / "marker_annotations.csv", index=False)
    for spot in sorted(picked):
        m = d["spots"] == spot
        meta = pd.DataFrame({
            "x": d.loc[m, "X:X"].to_numpy(), "y": d.loc[m, "Y:Y"].to_numpy(),
            "cell_type": d.loc[m, "ClusterName"].to_numpy(),
            "niche": d.loc[m, "neighborhood name"].to_numpy(),            # authors' cellular neighbourhoods
            "group": np.where(d.loc[m, "groups"].to_numpy() == 1, "CLR", "DII"),
            "patient": d.loc[m, "patients"].to_numpy(),
        })
        write("schurch_codex", f"spot{spot}", expr_all.loc[m], meta)


def imc_name(v: str) -> str:
    metal = re.search(r"([A-Z][a-z]?\d{2,3})Di", v)
    token = v.split()[-1]
    return f"{token}_{metal.group(1)}" if metal else v.replace(" ", "_")


def prep_jackson(cfg):
    a = ad.read_h5ad(RAW / cfg["file"])
    X = pd.DataFrame(dense(a.X), columns=[imc_name(v) for v in a.var_names])
    X = X.drop(columns=[c for c in X.columns if c == "Area" or c.startswith("Histone")])  # size / DNA channels
    meta = pd.DataFrame({"x": a.obsm["spatial"][:, 0], "y": a.obsm["spatial"][:, 1],
                         "cell_type": a.obs["cell type"].astype(str).to_numpy()})
    write("jackson_imc", "image1", transform_raw(X), meta)


def prep_hartmann(cfg):
    a = ad.read_h5ad(RAW / cfg["file"])
    X = pd.DataFrame(dense(a.X), columns=list(a.var_names))      # already scaled by the authors
    X = X.drop(columns=["H3"])                                          # nuclear channel
    for lib in sorted(a.obs["library_id"].unique()):
        m = (a.obs["library_id"] == lib).to_numpy()
        meta = pd.DataFrame({"x": a.obsm["spatial"][m, 0], "y": a.obsm["spatial"][m, 1],
                             "cell_type": a.obs.loc[m, "Cluster"].astype(str).to_numpy(),
                             "patient": a.obs.loc[m, "donor"].astype(str).to_numpy()})
        write("hartmann_mibi", str(lib), X.loc[m], meta)


if __name__ == "__main__":
    np.random.seed(SEED)
    prep_schurch(DATASETS["schurch_codex"])
    prep_jackson(DATASETS["jackson_imc"])
    prep_hartmann(DATASETS["hartmann_mibi"])
