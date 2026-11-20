"""Read/write the harmonised schema (see schema.md) as .h5ad or .npz, and make splits/graphs."""
from __future__ import annotations

import json
import os

import numpy as np

OBS_COLS = ["cell_id", "patient", "sample", "slide", "x_um", "y_um", "cell_type_coarse",
            "cell_type_fine", "label_source", "fold", "spatial_block"]


class Cells(dict):
    """dict of arrays with the schema keys; attribute access for convenience."""
    __getattr__ = dict.__getitem__

    @property
    def coords(self):
        return np.c_[self["x_um"], self["y_um"]].astype(float)

    @property
    def n(self):
        return len(self["cell_id"])


def load(path) -> Cells:
    if path.endswith(".npz"):
        z = np.load(path, allow_pickle=True)
        c = Cells({k: z[k] for k in z.files})
        c["pairing"] = json.loads(str(c.get("pairing", "{}")))
        return c
    import anndata as ad  # cluster only
    A = ad.read_h5ad(path)
    c = Cells({k: A.obs[k].to_numpy() for k in OBS_COLS if k in A.obs})
    for k in ("X_protein", "X_rna"):
        if k in A.obsm:
            c[k] = np.asarray(A.obsm[k])
    for k in ("protein_names", "gene_names", "protein_cognate_gene"):
        if k in A.uns:
            c[k] = np.asarray(A.uns[k])
    c["pairing"] = dict(A.uns.get("pairing", {}))
    return c


def save(c: dict, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if path.endswith(".npz"):
        d = {k: v for k, v in c.items() if k != "pairing"}
        np.savez_compressed(path, pairing=json.dumps(c.get("pairing", {})), **d)
        return
    import anndata as ad
    import pandas as pd
    obs = pd.DataFrame({k: c[k] for k in OBS_COLS if k in c}).set_index("cell_id", drop=False)
    obs.index = obs.index.astype(str)
    A = ad.AnnData(obs=obs)
    for k in ("X_protein", "X_rna"):
        if k in c:
            A.obsm[k] = np.asarray(c[k], dtype=np.float32)
    for k in ("protein_names", "gene_names", "protein_cognate_gene"):
        if k in c:
            A.uns[k] = list(map(str, c[k]))
    A.uns["pairing"] = c.get("pairing", {})
    A.write_h5ad(path, compression="gzip")


def validate(c: Cells):
    errs = []
    for k in ("cell_id", "patient", "x_um", "y_um"):
        if k not in c:
            errs.append(f"missing {k}")
    if "X_protein" not in c and "X_rna" not in c:
        errs.append("no X_protein or X_rna")
    if "X_protein" in c and "protein_names" in c and c["X_protein"].shape[1] != len(c["protein_names"]):
        errs.append("protein_names length mismatch")
    if c.get("pairing", {}).get("type") not in ("same_cell", "adjacent_section", "unpaired"):
        errs.append("pairing.type must be same_cell | adjacent_section | unpaired")
    if len(np.unique(c["patient"])) < 4:
        errs.append("fewer than 4 patients: patient-level UQ will be meaningless (report as descriptive)")
    return errs


def patient_folds(patient, k=5, seed=0):
    pats = np.unique(patient)
    rng = np.random.default_rng(seed)
    order = rng.permutation(pats)
    f = {p: i % k for i, p in enumerate(order)}
    return np.array([f[p] for p in patient])


def spatial_blocks(coords, slide, size_um=250.0):
    bx = np.floor(coords[:, 0] / size_um).astype(int)
    by = np.floor(coords[:, 1] / size_um).astype(int)
    keys = [f"{s}_{x}_{y}" for s, x, y in zip(slide, bx, by)]
    _, inv = np.unique(keys, return_inverse=True)
    return inv


def radius_edges(coords, group, radius=40.0):
    """Edge list (2, E) of a radius graph within each slide/patient group."""
    from sklearn.neighbors import NearestNeighbors
    src, dst = [], []
    for g in np.unique(group):
        idx = np.where(group == g)[0]
        nb = NearestNeighbors(radius=radius).fit(coords[idx]).radius_neighbors(coords[idx], return_distance=False)
        for i, js in enumerate(nb):
            js = js[js != i]
            src.extend(idx[js]); dst.extend([idx[i]] * len(js))
    return np.array([src, dst], dtype=np.int64)
