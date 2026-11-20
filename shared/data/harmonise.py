"""
Convert raw downloads into the harmonised schema (schema.md).

  python harmonise.py --registry ../../datasets/registry.yaml --dataset xenium_protein_rcc \
                      --raw /scratch/$ACCOUNT/$USER/raw --out /scratch/$ACCOUNT/$USER/harmonised

Loaders:
  xenium      10x Xenium outs/ (cell_feature_matrix.h5 + cells.parquet). Handles Xenium Protein,
              where the matrix holds both 'Gene Expression' and 'Protein Expression' features.
  cosmx       Bruker/NanoString flat files (*exprMat_file.csv + *metadata_file.csv); protein columns
              are matched by the regex given in the registry (release-specific; CHECK on download).
  table       one cell table (CSV/parquet) with columns named in the registry (CODEX/IMC/MIBI
              releases such as Schurch CRC, Jackson/Danenberg IMC, Keren MIBI, SPATCH CODEX).
Labels are attached from a separate table when the registry gives `labels:`.

Column names in vendor exports change between releases. Every loader prints what it found and
fails loudly if a required column is missing: check the printout before the first real run.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spt_io import patient_folds, save, spatial_blocks, validate  # noqa: E402

COGNATE = {  # marker -> encoding gene, extend per panel (datasets/marker_gene_map.csv overrides)
    "CD3": "CD3E", "CD3e": "CD3E", "CD4": "CD4", "CD8": "CD8A", "CD8a": "CD8A", "CD20": "MS4A1",
    "CD45": "PTPRC", "CD45RA": "PTPRC", "CD45RO": "PTPRC", "CD68": "CD68", "CD163": "CD163",
    "CD11c": "ITGAX", "CD11b": "ITGAM", "CD31": "PECAM1", "CD34": "CD34", "CD56": "NCAM1",
    "CD57": "B3GAT1", "FOXP3": "FOXP3", "Ki67": "MKI67", "Ki-67": "MKI67", "PD-1": "PDCD1",
    "PD1": "PDCD1", "PD-L1": "CD274", "PDL1": "CD274", "panCK": "", "Pan-Cytokeratin": "",
    "E-cadherin": "CDH1", "Vimentin": "VIM", "aSMA": "ACTA2", "SMA": "ACTA2", "HLA-DR": "HLA-DRA",
    "GZMB": "GZMB", "Granzyme B": "GZMB", "CD138": "SDC1", "CD15": "FUT4", "MPO": "MPO",
    "Collagen IV": "COL4A1", "Podoplanin": "PDPN", "CD38": "CD38", "CD19": "CD19", "CD14": "CD14",
    "CD16": "FCGR3A", "LAG3": "LAG3", "TIM3": "HAVCR2", "ICOS": "ICOS", "CD39": "ENTPD1",
    "pS6": "", "beta-catenin": "CTNNB1", "Histone H3": "", "DNA1": "", "DNA2": "",
}


def cognate(markers, override=None):
    m = dict(COGNATE)
    if override and os.path.exists(override):
        m.update(dict(pd.read_csv(override).astype(str).values))
    return [m.get(x, m.get(re.sub(r"[ _]", "", x), "")) for x in markers]


def arcsinh(x, cofactor=5.0):
    return np.arcsinh(np.asarray(x, float) / cofactor)


# ------------------------------------------------------------------ loaders
def load_xenium(spec, raw):
    import scanpy as sc
    cells_out = []
    for s in spec["samples"]:
        d = os.path.join(raw, s["path"])
        A = sc.read_10x_h5(os.path.join(d, "cell_feature_matrix.h5"), gex_only=False)
        ft = A.var["feature_types"].astype(str)
        print(f"  {s['path']}: feature types {ft.value_counts().to_dict()}")
        cp = os.path.join(d, "cells.parquet")
        cells = pd.read_parquet(cp) if os.path.exists(cp) else pd.read_csv(os.path.join(d, "cells.csv.gz"))
        cells = cells.set_index("cell_id").loc[A.obs_names]
        gex = (ft == "Gene Expression").values
        pro = (ft == "Protein Expression").values
        rec = dict(cell_id=np.array([f"{s['sample']}_{c}" for c in A.obs_names]),
                   patient=np.full(A.n_obs, s["patient"]), sample=np.full(A.n_obs, s["sample"]),
                   slide=np.full(A.n_obs, s.get("slide", s["sample"])),
                   x_um=cells["x_centroid"].values, y_um=cells["y_centroid"].values,
                   X_rna=A.X[:, gex].toarray(), gene_names=A.var_names[gex].values)
        if pro.any():
            rec["X_protein"] = arcsinh(A.X[:, pro].toarray(), spec.get("cofactor", 5.0))
            rec["protein_names"] = A.var_names[pro].values
        cells_out.append(rec)
    return _concat(cells_out)


def load_cosmx(spec, raw):
    out = []
    px = spec.get("um_per_px", 0.12028)   # CosMx global pixel size; CHECK in the release README
    for s in spec["samples"]:
        d = os.path.join(raw, s["path"])
        expr = pd.read_csv(glob.glob(os.path.join(d, "*exprMat_file.csv*"))[0])
        meta = pd.read_csv(glob.glob(os.path.join(d, "*metadata_file.csv*"))[0])
        key = ["fov", "cell_ID"]
        df = meta.merge(expr, on=key)
        df = df[df.cell_ID > 0]
        genes = [c for c in expr.columns if c not in key and not re.match(r"(Neg|SystemControl|Negative)", c)]
        prot_cols = [c for c in df.columns if re.match(spec.get("protein_regex", r"^Mean\."), c)]
        print(f"  {s['path']}: {len(df)} cells, {len(genes)} genes, {len(prot_cols)} protein columns")
        rec = dict(cell_id=np.array([f"{s['sample']}_{f}_{c}" for f, c in zip(df.fov, df.cell_ID)]),
                   patient=np.full(len(df), s["patient"]), sample=np.full(len(df), s["sample"]),
                   slide=np.full(len(df), s.get("slide", s["sample"])),
                   x_um=df["CenterX_global_px"].values * px, y_um=df["CenterY_global_px"].values * px,
                   X_rna=df[genes].values, gene_names=np.array(genes))
        if prot_cols:
            rec["X_protein"] = arcsinh(df[prot_cols].values, spec.get("cofactor", 5.0))
            rec["protein_names"] = np.array([re.sub(spec.get("protein_regex", r"^Mean\."), "", c) for c in prot_cols])
        out.append(rec)
    return _concat(out)


def load_table(spec, raw):
    f = os.path.join(raw, spec["path"])
    df = pd.read_parquet(f) if f.endswith(".parquet") else pd.read_csv(f)
    c = spec["columns"]
    missing = [v for v in c.values() if isinstance(v, str) and v not in df.columns]
    if missing:
        raise KeyError(f"columns not found in {f}: {missing}\navailable: {list(df.columns)[:60]}")
    markers = spec.get("markers") or [m for m in df.columns if m not in c.values()]
    rec = dict(cell_id=df[c["cell_id"]].astype(str).values if "cell_id" in c else np.arange(len(df)).astype(str),
               patient=df[c["patient"]].astype(str).values,
               sample=df[c.get("sample", c["patient"])].astype(str).values,
               slide=df[c.get("slide", c.get("sample", c["patient"]))].astype(str).values,
               x_um=df[c["x"]].values * spec.get("um_per_px", 1.0), y_um=df[c["y"]].values * spec.get("um_per_px", 1.0),
               X_protein=arcsinh(df[markers].values, spec.get("cofactor", 5.0)) if spec.get("transform", "arcsinh") == "arcsinh"
               else df[markers].values.astype(float),
               protein_names=np.array(markers))
    if "label_fine" in c:
        rec["cell_type_fine"] = df[c["label_fine"]].astype(str).values
    if "label_coarse" in c:
        rec["cell_type_coarse"] = df[c["label_coarse"]].astype(str).values
    return rec


def _concat(recs):
    keys = set.intersection(*[set(r) for r in recs])
    out = {}
    for k in keys:
        if k.endswith("_names"):
            if any(not np.array_equal(recs[0][k], r[k]) for r in recs):
                raise ValueError(f"{k} differ between samples: harmonise panels first")
            out[k] = recs[0][k]
        else:
            out[k] = np.concatenate([r[k] for r in recs])
    return out


LOADERS = {"xenium": load_xenium, "cosmx": load_cosmx, "table": load_table}


def main():
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    spec = yaml.safe_load(open(a.registry))["datasets"][a.dataset]
    print(f"[harmonise] {a.dataset} via {spec['loader']}")
    c = LOADERS[spec["loader"]](spec, a.raw)
    if "labels" in spec:
        lab = pd.read_csv(os.path.join(a.raw, spec["labels"]["path"])).set_index(spec["labels"]["cell_id"])
        for k_out, k_in in spec["labels"]["columns"].items():
            c[k_out] = lab[k_in].reindex(c["cell_id"]).fillna("NA").astype(str).values
    for k in ("cell_type_fine", "cell_type_coarse"):
        c.setdefault(k, np.full(len(c["cell_id"]), "NA"))
    c["label_source"] = np.full(len(c["cell_id"]), spec.get("label_source", "author"))
    if "protein_names" in c:
        c["protein_cognate_gene"] = np.array(cognate(c["protein_names"], spec.get("marker_gene_map")))
    c["pairing"] = spec.get("pairing", {"type": "unpaired"})
    coords = np.c_[c["x_um"], c["y_um"]]
    c["fold"] = patient_folds(c["patient"], k=spec.get("folds", 5))
    c["spatial_block"] = spatial_blocks(coords, c["slide"])
    errs = validate(c)
    for e in errs:
        print("  WARNING:", e)
    fn = os.path.join(a.out, f"{a.dataset}.h5ad")
    save(c, fn)
    print(f"  wrote {fn}: {len(c['cell_id'])} cells, {len(np.unique(c['patient']))} patients")


if __name__ == "__main__":
    main()
