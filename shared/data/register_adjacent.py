"""
Rigid registration of an adjacent protein section onto its RNA section (SPATCH-style data),
and the neighbourhood-RNA teacher input that Stage 2 uses for adjacent-section pairing.

Method: coarse grid search over rotation (+ optional flip) maximising the overlap of cell-density
images, then ICP on centroids of cells of the SAME coarse type (types are called separately in
each modality). The residual nearest-neighbour distance of matched same-type cells after
convergence is reported as the registration error; Stage 2 uses a radius >= 2x that error.

  python register_adjacent.py --rna rna.h5ad --protein prot.h5ad --out reg.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spt_io import load  # noqa: E402


def _rot(t):
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def density(xy, lo, hi, bins=64):
    H, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=bins, range=[[lo[0], hi[0]], [lo[1], hi[1]]])
    return H / (H.sum() + 1e-9)


def coarse_search(src, dst, n_angles=72, flips=(False, True)):
    lo, hi = np.minimum(src.min(0), dst.min(0)) - 200, np.maximum(src.max(0), dst.max(0)) + 200
    Hd = density(dst - dst.mean(0) + (lo + hi) / 2, lo, hi)
    best = (-1, None)
    for flip in flips:
        s = src * np.array([1, -1]) if flip else src
        for t in np.linspace(0, 2 * np.pi, n_angles, endpoint=False):
            r = (s - s.mean(0)) @ _rot(t).T + (lo + hi) / 2
            score = (density(r, lo, hi) * Hd).sum()
            if score > best[0]:
                best = (score, (flip, t))
    flip, t = best[1]
    F = np.diag([1, -1]) if flip else np.eye(2)
    R = _rot(t) @ F
    return R, dst.mean(0) - src.mean(0) @ R.T


def icp(src, dst, src_type, dst_type, R, b, iters=50, trim=0.8):
    for _ in range(iters):
        cur = src @ R.T + b
        P, Q = [], []
        for c in np.intersect1d(np.unique(src_type), np.unique(dst_type)):
            s, d = cur[src_type == c], dst[dst_type == c]
            if len(s) < 10 or len(d) < 10:
                continue
            dist, j = NearestNeighbors(n_neighbors=1).fit(d).kneighbors(s)
            keep = dist[:, 0] <= np.quantile(dist[:, 0], trim)
            P.append(src[src_type == c][keep]); Q.append(d[j[keep, 0]])
        P, Q = np.vstack(P), np.vstack(Q)
        mp, mq = P.mean(0), Q.mean(0)
        U, _, Vt = np.linalg.svd((P - mp).T @ (Q - mq))
        Rn = (U @ Vt).T
        if np.linalg.det(Rn) < 0 and np.linalg.det(R) > 0:
            Vt[-1] *= -1; Rn = (U @ Vt).T
        R, b = Rn, mq - mp @ Rn.T
    cur = src @ R.T + b
    err = np.median(NearestNeighbors(n_neighbors=1).fit(Q).kneighbors(P @ R.T + b)[0])
    return R, b, float(err)


def neighbourhood_rna(rna_xy, rna_X, query_xy, radius):
    """Mean RNA of RNA-section cells within `radius` of each registered protein cell (NaN if none)."""
    nb = NearestNeighbors(radius=radius).fit(rna_xy).radius_neighbors(query_xy, return_distance=False)
    out = np.full((len(query_xy), rna_X.shape[1]), np.nan, dtype=np.float32)
    for i, js in enumerate(nb):
        if len(js):
            out[i] = rna_X[js].mean(0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--type_col", default="cell_type_coarse")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R_, P_ = load(a.rna), load(a.protein)
    res = {}
    for slide in np.unique(P_["slide"]):
        pm = P_["slide"] == slide
        pat = P_["patient"][pm][0]
        rm = R_["patient"] == pat
        src, dst = P_.coords[pm], R_.coords[rm]
        R0, b0 = coarse_search(src, dst)
        R, b, err = icp(src, dst, P_[a.type_col][pm], R_[a.type_col][rm], R0, b0)
        res[str(slide)] = dict(patient=str(pat), R=R.tolist(), b=b.tolist(), error_um=err,
                               radius_um=float(max(30.0, 2 * err)))
        print(f"  {slide}: median matched distance {err:.1f} um")
    json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
