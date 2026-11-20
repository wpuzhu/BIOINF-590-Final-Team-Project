"""
Adjacent-section teacher input (SPATCH): for every protein (CODEX) cell, the mean log-normalised
RNA of the RNA-section cells within the registration radius of its registered position.

  python prepare_adjacent.py --rna $SCRATCH/harmonised/spatch_colon_xenium5k.h5ad \
      --protein $SCRATCH/harmonised/spatch_colon_codex.h5ad \
      --registration $SCRATCH/harmonised/spatch_colon_registration.json \
      --out $SCRATCH/harmonised/spatch_colon_codex_registered_rna.npy [--radius_scale 1.0]

`--radius_scale` > 1 is the adjacent-section analogue of the jitter ablation.
Cells with no RNA neighbour get the slide mean (and are counted in the printout).
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "shared", "data"))
from register_adjacent import neighbourhood_rna  # noqa: E402
from spt_io import load  # noqa: E402
from train_stage2 import lognorm  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--protein", required=True)
    ap.add_argument("--registration", required=True)
    ap.add_argument("--radius_scale", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R, P = load(a.rna), load(a.protein)
    reg = json.load(open(a.registration))
    X = lognorm(R["X_rna"]).astype(np.float32)
    out = np.zeros((P.n, X.shape[1]), np.float32)
    for slide, r in reg.items():
        pm = P["slide"] == slide
        rm = R["patient"] == r["patient"]
        q = P.coords[pm] @ np.array(r["R"]).T + np.array(r["b"])
        nb = neighbourhood_rna(R.coords[rm], X[rm], q, r["radius_um"] * a.radius_scale)
        miss = np.isnan(nb[:, 0])
        nb[miss] = X[rm].mean(0)
        out[pm] = nb
        print(f"  {slide}: radius {r['radius_um'] * a.radius_scale:.0f} um, {miss.mean():.1%} cells without RNA neighbours")
    np.save(a.out, out)


if __name__ == "__main__":
    main()
