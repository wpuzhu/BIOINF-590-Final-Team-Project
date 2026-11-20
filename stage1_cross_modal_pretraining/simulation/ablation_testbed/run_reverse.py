"""
Reverse direction of the Stage 1 testbed: pretrain on PROTEIN (Spatium-like), evaluate on
RNA and on protein. Together with run.py (RNA-pretrained, Nicheformer-like) this completes the
2 x 2 design:

                        target = RNA                  target = protein
  RNA-pretrained        in-modality  (run.py)         cross-modal  (run.py)
  protein-pretrained    cross-modal  (this file)      in-modality  (this file)

The protein corpus is written into the SAME shared vocabulary (marker j -> its cognate gene's
column), as a gene-token model would see it. On the RNA target only those P columns were ever
trained; the other genes are new to the model (as for Spatium reading 135 of 5,001 genes).

  python run_reverse.py --scenario coupled --seeds 5
"""
import argparse
import os
import time

import numpy as np
import pandas as pd

from model import finetune, init_weights, pretrain, protein_tokens, rna_tokens, shuffle_genes
from run import LABELS, TEST_P, per_patient_f1, sample_labels
from sim import N_TYPES, make_world, sample_protein, sample_rna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="coupled")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--P", type=int, default=16)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                                  "results", "simulation_ablation_testbed"))
    a = ap.parse_args()
    t0 = time.time()
    w = make_world(a.scenario, P=a.P, seed=11)          # same world as run.py
    tag = a.tag or a.scenario
    G, P = w.G, len(w.cognate)
    V = G + P
    cols = np.where(w.cognate >= 0, w.cognate, G + np.arange(P))
    rng = np.random.default_rng(2)

    # protein pretraining corpus (its own patients / batch)
    _, _, pc, ppat = sample_protein(w, 20000, rng)
    Tc = protein_tokens(pc, ppat, cols, V)
    ncomp = max(4, min(30, P // 2))
    pre = pretrain(Tc, pc, seed=0, n_comp=ncomp)
    pcs = shuffle_genes(pc, seed=0)
    noise = pretrain(protein_tokens(pcs, ppat, cols, V), pcs, seed=0, n_comp=ncomp)
    pre_w, noise_w = (pre[0], pre[1]), (noise[0], noise[1])
    print(f"[{tag}] protein-pretrained ({time.time()-t0:.0f}s)", flush=True)

    # targets
    y_p, _, prot, pat_p = sample_protein(w, 6000, np.random.default_rng(31))
    Xp = protein_tokens(prot, pat_p, cols, V)
    y_r, _, xr = sample_rna(w, 6000, np.random.default_rng(77))
    pat_r = np.random.default_rng(78).integers(0, 6, len(y_r))
    _, _, xcorp = sample_rna(w, 20000, np.random.default_rng(1))
    Xr = rna_tokens(xr, xcorp.mean(0), V)
    targets = {"protein_pre_protein_target": (Xp, y_p, pat_p), "protein_pre_RNA_target": (Xr, y_r, pat_r)}

    rows = []
    for seed in range(a.seeds):
        srng = np.random.default_rng(600 + seed)
        for mode, (X, y, pat) in targets.items():
            te = np.isin(pat, TEST_P)
            for k in LABELS:
                tr = sample_labels(y, pat, k, srng)
                for arm, base in (("scratch", None), ("noise_pre", noise_w), ("pretrained", pre_w)):
                    weights = init_weights(base, [], V, np.random.default_rng(seed + 7))
                    clf = finetune(X[tr], y[tr], weights, N_TYPES, seed)
                    for p_, f in per_patient_f1(y[te], clf.predict(X[te]), pat[te]).items():
                        rows.append(dict(scenario=tag, mode=mode, arm=arm, labels=k, seed=seed,
                                         test_patient=p_, macro_f1=f))
        print(f"[{tag}] seed {seed} done ({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, f"raw_reverse_{tag}.csv"), index=False)
    print(df.groupby(["mode", "arm", "labels"]).macro_f1.mean().unstack().round(3).to_string())


if __name__ == "__main__":
    main()
