"""
Run the H1 transfer-ablation grid for one scenario.

  python run.py --scenario realistic [--seeds 5] [--out DIR]

Writes <out>/raw_<scenario>.csv with one row per
  (mode: finetune | frozen, arm, labels_per_class, seed, test_patient) -> macro-F1
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score

from model import (column_map, embed, finetune, frozen_probe, init_weights, pretrain, protein_tokens,
                   random_weights, rna_tokens, shuffle_genes)
from sim import N_TYPES, make_world, sample_protein, sample_rna

FT_ARMS = ["scratch", "noise_pre", "pretrained", "A1_reinit", "A2_random_map",
           "A3_overlap_0.75", "A3_overlap_0.50", "A3_overlap_0.25"]
FROZEN_ARMS = ["random_features", "noise_pre", "pretrained", "A1_reinit", "A2_random_map"]
LABELS = [5, 20, 100]
TRAIN_P, TEST_P = [0, 1, 2], [3, 4, 5]


def per_patient_f1(y, yhat, patient):
    return {p: f1_score(y[patient == p], yhat[patient == p], average="macro",
                        labels=np.arange(N_TYPES), zero_division=0) for p in TEST_P}


def sample_labels(y, patient, k, rng):
    idx = []
    pool = np.where(np.isin(patient, TRAIN_P))[0]
    for c in range(N_TYPES):
        cand = pool[y[pool] == c]
        idx.extend(rng.choice(cand, size=min(k, len(cand)), replace=False))
    return np.array(idx)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="realistic")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--out", default="../../../08_results/simulation/h1_ablation")
    ap.add_argument("--rna_target", action="store_true", help="also run the in-modality RNA reference tasks")
    ap.add_argument("--P", type=int, default=16, help="protein panel size")
    ap.add_argument("--tag", default=None, help="name used in output file (default: scenario)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    w = make_world(a.scenario, P=a.P, seed=11)
    tag = a.tag or a.scenario
    G, P = w.G, len(w.cognate)
    V = G + P
    rng = np.random.default_rng(1)

    # ---- RNA pretraining corpus ("atlas": its own batch effect)
    _, _, x_corpus = sample_rna(w, 20000, rng)
    gene_mean = x_corpus.mean(0)
    T_corpus = np.hstack([rna_tokens(x_corpus, gene_mean, G), np.zeros((len(x_corpus), P))])
    pre = pretrain(T_corpus, x_corpus, seed=0)
    x_shuf = shuffle_genes(x_corpus, seed=0)
    T_shuf = np.hstack([rna_tokens(x_shuf, gene_mean, G), np.zeros((len(x_shuf), P))])
    noise = pretrain(T_shuf, x_shuf, seed=0)
    pre_w, noise_w = (pre[0], pre[1]), (noise[0], noise[1])
    print(f"[{a.scenario}] pretrained: final loss real {pre[2][-1]:.4f} vs shuffled-corpus {noise[2][-1]:.4f} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # ---- protein target data
    y, _, prot, patient = sample_protein(w, 6000, rng)
    test = np.isin(patient, TEST_P)

    rows = []
    for seed in range(a.seeds):
        srng = np.random.default_rng(100 + seed)
        maps = {arm: column_map(arm if not arm.startswith("A3") else "A3_overlap", w.cognate, G,
                                np.random.default_rng(1000 * seed + i),
                                q=float(arm.split("_")[-1]) if arm.startswith("A3") else 1.0)
                for i, arm in enumerate(FT_ARMS)}
        maps["random_features"] = maps["scratch"]
        tokens = {arm: protein_tokens(prot, patient, cols, V) for arm, (cols, _) in maps.items()}
        for k in LABELS:
            tr = sample_labels(y, patient, k, srng)
            # PCA baseline (frozen-mode reference, as in Justin's benchmark)
            Z = prot.copy()
            for p in np.unique(patient):
                m = patient == p
                Z[m] = (Z[m] - Z[m].mean(0)) / (Z[m].std(0) + 1e-9)
            pca = PCA(n_components=P).fit(Z[np.isin(patient, TRAIN_P)])
            lr = LogisticRegression(C=1.0, max_iter=3000).fit(pca.transform(Z[tr]), y[tr])
            for p_, f in per_patient_f1(y[test], lr.predict(pca.transform(Z[test])), patient[test]).items():
                rows.append(dict(mode="frozen", arm="PCA_probe", labels=k, seed=seed, test_patient=p_, macro_f1=f))
            # fine-tuned arms
            for arm in FT_ARMS:
                cols, reinit = maps[arm]
                X = tokens[arm]
                base = None if arm == "scratch" else (noise_w if arm == "noise_pre" else pre_w)
                weights = init_weights(base, reinit, V, np.random.default_rng(seed + 7))
                clf = finetune(X[tr], y[tr], weights, N_TYPES, seed)
                for p_, f in per_patient_f1(y[test], clf.predict(X[test]), patient[test]).items():
                    rows.append(dict(mode="finetune", arm=arm, labels=k, seed=seed, test_patient=p_, macro_f1=f))
            # frozen arms
            for arm in FROZEN_ARMS:
                cols, reinit = maps[arm]
                X = tokens[arm]
                if arm == "random_features":
                    weights = random_weights(V, np.random.default_rng(seed + 7))
                else:
                    base = noise_w if arm == "noise_pre" else pre_w
                    weights = init_weights(base, reinit, V, np.random.default_rng(seed + 7))
                yhat = frozen_probe(X[tr], y[tr], X[test], weights)
                for p_, f in per_patient_f1(y[test], yhat, patient[test]).items():
                    rows.append(dict(mode="frozen", arm=arm, labels=k, seed=seed, test_patient=p_, macro_f1=f))
        print(f"[{a.scenario}] seed {seed} done ({time.time()-t0:.0f}s)", flush=True)

    # ---- in-modality reference: RNA target from a new dataset (new batch), cells typed from RNA
    if a.rna_target:
        yr, _, xr = sample_rna(w, 6000, np.random.default_rng(77))
        pr = np.random.default_rng(78).integers(0, 6, len(yr))
        Xr = np.hstack([rna_tokens(xr, gene_mean, G), np.zeros((len(yr), P))])
        # RNA restricted to the cognate genes of the protein panel: same modality, protein-sized panel
        cg = w.cognate[w.cognate >= 0]
        rel = xr[:, cg] / (gene_mean[cg] + 1e-6)
        Xpanel = __import__("model").tokenize(np.where(xr[:, cg] > 0, rel, -np.inf), cg, V, L=len(cg))
        te = np.isin(pr, TEST_P)
        for seed in range(a.seeds):
            srng = np.random.default_rng(500 + seed)
            for k in LABELS:
                tr = sample_labels(yr, pr, k, srng)
                perm = np.random.default_rng(900 + seed).permutation(G)      # A2 for RNA: wrong gene->token map
                Xr_perm = Xr.copy()
                Xr_perm[:, :G] = Xr[:, :G][:, np.argsort(perm)]
                for arm, base, reinit in [("scratch", None, []), ("noise_pre", noise_w, []),
                                          ("pretrained", pre_w, []), ("A1_reinit", pre_w, list(range(G))),
                                          ("A2_random_map", pre_w, [])]:
                    weights = init_weights(base, reinit, V, np.random.default_rng(seed + 7))
                    if arm in ("A1_reinit", "A2_random_map"):
                        pairs = [("finetune_RNA_target", Xr_perm if arm == "A2_random_map" else Xr)]
                    else:
                        pairs = [("finetune_RNA_target", Xr), ("finetune_RNA_panel", Xpanel)]
                    for mode, XX in pairs:
                        clf = finetune(XX[tr], yr[tr], weights, N_TYPES, seed)
                        for p_, f in per_patient_f1(yr[te], clf.predict(XX[te]), pr[te]).items():
                            rows.append(dict(mode=mode, arm=arm, labels=k, seed=seed,
                                             test_patient=p_, macro_f1=f))
        print(f"[{a.scenario}] RNA reference done ({time.time()-t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.insert(0, "scenario", tag)
    df.to_csv(f"{a.out}/raw_{tag}.csv", index=False)
    print(df.groupby(["mode", "arm", "labels"]).macro_f1.mean().unstack().round(3).to_string())


if __name__ == "__main__":
    main()
