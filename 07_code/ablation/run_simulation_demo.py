"""
Self-check for xmodal_ablation.py on simulated paired spatial RNA + protein data with KNOWN
answers. Each simulated protein gets its information from a different source, so each
ablation should hurt a specific protein and leave the others alone:

  CD3     cell-intrinsic, read off its cognate gene CD3E          -> cognate ablation hurts
  Ki67    cell-intrinsic, cognate MKI67                           -> cognate ablation hurts
  HLADR   driven by a gene programme (CIITA + HLA genes), cognate weak -> group ablation hurts
  PDL1    in tumour cells, induced by IFNG from NEARBY T cells    -> neighbourhood helps;
                                                                     radius sweep peaks ~induction range;
                                                                     shuffling neighbours within type kills gain
  pS6     spatial gradient, no RNA signal (post-transcriptional)  -> nothing predicts it;
                                                                     residual stays spatially structured
  CD68    pure cell-type marker                                   -> cell type alone suffices

Usage:  python run_simulation_demo.py [--out results_sim]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from xmodal_ablation import (Features, PairedSpatial, context_controls, gene_ablation, pairing_degradation,
                             radius_sweep, region_holdout, residual_map, variance_ladder,
                             radius_graph)

TYPES = ["Tumor", "Tcell", "Macro", "Fibro"]
REGION_MIX = {  # composition by region (core | margin | stroma)
    "core":   [0.80, 0.05, 0.10, 0.05],
    "margin": [0.40, 0.30, 0.15, 0.15],
    "stroma": [0.10, 0.20, 0.20, 0.50],
}


def simulate(n_samples=4, n_cells=1500, n_genes=150, induction_um=30, seed=1):
    rng = np.random.default_rng(seed)
    named = ["CD3E", "MKI67", "CD274", "IFNG", "CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1",
             "CD68", "EPCAM", "COL1A1", "RPS6"]
    genes = named + [f"G{i:03d}" for i in range(n_genes - len(named))]
    gi = {g: i for i, g in enumerate(genes)}
    X, Y, C, T, S, R = [], [], [], [], [], []
    type_mean = rng.normal(0, 1, (len(TYPES), n_genes)) * (rng.random((len(TYPES), n_genes)) < 0.15)
    for s in range(n_samples):
        xy = rng.uniform(0, 1000, (n_cells, 2))
        region = np.where(xy[:, 0] < 333, "core", np.where(xy[:, 0] < 666, "margin", "stroma"))
        ct = np.array([rng.choice(4, p=REGION_MIX[r]) for r in region])
        x = type_mean[ct] + rng.normal(0, 0.6, (n_cells, n_genes))
        # marker genes
        x[:, gi["CD3E"]] += 2.0 * (ct == 1)
        x[:, gi["IFNG"]] += (1.5 + rng.normal(0, 0.8, n_cells)) * (ct == 1)
        x[:, gi["CD68"]] += 2.0 * (ct == 2)
        x[:, gi["EPCAM"]] += 2.0 * (ct == 0)
        x[:, gi["COL1A1"]] += 2.0 * (ct == 3)
        x[:, gi["MKI67"]] += rng.gamma(1.0, 0.8, n_cells) * (ct == 0)
        mhc = rng.normal(0, 1, n_cells) * (ct == 2) + 0.3 * rng.normal(0, 1, n_cells)
        for g in ["CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"]:
            x[:, gi[g]] += mhc + rng.normal(0, 0.5, n_cells)
        # neighbourhood IFNG within induction range (Gaussian kernel)
        d2 = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)
        K = np.exp(-d2 / (2 * induction_um ** 2)); np.fill_diagonal(K, 0)
        ifng_env = K @ (x[:, gi["IFNG"]] * (ct == 1)) / 3.0
        x[:, gi["CD274"]] += 0.15 * ifng_env * (ct == 0) + rng.normal(0, 0.5, n_cells)
        # proteins
        y = np.zeros((n_cells, 6))
        y[:, 0] = 0.9 * x[:, gi["CD3E"]] + rng.normal(0, 0.4, n_cells)                    # CD3
        y[:, 1] = 0.9 * x[:, gi["MKI67"]] + rng.normal(0, 0.4, n_cells)                   # Ki67
        y[:, 2] = 1.2 * mhc + 0.1 * x[:, gi["HLA-DRA"]] + rng.normal(0, 0.4, n_cells)     # HLADR
        y[:, 3] = (1.0 * ifng_env + rng.normal(0, 0.3, n_cells)) * (ct == 0) + rng.normal(0, 0.3, n_cells)  # PDL1
        y[:, 4] = (np.sin(xy[:, 1] / 160) * 1.2) * (ct == 0) + rng.normal(0, 0.4, n_cells)  # pS6
        y[:, 5] = 2.0 * (ct == 2) + rng.normal(0, 0.4, n_cells)                            # CD68
        X.append(x); Y.append(y); C.append(xy); T.append(np.array(TYPES)[ct])
        S.append(np.full(n_cells, f"S{s}")); R.append(region)
    Y = np.vstack(Y)
    Y = (Y - Y.mean(0)) / Y.std(0)
    return PairedSpatial(
        rna=np.vstack(X), protein=Y, coords=np.vstack(C), cell_type=np.concatenate(T),
        sample=np.concatenate(S), region=np.concatenate(R), genes=genes,
        proteins=["CD3", "Ki67", "HLADR", "PDL1", "pS6", "CD68"],
        cognate={"CD3": ["CD3E"], "Ki67": ["MKI67"], "HLADR": ["HLA-DRA"], "PDL1": ["CD274"],
                 "pS6": ["RPS6"], "CD68": ["CD68"]},
        gene_groups={"MHCII_programme": ["CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"],
                     "IFNG": ["IFNG"]})


def heat(df, col, title, ax, cmap="viridis", center=None):
    m = df.pivot_table(index="protein", columns="condition", values=col)
    vmax = np.nanmax(np.abs(m.values)) if center is not None else None
    im = ax.imshow(m.values, aspect="auto", cmap=cmap,
                   vmin=-vmax if center is not None else None, vmax=vmax)
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=8)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if (center is None and m.values[i, j] < np.nanmean(m.values)) else "black")
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.03)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="../../08_results/simulation/ablation")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    d = simulate()
    F = Features(d, n_pcs=30)
    print("cells:", len(d.sample), "genes:", len(d.genes), "proteins:", d.proteins)

    lad = variance_ladder(d, F); lad.to_csv(f"{a.out}/1_variance_ladder.csv", index=False)
    abl = gene_ablation(d, F, n_random=3); abl.to_csv(f"{a.out}/2_gene_ablation.csv", index=False)
    rad = radius_sweep(d, F); rad.to_csv(f"{a.out}/3_radius_sweep.csv", index=False)
    ctl = context_controls(d, F); ctl.to_csv(f"{a.out}/4_context_controls.csv", index=False)
    reg = region_holdout(d, F); reg.to_csv(f"{a.out}/5_region_holdout.csv", index=False)
    resid, _ = residual_map(d, F); resid.to_csv(f"{a.out}/6_residual_moran.csv", index=False)
    jit = pd.concat([pairing_degradation(d, F).assign(jitter="any_cell"),
                     pairing_degradation(d, F, same_type=True).assign(jitter="same_type")])
    jit.to_csv(f"{a.out}/8_pairing_degradation.csv", index=False)

    abl_s = abl.copy()
    abl_s.loc[abl_s.condition.str.startswith("drop_random"), "condition"] = "drop_random_mean"
    abl_s = abl_s.groupby(["protein", "condition"], as_index=False)["delta_vs_full"].mean()

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    heat(lad, "r2", "1. Variance ladder: held-out R2 as inputs are added", axes[0, 0])
    heat(abl_s, "delta_vs_full", "2. Gene ablation: change in R2 vs full input", axes[0, 1], cmap="RdBu", center=0)
    heat(ctl, "r2", "4. Context controls: is it real spatial arrangement?", axes[1, 0])
    for p, g in rad.groupby("protein"):
        axes[1, 1].plot([int(c.split("_")[1]) for c in g.condition], g.r2, marker="o", label=p)
    axes[1, 1].set_xlabel("neighbourhood radius (um)"); axes[1, 1].set_ylabel("held-out R2")
    axes[1, 1].set_title("3. Radius sweep", fontsize=9); axes[1, 1].legend(fontsize=7)
    fig.suptitle("Input-level ablation on simulated paired spatial RNA + protein (leave-one-sample-out)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{a.out}/ablation_overview.png", dpi=150)

    # console summary
    pd.set_option("display.width", 200); pd.set_option("display.max_columns", 20)
    print("\n1) variance ladder\n", lad.pivot_table(index="protein", columns="condition", values="r2").round(2))
    print("\n2) gene ablation (delta R2)\n", abl_s.pivot_table(index="protein", columns="condition", values="delta_vs_full").round(2))
    print("\n3) radius sweep\n", rad.pivot_table(index="protein", columns="condition", values="r2").round(2))
    print("\n4) context controls\n", ctl.pivot_table(index="protein", columns="condition", values="r2").round(2))
    print("\n5) region hold-out (pearson r)\n", reg.pivot_table(index=["protein"], columns=["heldout_region", "condition"], values="pearson").round(2))
    print("\n6) residual Moran's I in tumour cells\n", resid[resid.cell_type == "Tumor"].round(3))
    print("\n8) pairing degradation (R2)\n", jit.pivot_table(index=["jitter", "protein"], columns="condition", values="r2").round(2))


if __name__ == "__main__":
    main()
