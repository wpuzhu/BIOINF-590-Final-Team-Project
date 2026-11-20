"""
Stage 1 summary on the testbed: the 2 x 2 (pretraining modality x target modality) gains with
uncertainty, and cross-modal transfer efficiency.

  gain(A -> B)      = pretrained - scratch, fine-tuned on B with k labels/class
  generic(A -> B)   = noise_pre - scratch (a pretrained network with the biology removed)
  biology(A -> B)   = pretrained - noise_pre
  efficiency(A->B)  = gain(A -> B) / gain(B -> B)  (cross-modal gain as a share of in-modality gain)

Unit = test patient (3 per run) x seed (5). With 3 patients the smallest two-sided sign-flip
p is 0.25, so the testbed reports CIs and the pre-registered equivalence verdict; the real-data
pipeline has many more patients.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "shared", "uq"))
from uq_stats import hierarchical_bootstrap, paired_table, tost_equivalence  # noqa: E402

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "results", "simulation_ablation_testbed")
MARGIN = 0.02
CELLS = {("RNA", "RNA"): "finetune_RNA_target", ("RNA", "protein"): "finetune",
         ("protein", "protein"): "protein_pre_protein_target", ("protein", "RNA"): "protein_pre_RNA_target"}


def main():
    rows = []
    for sc in ("coupled", "realistic", "decoupled"):
        parts = [pd.read_csv(f) for f in (f"{D}/raw_{sc}.csv", f"{D}/raw_reverse_{sc}.csv") if os.path.exists(f)]
        if len(parts) < 2:
            continue
        df = pd.concat(parts).rename(columns={"test_patient": "patient"})
        df["unit"] = df["patient"].astype(str) + "_" + df["seed"].astype(str)
        for (src, tgt), mode in CELLS.items():
            for k in (5, 20, 100):
                g = df[(df["mode"] == mode) & (df.labels == k)]
                if g.empty:
                    continue
                for name, (a1, a0) in {"gain": ("pretrained", "scratch"), "generic": ("noise_pre", "scratch"),
                                       "biology": ("pretrained", "noise_pre")}.items():
                    if a1 not in g.arm.values or a0 not in g.arm.values:
                        continue
                    d = paired_table(g, "arm", "macro_f1", a1, a0, keys=("patient", "seed"))
                    r = hierarchical_bootstrap(d, B=4000)
                    r["tost_equivalent_0.02"] = tost_equivalence(d, MARGIN, B=4000)["equivalent"]
                    rows.append(dict(scenario=sc, pretrain=src, target=tgt, labels=k, contrast=name, **r))
    out = pd.DataFrame(rows)
    out.to_csv(f"{D}/stage1_2x2_gains.csv", index=False)
    g = out[out.contrast == "gain"].set_index(["scenario", "labels", "pretrain", "target"]).estimate
    eff = []
    for (sc, k, src, tgt), v in g.items():
        if src != tgt and (sc, k, tgt, tgt) in g.index:
            ref = g[(sc, k, tgt, tgt)]
            eff.append(dict(scenario=sc, labels=k, direction=f"{src}->{tgt}", cross_gain=v, in_modality_gain=ref,
                            efficiency=v / ref if ref > 0.01 else np.nan))
    eff = pd.DataFrame(eff)
    eff.to_csv(f"{D}/stage1_transfer_efficiency.csv", index=False)
    pd.set_option("display.width", 200)
    print(out[out.contrast == "gain"][["scenario", "pretrain", "target", "labels", "estimate", "ci_lo", "ci_hi",
                                        "tost_equivalent_0.02"]].round(3).to_string(index=False))
    print(eff.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
