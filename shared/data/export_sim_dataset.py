"""Write a simulated paired dataset (shared/sim) in the harmonised schema, for pipeline smoke tests.

  python export_sim_dataset.py --scenario transferable_protein --out /tmp/smoke/sim_tp.npz
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "sim")]
from paired_spatial_sim import SimConfig, simulate  # noqa: E402
from spt_io import patient_folds, save, spatial_blocks, validate  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--scenario", default="transferable_protein")
ap.add_argument("--cells", type=int, default=800)
ap.add_argument("--out", required=True)
a = ap.parse_args()
d = simulate(SimConfig(scenario=a.scenario, cells_per_patient=a.cells))
n = len(d.fine)
counts = np.random.default_rng(0).poisson(20 * np.exp(0.5 * d.rna))      # RNA as counts (~deep imaging panel)
c = dict(cell_id=np.arange(n).astype(str), patient=np.array([f"P{p}" for p in d.patient]),
         sample=np.array([f"P{p}" for p in d.patient]), slide=np.array([f"P{p}" for p in d.patient]),
         x_um=d.coords[:, 0], y_um=d.coords[:, 1],
         cell_type_coarse=np.array([f"T{t}" for t in d.coarse]), cell_type_fine=np.array([f"T{t // 2}_{t % 2}" for t in d.fine]),
         label_source=np.full(n, "simulated"), X_protein=d.protein.astype(np.float32), X_rna=counts.astype(np.float32),
         protein_names=np.array([f"M{j}" for j in range(d.protein.shape[1])]),
         gene_names=np.array([f"G{j}" for j in range(d.rna.shape[1])]),
         protein_cognate_gene=np.array([f"G{j}" for j in range(d.protein.shape[1])]),
         pairing={"type": "same_cell"})
c["fold"] = patient_folds(c["patient"], k=5)
c["spatial_block"] = spatial_blocks(np.c_[c["x_um"], c["y_um"]], c["slide"])
print("schema check:", validate(c) or "OK")
save(c, a.out)
print("wrote", a.out)
