"""
Expand a stage config into a jobs table: one line per SLURM array task.

  python make_jobs.py --stage 2 --out jobs_stage2.tsv [--settings default jitter_30 ...] [--datasets ...]
  sbatch --array=1-$(($(wc -l < jobs_stage2.tsv)-1))%40 slurm/stage2_array.sbatch jobs_stage2.tsv
Each array task reads line $SLURM_ARRAY_TASK_ID (1-based, after the header) and passes the
columns as --flags to the stage script.
"""
import argparse
import itertools
import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def stage1(cfg, a):
    g = cfg["grid"]
    for ds, spec in cfg["datasets"].items():
        if a.datasets and ds not in a.datasets:
            continue
        for mod in spec["modalities"]:
            for model in g["models"]:
                for arm in g["arms"]:
                    if model == "pca" and arm != "pretrained":
                        continue
                    for mode in g["modes"]:
                        if model == "pca" and mode != "frozen":
                            continue
                        qs = g["a3_q"] if arm == "A3_overlap" else [1.0]
                        for q, f, s in itertools.product(qs, range(cfg["folds"]), cfg["seeds"]):
                            yield dict(dataset=ds, modality=mod, model=model, arm=arm, q=q, mode=mode, fold=f, seed=s)


def stage2(cfg, a):
    sets = a.settings or list(cfg["settings"])
    for ds in cfg["datasets"]:
        if a.datasets and ds not in a.datasets:
            continue
        for st, f, s in itertools.product(sets, range(cfg["folds"]), cfg["seeds"]):
            yield dict(dataset=ds, setting=st, fold=f, seed=s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, choices=[1, 2], required=True)
    ap.add_argument("--config")
    ap.add_argument("--datasets", nargs="*")
    ap.add_argument("--settings", nargs="*")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    default = {1: "stage1_cross_modal_pretraining/pipeline/configs/stage1.yaml",
               2: "stage2_privileged_transfer/pipeline/configs/stage2.yaml"}[a.stage]
    cfg = yaml.safe_load(open(a.config or os.path.join(ROOT, default)))
    rows = list((stage1 if a.stage == 1 else stage2)(cfg, a))
    with open(a.out, "w") as fh:
        fh.write("\t".join(rows[0]) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r.values()) + "\n")
    print(f"{len(rows)} jobs -> {a.out}")


if __name__ == "__main__":
    main()
