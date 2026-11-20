# Running on UMich Great Lakes

Everything is config-driven: a YAML file describes the grid, `make_jobs.py` expands it into a jobs table, and one SLURM array task runs one row. Results go to `$SCRATCH`. Only the small aggregated tables go back into the repo.

## 0. One-time setup

```bash
ssh <uniqname>@greatlakes.arc-ts.umich.edu
git clone https://github.com/wpuzhu/BIOINF-590-Final-Team-Project.git ~/BIOINF-590-Final-Team-Project
cd ~/BIOINF-590-Final-Team-Project
vi cluster/env.sh                      # set ACCOUNT (your SLURM account from my.arc.umich.edu)
module load python3.11-anaconda/2024.02
conda env create -f cluster/environment.yml
# Spatium: install from github.com/ploughhh/Spatium following spatium.readthedocs.io
source cluster/env.sh
bash cluster/smoke/run_smoke.sh        # CPU end-to-end test on simulated data, ~3 min
```

`run_smoke.sh` exports simulated data in the harmonised schema and runs the data/split code of both stages (Stage 2 with the numpy backend, Stage 1 with the PCA arm). It then aggregates with the real UQ code. It must print `SMOKE TEST PASSED` before you spend GPU hours.

## 1. Data

```bash
bash datasets/download_scripts/download_zenodo_mendeley.sh $SCRATCH/raw      # Zenodo/Mendeley sets
# 10x / Bruker / SPATCH need a browser form: put them under $SCRATCH/raw/<name> (see datasets/README.md)
mkdir -p logs
sbatch --account=$ACCOUNT cluster/slurm/harmonise.sbatch                     # all registry entries
python shared/data/register_adjacent.py --rna $SCRATCH/harmonised/spatch_colon_xenium5k.h5ad \
   --protein $SCRATCH/harmonised/spatch_colon_codex.h5ad --out $SCRATCH/harmonised/spatch_colon_registration.json
python stage2_privileged_transfer/pipeline/prepare_adjacent.py --rna ... --protein ... --registration ... \
   --out $SCRATCH/harmonised/spatch_colon_codex_registered_rna.npy
```

Check every harmonise log: loaders print the feature types and columns they found and warn if a dataset has < 4 patients.

## 2. Stage 1 (cross-modal pretraining)

```bash
python cluster/make_jobs.py --stage 1 --out jobs_stage1.tsv          # all datasets x models x arms x modes x folds x seeds
N=$(($(wc -l < jobs_stage1.tsv) - 1))
JID=$(sbatch --parsable --account=$ACCOUNT --array=1-$N%40 cluster/slurm/stage1_array.sbatch jobs_stage1.tsv)
sbatch --account=$ACCOUNT --dependency=afterany:$JID cluster/slurm/aggregate.sbatch 1
```

Start with a pilot: `--datasets spatch_colon_codex xenium_protein_rcc`, and dry-run one row with `python stage1_cross_modal_pretraining/pipeline/run_stage1.py ... --dry_run`. The adapter lines marked `VERIFY` (Nicheformer context tokens, Spatium loading) must be checked against the installed package versions first.

## 3. Stage 2 (RNA as privileged information)

```bash
# (a) mechanism comparison on 'default'
python cluster/make_jobs.py --stage 2 --settings default --out jobs_s2_default.tsv
JID=$(sbatch --parsable --account=$ACCOUNT --array=1-$(($(wc -l < jobs_s2_default.tsv)-1))%40 \
      cluster/slurm/stage2_array.sbatch jobs_s2_default.tsv)
sbatch --account=$ACCOUNT --dependency=afterany:$JID cluster/slurm/aggregate.sbatch 2
# (b) choose the winner (largest RNA-specific contrast, mechanism - self_train), put it in the ablation
#     settings' `mechanisms:` list in configs/stage2.yaml, then run every other setting
python cluster/make_jobs.py --stage 2 --settings $(python -c "import yaml;print(' '.join(k for k in yaml.safe_load(open('stage2_privileged_transfer/pipeline/configs/stage2.yaml'))['settings'] if k!='default'))") --out jobs_s2_ablation.tsv
```

## Resources (starting points)

| Job | Partition | GPU | Memory | Time |
|---|---|---|---|---|
| harmonise | `standard` | — | 128 GB | 1–6 h per dataset |
| Stage 1 frozen embedding | `gpu` / `spgpu` | 1 | 64 GB | ≤ 1 h |
| Stage 1 fine-tune | `spgpu` (A40) | 1 | 64 GB | 2–4 h |
| Stage 2 one fold × seed, all mechanisms | `gpu` / `spgpu` | 1 | 64 GB | 1–4 h (CosMx ~800k cells) |
| aggregate | `standard` | — | 32 GB | < 1 h |

Throttle arrays with `%40` to stay within the account's GPU limit. Check partitions with `sinfo -s` and usage with `my_job_statistics <jobid>`.
