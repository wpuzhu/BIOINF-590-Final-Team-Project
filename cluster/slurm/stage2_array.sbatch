#!/bin/bash
#SBATCH --job-name=s2_privileged
#SBATCH --partition=gpu               # or spgpu (A40) / gpu_mig40; check `sinfo -s`
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Usage: sbatch --account=$ACCOUNT --array=1-N%40 cluster/slurm/stage2_array.sbatch jobs_stage2.tsv
set -euo pipefail
source "$REPO/cluster/env.sh"
JOBS=${1:?jobs tsv}
read -r -a H < <(head -1 "$JOBS")
read -r -a V < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$JOBS")
ARGS=(); for i in "${!H[@]}"; do ARGS+=("--${H[$i]}" "${V[$i]}"); done
cd "$REPO/stage2_privileged_transfer/pipeline"
python train_stage2.py "${ARGS[@]}" --out "$SCRATCH/results/stage2"
