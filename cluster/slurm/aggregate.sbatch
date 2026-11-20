#!/bin/bash
#SBATCH --job-name=aggregate_uq
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
# Usage: sbatch --account=$ACCOUNT --dependency=afterany:<array_job_id> cluster/slurm/aggregate.sbatch 2
set -euo pipefail
source "$REPO/cluster/env.sh"
python "$REPO/shared/uq/aggregate.py" --stage "${1:-2}" --results "$SCRATCH/results/stage${1:-2}" \
       --out "$REPO/stage${1:-2}_$( [ "${1:-2}" = 1 ] && echo cross_modal_pretraining || echo privileged_transfer )/results/real"
