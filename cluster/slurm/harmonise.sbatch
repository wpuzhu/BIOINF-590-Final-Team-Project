#!/bin/bash
#SBATCH --job-name=harmonise
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
# Usage: sbatch --account=$ACCOUNT cluster/slurm/harmonise.sbatch [dataset ...]   (default: all in registry)
set -euo pipefail
source "$REPO/cluster/env.sh"
cd "$REPO/shared/data"
DS=("$@")
if [ ${#DS[@]} -eq 0 ]; then mapfile -t DS < <(python -c "import yaml;print('\n'.join(yaml.safe_load(open('$REPO/datasets/registry.yaml'))['datasets']))"); fi
for d in "${DS[@]}"; do
  python harmonise.py --registry "$REPO/datasets/registry.yaml" --dataset "$d" --raw "$SCRATCH/raw" --out "$SCRATCH/harmonised" || echo "FAILED $d"
done
