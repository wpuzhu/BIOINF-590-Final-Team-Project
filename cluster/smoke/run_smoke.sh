#!/bin/bash
# End-to-end CPU smoke test of both pipelines on simulated data (no GPU, no torch, ~5 min).
#   bash cluster/smoke/run_smoke.sh
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd); T=/tmp/bioinf590_smoke; mkdir -p $T
for s in transferable_protein informative_only; do
  python "$REPO/shared/data/export_sim_dataset.py" --scenario $s --out $T/sim_$s.npz
done
python "$REPO/cluster/make_jobs.py" --stage 2 --config "$REPO/cluster/smoke/smoke_stage2.yaml" --out $T/jobs2.tsv
python "$REPO/cluster/make_jobs.py" --stage 1 --config "$REPO/cluster/smoke/smoke_stage1.yaml" --out $T/jobs1.tsv
cd "$REPO/stage2_privileged_transfer/pipeline"
python train_stage2.py --config "$REPO/cluster/smoke/smoke_stage2.yaml" --dataset sim_transferable_protein --dry_run
tail -n +2 $T/jobs2.tsv | while IFS=$'\t' read -r ds st f s; do
  python train_stage2.py --config "$REPO/cluster/smoke/smoke_stage2.yaml" --dataset $ds --setting $st --fold $f --seed $s \
         --backend numpy --out $T/results/stage2 > /dev/null
done
python "$REPO/shared/uq/aggregate.py" --stage 2 --results $T/results/stage2 --out $T/agg2
cd "$REPO/stage1_cross_modal_pretraining/pipeline"
tail -n +2 $T/jobs1.tsv | while IFS=$'\t' read -r ds mod model arm q mode f s; do
  python run_stage1.py --config "$REPO/cluster/smoke/smoke_stage1.yaml" --dataset $ds --modality $mod --model $model \
         --arm $arm --q $q --mode $mode --fold $f --seed $s --out $T/results/stage1 > /dev/null
done
python "$REPO/shared/uq/aggregate.py" --stage 1 --results $T/results/stage1 --out $T/agg1
echo "SMOKE TEST PASSED: tables in $T/agg1 and $T/agg2"
