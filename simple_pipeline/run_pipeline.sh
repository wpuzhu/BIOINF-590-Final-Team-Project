#!/usr/bin/env bash
# Simple small-batch zero-shot ctSVP pipeline: download -> prep -> detect (Py + R) -> aggregate -> validate -> report
# Usage:  bash run_pipeline.sh            (full run, ~30-40 min on a laptop; R methods are the slow part)
#         SKIP_R=1 bash run_pipeline.sh   (Python methods only, ~1 min)
set -euo pipefail
cd "$(dirname "$0")"
PY=.venv/bin/python

echo "== 0. environment"
[ -x "$PY" ] || { python3 -m venv --system-site-packages .venv && $PY -m pip install -q -r requirements.txt; }
if [ -z "${SKIP_R:-}" ] && [ ! -d .Rlib/CELINA ]; then bash scripts/install_r_methods.sh; fi

echo "== 1. download + preprocess"
bash scripts/download_data.sh
(cd src && ../$PY prep.py)

echo "== 2a. Python methods (observed + within-type permuted null)"
(cd src && ../$PY methods_py.py observed permuted)
if [ -z "${SKIP_R:-}" ]; then
  echo "== 2b. R methods from the benchmark (STANCE, CELINA, SPARK-X)"
  rm -f results/raw/r_*.csv
  bash src/run_r_parallel.sh
fi

echo "== 3-5. aggregate, validate, report"
(cd src && ../$PY aggregate.py observed permuted && ../$PY validate.py && ../$PY report.py)
echo "Done -> results/SUMMARY.md"
