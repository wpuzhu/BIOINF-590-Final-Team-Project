#!/usr/bin/env bash
# Run the R methods one image per process (observed + permuted), JOBS processes at a time.
set -euo pipefail
cd "$(dirname "$0")/.."
JOBS=${JOBS:-5}
export R_LIBS="$PWD/.Rlib" R_NCORES=1
ls data/processed/*/*__meta.csv | sed 's#.*/##; s#__meta.csv##' \
  | while read -r img; do for mode in observed permuted; do echo "$img $mode"; done; done \
  | xargs -P "$JOBS" -n 2 sh -c 'ONLY_IMAGE="$0" Rscript src/methods_r.R "$1" 2>&1 | grep -E "done|failed" || true'
