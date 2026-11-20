#!/usr/bin/env bash
# Fetch file listings / archives from Zenodo deposits.
# Usage: bash datasets/download_scripts/download_zenodo_mendeley.sh [OUT_DIR] [--list-only]
# On Great Lakes use OUT_DIR=/scratch/$SLURM_ACCOUNT/$USER/raw (never the repo).
# Files are large (several GB). Run with --list-only first to see what each record has.
set -euo pipefail
OUT=${1:-raw}/zenodo
[[ "${1:-}" == "--list-only" ]] && OUT=raw/zenodo
mkdir -p "$OUT"

# name:zenodo_record  (plain list so it works with macOS bash 3.2)
RECORDS="jackson2020_imc:3518284 danenberg2022_imc:5850952 dbitplus:16827046 xenium_comet_2026:21713660"
#        Jackson & Fischer 2020 IMC | Danenberg 2022 METABRIC IMC | DBiTplus (same-section RNA + CODEX)
#        | Duchini et al. 2026 Xenium + COMET + H&E same section (CC BY 4.0, 3.9 GB)

for entry in $RECORDS; do
  name=${entry%%:*}; rec=${entry##*:}
  echo "== $name (zenodo $rec)"
  urls=$(curl -s "https://zenodo.org/api/records/$rec" \
    | python3 -c 'import sys,json; [print(f["links"]["self"], f["key"], f["size"]) for f in json.load(sys.stdin)["files"]]')
  echo "$urls"
  if [[ "${*: -1}" != "--list-only" ]]; then
    mkdir -p "$OUT/$name"
    while read -r url key size; do
      curl -L -C - -o "$OUT/$name/$key" "$url"
    done <<< "$urls"
  fi
done

cat <<'MSG'

Mendeley Data deposits need a browser download (click "Download All"):
  Schurch 2020 CRC CODEX : https://data.mendeley.com/datasets/mpjzbtfgfr/1
Put the unzipped folders under $OUT_DIR/mendeley/<name>/
MSG
