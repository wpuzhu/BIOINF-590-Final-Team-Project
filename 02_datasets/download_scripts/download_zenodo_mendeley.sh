#!/usr/bin/env bash
# Fetch file listings / archives from Zenodo deposits.
# Usage (repo root): bash 02_datasets/download_scripts/download_zenodo_mendeley.sh [--list-only]
# Files are large (several GB). Run with --list-only first to see what each record has.
set -euo pipefail
OUT=02_datasets/raw/zenodo
mkdir -p "$OUT"

# name:zenodo_record  (plain list so it works with macOS bash 3.2)
RECORDS="jackson2020_imc:3518284 danenberg2022_imc:5850952 dbitplus:16827046"
#        Jackson & Fischer 2020 IMC | Danenberg 2022 METABRIC IMC | DBiTplus (same-section RNA + CODEX)

for entry in $RECORDS; do
  name=${entry%%:*}; rec=${entry##*:}
  echo "== $name (zenodo $rec)"
  urls=$(curl -s "https://zenodo.org/api/records/$rec" \
    | python3 -c 'import sys,json; [print(f["links"]["self"], f["key"], f["size"]) for f in json.load(sys.stdin)["files"]]')
  echo "$urls"
  if [[ "${1:-}" != "--list-only" ]]; then
    mkdir -p "$OUT/$name"
    while read -r url key size; do
      curl -L -C - -o "$OUT/$name/$key" "$url"
    done <<< "$urls"
  fi
done

cat <<'MSG'

Mendeley Data deposits need a browser download (click "Download All"):
  Schurch 2020 CRC CODEX : https://data.mendeley.com/datasets/mpjzbtfgfr/1
Put the unzipped folders under 02_datasets/raw/mendeley/<name>/
MSG
