#!/usr/bin/env bash
# Download the three annotated public datasets (~245 MB) and verify checksums.
set -euo pipefail
cd "$(dirname "$0")/../data" && mkdir -p raw && cd raw
get() {  # file url sha256
  if [ ! -f "$1" ] || [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" != "$3" ]; then
    echo "  downloading $1"; curl -sSL -o "$1" "$2"
  fi
  [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$3" ] && echo "  ok $1" || { echo "checksum mismatch: $1"; exit 1; }
}
get imc.h5ad https://exampledata.scverse.org/squidpy/imc.h5ad \
    950c44c785ea86c4262140b0229e0b4f77110a765c3b6874cdb5e0e52973c6fe
get mibitof.h5ad https://exampledata.scverse.org/squidpy/mibitof.h5ad \
    3f125c51695d78ed1c36d5485dc390ab400154d021f0c7715b89f8ee83978690
get schurch_CRC_codex.csv \
    "https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded" \
    416cc3926a7a900ce3b22a33be699535ca35ae85fca4585a8ba5dad6d7f3c677
