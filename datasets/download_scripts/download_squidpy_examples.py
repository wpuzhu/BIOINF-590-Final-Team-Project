"""Download Squidpy's built-in spatial proteomics examples as .h5ad.

Usage (repo root):  python datasets/download_scripts/download_squidpy_examples.py
"""
from pathlib import Path

import squidpy as sq

out = Path("datasets/raw")
out.mkdir(parents=True, exist_ok=True)

imc = sq.datasets.imc()          # IMC breast cancer (Jackson et al. 2020 subset)
imc.write_h5ad(out / "squidpy_imc.h5ad")

mibi = sq.datasets.mibitof()     # MIBI-TOF colorectal (Hartmann et al. 2021)
mibi.write_h5ad(out / "squidpy_mibitof.h5ad")

print(imc, mibi, sep="\n\n")
