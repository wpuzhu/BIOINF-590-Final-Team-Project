"""Pipeline settings. "Zero-shot" = every method runs with fixed default parameters:
no tuning, no training, no use of the validation annotations while detecting."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"

SEED = 590
MIN_CELLS_PER_TYPE = 40      # skip (image, cell type) groups smaller than this
MAX_CELLS_PER_TYPE = 1500    # subsample larger groups (keeps kernel methods O(n^2) small)
KNN = 10                     # spatial neighbours for Moran's I and niche composition
FDR = 0.05                   # BH level used for every "significant" call
LINEAGE_DELTA = 1.0          # z-score gap that makes a marker "lineage" for a cell type (validation only)
MAX_LINEAGE_PER_TYPE = 3
NATIVE_GAP = 1.5             # marker is 'native' to a type if within 1.5 SD of the best-expressing type (validation only)

# Small batch: a handful of images per dataset
DATASETS = {
    "schurch_codex": dict(
        file="schurch_CRC_codex.csv",
        url="https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/"
            "c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded",
        sha256="416cc3926a7a900ce3b22a33be699535ca35ae85fca4585a8ba5dad6d7f3c677",
        n_images_per_group=3,          # 3 CLR + 3 DII TMA spots
        drop_types=["dirt", "undefined"],
    ),
    "jackson_imc": dict(
        file="imc.h5ad",
        url="https://exampledata.scverse.org/squidpy/imc.h5ad",
        sha256="950c44c785ea86c4262140b0229e0b4f77110a765c3b6874cdb5e0e52973c6fe",
    ),
    "hartmann_mibi": dict(
        file="mibitof.h5ad",
        url="https://exampledata.scverse.org/squidpy/mibitof.h5ad",
        sha256="3f125c51695d78ed1c36d5485dc390ab400154d021f0c7715b89f8ee83978690",
    ),
}

# Coarse harmonisation, used only to compare findings across datasets
COARSE_TYPE = [
    ("tumor|epithel|ck", "tumor/epithelial"),
    ("cd8|cytotoxic", "CD8 T"),
    ("treg", "CD4 T"),
    ("cd4|t cell|cd3", "CD4 T"),
    ("macroph|myeloid|monocyt|cd68|cd163|dc", "myeloid"),
    ("b cell|plasma", "B/plasma"),
    ("stroma|fibro|smooth|elongated|vimentin", "stromal"),
    ("endothel|vascul|lymphatic", "endothelial"),
    ("granulocyte|nk", "other immune"),
]
