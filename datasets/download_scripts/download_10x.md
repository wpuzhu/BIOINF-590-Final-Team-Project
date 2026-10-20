# Downloading 10x Genomics datasets

10x asks for a short form before showing download links. After you submit it, each dataset page lists `curl`/`wget` commands. Copy them into your terminal from `datasets/raw/10x/<dataset>/`.

## Visium CytAssist Gene + Protein (spatial protein **counts**)
Human tonsil, FFPE, 35-antibody panel (CC BY 4.0):
- https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard
- https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard

Files needed: `*_filtered_feature_bc_matrix.h5` (contains both "Gene Expression" and "Antibody Capture" features) and `*_spatial.tar.gz` (coordinates and images).

Load in Python:
```python
import scanpy as sc
ad = sc.read_visium("datasets/raw/10x/tonsil", count_file="filtered_feature_bc_matrix.h5", gex_only=False)
prot = ad[:, ad.var["feature_types"] == "Antibody Capture"].copy()
rna  = ad[:, ad.var["feature_types"] == "Gene Expression"].copy()
```
Load in R:
```r
x <- Seurat::Read10X_h5("filtered_feature_bc_matrix.h5")   # list: `Gene Expression`, `Antibody Capture`
```

## Xenium (source data of the Li et al. simulation)
- Human breast: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
- Others (ovarian, lymph node, lung, pancreas, mouse colon): https://www.10xgenomics.com/datasets (filter Product = Xenium)

## Visium HD (ctSVG paper)
See links in `datasets/README.md` §D16.
