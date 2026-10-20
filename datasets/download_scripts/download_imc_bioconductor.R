# Download public IMC / MIBI datasets via Bioconductor.
# Run from the repo root:  Rscript datasets/download_scripts/download_imc_bioconductor.R
# Output: datasets/raw/*.rds  (SingleCellExperiment / SpatialExperiment objects)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("imcdatasets", "Statial", "SpatialExperiment"), update = FALSE)

library(imcdatasets)
library(Statial)

out <- "datasets/raw"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# imcdatasets returns cells / images / masks; we save only the single-cell objects.
# (Pass data_type = "images" or "masks" to get pixel data.)
sets <- list(
  JacksonFischer_2020 = function() JacksonFischer_2020_BreastCancer(data_type = "spe"),
  Damond_2019         = function() Damond_2019_Pancreas(data_type = "spe"),
  HochSchulz_2022     = function() HochSchulz_2022_Melanoma(data_type = "spe"),
  IMMUcan_2022        = function() IMMUcan_2022_CancerExample(data_type = "spe")
)
for (nm in names(sets)) {
  message("Downloading ", nm, " ...")
  obj <- tryCatch(sets[[nm]](), error = function(e) { message("  failed: ", e$message); NULL })
  if (!is.null(obj)) saveRDS(obj, file.path(out, paste0(nm, "_spe.rds")))
}

# Keren et al. 2018 TNBC MIBI-TOF subset (10 images, 36 proteins, with survival)
data("kerenSCE", package = "Statial")
saveRDS(kerenSCE, file.path(out, "Keren_2018_kerenSCE.rds"))

message("Done. See ?imcdatasets::listDatasets for all available datasets.")
