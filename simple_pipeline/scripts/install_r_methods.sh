#!/usr/bin/env bash
# Install SPARK, STANCE and CELINA into a project-local library (simple_pipeline/.Rlib).
# Nothing is written to your global R library.
# Notes for macOS without the gfortran toolchain: an empty FLIBS (via .Rmakevars) lets SPARK link,
# and gaston must be the CRAN binary (a source build fails to link RcppParallel).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p .Rlib && printf 'FLIBS=\n' > .Rmakevars
export R_LIBS="$PWD/.Rlib" R_MAKEVARS_USER="$PWD/.Rmakevars"
Rscript -e '
lib <- normalizePath(".Rlib"); options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
for (p in c("remotes", "BiocManager")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib)
bin <- if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") "binary" else getOption("pkgType")
install.packages(c("gaston", "RcppParallel"), lib = lib, type = bin)
BiocManager::install(c("scater", "SingleCellExperiment"), lib = lib, update = FALSE, ask = FALSE)
remotes::install_github("xzhoulab/SPARK", lib = lib, dependencies = TRUE, upgrade = "never")
remotes::install_github("Cui-STT-Lab/STANCE", lib = lib, dependencies = TRUE, upgrade = "never")
remotes::install_github("pekjoonwu/CELINA", lib = lib, dependencies = TRUE, upgrade = "never")
print(sapply(c("SPARK", "STANCE", "CELINA"), requireNamespace, quietly = TRUE))'
