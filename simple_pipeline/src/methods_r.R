# Step 2b: methods from the anchor benchmark (zero-shot, package defaults), run in R.
#   ct_STANCE      STANCE two-stage test, stage-2 per-cell-type p (Su et al. 2025); normalized = TRUE path;
#                  run on a stratified ~600-cell subsample per image (runtime)
#   ct_CELINA      CELINA (Shang, Wu, Zhou 2025); its count-based gene QC is bypassed (every marker
#                  is tested for every cell type). The benchmark found that QC drops true signals.
#   ct_sparkx_R    official SPARK-X (xzhoulab/SPARK) within each cell type (cross-checks the Python port)
# Single-cell data -> the cell-type "composition" matrix is one-hot.
# Usage: [ONLY_IMAGE=<image>] Rscript methods_r.R [observed|permuted]   (R_LIBS -> ../.Rlib; see run_pipeline.sh)

suppressPackageStartupMessages({ library(SPARK); library(STANCE); library(CELINA) })
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "observed"
MIN_CELLS <- 40; MAX_CELLS <- 1500; SEED <- 590; NCORES <- as.integer(Sys.getenv("R_NCORES", "4"))
STANCE_CELLS <- 600   # STANCE stage 2 is ~O(n^3) single-core on macOS: stratified small batch per image
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
proc <- file.path(here, "..", "data", "processed"); out <- file.path(here, "..", "results", "raw")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

rows <- list(); timing <- list()
add <- function(dataset, image, ct, feats, method, p, n) {
  rows[[length(rows) + 1]] <<- data.frame(dataset = dataset, image = image, cell_type = ct,
                                          feature = feats, method = method, stat = NA_real_,
                                          p_value = as.numeric(p), n_cells = n)
}

only <- Sys.getenv("ONLY_IMAGE")                    # optional: restrict to one image (parallel runs)
for (meta_f in sort(Sys.glob(file.path(proc, "*", "*__meta.csv")))) {
  if (nzchar(only) && !grepl(only, meta_f, fixed = TRUE)) next
  dataset <- basename(dirname(meta_f)); image <- sub("__meta.csv", "", basename(meta_f))
  meta <- read.csv(meta_f, check.names = FALSE)
  expr <- as.matrix(read.csv(sub("__meta", "__expr", meta_f), row.names = 1, check.names = FALSE))

  # keep cell types with >= MIN_CELLS, cap each at MAX_CELLS (same rule as the Python step)
  tab <- table(meta$cell_type); types <- sort(names(tab)[tab >= MIN_CELLS])
  keep <- unlist(lapply(types, function(t) {
    i <- which(meta$cell_type == t); if (length(i) > MAX_CELLS) sort(sample(i, MAX_CELLS)) else i }))
  meta <- meta[keep, ]; expr <- expr[keep, , drop = FALSE]
  pos <- as.matrix(meta[, c("x", "y")]); rownames(pos) <- meta$cell_id
  if (mode == "permuted") for (t in types) {                     # shuffle locations within type
    i <- which(meta$cell_type == t); pos[i, ] <- pos[sample(i), , drop = FALSE] }
  rownames(expr) <- meta$cell_id
  counts <- t(expr)                                              # markers x cells
  counts <- counts[apply(counts, 1, sd) > 0, , drop = FALSE]
  onehot <- sapply(types, function(t) as.numeric(meta$cell_type == t))
  rownames(onehot) <- meta$cell_id
  n <- table(meta$cell_type)

  # --- STANCE (stratified subsample: each type keeps max(MIN_CELLS, its share of STANCE_CELLS))
  t0 <- Sys.time()
  sidx <- sort(unlist(lapply(types, function(t) {
    i <- which(meta$cell_type == t)
    k <- min(length(i), max(MIN_CELLS, round(STANCE_CELLS * length(i) / nrow(meta))))
    if (k < length(i)) sample(i, k) else i })))
  st <- try({
    o <- creatSTANCEobject(counts = counts[, sidx], pos = pos[sidx, ], prop = onehot[sidx, ], covariates = NULL)
    capture.output(o <- data_preprocess(o, normalized = TRUE))
    capture.output(o <- build_kernelMatrix(o))
    capture.output(o <- runTest1(o, correction = FALSE, pv.adjust = "BY"))
    capture.output(o <- runTest2(o, Genes_to_test = rownames(counts), Cell_types_to_test = NULL,
                                 correction = FALSE, ncores = 1))   # mclapply forking fails on macOS
    o }, silent = TRUE)
  if (!inherits(st, "try-error")) {
    for (t in names(st@Test_2)) {
      r <- st@Test_2[[t]]
      pcol <- grep("^p_value$|p_value", colnames(r), value = TRUE)[1]
      add(dataset, image, t, rownames(r), "ct_STANCE", r[[pcol]], sum(meta$cell_type[sidx] == t))
    }
  } else message("STANCE failed on ", image, ": ", st)
  timing[[length(timing) + 1]] <- data.frame(image = image, method = "ct_STANCE",
                                             seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))

  # --- CELINA (QC bypassed: every marker tested for every type)
  t0 <- Sys.time()
  ce <- try({
    o <- Create_Celina_Object(celltype_mat = t(onehot), gene_expression_mat = counts,
                              location = pos, covariates = NULL, project = image)
    o@genes_list <- setNames(rep(list(rownames(counts)), length(types)), types)
    capture.output(o <- Calculate_Kernel(o, approximation = FALSE))
    capture.output(o <- Testing_interaction_all(o, celltype_to_test = types, num_cores = NCORES))
    o }, silent = TRUE)
  if (!inherits(ce, "try-error")) {
    for (t in names(ce@result)) add(dataset, image, t, rownames(ce@result[[t]]), "ct_CELINA",
                                    ce@result[[t]][["CombinedPvals"]], n[[t]])
  } else message("CELINA failed on ", image, ": ", ce)
  timing[[length(timing) + 1]] <- data.frame(image = image, method = "ct_CELINA",
                                             seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))

  # --- SPARK-X within each cell type
  t0 <- Sys.time()
  for (t in types) {
    i <- which(meta$cell_type == t)
    y <- counts[, i, drop = FALSE]; y <- y - apply(y, 1, min)          # SPARK-X expects non-negative input
    y <- y[apply(y, 1, sd) > 0, , drop = FALSE]
    sx <- try(sparkx(count_in = y, locus_in = pos[i, ], numCores = 1, option = "mixture",
                     verbose = FALSE), silent = TRUE)
    if (!inherits(sx, "try-error"))
      add(dataset, image, t, rownames(sx$res_mtest), "ct_sparkx_R", sx$res_mtest$combinedPval, length(i))
  }
  timing[[length(timing) + 1]] <- data.frame(image = image, method = "ct_sparkx_R",
                                             seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  message(sprintf("  [%s] %s/%s done", mode, dataset, image))
}

tag <- if (nzchar(only)) paste0("__", only) else ""          # one file per image when run in parallel
write.csv(do.call(rbind, rows), file.path(out, paste0("r_", mode, tag, ".csv")), row.names = FALSE)
write.csv(do.call(rbind, timing), file.path(out, paste0("r_timing_", mode, tag, ".csv")), row.names = FALSE)
