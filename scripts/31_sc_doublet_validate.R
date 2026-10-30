# =============================================================================
# 31_sc_doublet_validate.R
# Independent doublet validation with scDblFinder, per sample, on raw counts.
# Runs on the QC-passing subset (proposed MAD gate) but writes NO filtered object.
#
# Modes:
#   --mode=sample --sample=GSM8848306
#   --mode=finalize
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}
mode <- arg_value("--mode", "sample")
sample_arg <- arg_value("--sample", NA_character_)

proj_dir <- "E:/another data file/GEObulk analysis project"
out_dir <- file.path(proj_dir, "results", "Singlecel")
obj_dir <- file.path(out_dir, "objects")
table_dir <- file.path(out_dir, "tables")
dbl_dir <- file.path(table_dir, "doublet_validation")
log_dir <- file.path(out_dir, "logs")
for (d in c(dbl_dir, log_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

manifest_path <- file.path(proj_dir, "clean_data", "GSE292067_sample_manifest.csv")
manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)

if (identical(mode, "sample")) {
  if (is.na(sample_arg)) stop("--sample is required")
  hit <- which(manifest$sample == sample_arg)
  if (length(hit) != 1L) stop("sample not found: ", sample_arg)
  row <- manifest[hit, , drop = FALSE]

  obj <- readRDS(file.path(obj_dir, paste0(sample_arg, "_unfiltered_qc.rds")))
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "CsparseMatrix")

  n_total <- ncol(obj)
  qc_pass <- !obj$qc_candidate_remove
  counts_qc <- counts[, qc_pass, drop = FALSE]
  barcodes_qc <- colnames(counts_qc)
  if (ncol(counts_qc) < 200) stop("Too few QC-passing cells for ", sample_arg)

  # Physical loading rate is a property of the capture, so derive it from the
  # total recovered nuclei, not from the post-QC subset.
  dbr <- 0.008 * n_total / 1000

  sce <- SingleCellExperiment(assays = list(counts = counts_qc))
  colData(sce)$sample_id <- sample_arg
  sce <- scDblFinder(sce, dbr = dbr, verbose = FALSE)

  calls <- data.frame(
    sample = sample_arg,
    barcode_raw = barcodes_qc,
    scDblFinder_score = as.numeric(sce$scDblFinder.score),
    scDblFinder_class = as.character(sce$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  write.csv(calls, file.path(dbl_dir, paste0(sample_arg, "_scDblFinder.csv")), row.names = FALSE)

  # --- agreement with the source-supplied predictions, on the same cells ---
  src_idx <- match(barcodes_qc, obj$barcode_raw)
  src_dbl <- as.character(obj$doublet_prediction_source)[src_idx] == "doublet"
  new_dbl <- calls$scDblFinder_class == "doublet"
  tab <- table(source = src_dbl, scDblFinder = new_dbl)
  tp <- sum(src_dbl & new_dbl)
  jaccard <- if (sum(src_dbl | new_dbl) > 0) tp / sum(src_dbl | new_dbl) else NA_real_

  summ <- data.frame(
    sample = sample_arg,
    n_total = n_total,
    n_qc_pass = ncol(counts_qc),
    dbr_used = dbr,
    scDblFinder_version = as.character(packageVersion("scDblFinder")),
    n_scDblFinder_doublet = sum(new_dbl),
    pct_scDblFinder_doublet = 100 * sum(new_dbl) / ncol(counts_qc),
    n_source_doublet_on_qc_cells = sum(src_dbl),
    pct_source_doublet_on_qc_cells = 100 * sum(src_dbl) / ncol(counts_qc),
    concordant_doublet = tp,
    discordant_scDblFinder_only = sum(!src_dbl & new_dbl),
    discordant_source_only = sum(src_dbl & !new_dbl),
    both_singlet = sum(!src_dbl & !new_dbl),
    jaccard = jaccard,
    agreement_pct = 100 * (tp + sum(!src_dbl & !new_dbl)) / length(src_dbl),
    stringsAsFactors = FALSE
  )
  write.csv(summ, file.path(dbl_dir, paste0(sample_arg, "_doublet_summary.csv")), row.names = FALSE)

  # per-cell call table for the whole object (QC-failing cells marked not_evaluated)
  all_calls <- data.frame(
    cell = colnames(obj),
    barcode_raw = obj$barcode_raw,
    sample = sample_arg,
    qc_pass = as.logical(qc_pass),
    scDblFinder_class = "not_evaluated",
    scDblFinder_score = NA_real_,
    stringsAsFactors = FALSE
  )
  m <- match(all_calls$barcode_raw[qc_pass], calls$barcode_raw)
  all_calls$scDblFinder_class[qc_pass] <- calls$scDblFinder_class[m]
  all_calls$scDblFinder_score[qc_pass] <- calls$scDblFinder_score[m]
  con <- gzfile(file.path(dbl_dir, paste0(sample_arg, "_doublet_calls.csv.gz")), "wt")
  write.csv(all_calls, con, row.names = FALSE)
  close(con)

  writeLines(c(
    paste0("sample=", sample_arg),
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("n_total=", n_total),
    paste0("n_qc_pass=", ncol(counts_qc)),
    paste0("dbr=", signif(dbr, 5)),
    paste0("scDblFinder_doublets=", sum(new_dbl)),
    paste0("source_doublets_on_same_cells=", sum(src_dbl)),
    paste0("jaccard=", signif(jaccard, 4)),
    "object_modified=FALSE"
  ), file.path(log_dir, paste0(sample_arg, "_scDblFinder.log")), useBytes = TRUE)

  rm(obj, counts, counts_qc, sce); invisible(gc())
  cat("done", sample_arg, "\n")
}

if (identical(mode, "finalize")) {
  files <- list.files(dbl_dir, pattern = "_doublet_summary[.]csv$", full.names = TRUE)
  if (length(files) != nrow(manifest)) stop("expected ", nrow(manifest), " summaries, got ", length(files))
  s <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  s <- s[match(manifest$sample, s$sample), , drop = FALSE]
  s$condition <- manifest$condition[match(s$sample, manifest$sample)]
  write.csv(s, file.path(table_dir, "doublet_validation_summary.csv"), row.names = FALSE)
  cat("samples:", nrow(s), "\n")
  cat("total QC-passing cells:", sum(s$n_qc_pass), "\n")
  cat("scDblFinder doublets:", sum(s$n_scDblFinder_doublet),
      sprintf("(%.2f%%)\n", 100 * sum(s$n_scDblFinder_doublet) / sum(s$n_qc_pass)))
  cat("source doublets on same cells:", sum(s$n_source_doublet_on_qc_cells),
      sprintf("(%.2f%%)\n", 100 * sum(s$n_source_doublet_on_qc_cells) / sum(s$n_qc_pass)))
  cat("median jaccard:", sprintf("%.3f\n", median(s$jaccard, na.rm = TRUE)))
  cat("median agreement %:", sprintf("%.2f\n", median(s$agreement_pct)))
}
