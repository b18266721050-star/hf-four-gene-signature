# =============================================================================
# 33_sc_doublet_doubletfinder.R
# Independent doublet validation with DoubletFinder (Seurat-native), because
# SingleCellExperiment/scDblFinder crash R on this installation.
# Per sample, on raw counts, QC-passing subset only; writes no filtered object.
#
# Modes: --mode=sample --sample=GSMxxxx | --mode=finalize
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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

prog_file <- file.path(log_dir, "doublet_progress.log")
prog <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n",
                          file = prog_file, append = TRUE)

manifest <- read.csv(file.path(proj_dir, "clean_data", "GSE292067_sample_manifest.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(DoubletFinder))

if (identical(mode, "sample")) {
  s <- sample_arg
  row <- manifest[manifest$sample == s, , drop = FALSE]
  prog("START", s)

  obj <- readRDS(file.path(obj_dir, paste0(s, "_unfiltered_qc.rds")))
  prog("loaded", s, ncol(obj), "cells")

  n_total <- ncol(obj)
  keep <- !obj$qc_candidate_remove
  sub <- subset(obj, cells = colnames(obj)[keep])
  rm(obj); invisible(gc())
  prog("subset", s, ncol(sub))

  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.5, verbose = FALSE)
  prog("preprocessed", s, ncol(sub))

  sweep <- paramSweep(sub, PCs = 1:20, sct = FALSE)
  bc <- summarizeSweep(sweep)
  pk <- find.pK(bc)
  pk$BCmetric <- as.numeric(pk$BCmetric)
  best_pk <- as.numeric(as.character(pk$pK[which.max(pk$BCmetric)]))
  prog("pK", s, best_pk)

  rate <- 0.008 * n_total / 1000
  nExp <- round(rate * ncol(sub))
  homotypic <- modelHomotypic(sub$seurat_clusters)
  nExp_adj <- round(nExp * (1 - homotypic))
  if (nExp_adj < 1) nExp_adj <- 1
  prog("nExp", s, nExp, "adj", nExp_adj, "homotypic", round(homotypic, 4))

  sub <- doubletFinder(sub, PCs = 1:20, pN = 0.25, pK = best_pk,
                       nExp = nExp_adj, sct = FALSE, verbose = FALSE)

  df_col <- grep("^DF.classifications", colnames(sub@meta.data), value = TRUE)
  if (length(df_col) != 1L) stop("cannot locate DF classification column for ", s)
  cls <- as.character(sub@meta.data[[df_col]])
  is_dbl <- cls == "Doublet"

  calls <- data.frame(
    sample = s,
    barcode_raw = sub$barcode_raw,
    DoubletFinder_class = cls,
    DoubletFinder_pANN = as.numeric(sub$pANN_0.25)[seq_along(cls)],
    stringsAsFactors = FALSE
  )
  # pANN column name carries pK/nExp; grab by pattern instead
  ann_col <- grep("^pANN_", colnames(sub@meta.data), value = TRUE)
  if (length(ann_col) == 1L) calls$DoubletFinder_pANN <- as.numeric(sub@meta.data[[ann_col]])
  write.csv(calls, file.path(dbl_dir, paste0(s, "_DoubletFinder.csv")), row.names = FALSE)

  src_idx <- match(calls$barcode_raw, read.csv(gzfile(file.path(table_dir, "cell_qc",
              paste0(s, "_cell_qc.csv.gz"))), stringsAsFactors = FALSE)$barcode_raw)
  qc <- read.csv(gzfile(file.path(table_dir, "cell_qc", paste0(s, "_cell_qc.csv.gz"))),
                 stringsAsFactors = FALSE, check.names = FALSE)
  src <- qc$doublet_prediction_source[src_idx] == "doublet"
  tp <- sum(src & is_dbl)
  jac <- if (sum(src | is_dbl) > 0) tp / sum(src | is_dbl) else NA_real_

  summ <- data.frame(
    sample = s, n_total = n_total, n_evaluated = ncol(sub),
    pK = best_pk, expected_rate = rate, nExp_raw = nExp, nExp_adj = nExp_adj,
    homotypic_fraction = homotypic,
    n_DoubletFinder = sum(is_dbl), pct_DoubletFinder = 100 * mean(is_dbl),
    n_source_doublet = sum(src), pct_source_doublet = 100 * mean(src),
    concordant = tp,
    DF_only = sum(!src & is_dbl), source_only = sum(src & !is_dbl),
    jaccard = jac,
    agreement_pct = 100 * mean(src == is_dbl),
    stringsAsFactors = FALSE
  )
  write.csv(summ, file.path(dbl_dir, paste0(s, "_DoubletFinder_summary.csv")), row.names = FALSE)
  prog("DONE", s, sum(is_dbl), "doublets")
  rm(sub); invisible(gc())
}

if (identical(mode, "finalize")) {
  files <- list.files(dbl_dir, pattern = "_DoubletFinder_summary[.]csv$", full.names = TRUE)
  if (!length(files)) stop("no summaries yet")
  s <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  s <- s[match(manifest$sample, s$sample), , drop = FALSE]
  s$condition <- manifest$condition[match(s$sample, manifest$sample)]
  write.csv(s, file.path(table_dir, "doublet_validation_summary.csv"), row.names = FALSE)
  cat("samples:", sum(!is.na(s$sample)), "/", nrow(manifest), "\n")
  cat("evaluated cells:", sum(s$n_evaluated, na.rm = TRUE), "\n")
  cat(sprintf("DoubletFinder doublets: %d (%.2f%%)\n",
              sum(s$n_DoubletFinder, na.rm = TRUE),
              100 * sum(s$n_DoubletFinder, na.rm = TRUE) / sum(s$n_evaluated, na.rm = TRUE)))
  cat(sprintf("source doublets on same cells: %d (%.2f%%)\n",
              sum(s$n_source_doublet, na.rm = TRUE),
              100 * sum(s$n_source_doublet, na.rm = TRUE) / sum(s$n_evaluated, na.rm = TRUE)))
  cat(sprintf("median jaccard: %.3f | median agreement: %.2f%%\n",
              median(s$jaccard, na.rm = TRUE), median(s$agreement_pct, na.rm = TRUE)))
}
