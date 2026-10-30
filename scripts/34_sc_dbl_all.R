# =============================================================================
# 34_sc_dbl_all.R
# Independent doublet validation for ALL samples in ONE R process.
#
# Why a single process: this machine only has ~5 GB free RAM. Launching several
# R sessions that each load Seurat + a 36k x 9k sparse matrix exhausts memory and
# the processes die with a Windows access violation (0xc0000005). One process,
# one sample at a time, is the stable configuration.
#
# Logic is identical to 31_sc_doublet_validate.R --mode=sample:
#   scDblFinder on the QC-passing raw-count subset, dbr from total recovered
#   nuclei, agreement computed against the source-supplied predictions.
#
# Resumable: samples whose _doublet_summary.csv already exists are skipped.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

proj_dir <- "E:/another data file/GEObulk analysis project"
out_dir  <- file.path(proj_dir, "results", "Singlecel")
obj_dir  <- file.path(out_dir, "objects")
table_dir <- file.path(out_dir, "tables")
dbl_dir  <- file.path(table_dir, "doublet_validation")
log_dir  <- file.path(out_dir, "logs")
for (d in c(dbl_dir, log_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

manifest <- read.csv(file.path(proj_dir, "clean_data", "GSE292067_sample_manifest.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)

# Optional batching controls (used when a single R process must stay short):
#   SC_DBL_ONLY=GSM8848306   -> process just this sample
#   SC_DBL_MAX=3             -> process at most N pending samples
n_all <- nrow(manifest)
full_manifest <- manifest
only <- Sys.getenv("SC_DBL_ONLY", "")
if (nzchar(only)) manifest <- manifest[manifest$sample == only, , drop = FALSE]
max_n <- suppressWarnings(as.integer(Sys.getenv("SC_DBL_MAX", "")))
if (is.finite(max_n) && max_n > 0L) {
  done <- file.exists(file.path(table_dir, "doublet_validation",
                                paste0(manifest$sample, "_doublet_summary.csv")))
  manifest <- manifest[order(done), , drop = FALSE]
  manifest <- manifest[seq_len(min(max_n, nrow(manifest))), , drop = FALSE]
}

# Optional sharding so several R processes can share the work without overlap:
#   Rscript 34_sc_dbl_all.R 2/3    -> process every 3rd sample starting at 2
shard <- Sys.getenv("SC_DBL_SHARD", "")
cli <- commandArgs(trailingOnly = TRUE)
if (!length(shard) && length(cli) && grepl("^[0-9]+/[0-9]+$", cli[1])) shard <- cli[1]
if (nzchar(shard)) {
  parts <- as.integer(strsplit(shard, "/")[[1]])
  k <- parts[1]; n <- parts[2]
  if (is.na(k) || is.na(n) || k < 1L || k > n) stop("bad shard spec: ", shard)
  manifest <- manifest[seq.int(k, nrow(manifest), by = n), , drop = FALSE]
}

prog <- file.path(log_dir, "doublet_validate_progress.log")
say <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = prog, append = TRUE)
}
say("scDblFinder ", as.character(packageVersion("scDblFinder")),
    " | R ", as.character(getRversion()), " | samples: ", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  smp <- manifest$sample[i]
  summ_file <- file.path(dbl_dir, paste0(smp, "_doublet_summary.csv"))
  if (file.exists(summ_file)) {
    say(sprintf("[%2d/%d] %s  SKIP (already done)", i, nrow(manifest), smp))
    next
  }
  ok <- tryCatch({
    t0 <- Sys.time()
    obj <- readRDS(file.path(obj_dir, paste0(smp, "_unfiltered_qc.rds")))
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "CsparseMatrix")

    n_total <- ncol(obj)
    qc_pass <- !obj$qc_candidate_remove
    counts_qc <- counts[, qc_pass, drop = FALSE]
    # colnames(obj) carry a "<sample>_" prefix; barcode_raw holds the original
    # 10x barcode ("AAAC...-1"). They are positionally aligned, so take the
    # barcode column directly instead of matching cell names against barcodes.
    barcodes_qc <- obj$barcode_raw[qc_pass]
    if (anyNA(barcodes_qc)) stop("NA in barcode_raw for QC-passing cells")
    if (ncol(counts_qc) < 200) stop("too few QC-passing cells: ", ncol(counts_qc))

    # Loading rate is a property of the capture -> use total recovered nuclei.
    dbr <- 0.008 * n_total / 1000

    set.seed(1234)
    sce <- SingleCellExperiment(assays = list(counts = counts_qc))
    colData(sce)$sample_id <- smp
    sce <- scDblFinder(sce, dbr = dbr, verbose = FALSE)

    calls <- data.frame(
      sample = smp,
      barcode_raw = barcodes_qc,
      scDblFinder_score = as.numeric(sce$scDblFinder.score),
      scDblFinder_class = as.character(sce$scDblFinder.class),
      stringsAsFactors = FALSE
    )
    write.csv(calls, file.path(dbl_dir, paste0(smp, "_scDblFinder.csv")), row.names = FALSE)

    src_idx <- match(barcodes_qc, obj$barcode_raw)
    src_dbl <- as.character(obj$doublet_prediction_source)[src_idx] == "doublet"
    new_dbl <- calls$scDblFinder_class == "doublet"
    tp <- sum(src_dbl & new_dbl)
    jaccard <- if (sum(src_dbl | new_dbl) > 0) tp / sum(src_dbl | new_dbl) else NA_real_

    summ <- data.frame(
      sample = smp,
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
    write.csv(summ, file.path(dbl_dir, paste0(smp, "_doublet_summary.csv")), row.names = FALSE)

    all_calls <- data.frame(
      cell = colnames(obj),
      barcode_raw = obj$barcode_raw,
      sample = smp,
      qc_pass = as.logical(qc_pass),
      scDblFinder_class = "not_evaluated",
      scDblFinder_score = NA_real_,
      stringsAsFactors = FALSE
    )
    m <- match(all_calls$barcode_raw[qc_pass], calls$barcode_raw)
    all_calls$scDblFinder_class[qc_pass] <- calls$scDblFinder_class[m]
    all_calls$scDblFinder_score[qc_pass] <- calls$scDblFinder_score[m]
    con <- gzfile(file.path(dbl_dir, paste0(smp, "_doublet_calls.csv.gz")), "wt")
    write.csv(all_calls, con, row.names = FALSE)
    close(con)

    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    say(sprintf("[%2d/%d] %s  n=%d (qc-pass %d)  scDbl=%d (%.2f%%)  src=%d (%.2f%%)  Jaccard=%.3f  %.0fs",
                i, nrow(manifest), smp, n_total, ncol(counts_qc),
                sum(new_dbl), 100 * sum(new_dbl) / ncol(counts_qc),
                sum(src_dbl), 100 * sum(src_dbl) / ncol(counts_qc),
                jaccard, dt))

    rm(obj, counts, counts_qc, sce, calls); invisible(gc())
    TRUE
  }, error = function(e) {
    say(sprintf("[%2d/%d] %s  ERROR: %s", i, nrow(manifest), smp, conditionMessage(e)))
    FALSE
  })
  if (!ok) next
}

# ---- combined summary ----
files <- list.files(dbl_dir, pattern = "_doublet_summary[.]csv$", full.names = TRUE)
say("per-sample summaries found: ", length(files), " / ", n_all)
if (length(files) == n_all) {
  s <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  s <- s[match(full_manifest$sample, s$sample), , drop = FALSE]
  s$condition <- full_manifest$condition[match(s$sample, full_manifest$sample)]
  write.csv(s, file.path(table_dir, "doublet_validation_summary.csv"), row.names = FALSE)
  say("COMBINED rows: ", nrow(s))
  say(sprintf("QC-passing cells: %d", sum(s$n_qc_pass)))
  say(sprintf("scDblFinder doublets: %d (%.2f%%)", sum(s$n_scDblFinder_doublet),
              100 * sum(s$n_scDblFinder_doublet) / sum(s$n_qc_pass)))
  say(sprintf("source doublets (same cells): %d (%.2f%%)", sum(s$n_source_doublet_on_qc_cells),
              100 * sum(s$n_source_doublet_on_qc_cells) / sum(s$n_qc_pass)))
  say(sprintf("median Jaccard: %.3f | median agreement: %.2f%%",
              median(s$jaccard, na.rm = TRUE), median(s$agreement_pct)))
} else {
  say("NOT all samples finished; combined table not written yet.")
}
say("END")
