# =============================================================================
# 40_sc_stage2_filter.R
# Stage 2 - Step 1: apply the agreed per-sample thresholds + doublet removal,
# keep RAW COUNTS only, and write one filtered object per sample.
#
# Filtering rule (final, agreed 2026-10-28):
#   QC gate : nFeature / nCount / percent.mt within the per-sample 3-MAD bounds
#             clamped to biological floor/ceiling (from qc_threshold_final_per_sample.csv)
#   Doublets: UNION of scDblFinder calls and the source-supplied predictions
#             (19,149 + 2,534 = 21,683 cells over the whole cohort)
#
# Resumable: samples whose filtered object already exists are skipped.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
set.seed(20261027)

proj <- "E:/another data file/GEObulk analysis project"
sdir <- file.path(proj, "results", "Singlecel")
obj_in <- file.path(sdir, "objects")
obj_out <- file.path(sdir, "objects", "02_filtered")
tab_dir <- file.path(sdir, "tables")
dbl_dir <- file.path(tab_dir, "doublet_validation")
log_dir <- file.path(sdir, "logs")
for (d in c(obj_out, log_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file(file.path(log_dir, "40_filter_progress.log"), "wt")
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf); flush(lgf)
}

thr <- read.csv(file.path(tab_dir, "qc_threshold_final_per_sample.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
lg("threshold rows: ", nrow(thr), " | rule: 3-MAD clamped | doublet: union(scDblFinder, source)")

rows <- list()
for (i in seq_len(nrow(thr))) {
  smp <- thr$sample[i]
  fout <- file.path(obj_out, paste0(smp, "_filtered.rds"))
  if (file.exists(fout)) { lg(sprintf("[%2d/%d] %s SKIP (done)", i, nrow(thr), smp)); next }

  t0 <- Sys.time()
  ok <- tryCatch({
    obj <- readRDS(file.path(obj_in, paste0(smp, "_unfiltered_qc.rds")))
    n_before <- ncol(obj)

    ## ---- QC gate recomputed from the FINAL clamped thresholds ----
    keep_qc <- obj$nFeature_RNA >= thr$nFeature_min[i] &
               obj$nFeature_RNA <= thr$nFeature_max[i] &
               obj$nCount_RNA   >= thr$nCount_min[i]   &
               obj$nCount_RNA   <= thr$nCount_max[i]   &
               obj$percent.mt   <= thr$percent_mt_max[i]
    keep_qc[is.na(keep_qc)] <- FALSE
    n_qc_rm <- sum(!keep_qc)
    if (!is.na(thr$qc_flagged[i]) && n_qc_rm != thr$qc_flagged[i]) {
      lg(sprintf("  !! %s QC mismatch: recomputed %d vs table %d", smp, n_qc_rm, thr$qc_flagged[i]))
    }

    ## ---- doublets: union of scDblFinder and source predictions ----
    calls <- read.csv(gzfile(file.path(dbl_dir, paste0(smp, "_doublet_calls.csv.gz"))),
                      stringsAsFactors = FALSE)
    m <- match(obj$barcode_raw, calls$barcode_raw)
    if (anyNA(m)) stop("barcode mapping failed for ", smp)
    sc_dbl <- calls$scDblFinder_class[m] == "doublet"
    sc_dbl[is.na(sc_dbl)] <- FALSE
    src_dbl <- as.character(obj$doublet_prediction_source) == "doublet"
    src_dbl[is.na(src_dbl)] <- FALSE
    is_dbl <- sc_dbl | src_dbl

    keep <- keep_qc & !is_dbl

    ## ---- subset, keep counts only ----
    so <- subset(obj, cells = colnames(obj)[keep])
    lyr <- Layers(so[["RNA"]])
    if (!identical(lyr, "counts")) {
      for (l in setdiff(lyr, "counts")) LayerData(so[["RNA"]], layer = l) <- NULL
    }
    stopifnot(identical(Layers(so[["RNA"]]), "counts"))
    so$filter_status <- "retained"
    so <- DietSeurat(so, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
    saveRDS(so, fout)

    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    lg(sprintf("[%2d/%d] %s  %d -> %d (QC -%d, dbl -%d, both %d)  %.0fs",
               i, nrow(thr), smp, n_before, ncol(so), n_qc_rm,
               sum(!keep_qc & is_dbl), sum(!keep_qc & is_dbl), dt))

    rm(obj, so, calls); invisible(gc())
    data.frame(sample = smp, condition = thr$condition[i],
               n_before = n_before,
               n_removed_qc_only = sum(!keep_qc & !is_dbl),
               n_removed_doublet_only = sum(keep_qc & is_dbl),
               n_removed_both = sum(!keep_qc & is_dbl),
               n_retained = sum(keep),
               pct_retained = 100 * sum(keep) / n_before,
               scDblFinder_doublet = sum(sc_dbl),
               source_doublet = sum(src_dbl),
               union_doublet = sum(is_dbl),
               stringsAsFactors = FALSE)
  }, error = function(e) { lg(sprintf("[%2d/%d] %s ERROR: %s", i, nrow(thr), smp, conditionMessage(e))); NULL })
  if (!is.null(ok)) rows[[length(rows) + 1]] <- ok
}

if (length(rows)) {
  s <- do.call(rbind, rows)
  old <- file.path(tab_dir, "filtering_summary.csv")
  if (file.exists(old)) {
    prev <- read.csv(old, stringsAsFactors = FALSE)
    s <- rbind(prev[!prev$sample %in% s$sample, , drop = FALSE], s)
  }
  s <- s[match(thr$sample, s$sample), , drop = FALSE]
  write.csv(s, old, row.names = FALSE)
  lg(sprintf("SUMMARY: %d samples | before %d -> retained %d (%.2f%%)",
             nrow(s), sum(s$n_before), sum(s$n_retained, na.rm = TRUE),
             100 * sum(s$n_retained, na.rm = TRUE) / sum(s$n_before)))
}
lg("END")
close(lgf)
