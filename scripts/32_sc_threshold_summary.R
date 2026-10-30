# =============================================================================
# 32_sc_threshold_summary.R  (v2)
# Final QC threshold recommendation, per sample, with a single-sample sanity
# comparison against uniform rules. Base R only; nothing is filtered on disk.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

proj_dir <- "E:/another data file/GEObulk analysis project"
out_dir <- file.path(proj_dir, "results", "Singlecel")
cell_qc_dir <- file.path(out_dir, "tables", "cell_qc")
table_dir <- file.path(out_dir, "tables")
dbl_dir <- file.path(table_dir, "doublet_validation")

manifest <- read.csv(file.path(proj_dir, "clean_data", "GSE292067_sample_manifest.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
thr <- read.csv(file.path(table_dir, "qc_threshold_recommendations.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
thr <- thr[match(manifest$sample, thr$sample), ]

# biological floor / ceiling used to keep per-sample MAD estimates sane
FLOOR_NF <- 200; CEIL_NF <- 10000
FLOOR_NC <- 300; CEIL_NC <- 50000
CEIL_MT  <- 3

acc <- list()
for (i in seq_len(nrow(manifest))) {
  s <- manifest$sample[i]
  d <- read.csv(gzfile(file.path(cell_qc_dir, paste0(s, "_cell_qc.csv.gz"))),
                check.names = FALSE, stringsAsFactors = FALSE)
  t <- thr[i, ]
  n <- nrow(d)

  nf_lo <- max(t$nFeature_lower, FLOOR_NF)
  nf_hi <- min(t$nFeature_upper, CEIL_NF)
  nc_lo <- max(t$nCount_lower, FLOOR_NC)
  nc_hi <- min(t$nCount_upper, CEIL_NC)
  mt_hi <- min(t$percent_mt_upper, CEIL_MT)

  qc_H <- d$nFeature_RNA < nf_lo | d$nFeature_RNA > nf_hi |
          d$nCount_RNA < nc_lo | d$nCount_RNA > nc_hi | d$percent.mt > mt_hi

  # uniform alternatives
  qc_U1 <- d$nFeature_RNA < 200 | d$nFeature_RNA > 8000 |
           d$nCount_RNA < 300 | d$nCount_RNA > 40000 | d$percent.mt >= 5
  qc_U2 <- d$nFeature_RNA < 200 | d$nFeature_RNA > 6000 |
           d$nCount_RNA < 300 | d$nCount_RNA > 30000 | d$percent.mt >= 3

  src <- d$doublet_prediction_source == "doublet"

  acc[[i]] <- data.frame(
    sample = s,
    condition = manifest$condition[i],
    n_cells = n,
    nFeature_min = round(nf_lo), nFeature_max = round(nf_hi),
    nCount_min = round(nc_lo), nCount_max = round(nc_hi),
    percent_mt_max = round(mt_hi, 2),
    qc_flagged = sum(qc_H),
    qc_flagged_pct = round(100 * mean(qc_H), 2),
    source_doublets = sum(src),
    total_removed = sum(qc_H | src),
    total_removed_pct = round(100 * mean(qc_H | src), 2),
    retain = sum(!(qc_H | src)),
    retain_pct = round(100 * mean(!(qc_H | src)), 2),
    U1_qc = sum(qc_U1), U2_qc = sum(qc_U2),
    stringsAsFactors = FALSE
  )
  rm(d); invisible(gc())
}
summ <- do.call(rbind, acc)
write.csv(summ, file.path(table_dir, "qc_threshold_final_per_sample.csv"), row.names = FALSE)

tot <- sum(summ$n_cells)
cat("total cells:", tot, "\n\n")
cat(sprintf("H  per-sample (floored/ceiled): QC %d (%.2f%%) | with source doublets %d (%.2f%%) | retain %d (%.2f%%)\n",
            sum(summ$qc_flagged), 100 * sum(summ$qc_flagged) / tot,
            sum(summ$total_removed), 100 * sum(summ$total_removed) / tot,
            sum(summ$retain), 100 * sum(summ$retain) / tot))
cat(sprintf("U1 uniform 200-8000 / 300-40000 / mt<5 : QC %d (%.2f%%)\n",
            sum(summ$U1_qc), 100 * sum(summ$U1_qc) / tot))
cat(sprintf("U2 uniform 200-6000 / 300-30000 / mt<3 : QC %d (%.2f%%)\n",
            sum(summ$U2_qc), 100 * sum(summ$U2_qc) / tot))
cat("\n")
print(summ[, c("sample", "nFeature_min", "nFeature_max", "nCount_min", "nCount_max",
               "percent_mt_max", "qc_flagged", "source_doublets", "total_removed",
               "retain", "retain_pct")], row.names = FALSE)
