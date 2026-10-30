###############################################################################
## 72_global_shift_check.R
## GSE57338 -- is the "everything down-regulated" pattern a normalisation /
##             RNA-quality artefact, or a real group-level shift?
##
## Checks
##   1. per-sample global summaries (median, mean, IQR, 90th percentile)
##   2. housekeeping genes (GAPDH, ACTB, B2M, PPIA, TBP, RPLP0, HPRT1, TUBB)
##   3. correlation of the per-sample median with PC1 / PC2
##   4. mean-vs-log2FC relationship for all genes (MA-style check):
##      if low-expressed genes are systematically down, the shift is global
##   5. the 70 WGCNA-DEG candidate genes: mean expression in cases vs controls
##
## Outputs (results/Batch/)
##   01_tables/T15_global_intensity_by_group.csv
##   01_tables/T16_housekeeping_by_group.csv
##   01_tables/T17_ma_shift_by_expression_bin.csv
##   02_figures/F09_global_intensity_and_housekeeping.{pdf,png}
##
## Run:  Rscript scripts/72_global_shift_check.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma) })

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
clean_dir <- file.path(proj_dir, "clean_data")
res_dir   <- file.path(proj_dir, "results")
d_tab <- file.path(res_dir, "Batch", "01_tables")
d_fig <- file.path(res_dir, "Batch", "02_figures")
d_log <- file.path(res_dir, "Batch", "03_logs")

run_log <- character(0)
log_msg <- function(...) { t <- paste0(...); run_log <<- c(run_log, t); cat(t, "\n", sep = "") }
write_csv <- function(df, p) { utils::write.csv(df, p, row.names = FALSE, fileEncoding = "UTF-8", na = ""); invisible(p) }

## ---------------------------------------------------------------------------
log_msg("[1] loading ...")
expr_dt <- fread(file.path(clean_dir, "GSE57338_clean_expression_matrix.csv"),
                 data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())
bmeta <- fread(file.path(d_tab, "T01_sample_batch_metadata.csv"))
bmeta <- bmeta[match(colnames(expr), sample)]
bmeta[, group2 := factor(group2, levels = c("Control", "Disease"))]
bmeta[, HF := as.integer(group2 == "Disease")]
log_msg("  ", nrow(expr), " genes x ", ncol(expr), " samples")

## ---------------------------------------------------------------------------
log_msg("")
log_msg("[2] per-sample global summaries ...")

gsum <- data.table(
  sample = colnames(expr),
  median_expr = apply(expr, 2, median),
  mean_expr   = apply(expr, 2, mean),
  iqr_expr    = apply(expr, 2, IQR),
  p90_expr    = apply(expr, 2, quantile, 0.9),
  n_detected  = apply(expr, 2, function(x) sum(x > 6))
)
gsum <- merge(gsum, bmeta[, .(sample, group, group2, batch_order)], by = "sample")

summ_rows <- rbindlist(lapply(c("median_expr", "mean_expr", "iqr_expr", "p90_expr", "n_detected"),
  function(v) {
    ct <- gsum[group2 == "Control"][[v]]; ds <- gsum[group2 == "Disease"][[v]]
    w  <- wilcox.test(ds, ct)
    data.table(
      metric = v,
      control_mean = round(mean(ct), 3), control_sd = round(sd(ct), 3),
      disease_mean = round(mean(ds), 3), disease_sd = round(sd(ds), 3),
      difference = round(mean(ds) - mean(ct), 3),
      pct_difference = round(100 * (mean(ds) - mean(ct)) / mean(ct), 2),
      p_wilcox = signif(w$p.value, 4)
    )
  }))
write_csv(summ_rows, file.path(d_tab, "T15_global_intensity_by_group.csv"))
log_msg("  global summaries (disease - control):")
for (i in seq_len(nrow(summ_rows))) {
  log_msg(sprintf("    %-12s  %.3f vs %.3f  (diff %+.3f, %+.1f%%, p=%.3g)",
                  summ_rows$metric[i], summ_rows$disease_mean[i], summ_rows$control_mean[i],
                  summ_rows$difference[i], summ_rows$pct_difference[i], summ_rows$p_wilcox[i]))
}

## ---------------------------------------------------------------------------
log_msg("")
log_msg("[3] housekeeping genes ...")

hk <- c("GAPDH", "ACTB", "B2M", "PPIA", "TBP", "RPLP0", "HPRT1", "TUBB", "EEF1A1", "RPL13A")
hk <- hk[hk %in% rownames(expr)]
log_msg("  found: ", paste(hk, collapse = ", "))

hk_rows <- rbindlist(lapply(hk, function(g) {
  ct <- expr[g, bmeta$group2 == "Control"]; ds <- expr[g, bmeta$group2 == "Disease"]
  data.table(gene = g,
             control_mean = round(mean(ct), 3), control_sd = round(sd(ct), 3),
             disease_mean = round(mean(ds), 3), disease_sd = round(sd(ds), 3),
             diff = round(mean(ds) - mean(ct), 3),
             p_wilcox = signif(wilcox.test(ds, ct)$p.value, 4))
}))
hk_rows[, direction := ifelse(diff < 0, "lower in disease", "higher in disease")]
write_csv(hk_rows, file.path(d_tab, "T16_housekeeping_by_group.csv"))
log_msg("  housekeeping: ", paste(sprintf("%s %+.2f", hk_rows$gene, hk_rows$diff), collapse = ", "))

## ---------------------------------------------------------------------------
log_msg("")
log_msg("[4] PC1 vs global intensity, and MA-style shift by expression bin ...")

gv <- apply(expr, 1, var)
top <- names(sort(gv, decreasing = TRUE))[1:2000]
pca <- prcomp(t(expr[top, , drop = FALSE]), center = TRUE, scale. = FALSE)

cor_tab <- data.table(
  metric = c("median_expr", "mean_expr", "iqr_expr"),
  spearman_vs_PC1 = sapply(c("median_expr", "mean_expr", "iqr_expr"),
                           function(v) round(cor(gsum[[v]], pca$x[, 1], method = "spearman"), 3)),
  spearman_vs_PC2 = sapply(c("median_expr", "mean_expr", "iqr_expr"),
                           function(v) round(cor(gsum[[v]], pca$x[, 2], method = "spearman"), 3)),
  eta2_group = sapply(c("median_expr", "mean_expr", "iqr_expr"), function(v) {
    ss <- summary(aov(gsum[[v]] ~ gsum$group2))[[1]][["Sum Sq"]]; round(ss[1] / sum(ss), 4)
  })
)
write_csv(cor_tab, file.path(d_tab, "T15b_global_intensity_vs_pc.csv"))
log_msg("  PC1 vs median intensity: rho = ", cor_tab$spearman_vs_PC1[1],
        "; eta2 group of median intensity = ", cor_tab$eta2_group[1])

## log2FC vs mean expression, in expression deciles
lf <- fread(file.path(res_dir, "DEG", "GSE57338_DEGall.csv"))
lf <- lf[match(rownames(expr), GeneSymbol)]
lf <- lf[!is.na(GeneSymbol)]
lf[, bin := cut(AveExpr, breaks = quantile(AveExpr, seq(0, 1, 0.1)), include.lowest = TRUE)]
ma <- lf[, .(n_genes = .N,
             mean_AveExpr = round(mean(AveExpr), 2),
             median_log2FC = round(median(log2FoldChange), 3),
             mean_log2FC = round(mean(log2FoldChange), 3),
             pct_up = round(100 * mean(log2FoldChange > 0), 1),
             pct_sig_down = round(100 * mean(padj < 0.05 & log2FoldChange < -0.585), 1)),
         by = bin][order(mean_AveExpr)]
write_csv(ma, file.path(d_tab, "T17_ma_shift_by_expression_bin.csv"))
log_msg("  log2FC by expression decile (median): ",
        paste(sprintf("%.1f:%+.2f", ma$mean_AveExpr, ma$median_log2FC), collapse = "  "))

## ---------------------------------------------------------------------------
log_msg("")
log_msg("[5] 70 candidate genes: expression level in cases vs controls ...")

cand <- fread(file.path(res_dir, "HubGene", "GSE57338_candidate_genes.csv"))
cand <- cand[GeneSymbol %in% rownames(expr)]
cand[, mean_control := sapply(GeneSymbol, function(g) mean(expr[g, bmeta$group2 == "Control"]))]
cand[, mean_disease := sapply(GeneSymbol, function(g) mean(expr[g, bmeta$group2 == "Disease"]))]
cand[, detected_pct := sapply(GeneSymbol, function(g) round(100 * mean(expr[g, ] > 6), 1))]
write_csv(cand[, .(GeneSymbol, log2FoldChange, AveExpr, padj, mean_control, mean_disease,
                   detected_pct, Regulation)],
          file.path(d_tab, "T18_candidate70_expression_level.csv"))
log_msg("  candidate genes: mean AveExpr = ", round(mean(cand$AveExpr), 2),
        ", median detected% = ", median(cand$detected_pct))

## ---------------------------------------------------------------------------
log_msg("")
log_msg("[6] figure ...")

save_both <- function(stem, e, width = 10, height = 4.6) {
  ee <- substitute(e); env <- parent.frame()
  grDevices::pdf(file.path(d_fig, paste0(stem, ".pdf")), width = width, height = height)
  eval(ee, env); grDevices::dev.off()
  grDevices::png(file.path(d_fig, paste0(stem, ".png")), width = width, height = height,
                 units = "in", res = 300)
  eval(ee, env); grDevices::dev.off()
}
col2 <- c("#9ECAE1", "#FC9272")

save_both("F09_global_intensity_and_housekeeping", {
  op <- par(mfrow = c(1, 3), mar = c(3.2, 4, 3, 1), cex = 0.8)
  boxplot(median_expr ~ group2, data = gsum, col = col2, outline = FALSE,
          main = "Per-sample median expression", ylab = "median log2 intensity", xlab = "")
  boxplot(iqr_expr ~ group2, data = gsum, col = col2, outline = FALSE,
          main = "Per-sample IQR", ylab = "IQR of log2 intensity", xlab = "")
  matplot(hk_rows[, .(control_mean, disease_mean)], type = "b", pch = c(19, 17), lty = 1,
          col = c("#2C7FB8", "#D7301F"), xaxt = "n", ylab = "mean log2 expression",
          xlab = "", main = "Housekeeping genes")
  axis(1, at = seq_len(nrow(hk_rows)), labels = hk_rows$gene, las = 2, cex.axis = 0.7)
  legend("bottomleft", pch = c(19, 17), lty = 1, col = c("#2C7FB8", "#D7301F"),
         legend = c("control", "disease"), bty = "n", cex = 0.75)
  par(op)
}, width = 11, height = 4.8)

writeLines(run_log, file.path(d_log, "72_run_log.txt"), useBytes = TRUE)
con <- file(file.path(d_log, "sessionInfo_72.txt"), open = "w"); writeLines(capture.output(sessionInfo()), con); close(con)
log_msg("")
log_msg("Done.")
