###############################################################################
## 139_power_and_compartments.R
##
## Two audit items in one pass.
##
## R6  What difference could the external cohorts actually detect?
##     Minimum detectable AUC difference at 80% power (two-sided 5%) for each
##     external cohort against the internal split, using the Hanley-McNeil
##     variance of the AUC.
##
## R11 Which deconvolution compartments may be used as evidence?
##     A curated table that keeps only compartments whose marker genes are
##     plausibly specific to the tissue (macrophages, fibroblasts/CAFs,
##     endothelium, myocytes) and flags the rest as signature cross-reactivity
##     (neurons, chondrocytes, melanocytes, keratinocytes, skeletal muscle).
##
## Outputs
##   results/ExternalValidation/01_tables/T17_external_power.csv
##   results/Composition/01_tables/T25_compartments_curated.csv
##   results/Composition/02_figures/F11_curated_compartments.{pdf,png}
##   logs
##
## Run: Rscript scripts/139_power_and_compartments.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

proj <- "E:/another data file/GEObulk analysis project"
run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = ""); flush.console()
}
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

## ===========================================================================
## R6  minimum detectable difference
## ===========================================================================
se_auc <- function(a, n1, n0) {
  q1 <- a / (2 - a)
  q2 <- 2 * a^2 / (1 + a)
  sqrt((a * (1 - a) + (n1 - 1) * (q1 - a^2) + (n0 - 1) * (q2 - a^2)) / (n1 * n0))
}

cohorts <- data.table(
  cohort = c("Internal split (reference)", "GSE5406", "GSE116250", "GSE55296"),
  n_case = c(53, 194, 50, 23),
  n_ctrl = c(41, 16, 12, 10),
  auc    = c(0.9926, 0.9191, 0.9267, 0.8435)
)
cohorts[, se := se_auc(auc, n_case, n_ctrl)]
ref <- cohorts[1]
cohorts[, diff_vs_internal := round(ref$auc - auc, 4)]
cohorts[, se_diff := sqrt(se^2 + ref$se^2)]
cohorts[, mdd_80 := round(2.80 * se_diff, 4)]          # 1.96 + 0.84
cohorts[, power_observed := round(pnorm(abs(diff_vs_internal) / se_diff - 1.96), 3)]
cohorts[, note := c("reference",
                    "16 controls: only a large drop is detectable",
                    "12 controls: only a large drop is detectable",
                    "10 controls: only a large drop is detectable")]
write_csv(cohorts[, .(cohort, n_case, n_ctrl, auc, se = round(se, 4),
                      diff_vs_internal, mdd_80, power_observed, note)],
          file.path(proj, "results", "ExternalValidation", "01_tables",
                    "T17_external_power.csv"))
for (i in seq_len(nrow(cohorts))) {
  log_msg(sprintf("  %-28s n=%3d/%2d  AUC=%.3f  observed drop=%.3f  minimum detectable (80%% power)=%.3f",
                  cohorts$cohort[i], cohorts$n_case[i], cohorts$n_ctrl[i], cohorts$auc[i],
                  abs(cohorts$diff_vs_internal[i]), cohorts$mdd_80[i]))
}

## ===========================================================================
## R11  curated compartments
## ===========================================================================
log_msg("")
log_msg("[R11] curated compartment table")
t07 <- data.table::fread(file.path(proj, "results", "Composition", "01_tables",
                                   "T07_single_feature_auc.csv"))
keep_x <- c("Macrophages", "Macrophages_M1", "Macrophages_M2", "Fibroblasts",
            "CAFs", "Endothelial_cells", "mv_Endothelial_cells", "Myocytes")
keep_e <- c("Macrophages", "CAFs", "Endothelial")

t07[, keep := (method == "xCell" & feature %in% keep_x) |
             (method == "EPIC"  & feature %in% keep_e)]
t07[, evidence := ifelse(keep, "marker genes plausibly tissue-specific",
                         "signature cross-reactivity - do not interpret as cell abundance")]
cur <- t07[keep == TRUE][order(method, -auc_dist_from_half)]
write_csv(cur[, .(method, feature, auc, delta, fdr, direction = ifelse(delta > 0,
                                                                       "higher in disease", "lower in disease"))],
          file.path(proj, "results", "Composition", "01_tables", "T25_compartments_curated.csv"))

excluded <- t07[keep == FALSE][order(-auc_dist_from_half)][1:6]
log_msg("  retained compartments: ", paste(cur$feature, collapse = ", "))
log_msg("  excluded examples: ", paste(excluded$feature, collapse = ", "))

## concordance of the two deconvolution methods for macrophages
t04 <- data.table::fread(file.path(proj, "results", "Composition", "01_tables",
                                   "T04_composition_merged.csv"))
rho <- cor(t04$xcell_Macrophages, t04$epic_Macrophages, method = "spearman")
log_msg(sprintf("  xCell vs EPIC macrophage score, Spearman rho = %.3f", rho))
log_msg(sprintf("  AUC: xCell Macrophages %.3f ; EPIC Macrophages %.3f",
                t07[method == "xCell" & feature == "Macrophages"]$auc,
                t07[method == "EPIC" & feature == "Macrophages"]$auc))

## figure: curated compartments only
cur[, label := factor(paste(method, feature), levels = rev(paste(method, feature)))]
p <- ggplot(cur, aes(auc, label, colour = method)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = 0.5, xend = auc, yend = label), colour = "grey85") +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(xCell = "#2B6CB0", EPIC = "#C53030"), name = NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "AUC, disease versus non-failing (values below 0.5 = lower in disease)",
       y = NULL,
       title = "Only compartments with plausibly specific markers",
       subtitle = "Compartments whose signatures cross-react with non-cardiac tissues are not shown") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
d_fig <- file.path(proj, "results", "Composition", "02_figures")
ggsave(file.path(d_fig, "F11_curated_compartments.pdf"), p, width = 6.6, height = 4.2)
ggsave(file.path(d_fig, "F11_curated_compartments.png"), p, width = 6.6, height = 4.2, dpi = 600)

log_msg("Done.")
writeLines(run_log, file.path(proj, "results", "Composition", "03_logs", "139_run_log.txt"),
           useBytes = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(proj, "results", "Composition", "03_logs", "sessionInfo_139.txt"))
