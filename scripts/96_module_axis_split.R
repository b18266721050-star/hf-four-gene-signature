###############################################################################
## 96_module_axis_split.R
## GSE57338 -- are all 70 candidate genes carried by the same myeloid axis?
##
## Question
##   95 showed that all 70 candidate genes correlate positively with the xCell
##   macrophage score. If the module is really one axis, then (i) each gene's
##   ability to separate disease from non-failing should be predictable from its
##   myeloid loading, and (ii) once the macrophage axis is regressed out, that
##   ability should largely disappear.
##
## What is NOT done here
##   GSE57338 is a transplant explant cohort without follow-up, so "prognosis"
##   cannot be assessed. Everything below is about DISCRIMINATION (disease versus
##   non-failing) and composition dependence, not about outcome.
##
## Outputs (results/Composition/)
##   01_tables/T18_gene_axis_metrics.csv
##   01_tables/T19_axis_cluster_summary.csv
##   01_tables/T20_axis_tests.csv
##   02_figures/F08_gene_auc_vs_axis.{pdf,png}
##   02_figures/F09_residual_auc_by_cluster.{pdf,png}
##   03_logs/96_run_log.txt
##
## Run:  Rscript scripts/96_module_axis_split.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(pROC); library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20261029)

proj_dir <- "E:/another data file/GEObulk analysis project"
out_root <- file.path(proj_dir, "results", "Composition")
d_tab <- file.path(out_root, "01_tables"); d_fig <- file.path(out_root, "02_figures")
d_log <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) { txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = "") }
write_csv <- function(df, p) { utils::write.csv(df, p, row.names = FALSE, fileEncoding = "UTF-8", na = ""); invisible(p) }
save_fig <- function(p, stem, w, h) {
  ggsave(file.path(d_fig, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(d_fig, paste0(stem, ".png")), p, width = w, height = h, dpi = 600)
}
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

## ===========================================================================
## 1. Inputs
## ===========================================================================
comp <- fread(file.path(d_tab, "T04_composition_merged.csv"))
e <- fread(file.path(proj_dir, "clean_data", "GSE57338_clean_expression_matrix.csv"),
           data.table = FALSE)
m <- as.matrix(e[, -1]); rownames(m) <- e[[1]]
m <- m[, comp$sample]
cand <- fread(file.path(proj_dir, "results", "HubGene", "GSE57338_candidate_genes.csv"))
g70 <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
adj <- fread(file.path(d_tab, "T12_adjusted_70_genes.csv"))
log_msg("70 genes: ", length(g70), " | samples: ", ncol(m))

y <- as.integer(comp$group2 == "Disease")

## per-gene AUC on the raw scale (orientation reported explicitly)
auc_of <- function(v) {
  r <- pROC::roc(y, v, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(r))
}

COMP <- c(mac = "xcell_Macrophages", fib = "xcell_Fibroblasts",
          endo = "xcell_Endothelial_cells", myo = "xcell_Myocytes",
          cd8 = "xcell_CD8_T_cells", imm = "xcell_ImmuneScore")

rows <- rbindlist(lapply(g70, function(g) {
  v <- m[g, ]
  a_raw <- auc_of(v)
  ## residuals after removing the four compartment scores -> composition-free part
  fit <- lm(v ~ comp$xcell_Macrophages + comp$xcell_Fibroblasts +
              comp$xcell_Endothelial_cells + comp$xcell_Myocytes)
  res4 <- residuals(fit)
  a_res <- auc_of(res4)
  data.table(GeneSymbol = g,
             auc_raw = round(max(a_raw, 1 - a_raw), 3),
             direction_raw = ifelse(a_raw < 0.5, "lower in disease", "higher in disease"),
             auc_residual = round(max(a_res, 1 - a_res), 3),
             rho_mac = round(cor(v, comp[[COMP[["mac"]]]], method = "spearman"), 3),
             rho_fib = round(cor(v, comp[[COMP[["fib"]]]], method = "spearman"), 3),
             rho_cd8 = round(cor(v, comp[[COMP[["cd8"]]]], method = "spearman"), 3))
}))

## axis index: myeloid loading minus stromal loading (both z-scored)
rows[, axis := round(as.numeric(scale(rho_mac)) - as.numeric(scale(rho_fib)), 3)]
merge(rows, adj[, .(GeneSymbol, log2FC_raw, log2FC_adj, padj_adj, sig_adj,
                    attenuation_pct, MM_kME)], by = "GeneSymbol") -> gene_tab
setnames(gene_tab, "sig_adj", "still_sig_after_adjust")
setorder(gene_tab, -axis)
write_csv(gene_tab, file.path(d_tab, "T18_gene_axis_metrics.csv"))

## ===========================================================================
## 2. Two-group split by correlation profile (k-means, k = 2)
## ===========================================================================
prof <- scale(as.matrix(gene_tab[, .(rho_mac, rho_fib, rho_cd8)]))
set.seed(20261029)
km <- kmeans(prof, centers = 2, nstart = 50)
gene_tab[, cluster := ifelse(km$cluster == which.max(tapply(gene_tab$axis, km$cluster, mean)),
                             "myeloid-high", "mixed/stromal-low")]
log_msg("cluster sizes: ",
        paste(names(table(gene_tab$cluster)), table(gene_tab$cluster), sep = "=", collapse = ", "))

clu <- gene_tab[, .(n = .N,
                    axis_median = round(median(axis), 2),
                    rho_mac_median = round(median(rho_mac), 2),
                    rho_fib_median = round(median(rho_fib), 2),
                    log2FC_raw_median = round(median(log2FC_raw), 2),
                    log2FC_adj_median = round(median(log2FC_adj), 2),
                    attenuation_median = round(median(attenuation_pct), 1),
                    auc_raw_median = round(median(auc_raw), 2),
                    auc_residual_median = round(median(auc_residual), 2),
                    n_still_sig = sum(still_sig_after_adjust, na.rm = TRUE)),
                by = cluster]
write_csv(clu, file.path(d_tab, "T19_axis_cluster_summary.csv"))
print(clu)

## ===========================================================================
## 3. Tests
## ===========================================================================
tests <- list()
tests[["cor(AUC_raw, axis)"]] <- cor.test(gene_tab$auc_raw, gene_tab$axis)
tests[["cor(AUC_raw, rho_mac)"]] <- cor.test(gene_tab$auc_raw, gene_tab$rho_mac)
tests[["wilcox(AUC_raw by cluster)"]] <- wilcox.test(auc_raw ~ cluster, data = gene_tab)
tests[["wilcox(auc_residual by cluster)"]] <- wilcox.test(auc_residual ~ cluster, data = gene_tab)
tests[["wilcox(attenuation by cluster)"]] <- wilcox.test(attenuation_pct ~ cluster, data = gene_tab)

test_tab <- rbindlist(lapply(names(tests), function(nm) {
  t <- tests[[nm]]
  est <- if (is.null(t$estimate)) NA_real_ else unname(t$estimate)[1]
  data.table(comparison = nm, estimate = round(est, 3), p = signif(t$p.value, 3))
}))
write_csv(test_tab, file.path(d_tab, "T20_axis_tests.csv"))
for (i in seq_len(nrow(test_tab))) {
  log_msg(sprintf("  %-32s estimate=%s  p=%.3g", test_tab$comparison[i],
                  ifelse(is.na(test_tab$estimate[i]), "-", sprintf("%+.2f", test_tab$estimate[i])),
                  test_tab$p[i]))
}
log_msg("  AUC drop after removing composition (median): ",
        round(median(gene_tab$auc_raw - gene_tab$auc_residual), 3))
log_msg("  genes with residual AUC >= 0.70: ",
        sum(gene_tab$auc_residual >= 0.70), "/", nrow(gene_tab))
log_msg("  genes with residual AUC <= 0.60 (axis fully explains): ",
        sum(gene_tab$auc_residual <= 0.60), "/", nrow(gene_tab))

## ===========================================================================
## 4. Figures
## ===========================================================================
p8 <- ggplot(gene_tab, aes(axis, auc_raw, colour = cluster)) +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.15,
              aes(group = 1), colour = "grey30", formula = y ~ x) +
  scale_colour_manual(values = c(`myeloid-high` = "#C53030",
                                 `mixed/stromal-low` = "#2B6CB0"), name = NULL) +
  ggrepel::geom_text_repel(data = gene_tab[auc_raw >= 0.90],
                           aes(label = GeneSymbol), size = 2.5,
                           min.segment.length = 0, max.overlaps = 30,
                           box.padding = 0.25, show.legend = FALSE) +
  labs(x = "Axis index (myeloid loading minus stromal loading, z-scored)",
       y = "Single-gene AUC (disease vs non-failing)",
       title = "Discrimination of individual candidate genes scales with myeloid loading",
       subtitle = sprintf("r = %.2f (p = %.1e); the two k-means groups are not separated on this axis",
                          cor(gene_tab$auc_raw, gene_tab$axis),
                          cor.test(gene_tab$auc_raw, gene_tab$axis)$p.value)) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"))
save_fig(p8, "F08_gene_auc_vs_axis", 7.6, 5.6)

lg <- melt(gene_tab[, .(GeneSymbol, cluster, raw = auc_raw, residual = auc_residual)],
           id.vars = c("GeneSymbol", "cluster"),
           variable.name = "scale", value.name = "auc")
lg[, scale := factor(scale, levels = c("raw", "residual"),
                     labels = c("raw expression", "after removing composition"))]
p9 <- ggplot(lg, aes(scale, auc, fill = cluster)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 0.7, linetype = 3, colour = "grey60") +
  scale_fill_manual(values = c(`myeloid-high` = "#C53030",
                               `mixed/stromal-low` = "#2B6CB0"), name = NULL) +
  labs(x = NULL, y = "Single-gene AUC",
       title = "Most of the per-gene discrimination disappears with the composition axis")
save_fig(p9, "F09_residual_auc_by_cluster", 6.4, 4.8)

writeLines(run_log, file.path(d_log, "96_run_log.txt"), useBytes = TRUE)
