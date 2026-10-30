###############################################################################
## 93_composition_sc_crosscheck.R
## Independent cross-check of the bulk composition finding in the single-cell
## cohort (GSE292067), using nucleus counts per donor x cell type that were
## already exported by 65_sc_fourgene_plots_donor.R.
##
## What is checked
##   Bulk deconvolution says failing myocardium has FEWER macrophages and MORE
##   fibroblasts than non-failing myocardium. If that is a real tissue difference
##   it should be visible as a per-donor proportion shift in an independent
##   single-nucleus cohort too. The two cohorts differ in aetiology (bulk:
##   DCM/ICM; single-cell: NICM/DoxCM), so this is a directional cross-check,
##   not a replication.
##
## Caveat
##   n_nuclei per donor is affected by dissociation and by how many nuclei were
##   captured per sample; proportions are therefore compositional and only
##   interpretable as relative shifts.
##
## Inputs
##   results/Singlecel/fourgene/tables/65_pseudobulk_sample_x_celltype.csv
##   results/Singlecel/fourgene/tables/64_sample_group.csv
##
## Outputs (results/Composition/)
##   01_tables/T15_sc_donor_composition.csv
##   01_tables/T16_sc_group_comparison.csv
##   02_figures/F07_sc_composition_by_group.{pdf,png}
##   03_logs/93_run_log.txt
##
## Run:  Rscript scripts/93_composition_sc_crosscheck.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

options(stringsAsFactors = FALSE)
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

base <- file.path(proj_dir, "results", "Singlecel", "fourgene", "tables")
pb  <- fread(file.path(base, "65_pseudobulk_sample_x_celltype.csv"))
grp <- fread(file.path(base, "64_sample_group.csv"))

## one row per sample x cell type (n_nuclei is repeated for every gene)
ct <- unique(pb[, .(sample, group, celltype, n_nuclei)])
wide <- dcast(ct, sample ~ celltype, value.var = "n_nuclei", fill = 0)
tot <- rowSums(wide[, -1])
prop <- cbind(sample = wide$sample, as.data.table(wide[, -1] / tot))
prop <- merge(grp[, .(sample, group, condition)], prop, by = "sample")
log_msg("cell types: ", paste(setdiff(names(wide), "sample"), collapse = ", "))
log_msg("donors: ", nrow(prop), " (", paste(names(table(prop$group)), table(prop$group),
                                          sep = "=", collapse = ", "), ")")

write_csv(prop, file.path(d_tab, "T15_sc_donor_composition.csv"))

## group comparison at donor level (Wilcoxon; only 25 donors, so report effect size)
cts <- setdiff(names(wide), "sample")
cmp <- rbindlist(lapply(cts, function(cc) {
  x <- prop[group != "NF"][[cc]]; y <- prop[group == "NF"][[cc]]
  u <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  data.table(celltype = cc,
             median_NF_pct = round(100 * median(y), 2),
             median_disease_pct = round(100 * median(x), 2),
             delta = round(2 * as.numeric(u$statistic) / (length(x) * length(y)) - 1, 3),
             p = signif(u$p.value, 3))
}))
cmp[, fdr := signif(p.adjust(p, "BH"), 3)]
setorder(cmp, p)
write_csv(cmp, file.path(d_tab, "T16_sc_group_comparison.csv"))
for (i in seq_len(nrow(cmp))) {
  log_msg(sprintf("  %-22s NF=%5.2f%%  disease=%5.2f%%  delta=%+.2f  p=%.3f  FDR=%.3f",
                  cmp$celltype[i], cmp$median_NF_pct[i], cmp$median_disease_pct[i],
                  cmp$delta[i], cmp$p[i], cmp$fdr[i]))
}

## figure: proportions of the four compartments of interest
keep <- intersect(c("Cardiomyocyte", "Fibroblast", "Endothelial cell", "Macrophage",
                    "Pericyte", "Smooth muscle cell", "T cell", "B cell"),
                  names(prop))
pl <- lapply(keep, function(cc) {
  dt <- prop[, .(group, value = 100 * get(cc))]
  ggplot(dt, aes(group, value, fill = group)) +
    geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.6) +
    scale_fill_manual(values = c(NF = "#4C78A8", NICM = "#E45756", DoxCM = "#F2A93B"),
                      guide = "none") +
    labs(x = NULL, y = "% of nuclei", title = cc,
         subtitle = sprintf("p = %.3f", cmp[celltype == cc]$p[1])) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8, colour = "grey30"))
})
if (length(pl) > 0) {
  p <- patchwork::wrap_plots(pl, ncol = 4)
  save_fig(p, "F07_sc_composition_by_group", 10, 5)
}

writeLines(run_log, file.path(d_log, "93_run_log.txt"), useBytes = TRUE)
log_msg("Done.")
