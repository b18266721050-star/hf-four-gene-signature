###############################################################################
# 05_immune_report.R
# Session info, parameter record, consolidated result tables
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table) })

ROOT   <- "E:/another data file/GEObulk analysis project"
OUT    <- file.path(ROOT, "results", "Immune")
DIR_LOG<- file.path(OUT, "00_logs")
DIR_T  <- file.path(OUT, "04_stats")
DIR_S  <- file.path(OUT, "02_ssgsea")

pkgs <- c("GSVA", "GSEABase", "Biobase", "data.table", "ggplot2",
          "ComplexHeatmap", "circlize", "RColorBrewer", "car", "matrixStats",
          "ggpubr", "e1071", "msigdbr")
ver <- data.frame(
  package = pkgs,
  installed = sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE)),
  version = sapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)),
                                              error = function(e) "not installed")),
  stringsAsFactors = FALSE
)
write.csv(ver, file.path(DIR_LOG, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo.txt"), useBytes = TRUE)

st  <- readRDS(file.path(DIR_T, "immune_stats.rds"))
kw  <- st$kw; pw <- st$pw; ct <- st$cor_tab
obj <- readRDS(file.path(OUT, "01_qc", "qc_input.rds"))
si  <- obj$si

## consolidated significant table
sig_cells <- kw[FDR < 0.05][order(p_value)]
tab <- rbindlist(lapply(sig_cells$cell_type, function(cl) {
  k  <- kw[cell_type == cl]
  d1 <- pw[cell_type == cl & contrast == "NF_vs_DCM"]
  d2 <- pw[cell_type == cl & contrast == "NF_vs_ICM"]
  d3 <- pw[cell_type == cl & contrast == "DCM_vs_ICM"]
  data.table(
    cell_type = cl,
    KW_chi2 = round(k$kw_chi2, 3),
    KW_P = signif(k$p_value, 4),
    KW_FDR = signif(k$FDR, 4),
    epsilon2 = round(k$epsilon2, 3),
    median_NF = round(d1$median_1, 4),
    median_DCM = round(d1$median_2, 4),
    median_ICM = round(d2$median_2, 4),
    diff_NF_DCM = round(d1$median_diff, 4),
    FDR_NF_DCM = signif(d1$FDR, 4),
    cliffs_delta_NF_DCM = round(d1$cliffs_delta, 3),
    direction_NF_DCM = d1$direction,
    diff_NF_ICM = round(d2$median_diff, 4),
    FDR_NF_ICM = signif(d2$FDR, 4),
    cliffs_delta_NF_ICM = round(d2$cliffs_delta, 3),
    direction_NF_ICM = d2$direction,
    FDR_DCM_ICM = signif(d3$FDR, 4),
    direction_DCM_ICM = d3$direction
  )
}))
write.csv(tab, file.path(DIR_T, "immune_significant_cells_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

sig_cor <- ct[FDR < 0.05][order(gene, -abs(rho))]
write.csv(sig_cor, file.path(DIR_T, "immune_keygene_significant_correlations.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

params <- c(
  "GSE57338 immune infiltration - analysis parameters",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "CST"),
  "",
  "[input]",
  "  expression matrix : clean_data/GSE57338_clean_expression_matrix.csv",
  "  sample info       : clean_data/GSE57338_sample_info.csv (group column)",
  "  key genes         : AQP3, GGT5, ISLR, LAD1",
  "  genes x samples   : 18726 x 313",
  "  scale             : log2 (RMA), continuous, no missing values",
  "  groups            : NF=136, DCM=82, ICM=95",
  "",
  "[method]",
  "  algorithm         : ssGSEA via GSVA::gsva(gsvaParam(...))",
  "  GSVA version      : 2.6.6",
  "  kcdf              : Gaussian",
  "  tau               : 1",
  "  maxDiff           : TRUE",
  "  absRanking        : FALSE",
  "  minSize / maxSize : 1 / Inf",
  "  sparse            : FALSE",
  "  checkNA           : yes",
  "  filterRows        : TRUE",
  "  random seed       : 20261027",
  "",
  "[reference gene set]",
  "  source            : Charoentong P, et al. Cell Reports 2017;18:248-262, Table S6",
  "  content           : 28 human immune cell types, 782 marker genes",
  "  matched in matrix : 705 / 782 (90.15%)",
  "  gmt               : 05_reference/Charoentong2017_28immune_celltypes_matched.gmt",
  "",
  "[statistics]",
  "  normality         : Shapiro-Wilk per cell type x group",
  "  variance          : Levene test (median-centred)",
  "  global test       : Kruskal-Wallis rank-sum, effect size epsilon^2",
  "  post-hoc          : pairwise Wilcoxon rank-sum (normal approx., continuity correction)",
  "  effect size       : Cliff's delta",
  "  multiplicity      : Benjamini-Hochberg FDR (within contrast, 28 tests)",
  "  correlation       : Spearman, BH FDR within gene",
  "  alpha             : 0.05",
  "",
  "[software]",
  paste0("  ", ver$package, " ", ver$version)
)
writeLines(params, file.path(DIR_LOG, "analysis_parameters.txt"), useBytes = TRUE)

cat("done\n")
