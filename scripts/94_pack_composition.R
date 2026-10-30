###############################################################################
## 94_pack_composition.R
## Pack the cell-composition results into one Excel workbook + figure legends.
## NB: this script is deliberately pure ASCII. R on this machine runs in a
## non-UTF-8 locale and refuses to parse Chinese literals; the Chinese readme
## (00_readme_cn.txt) is therefore written outside R.
##
## Output
##   results/Composition/GSE57338_cell_composition.xlsx
##   results/Composition/02_figures/figure_legends_composition.txt
##
## Run:  Rscript scripts/94_pack_composition.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })

proj_dir <- "E:/another data file/GEObulk analysis project"
out_root <- file.path(proj_dir, "results", "Composition")
d_tab <- file.path(out_root, "01_tables")
d_fig <- file.path(out_root, "02_figures")
xlsx_out <- file.path(out_root, "GSE57338_cell_composition.xlsx")

read_utf8 <- function(p) {
  raw <- readBin(p, "raw", file.info(p)$size)
  txt <- iconv(rawToChar(raw), from = "UTF-8", to = "UTF-8", sub = "?")
  strsplit(txt, "\n", fixed = TRUE)[[1]]
}

read_tab <- function(p) {
  d <- tryCatch(data.table::fread(p, encoding = "UTF-8", showProgress = FALSE),
                error = function(e) NULL)
  if (is.null(d) || ncol(d) == 1) {
    txt <- read_utf8(p)
    d <- data.table::fread(text = paste(txt, collapse = "\n"), showProgress = FALSE)
  }
  as.data.frame(d)
}

sheet_map <- data.table(
  file = c("T05_xcell_group_comparison.csv",
           "T06_epic_group_comparison.csv",
           "T07_single_feature_auc.csv",
           "T08_composition_vs_genes_correlation.csv",
           "T09_composition_classifier_auc.csv",
           "T11_adjusted_deg_summary.csv",
           "T12_adjusted_70_genes.csv",
           "T13_adjusted_key_genes.csv",
           "T14_incremental_value_auc.csv",
           "T15_sc_donor_composition.csv",
           "T16_sc_group_comparison.csv",
           "T10_composition_pca.csv",
           "T04_composition_merged.csv",
           "T02_xcell_scores.csv",
           "T03_epic_fractions.csv",
           "T17_candidate70_vs_composition.csv",
           "T18_gene_axis_metrics.csv",
           "T19_axis_cluster_summary.csv",
           "T20_axis_tests.csv"),
  sheet = c("02_xcell_group_comparison",
            "03_epic_group_comparison",
            "04_single_feature_auc",
            "05_composition_vs_genes",
            "06_composition_classifier",
            "07_adjusted_deg_summary",
            "08_adjusted_70_genes",
            "09_adjusted_key_genes",
            "10_incremental_auc",
            "11_sc_donor_composition",
            "12_sc_group_comparison",
            "13_composition_pca",
            "14_merged_scores",
            "15_xcell_scores_raw",
            "16_epic_fractions_raw",
            "17_candidate70_vs_composition",
            "18_gene_axis_metrics",
            "19_axis_cluster_summary",
            "20_axis_tests")
)

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 10, fontName = "Calibri")

readme <- data.frame(
  item = c("Dataset", "Date", "Deconvolution", "Input", "Samples",
           "Finding 1", "Finding 2", "Finding 3", "Caveat", "Scripts"),
  detail = c("GSE57338 (313 samples: NF 136 / DCM 82 / ICM 95)",
             "2026-10-29",
             "xCell 1.1.0 (ssGSEA, primary) + EPIC (absolute fractions, sensitivity) + single-nucleus donor proportions (cross-check)",
             "clean_data/GSE57338_clean_expression_matrix.csv (RMA log2, 18,726 genes)",
             "bulk 313; single-nucleus 25 donors (NF 8 / NICM 6 / DoxCM 11)",
             "Macrophage and endothelial scores are lower, fibroblast/CAF scores higher, in failing myocardium (42/67 xCell features FDR<0.05)",
             "Composition alone discriminates disease from non-failing at validation AUC 0.944 (6 features) to 0.979 (67 features); the four-gene model gives 0.993",
             "After adjustment for four compartment scores only 22/70 candidate genes remain significant and PPI hubs attenuate by 40-77%, whereas the four model genes attenuate by only 15-33% and stay highly significant",
             "Macrophage nuclei are NOT reduced in the independent single-nucleus cohort (7.7% vs 9.3%, p=0.47), so the lower bulk macrophage score cannot be read as fewer macrophages; see the report",
             "scripts/90_composition_deconv.R -> 91 -> 92 -> 93 -> 94")
)
addWorksheet(wb, "01_readme")
writeData(wb, "01_readme", readme)
setColWidths(wb, "01_readme", cols = 1:2, widths = c(14, 120))
addStyle(wb, "01_readme", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE)
setRowHeights(wb, "01_readme", rows = 2:(nrow(readme) + 1), heights = 30)

for (i in seq_len(nrow(sheet_map))) {
  p <- file.path(d_tab, sheet_map$file[i])
  if (!file.exists(p)) next
  d <- read_tab(p)
  addWorksheet(wb, sheet_map$sheet[i])
  writeData(wb, sheet_map$sheet[i], d)
  freezePane(wb, sheet_map$sheet[i], firstRow = TRUE)
  setColWidths(wb, sheet_map$sheet[i], cols = seq_len(ncol(d)),
               widths = pmin(pmax(nchar(names(d)) + 3, 11), 42))
}
saveWorkbook(wb, xlsx_out, overwrite = TRUE)
cat("xlsx written:", xlsx_out, "\n")

leg <- c(
  "F01_key_compartment_boxplots. xCell compartment scores (spillover-compensated ssGSEA, relative scale) in non-failing (NF), dilated (DCM) and ischaemic (ICM) cardiomyopathy. Boxes show the median and interquartile range; whiskers extend to 1.5x IQR. Two-sided Wilcoxon rank-sum p values (disease versus NF) are given above each panel.",
  "F02_composition_gene_correlation. Spearman correlations between compartment scores (xCell rows, plus two EPIC compartments) and the four panel genes, the seven PPI hub genes and the frozen four-gene linear predictor (lp_M2). Only |rho| >= 0.30 is labelled. Deconvolution scores are estimated from marker transcripts, so part of this correlation is intrinsic.",
  "F03_composition_pca. Principal component analysis of the 67 xCell features across 313 samples, coloured by diagnosis. eta-squared of diagnosis on PC1 is given; PC2 separates disease from non-failing myocardium far better than PC1.",
  "F04_single_feature_auc. Single-feature discrimination (disease versus non-failing) for the strongest xCell and EPIC features. Values are oriented so that larger means higher in disease; the direction is stated in T05/T06.",
  "F05_log2fc_before_after_adjustment. log2 fold change of the 269 differentially expressed genes before and after adding four compartment scores (macrophages, fibroblasts, endothelial cells, myocytes) as covariates. Red points are the 70 WGCNA-intersect candidate genes; the dashed line is identity and the blue line is the fitted slope.",
  "F06_incremental_auc. Validation AUC of a composition-only model, the four-gene model and their combination, all fitted on the training split (n = 219) and evaluated on the same internal validation split (n = 94). Error bars are DeLong 95% confidence intervals.",
  "F07_sc_composition_by_group. Per-donor nucleus proportions in the independent single-nucleus cohort (GSE292067; NF 8, NICM 6, DoxCM 11). Two-sided Wilcoxon p values; macrophages and fibroblasts do not differ, whereas pericytes fall and T/NK cells rise.",
  "F08_gene_auc_vs_axis. For each of the 70 candidate genes, single-gene discrimination of disease versus non-failing is plotted against the axis index (z-scored macrophage loading minus z-scored fibroblast loading). The fitted line has r = 0.75 (p = 5.1e-14). Colours are k-means groups computed from the correlation profiles; the two groups overlap completely on this axis, so the module does not split into a myeloid and a stromal subpopulation.",
  "F09_residual_auc_by_cluster. Single-gene AUC before and after regressing out the four compartment scores. Raw AUC differs between the k-means groups (median 0.89 vs 0.85, p = 0.0014), but after removal of the composition axis both collapse to 0.63 and 0.62 (p = 0.78), i.e. the per-gene discrimination is carried by the same composition axis in both groups."
)
writeLines(leg, file.path(d_fig, "figure_legends_composition.txt"), useBytes = TRUE)
cat("done\n")
