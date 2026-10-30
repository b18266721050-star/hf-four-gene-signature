###############################################################################
## 73_pack_batch_results.R
## Pack every table from results/Batch/01_tables into one Excel workbook,
## with the main before/after comparison as the first sheet.
##
## Output: results/Batch/GSE57338_batch_diagnostics_and_combat.xlsx
## Run:    Rscript scripts/73_pack_batch_results.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
out_root <- file.path(proj_dir, "results", "Batch")
d_tab    <- file.path(out_root, "01_tables")
xlsx_out <- file.path(out_root, "GSE57338_batch_diagnostics_and_combat.xlsx")

## byte-safe UTF-8 reader: file(encoding="UTF-8") fails on this Windows locale
read_utf8 <- function(p) {
  raw <- readBin(p, "raw", file.info(p)$size)
  txt <- iconv(rawToChar(raw), from = "UTF-8", to = "UTF-8", sub = "?")
  strsplit(txt, "
", fixed = TRUE)[[1]]
}

## sheet order matters: the main table first
sheet_map <- data.table(
  file = c("T14_comparison_table_main.csv",
           "T03_batch_association_summary.csv",
           "T04_variance_decomposition.csv",
           "T05_batch_predicts_group_auc.csv",
           "T07_id_axis_vs_order.csv",
           "T08_deg_summary_before_after.csv",
           "T09_key_gene_log2fc_before_after.csv",
           "T10_candidate70_direction_before_after.csv",
           "T11_model_auc_before_after.csv",
           "T12_cv_scheme_comparison.csv",
           "T13_pca_variance_before_after.csv",
           "T15_global_intensity_by_group.csv",
           "T16_housekeeping_by_group.csv",
           "T17_ma_shift_by_expression_bin.csv",
           "T06_surrogate_variables.csv"),
  sheet = c("01_main_comparison",
            "02_batch_association",
            "03_variance_decomposition",
            "04_batch_predicts_group",
            "05_id_axis_vs_order",
            "06_deg_summary",
            "07_key_genes",
            "08_candidate70",
            "09_model_auc",
            "10_cv_schemes",
            "11_pca_before_after",
            "12_global_intensity",
            "13_housekeeping",
            "14_ma_shift",
            "15_surrogate_variables")
)

wb <- createWorkbook()

## --- READ ME sheet (Chinese text read from a UTF-8 file) ---
addWorksheet(wb, "00_READ_ME")
readme <- read_utf8(file.path(out_root, "00_readme_cn.txt"))
writeData(wb, "00_READ_ME", data.frame(text = readme), colNames = FALSE)
setColWidths(wb, "00_READ_ME", cols = 1, widths = 110)
addStyle(wb, "00_READ_ME", createStyle(wrapText = FALSE, valign = "top"),
         rows = 1:length(readme), cols = 1, gridExpand = TRUE)

## --- data sheets ---
for (i in seq_len(nrow(sheet_map))) {
  f <- file.path(d_tab, sheet_map$file[i])
  if (!file.exists(f)) { message("skip (missing): ", sheet_map$file[i]); next }
  d <- fread(f, encoding = "UTF-8")
  addWorksheet(wb, sheet_map$sheet[i])
  writeData(wb, sheet_map$sheet[i], d, keepNA = FALSE)
  freezePane(wb, sheet_map$sheet[i], firstRow = TRUE)
  setColWidths(wb, sheet_map$sheet[i], cols = seq_len(ncol(d)),
               widths = "auto")
  addStyle(wb, sheet_map$sheet[i],
           createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom"),
           rows = 1, cols = seq_len(ncol(d)), gridExpand = TRUE)
}

## --- highlight the key rows of the main table ---
key_metrics <- c("Significant DE genes (n)", "M2 internal-validation AUC",
                 "Validation AUC: batch block alone", "Batch-aware CV AUC")
d_main <- fread(file.path(d_tab, "T14_comparison_table_main.csv"), encoding = "UTF-8")
rows <- which(d_main$metric %in% key_metrics) + 1
if (length(rows) > 0) {
  addStyle(wb, "01_main_comparison",
           createStyle(fgFill = "#FFF2CC", textDecoration = "bold"),
           rows = rows, cols = 1:ncol(d_main), gridExpand = TRUE, stack = TRUE)
}

saveWorkbook(wb, xlsx_out, overwrite = TRUE)
cat("written:", xlsx_out, "\n")
cat("sheets:", paste(names(wb), collapse = " | "), "\n")
