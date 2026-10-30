# =============================================================================
# 82_pack_ext.R —— 外部验证结果打包（xlsx + 主汇总表）
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })

root <- "E:/another data file/GEObulk analysis project"
outd <- file.path(root, "results", "ExternalValidation")
tabd <- file.path(outd, "01_tables"); datd <- file.path(outd, "04_data")

rd <- function(f) fread(file.path(tabd, f), data.table = FALSE)
auc <- rd("T02_external_auc.csv"); dir <- rd("T03_direction_external.csv")
dirsum <- rd("T04_direction_summary.csv"); cal <- rd("T05_calibration_brier.csv")
sc <- rd("T06_scale_diagnostics.csv"); dl <- rd("T07_delong_internal_vs_external.csv")
int <- rd("T01_internal_auc_recomputed.csv")
boot <- rd("T09_bootstrap_auc.csv"); base <- rd("T10_cohort_baseline.csv")
dca <- rbindlist(lapply(c("GSE5406", "GSE116250", "GSE55296"), function(cn)
  fread(file.path(datd, paste0("dca_", cn, ".csv")))))
caldisc <- fread(file.path(datd, "calibration_deciles_GSE57338_internal.csv"))

# ---- 主汇总：每个队列一行
coh <- c("GSE5406", "GSE116250", "GSE55296")
get <- function(cn, subset, ver) {
  v <- auc$auc[auc$cohort == cn & auc$subset == subset & auc$version == ver]
  if (!length(v)) NA_real_ else v[1]
}
master <- data.frame(
  cohort = coh,
  platform = c("Affymetrix HG-U133A (RMA log2)", "Illumina HiSeq 2500 (log2 RPKM)",
               "RNA-seq normalized counts (log2)"),
  n_total = sapply(coh, function(cn) auc$n[auc$cohort == cn & auc$subset == "all disease vs NF" & auc$version == "V2"]),
  n_case  = sapply(coh, function(cn) auc$n_case[auc$cohort == cn & auc$subset == "all disease vs NF" & auc$version == "V2"]),
  n_ctrl  = sapply(coh, function(cn) auc$n_ctrl[auc$cohort == cn & auc$subset == "all disease vs NF" & auc$version == "V2"]),
  auc_V2 = sapply(coh, get, subset = "all disease vs NF", ver = "V2"),
  lo_V2  = sapply(coh, function(cn) auc$low[auc$cohort == cn & auc$subset == "all disease vs NF" & auc$version == "V2"]),
  hi_V2  = sapply(coh, function(cn) auc$high[auc$cohort == cn & auc$subset == "all disease vs NF" & auc$version == "V2"]),
  auc_DCM = sapply(coh, get, subset = "DCM vs NF", ver = "V2"),
  auc_ICM = sapply(coh, get, subset = "ICM vs NF", ver = "V2"),
  auc_V1_raw = sapply(coh, get, subset = "all disease vs NF", ver = "V1"),
  auc_V3_within = sapply(coh, get, subset = "all disease vs NF", ver = "V3"),
  auc_dropLAD1 = sapply(coh, get, subset = "score: 3 genes (drop LAD1)", ver = "V2"),
  auc_dropISLR = sapply(coh, get, subset = "score: 3 genes (drop ISLR)", ver = "V2"),
  auc_ISLR_only = sapply(coh, get, subset = "single gene ISLR", ver = "V2"),
  cal_slope = cal$cal_slope[match(coh, cal$cohort)],
  brier = cal$brier[match(coh, cal$cohort)],
  boot_low = boot$boot_low[match(coh, boot$cohort)],
  boot_high = boot$boot_high[match(coh, boot$cohort)],
  direction_consistent = c("6/6 (all genes, both aetiologies)",
                           "5/6 (GGT5 in ICM reversed, p=0.85)",
                           "5/6 (LAD1 in DCM reversed, p=0.60)"),
  stringsAsFactors = FALSE)
master <- master[, c("cohort","platform","n_total","n_case","n_ctrl","auc_V2","lo_V2","hi_V2",
                     "boot_low","boot_high","auc_DCM","auc_ICM","auc_V1_raw","auc_V3_within",
                     "auc_dropLAD1","auc_dropISLR","auc_ISLR_only","cal_slope","brier",
                     "direction_consistent")]
write.csv(master, file.path(tabd, "T08_master_summary.csv"), row.names = FALSE)

smp <- list(
  GSE5406  = fread(file.path(outd, "01_data", "GSE5406_sample_info.csv"), data.table = FALSE),
  GSE116250 = fread(file.path(outd, "01_data", "GSE116250_sample_info.csv"), data.table = FALSE),
  GSE55296 = fread(file.path(outd, "01_data", "GSE55296_sample_info.csv"), data.table = FALSE))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")
add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("00_SUMMARY", master)
add_sheet("01_AUC_allversions", auc)
add_sheet("02_direction_bycohort", dir)
add_sheet("03_direction_summary", dirsum)
add_sheet("04_calibration", cal)
add_sheet("05_calibration_deciles", rbind(caldisc, rbindlist(lapply(
  c("GSE5406", "GSE116250", "GSE55296"), function(cn)
    fread(file.path(datd, paste0("calibration_deciles_", cn, ".csv")))))))
add_sheet("06_dca", dca)
add_sheet("07_bootstrap_auc", boot)
add_sheet("08_cohort_baseline", base)
add_sheet("09_scale_diagnostics", sc)
add_sheet("10_delong_vs_internal", dl)
add_sheet("11_internal_auc", int)
add_sheet("12_samples_GSE5406", smp$GSE5406)
add_sheet("13_samples_GSE116250", smp$GSE116250)
add_sheet("14_samples_GSE55296", smp$GSE55296)

notes <- data.frame(note = c(
 "External validation of the frozen 4-gene model (M2) of GSE57338.",
 "Model: LP = 0.0185 - 1.4476*z(AQP3) - 1.4385*z(GGT5) + 2.1951*z(ISLR) - 1.1838*z(LAD1).",
 "z-score parameters (mean/SD) are those of the training split (T06_zscore_parameters_train_only.csv); no refitting.",
 "V2 (primary): each gene's external values were quantile-mapped to the training distribution, then z-scored with the frozen parameters.",
 "V1 (diagnostic): frozen parameters applied to raw external values without any scale alignment.",
 "V3 (sensitivity): within-cohort z-standardisation, i.e. the frozen weights applied to cohort-specific SD units.",
 "Direction: limma with group as the only factor (GSE5406 lacks age/sex; GSE116250 has them but was kept unadjusted for comparability).",
 "GSE116250: LAD1 contained the placeholder value 999999999 in two non-failing samples (GSM3219560, GSM3219561); set to NA, so the 4-gene model uses n = 62.",
 "GSE55296: three count columns (G1, G105, G110) are absent from the submitter readme and could not be assigned a disease state; 33 of 36 samples used.",
 "Calibration intercept/slope are estimated by logistic recalibration of the frozen linear predictor; slopes < 1 indicate overfitting of the discovery model.",
 "Bootstrap confidence intervals (sheet 07) are stratified by outcome, 2000 resamples, percentile method; reported alongside DeLong intervals because DeLong is anti-conservative with few controls.",
 "Decision curve analysis (sheet 06) uses the frozen model only; no clinical variables were available in the external cohorts, so the net benefit refers to a gene-only strategy.",
 "Sheet 08 lists baseline characteristics where the GEO record provides them; age and sex are absent for GSE5406 and GSE55296.",
 "Brier scores are inflated by case-enriched designs (92%, 81%, 70% cases) and should not be compared across cohorts without prevalence adjustment.",
 "Sheets 15-18 report the sample-overlap audit between GSE5406 and the discovery cohort: individual-specific residuals were correlated across 11,808 shared genes after removing group means and 10 within-cohort PCs; same-platform self-splits provide the ceiling attainable without any duplication (max r = 0.368), which the tested pair did not reach (max r = 0.283)."),
 stringsAsFactors = FALSE)
add_sheet("99_README", notes)

# 若样本重叠审计（83）已跑过，一并打包
ov_files <- c("T11_overlap_top_pairs.csv", "T12_overlap_summary.csv",
              "T13_overlap_group_concordance.csv", "T14_overlap_mutual_best.csv")
ov_names <- c("15_overlap_top_pairs", "16_overlap_summary",
              "17_overlap_concordance", "18_overlap_mutual_best")
for (i in seq_along(ov_files)) {
  fp <- file.path(tabd, ov_files[i])
  if (file.exists(fp)) add_sheet(ov_names[i], fread(fp, data.table = FALSE))
}

f <- file.path(outd, "GSE57338_external_validation.xlsx")
saveWorkbook(wb, f, overwrite = TRUE)
cat("written:", f, "\n")
print(master[, c("cohort","n_total","n_ctrl","auc_V2","lo_V2","hi_V2","auc_DCM","auc_ICM","cal_slope")])
