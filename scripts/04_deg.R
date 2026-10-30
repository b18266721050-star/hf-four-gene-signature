# =============================================================================
# 04_deg.R  ——  GSE57338 差异表达分析（疾病 vs 对照，limma）
#
# 输入（clean_data/）：
#   GSE57338_clean_expression_matrix.csv   log2 表达矩阵（18,726 基因 × 313 样本）
#   GSE57338_sample_info.csv               样本信息（含 group: NF / DCM / ICM）
#
# 输出（results/DEG/）：
#   GSE57338_DEGall.csv                    全部基因的差异分析结果（未筛选）
#   GSE57338_DEG_significant.csv           显著差异基因（|log2FC|>0.585 且 padj<0.05）
#   GSE57338_DEG_up.csv / _down.csv        显著基因按上调/下调拆分
#   GSE57338_DEG_summary.txt               分析参数与结果汇总（机器可读的 txt）
#   补充分析：
#     GSE57338_DEGall_DCM_vs_Control.csv   DCM 亚型 vs 对照
#     GSE57338_DEGall_ICM_vs_Control.csv   ICM 亚型 vs 对照
#     GSE57338_DEGall_noCovariate.csv      敏感性分析：仅 ~group（不校正性别/年龄）
#     GSE57338_DEGall_classicEB.csv        敏感性分析：经典 eBayes（无趋势/稳健估计）
# 输出（report/）：
#   sessionInfo_04_deg.txt
#
# 用法（在项目根目录下）：
#   Rscript scripts/04_deg.R
#   或在 RStudio 中把工作目录设到项目根后 source("scripts/04_deg.R")
#
# 注意：脚本内的运行信息（console / txt）统一用英文，是为了避免 Windows 中文
#       环境（GBK / C locale）下 source 文件的字符串被误读成乱码；中文说明写在
#       注释里（注释不参与解析，不受编码影响）。
#
# -----------------------------------------------------------------------------
# 方法学说明（为什么用 limma 而不是 DESeq2 / edgeR）
# -----------------------------------------------------------------------------
#   1) 本数据集是 Affymetrix HuGene-1.1-ST 芯片（GPL11532），clean_data 里是
#      RMA 归一化后的 **log2 连续表达值**（range 1.61~13.05），不是整数 counts。
#      DESeq2 / edgeR 的负二项模型要求原始 count，强行套用属于方法误用。
#   2) limma 的线性模型 + 经验贝叶斯方差收缩正是为这类（近似正态的 log 强度）
#      数据设计的，样本量大（313）时稳健性很好。
#   3) 采用 eBayes(trend = TRUE, robust = TRUE)：
#        trend  —— 允许基因方差随平均表达强度变化（芯片数据普遍存在的
#                  强度-方差依赖），比假设方差齐性更贴近数据；
#        robust —— 稳健经验贝叶斯，降低少数离群基因对先验方差的拉偏。
#      另跑经典 eBayes(fit) 作为敏感性分析，两者显著基因重叠率写入汇总。
#   4) 设计矩阵纳入 sex 与 age 作为协变量：本队列性别分布与分组强相关
#      （对照 female 63/136，DCM 19/82，ICM 14/95），不校正会引入性别假阳性；
#      年龄中位数在组间亦有差异（对照 52、DCM 53、ICM 60）。协变量是"加进模型"，
#      不是先做批次校正再检验（避免 Nygaard 2016 指出的"校正后重复检验"错误）。
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 0. 路径：优先用当前工作目录（项目根），否则用脚本所在目录的上一级
## ---------------------------------------------------------------------------
get_this_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  return(getwd())
}

wd <- getwd()
proj_dir <- if (dir.exists(file.path(wd, "clean_data"))) normalizePath(wd) else dirname(get_this_dir())
clean_dir <- file.path(proj_dir, "clean_data")
res_dir   <- file.path(proj_dir, "results", "DEG")
rep_dir   <- file.path(proj_dir, "report")

for (d in c(res_dir, rep_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(clean_dir)) stop("clean_data/ not found. Run 02_clean.R first, or run this script from the project root.")

## ---------------------------------------------------------------------------
## 1. 参数
## ---------------------------------------------------------------------------
ACCESSION <- "GSE57338"

# 分组定义：对照 = NF（non-failing，非衰竭心脏）；疾病 = DCM + ICM（心肌病）
CTRL_LABEL   <- "NF"
DISEASE_GRPS <- c("DCM", "ICM")

# 显著性阈值：|log2FC| > 0.585 约等于 FC > 1.5；padj 用 BH 校正
LFC_CUT  <- 0.585
PADJ_CUT <- 0.05

# 若主阈值筛不到基因，按此序列逐级放宽（脚本会自动记录实际使用的阈值）
LFC_FALLBACK  <- c(0.585, 0.500, 0.400, 0.300)
PADJ_FALLBACK <- c(0.05, 0.10)

# limma 参数
USE_TREND  <- TRUE   # 方差随表达强度变化
USE_ROBUST <- TRUE   # 稳健经验贝叶斯
USE_COVARS <- TRUE   # 模型中纳入 sex + age 协变量
ADJ_METHOD <- "BH"   # 多重检验校正

## ---------------------------------------------------------------------------
## 2. 小工具
## ---------------------------------------------------------------------------
msg <- function(...) cat(..., "\n", sep = "")
hr  <- function(ch = "-", n = 72) cat(strrep(ch, n), "\n", sep = "")

# 统一以 UTF-8 写文本
write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
}

## ---------------------------------------------------------------------------
## 3. 读入数据并做严格对齐校验
## ---------------------------------------------------------------------------
msg("Reading expression matrix and sample metadata ...")
expr_dt <- data.table::fread(file.path(clean_dir, paste0(ACCESSION, "_clean_expression_matrix.csv")),
                             data.table = FALSE)
gene_id <- expr_dt[[1]]
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
rownames(expr) <- gene_id
rm(expr_dt); invisible(gc())

meta <- read.csv(file.path(clean_dir, paste0(ACCESSION, "_sample_info.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)

# 3.1 维度与唯一性
stopifnot(!anyDuplicated(gene_id))
stopifnot(ncol(expr) == nrow(meta))

# 3.2 列顺序必须与 meta 行顺序完全一致（否则样本会错配 —— 最致命的一类 bug）
if (!identical(colnames(expr), meta$sample)) {
  if (setequal(colnames(expr), meta$sample)) {
    meta <- meta[match(colnames(expr), meta$sample), , drop = FALSE]
    warning("Sample order differed between matrix and metadata; re-ordered by column names.")
  } else {
    stop("Column names of the matrix do not match sample_info$sample.")
  }
}

# 3.3 表达值范围与缺失
if (anyNA(expr)) stop("Expression matrix contains NA. Fix it in 02_clean.R first.")
expr_range <- range(expr)
if (max(expr, na.rm = TRUE) > 30) warning("Max expression > 30; confirm the data are log2-transformed.")

msg(sprintf("Matrix: %d genes x %d samples; expression range %.3f - %.3f (log2)",
            nrow(expr), ncol(expr), expr_range[1], expr_range[2]))

## ---------------------------------------------------------------------------
## 4. 构建分组与协变量
## ---------------------------------------------------------------------------
# 4.1 二分类主分组：Control（NF） vs Disease（DCM + ICM）
#     因子水平显式指定，保证 limma 的对比方向是 "Disease - Control"（疾病组升高为正）
meta$group2 <- factor(ifelse(meta$group == CTRL_LABEL, "Control", "Disease"),
                      levels = c("Control", "Disease"))

# 4.2 三分类亚型分组（仅用于补充分析，不改变主分析）
meta$subtype <- factor(ifelse(meta$group == CTRL_LABEL, "Control", meta$group),
                       levels = c("Control", "DCM", "ICM"))

# 4.3 性别协变量
meta$sex <- factor(meta$gender, levels = c("female", "male"))

# 4.4 年龄协变量：0 岁在本数据中属录入缺失（仅 1 例，DCM 组），置 NA 后用中位数填补
meta$age_num <- suppressWarnings(as.numeric(meta$age))
n_age_zero <- sum(meta$age_num == 0, na.rm = TRUE)
meta$age_num[!is.na(meta$age_num) & meta$age_num == 0] <- NA
n_age_na <- sum(is.na(meta$age_num))
if (n_age_na > 0) {
  meta$age_num[is.na(meta$age_num)] <- median(meta$age_num, na.rm = TRUE)
}

msg("")
msg("Group sizes (Control = NF, Disease = DCM + ICM):")
print(table(meta$group2))
msg("")
msg("Sex x group:")
print(table(meta$sex, meta$group2))
msg("")
med_age <- tapply(meta$age_num, meta$group, median)
age_desc <- paste(sprintf("%s=%.0f", names(med_age), as.numeric(med_age)), collapse = ", ")
msg(sprintf("Median age by subgroup: %s; age==0 set to NA: %d; imputed with median: %d",
            age_desc, n_age_zero, n_age_na))

## ---------------------------------------------------------------------------
## 5. limma 拟合与差异检验
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(limma))

#' 拟合 limma 模型并返回指定系数的 topTable
#'
#' @param expr_mat  基因 × 样本 的 log2 表达矩阵
#' @param meta_df   样本表（行顺序与 expr_mat 列顺序一致）
#' @param form      model.matrix 用的公式
#' @param coef_name 要提取的系数名（必须存在于 design 的列名中）
#' @param trend     是否允许方差-表达强度趋势
#' @param robust    是否使用稳健经验贝叶斯
#' @param adj       多重检验校正方法
fit_limma <- function(expr_mat, meta_df, form, coef_name,
                      trend = TRUE, robust = TRUE, adj = "BH") {
  design <- model.matrix(form, data = meta_df)
  if (!coef_name %in% colnames(design)) {
    stop("Coefficient ", coef_name, " not in design. Available: ",
         paste(colnames(design), collapse = ", "))
  }
  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit, trend = trend, robust = robust)
  topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = adj)
}

#' topTable → 统一命名的结果表（满足"基因名 / log2FoldChange / P / padj"的最低要求）
#' 注意：log2FC 先四舍五入到 4 位小数，再做显著性判定，
#'       保证 csv 里"按该列筛选"能精确复现显著基因表（避免 0.5850 这类边界不一致）
tidy_tt <- function(tt) {
  data.frame(
    GeneSymbol       = rownames(tt),
    log2FoldChange   = round(tt$logFC, 4),
    AveExpr          = round(tt$AveExpr, 4),
    t                = round(tt$t, 4),
    P.Value          = tt$P.Value,
    padj             = tt$adj.P.Val,
    B                = round(tt$B, 4),
    stringsAsFactors = FALSE
  )
}

#' 按阈值标记显著性与上/下调方向
flag_sig <- function(df, lfc_cut, padj_cut) {
  sig <- !is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) > lfc_cut
  df$Regulation <- ifelse(sig & df$log2FoldChange > 0, "Up",
                          ifelse(sig & df$log2FoldChange < 0, "Down", "NS"))
  df$Significant <- sig
  df
}

# 5.1 主模型：Disease vs Control，校正性别与年龄
form_main <- if (USE_COVARS) ~ group2 + sex + age_num else ~ group2
msg("")
hr("=")
msg("Primary model: ", paste(deparse(form_main), collapse = ""))
msg("  Contrast: Disease - Control (log2FC > 0 means up-regulated in disease)")

tt_main  <- fit_limma(expr, meta, form_main, "group2Disease",
                      trend = USE_TREND, robust = USE_ROBUST, adj = ADJ_METHOD)
res_main <- tidy_tt(tt_main)

# 5.2 显著性筛选（含阈值的自动放宽逻辑）
lfc_used  <- LFC_CUT
padj_used <- PADJ_CUT
res_main  <- flag_sig(res_main, lfc_used, padj_used)
n_sig     <- sum(res_main$Significant)

if (n_sig == 0) {
  msg("")
  msg("[NOTE] No gene passed |log2FC| > ", LFC_CUT, " & padj < ", PADJ_CUT,
      "; relaxing thresholds stepwise ...")
  for (p in PADJ_FALLBACK) {
    for (l in LFC_FALLBACK[-1]) {
      tmp <- flag_sig(res_main, l, p)
      if (sum(tmp$Significant) > 0) {
        res_main <- tmp; lfc_used <- l; padj_used <- p; n_sig <- sum(tmp$Significant)
        msg(sprintf("  relaxed to |log2FC| > %.3f and padj < %.2f -> %d genes", l, p, n_sig))
        break
      }
    }
    if (n_sig > 0) break
  }
}

sig_main <- res_main[res_main$Significant, , drop = FALSE]
sig_main <- sig_main[order(sig_main$padj, -abs(sig_main$log2FoldChange)), , drop = FALSE]
n_up     <- sum(sig_main$Regulation == "Up")
n_down   <- sum(sig_main$Regulation == "Down")

msg("")
msg(sprintf("Genes tested: %d | significant DE genes: %d (Up %d / Down %d) | cutoffs: |log2FC| > %.3f, padj < %.2f",
            nrow(res_main), n_sig, n_up, n_down, lfc_used, padj_used))

## ---------------------------------------------------------------------------
## 6. 补充分析与敏感性检验
## ---------------------------------------------------------------------------
msg("")
msg("Supplementary: DCM / ICM subtype vs Control ...")
form_sub <- if (USE_COVARS) ~ subtype + sex + age_num else ~ subtype
res_dcm <- flag_sig(tidy_tt(fit_limma(expr, meta, form_sub, "subtypeDCM",
                                      USE_TREND, USE_ROBUST, ADJ_METHOD)),
                    lfc_used, padj_used)
res_icm <- flag_sig(tidy_tt(fit_limma(expr, meta, form_sub, "subtypeICM",
                                      USE_TREND, USE_ROBUST, ADJ_METHOD)),
                    lfc_used, padj_used)

msg("Sensitivity 1: no covariates (~ group2) ...")
res_nocov <- flag_sig(tidy_tt(fit_limma(expr, meta, ~ group2, "group2Disease",
                                        USE_TREND, USE_ROBUST, ADJ_METHOD)),
                      lfc_used, padj_used)

msg("Sensitivity 2: classic eBayes (trend = FALSE, robust = FALSE) ...")
res_classic <- flag_sig(tidy_tt(fit_limma(expr, meta, form_main, "group2Disease",
                                          trend = FALSE, robust = FALSE, adj = ADJ_METHOD)),
                        lfc_used, padj_used)

# 与主分析显著基因的重叠率（以主分析为分母）
overlap_rate <- function(a, b) {
  sa <- a$GeneSymbol[a$Significant]; sb <- b$GeneSymbol[b$Significant]
  sprintf("%d / %d (%.1f%%)", length(intersect(sa, sb)), length(sa),
          100 * length(intersect(sa, sb)) / max(1, length(sa)))
}

## ---------------------------------------------------------------------------
## 7. 输出结果表
## ---------------------------------------------------------------------------
# 只对效应量类数值四舍五入；P 值与 padj 保留全精度，便于用户自行改阈值
round_cols <- function(df) {
  for (cc in c("log2FoldChange", "AveExpr", "t", "B")) {
    if (cc %in% names(df)) df[[cc]] <- round(df[[cc]], 4)
  }
  df
}

f_all  <- file.path(res_dir, paste0(ACCESSION, "_DEGall.csv"))
f_sig  <- file.path(res_dir, paste0(ACCESSION, "_DEG_significant.csv"))
f_up   <- file.path(res_dir, paste0(ACCESSION, "_DEG_up.csv"))
f_down <- file.path(res_dir, paste0(ACCESSION, "_DEG_down.csv"))
f_sum  <- file.path(res_dir, paste0(ACCESSION, "_DEG_summary.txt"))

write.csv(round_cols(res_main), f_all, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(sig_main), f_sig, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(sig_main[sig_main$Regulation == "Up", , drop = FALSE]),
          f_up, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(sig_main[sig_main$Regulation == "Down", , drop = FALSE]),
          f_down, row.names = FALSE, fileEncoding = "UTF-8")

write.csv(round_cols(res_dcm), file.path(res_dir, paste0(ACCESSION, "_DEGall_DCM_vs_Control.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(res_icm), file.path(res_dir, paste0(ACCESSION, "_DEGall_ICM_vs_Control.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(res_nocov), file.path(res_dir, paste0(ACCESSION, "_DEGall_noCovariate.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(round_cols(res_classic), file.path(res_dir, paste0(ACCESSION, "_DEGall_classicEB.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------------------------------------------------------
## 8. 汇总 txt（参数 + 结果 + 方法段落）
## ---------------------------------------------------------------------------
top15 <- head(sig_main, 15)
top_lines <- c(
  sprintf("%-12s %10s %10s %12s %12s", "GeneSymbol", "log2FC", "AveExpr", "P.Value", "padj"),
  sprintf("%-12s %10s %10s %12s %12s", strrep("-", 12), strrep("-", 10),
          strrep("-", 10), strrep("-", 12), strrep("-", 12)),
  sprintf("%-12s %10.3f %10.3f %12.3e %12.3e",
          top15$GeneSymbol, top15$log2FoldChange, top15$AveExpr,
          top15$P.Value, top15$padj)
)

summary_lines <- c(
  "================================================================",
  paste0(ACCESSION, " differential expression summary (limma)"),
  "================================================================",
  paste0("Generated    : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " (CST)"),
  paste0("Dataset      : ", ACCESSION, " / GPL11532 (Affymetrix HuGene-1.1-ST)"),
  "Data type    : microarray, RMA-normalised log2 (continuous, not counts)",
  "",
  "-------------------------- Groups --------------------------",
  sprintf("Control      : %s, n = %d", CTRL_LABEL, sum(meta$group2 == "Control")),
  sprintf("Disease      : %s, n = %d", paste(DISEASE_GRPS, collapse = "+"),
          sum(meta$group2 == "Disease")),
  sprintf("  DCM n = %d, ICM n = %d", sum(meta$group == "DCM"), sum(meta$group == "ICM")),
  "Contrast     : Disease - Control (log2FC > 0 = higher in disease)",
  "",
  "-------------------------- Method --------------------------",
  "Software     : R + limma (lmFit -> eBayes)",
  paste0("Design       : ", paste(deparse(form_main), collapse = ""),
         if (USE_COVARS) "  (sex and age as covariates)" else "  (no covariates)"),
  sprintf("eBayes       : trend = %s, robust = %s", USE_TREND, USE_ROBUST),
  sprintf("P adjustment : %s (Benjamini-Hochberg)", ADJ_METHOD),
  "Model        : per-gene linear model on log2 expression, empirical Bayes variance shrinkage",
  "",
  "-------------------------- Cutoffs -------------------------",
  sprintf("Significance : |log2FC| > %.3f (FC > %.2f) AND padj < %.2f",
          lfc_used, 2^lfc_used, padj_used),
  sprintf("Relaxed?     : %s", ifelse(lfc_used == LFC_CUT && padj_used == PADJ_CUT,
                                     "no (primary cutoffs were sufficient)", "yes")),
  sprintf("Genes tested : %d", nrow(res_main)),
  "",
  "-------------------------- Results -------------------------",
  sprintf("Significant  : %d", n_sig),
  sprintf("  Up         : %d", n_up),
  sprintf("  Down       : %d", n_down),
  sprintf("padj < 0.05  : %d", sum(res_main$padj < 0.05, na.rm = TRUE)),
  sprintf("P.Value<0.01 : %d", sum(res_main$P.Value < 0.01, na.rm = TRUE)),
  sprintf("|log2FC|>1   : %d", sum(abs(res_main$log2FoldChange) > 1 &
                                   res_main$padj < 0.05, na.rm = TRUE)),
  sprintf("Smallest padj: %.3e (%s)", min(res_main$padj, na.rm = TRUE),
          res_main$GeneSymbol[which.min(res_main$padj)]),
  "",
  "Top 15 genes by padj:",
  paste0("  ", top_lines),
  "",
  "---------------- Supplementary / sensitivity ----------------",
  sprintf("DCM vs Control : %d significant (Up %d / Down %d)",
          sum(res_dcm$Significant), sum(res_dcm$Regulation == "Up"),
          sum(res_dcm$Regulation == "Down")),
  sprintf("ICM vs Control : %d significant (Up %d / Down %d)",
          sum(res_icm$Significant), sum(res_icm$Regulation == "Up"),
          sum(res_icm$Regulation == "Down")),
  sprintf("No covariates  : %d significant, overlap with primary %s",
          sum(res_nocov$Significant), overlap_rate(res_main, res_nocov)),
  sprintf("Classic eBayes : %d significant, overlap with primary %s",
          sum(res_classic$Significant), overlap_rate(res_main, res_classic)),
  sprintf("DCM/ICM shared : %d", length(intersect(
    res_dcm$GeneSymbol[res_dcm$Significant], res_icm$GeneSymbol[res_icm$Significant]))),
  "",
  "-------------------------- Files ---------------------------",
  "GSE57338_DEGall.csv              all genes (GeneSymbol/log2FC/AveExpr/t/P.Value/padj/B)",
  "GSE57338_DEG_significant.csv     significant genes (same columns + Regulation)",
  "GSE57338_DEG_up.csv              significant up-regulated genes",
  "GSE57338_DEG_down.csv            significant down-regulated genes",
  "GSE57338_DEGall_DCM_vs_Control.csv / _ICM_vs_Control.csv   subtype analyses",
  "GSE57338_DEGall_noCovariate.csv  sensitivity analysis (no covariates)",
  "GSE57338_DEGall_classicEB.csv    sensitivity analysis (classic eBayes)",
  "GSE57338_volcano.pdf / .png      volcano plot",
  "GSE57338_heatmap.pdf / .png      top-50 DE gene heatmap",
  "GSE57338_QC_diagnostics.pdf/.png PCA + p-value histogram",
  "",
  "================================================================"
)
write_utf8(summary_lines, f_sum)

## ---------------------------------------------------------------------------
## 9. 会话信息（复现必需）
## ---------------------------------------------------------------------------
si_file <- file.path(rep_dir, "sessionInfo_04_deg.txt")
con <- file(si_file, "w", encoding = "UTF-8")
writeLines(capture.output(sessionInfo()), con)
close(con)

msg("")
hr("=")
msg("Done. Output directory: ", res_dir)
msg("  DEGall   : ", basename(f_all))
msg("  Sig genes: ", basename(f_sig), "  (n = ", n_sig, ")")
msg("  Summary  : ", basename(f_sum))
hr("=")
