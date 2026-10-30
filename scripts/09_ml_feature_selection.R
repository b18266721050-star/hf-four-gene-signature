###############################################################################
## 09_ml_feature_selection.R
## GSE57338 —— LASSO / SVM-RFE / 随机森林 三种方法筛选特征基因
##
## 输入：
##   clean_data/GSE57338_clean_expression_matrix.csv  (基因 × 样本, log2)
##   clean_data/GSE57338_sample_info.csv              (样本表, 用 group 列)
##   clean_data/GSE57338candidate_genes.csv           (候选基因, 若缺失则按规则重建)
## 输出：
##   results/ML/**  (CSV / PDF / RDS) + report/sessionInfo_09_ml.txt
##
## 运行：Rscript scripts/09_ml_feature_selection.R
## 说明：脚本可从头独立跑通；所有随机过程使用同一个固定种子 SEED。
###############################################################################

rm(list = ls())

## ===========================================================================
## 0. 编码守卫 + 全局参数（换机器只需改 proj_dir 一处）
## ===========================================================================
## 某些 shell / 调度器会传入 LC_ALL=C.UTF-8，使 R 落入 ASCII locale（l10n_info()$UTF-8
## 为 FALSE），此时中文无法写入文件（writeLines 报 "invalid char string"、CSV 中文变 ?）。
## 这里显式切回系统 UTF-8 locale，保证中文注释、日志与表格都能正确落盘。
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("Chinese (Simplified)_China.utf8", "Chinese", "C.UTF-8", "en_US.UTF-8")) {
    ok <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error = function(e) ""))
    if (nzchar(ok) && isTRUE(l10n_info()$`UTF-8`)) break
  }
}
if (!isTRUE(l10n_info()$`UTF-8`)) {
  warning("当前 R 会话不是 UTF-8 locale，中文可能无法正确写入文件；建议在 UTF-8 终端下运行")
}

SEED <- 20261027L            # 固定随机种子（本次运行日期，全脚本唯一随机源）
N_FOLD <- 5L                 # 交叉验证折数（三种方法统一 5 折，便于横向比较）
CASE_LABEL <- "Disease"      # 阳性类标签
CTRL_LABEL <- "Control"      # 阴性类标签
CASE_CODE <- 1L              # 阳性编码
CTRL_CODE <- 0L              # 阴性编码
LFC_CUT <- 0.585             # 与 04_deg.R 完全一致的效应量阈值 (FC>1.5)
PADJ_CUT <- 0.05             # 与 04_deg.R 完全一致的显著性阈值
RF_NTREE <- 500L             # 随机森林树数
RF_TOPN <- 30L               # 随机森林重要性 Top-N（补充规则）
SVM_COST <- 1                # 线性 SVM 的 cost（固定，不做调参，避免额外随机性）
SVM_DROP_PROP <- 0.10        # SVM-RFE 每轮剔除比例
SVM_MIN_FEAT <- 2L           # SVM-RFE 迭代终止的最少特征数
LASSO_MEASURE <- "deviance"  # 交叉验证准则（二分类 deviance）

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()   # 换机器时回退到当前目录

in_expr  <- file.path(proj_dir, "clean_data", "GSE57338_clean_expression_matrix.csv")
in_meta  <- file.path(proj_dir, "clean_data", "GSE57338_sample_info.csv")
in_cand  <- file.path(proj_dir, "clean_data", "GSE57338candidate_genes.csv")
in_dcm   <- file.path(proj_dir, "results", "DEG", "GSE57338_DEGall_DCM_vs_Control.csv")
in_icm   <- file.path(proj_dir, "results", "DEG", "GSE57338_DEGall_ICM_vs_Control.csv")

out_root <- file.path(proj_dir, "results", "ML")
d_mat    <- file.path(out_root, "01_matrix")
d_lasso  <- file.path(out_root, "02_lasso")
d_svm    <- file.path(out_root, "03_svmrfe")
d_rf     <- file.path(out_root, "04_randomforest")
d_cmp    <- file.path(out_root, "05_comparison")
d_rep    <- file.path(proj_dir, "report")
for (d in c(out_root, d_mat, d_lasso, d_svm, d_rf, d_cmp, d_rep)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- 运行日志：屏幕 + 文件双写（二进制连接写 UTF-8 原始字节，避免转码失败）----
log_con <- file(file.path(out_root, "GSE57338_ML_run_log.txt"), open = "wb")
msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  writeLines(txt, log_con, useBytes = TRUE)
  flush(log_con)
}

## ---- 统一的 CSV 写出函数：固定 UTF-8，保证中文表头/字段可读 ----
## 注意：fileEncoding 必须是形参（而非写死在内部再接收 ...），
##       否则调用处再传 fileEncoding 会触发 "matched by multiple actual arguments"
wcsv <- function(x, file, row.names = FALSE, fileEncoding = "UTF-8", ...) {
  utils::write.csv(x, file, row.names = row.names, fileEncoding = fileEncoding, ...)
}

## ---- 异常登记表：任何"不满足预期"的情况都必须落到这里，禁止静默 -------------
issues <- data.frame(step = character(), level = character(), item = character(),
                     detail = character(), stringsAsFactors = FALSE)
add_issue <- function(step, level, item, detail) {
  issues <<- rbind(issues, data.frame(step = step, level = level, item = item,
                                      detail = detail, stringsAsFactors = FALSE))
  msg(sprintf("[%s] %s | %s | %s", level, step, item, detail))
}

## ---- 统一随机种子重置（显式指定三类采样器，保证跨 R 版本可复现）------------
reset_seed <- function() {
  set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
}

msg("==================================================================")
msg("GSE57338 特征基因筛选：LASSO / SVM-RFE / Random Forest")
msg("开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
msg("随机种子 SEED = ", SEED, "；交叉验证折数 = ", N_FOLD)
msg("==================================================================")

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(ggplot2)
  library(dplyr)
})
has_venn  <- requireNamespace("VennDiagram", quietly = TRUE)
has_upset <- requireNamespace("UpSetR", quietly = TRUE)
if (!has_venn)  add_issue("环境", "WARN", "VennDiagram", "包不可用，Venn 图将跳过")
if (!has_upset) add_issue("环境", "WARN", "UpSetR", "包不可用，UpSet 图将跳过")

###############################################################################
## 1. 读取数据并进行严格核对
###############################################################################
msg("\n---- 1. 数据读取与核对 ----")

expr <- read.csv(in_expr, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(in_meta, check.names = FALSE, stringsAsFactors = FALSE)

if (!"GeneSymbol" %in% colnames(expr)) stop("表达矩阵缺少 GeneSymbol 列")
if (!"group" %in% colnames(meta))      stop("样本表缺少 group 列")
if (!"sample" %in% colnames(meta))     stop("样本表缺少 sample 列")

## 1.1 维度
msg(sprintf("表达矩阵：%d 基因 × %d 样本（含 GeneSymbol 列）", nrow(expr), ncol(expr) - 1))
msg(sprintf("样本表  ：%d 行 × %d 列", nrow(meta), ncol(meta)))

## 1.2 基因名唯一性与完整性
if (anyDuplicated(expr$GeneSymbol) > 0) {
  add_issue("数据核对", "ERROR", "GeneSymbol 重复",
            sprintf("%d 个重复基因名，后续按行号去重保留第一条", sum(duplicated(expr$GeneSymbol))))
  expr <- expr[!duplicated(expr$GeneSymbol), , drop = FALSE]
}
if (any(is.na(expr$GeneSymbol) | expr$GeneSymbol == "")) {
  add_issue("数据核对", "ERROR", "GeneSymbol 空值",
            sprintf("%d 行基因名为空，已剔除", sum(is.na(expr$GeneSymbol) | expr$GeneSymbol == "")))
  expr <- expr[!(is.na(expr$GeneSymbol) | expr$GeneSymbol == ""), , drop = FALSE]
}

## 1.3 样本名一一对应（顺序也必须一致，否则标签会错配）
mat_samples <- colnames(expr)[-1]
if (!identical(mat_samples, meta$sample)) {
  only_mat <- setdiff(mat_samples, meta$sample)
  only_meta <- setdiff(meta$sample, mat_samples)
  add_issue("数据核对", "ERROR", "样本名/顺序不一致",
            sprintf("仅矩阵有 %d 个，仅样本表有 %d 个，是否同集合=%s",
                    length(only_mat), length(only_meta), setequal(mat_samples, meta$sample)))
  if (!setequal(mat_samples, meta$sample)) {
    stop("表达矩阵列与样本表 sample 不是同一集合，无法安全配对，脚本终止")
  }
  ## 集合相同仅顺序不同：按矩阵列序重排样本表（安全修正，并留痕）
  meta <- meta[match(mat_samples, meta$sample), , drop = FALSE]
  add_issue("数据核对", "FIX", "样本顺序",
            "集合一致但顺序不同，已按表达矩阵列序重排样本表")
}
if (anyDuplicated(meta$sample) > 0) {
  add_issue("数据核对", "ERROR", "样本 ID 重复", "样本表存在重复 sample，脚本终止")
  stop("样本表 sample 列存在重复，无法确定唯一标签")
}
msg("样本名对齐检查：通过（矩阵列序 == 样本表行序）")

## 1.4 分组结构
grp_tab <- table(meta$group)
msg("group 列分布：")
for (g in names(grp_tab)) msg(sprintf("  %-6s %4d", g, grp_tab[[g]]))

nf_groups <- names(grp_tab)[grepl("^(NF|Control|control|Normal|normal)", names(grp_tab))]
dis_groups <- setdiff(names(grp_tab), nf_groups)
if (length(nf_groups) == 0) stop("group 列中找不到对照组（NF/Control/Normal）")
if (length(dis_groups) == 0) stop("group 列中找不到疾病组")

## 1.5 二分类编码（明确定义，写入标签表）
meta$group_binary <- ifelse(meta$group %in% nf_groups, CTRL_LABEL, CASE_LABEL)
meta$y <- ifelse(meta$group_binary == CASE_LABEL, CASE_CODE, CTRL_CODE)
msg(sprintf("二分类编码：Case = %s (code %d, 组: %s)；Control = %s (code %d, 组: %s)",
            CASE_LABEL, CASE_CODE, paste(dis_groups, collapse = "+"),
            CTRL_LABEL, CTRL_CODE, paste(nf_groups, collapse = "+")))
msg(sprintf("样本数：Case %d / Control %d（不平衡比 %.2f : 1）",
            sum(meta$y == CASE_CODE), sum(meta$y == CTRL_CODE),
            sum(meta$y == CASE_CODE) / sum(meta$y == CTRL_CODE)))

if (any(table(meta$y) < N_FOLD)) {
  add_issue("数据核对", "ERROR", "某类样本数 < 折数",
            "分层交叉验证无法保证每折含两类样本")
}

## 1.6 数值完整性
num_mat <- as.matrix(expr[, -1, drop = FALSE])
storage.mode(num_mat) <- "double"
if (anyNA(num_mat)) {
  add_issue("数据核对", "ERROR", "表达矩阵含 NA",
            sprintf("%d 个缺失值，按基因中位数填补", sum(is.na(num_mat))))
  for (j in seq_len(ncol(num_mat))) {
    idx <- is.na(num_mat[, j])
    if (any(idx)) num_mat[idx, j] <- median(num_mat[, j], na.rm = TRUE)
  }
}
msg(sprintf("表达值范围：%.3f ~ %.3f（log2 尺度，无缺失）", min(num_mat), max(num_mat)))

###############################################################################
## 2. 候选基因：读取 or 按明确规则重建
###############################################################################
msg("\n---- 2. 候选基因集合 ----")

cand_source <- "原始文件"
if (file.exists(in_cand)) {
  cand_raw <- read.csv(in_cand, check.names = FALSE, stringsAsFactors = FALSE)
  ## 兼容多种列名/单列无表头格式
  cand_col <- intersect(c("GeneSymbol", "gene", "Gene", "symbol", "SYMBOL"), colnames(cand_raw))
  if (length(cand_col) > 0) {
    cand_genes <- unique(trimws(as.character(cand_raw[[cand_col[1]]])))
  } else {
    cand_genes <- unique(trimws(as.character(cand_raw[[1]])))
    add_issue("候选基因", "WARN", "候选文件列名",
              sprintf("未识别到标准基因列，改用第 1 列：%s", colnames(cand_raw)[1]))
  }
  cand_genes <- cand_genes[!is.na(cand_genes) & cand_genes != ""]
  msg(sprintf("读取候选基因文件：%s（%d 个唯一基因）", basename(in_cand), length(cand_genes)))
  ## 若该文件本身是本脚本早前重建的产物，必须在本次运行中继续留痕：
  ## 否则后续重跑只看得到"无异常"，会掩盖"原始候选清单缺失"这一事实。
  if ("source" %in% colnames(cand_raw) && any(grepl("rebuilt", cand_raw$source))) {
    add_issue("候选基因", "NOTE", "候选基因为脚本重建产物",
              "该文件由本脚本在缺少原始候选清单时生成（DCM/ICM 亚组显著基因交集），并非上游原始文件；如另有来源请替换后重跑")
  }
} else {
  ## ---------------------------------------------------------------------
  ## 【重要】原始 candidate_genes.csv 不存在。
  ## 不做任何"凑数"或放宽阈值：严格沿用 04_deg.R 的判定标准
  ##   padj < 0.05 且 |log2FC| > 0.585
  ## 取 DCM vs Control 与 ICM vs Control 两个亚组各自显著基因的交集，
  ## 作为"交集基因"候选集，并把重建结果写回候选文件，保证下游可溯源。
  ## ---------------------------------------------------------------------
  add_issue("候选基因", "WARN", "candidate_genes.csv 缺失",
            paste0("未找到 ", in_cand, "；按 04_deg.R 同阈值重建（DCM/ICM 亚组显著基因交集）"))
  if (!file.exists(in_dcm) || !file.exists(in_icm)) {
    stop("候选基因文件缺失，且 DEG 亚组结果也不存在，无法重建候选集")
  }
  deg_dcm <- read.csv(in_dcm, stringsAsFactors = FALSE)
  deg_icm <- read.csv(in_icm, stringsAsFactors = FALSE)
  sig_of <- function(df) {
    df$GeneSymbol[!is.na(df$padj) & df$padj < PADJ_CUT & abs(df$log2FoldChange) > LFC_CUT]
  }
  s_dcm <- unique(sig_of(deg_dcm))
  s_icm <- unique(sig_of(deg_icm))
  cand_genes <- sort(intersect(s_dcm, s_icm))
  msg(sprintf("DCM_vs_Control 显著基因 %d 个；ICM_vs_Control 显著基因 %d 个；交集 %d 个",
              length(s_dcm), length(s_icm), length(cand_genes)))
  if (length(cand_genes) == 0) {
    add_issue("候选基因", "ERROR", "重建交集为空",
              "两亚组显著基因无交集，脚本终止（不放宽阈值）")
    stop("候选基因交集为空")
  }
  rebuild <- data.frame(
    GeneSymbol = cand_genes,
    source = "rebuilt_from_DEG",
    rule = sprintf("padj<%.2f & |log2FC|>%.3f; intersect(DCM_vs_Control, ICM_vs_Control)", PADJ_CUT, LFC_CUT),
    n_DCM_sig = length(s_dcm),
    n_ICM_sig = length(s_icm),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  wcsv(rebuild, in_cand, row.names = FALSE, fileEncoding = "UTF-8")
  cand_source <- "重建（DEG 亚组交集）"
  msg(sprintf("已写入重建的候选基因文件：%s", in_cand))
}

## 2.1 候选基因与表达矩阵的匹配核查
cand_in_expr <- intersect(cand_genes, expr$GeneSymbol)
cand_missing <- setdiff(cand_genes, expr$GeneSymbol)
msg(sprintf("候选基因 %d 个，其中 %d 个在表达矩阵中匹配到；缺失 %d 个",
            length(cand_genes), length(cand_in_expr), length(cand_missing)))
if (length(cand_missing) > 0) {
  add_issue("候选基因", "ERROR", "候选基因缺失",
            sprintf("%d 个候选基因不在表达矩阵：%s", length(cand_missing),
                    paste(head(cand_missing, 20), collapse = ", ")))
  wcsv(data.frame(GeneSymbol = cand_missing),
            file.path(d_mat, "GSE57338_candidate_genes_missing.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}
if (length(cand_in_expr) < 3) {
  stop(sprintf("可用候选基因仅 %d 个，无法进行特征筛选", length(cand_in_expr)))
}
if (length(cand_missing) > 0) {
  writeLines(sprintf("候选基因共 %d 个，表达矩阵缺失 %d 个（未做任何补全，直接剔除并在报告中列明）",
                     length(cand_genes), length(cand_missing)),
             file.path(d_mat, "GSE57338_candidate_genes_missing_note.txt"))
}

###############################################################################
## 3. 构建"样本 × 基因"的机器学习矩阵
###############################################################################
msg("\n---- 3. 构建机器学习矩阵 ----")

idx_row <- match(cand_in_expr, expr$GeneSymbol)
X_raw <- t(num_mat[idx_row, , drop = FALSE])          # 转置：行=样本，列=基因
colnames(X_raw) <- cand_in_expr
rownames(X_raw) <- mat_samples

X_df <- as.data.frame(X_raw, check.names = FALSE)
X_df <- cbind(sample = rownames(X_raw),
              group = meta$group,
              group_binary = meta$group_binary,
              y = meta$y,
              X_df, stringsAsFactors = FALSE)
rownames(X_df) <- NULL

msg(sprintf("机器学习矩阵：%d 样本 × %d 基因（样本 × 基因）", nrow(X_raw), ncol(X_raw)))
msg(sprintf("列名前 5 个：%s", paste(head(colnames(X_raw), 5), collapse = ", ")))

## 3.1 保存矩阵与标签（原始 log2 值，未做任何中心化）
wcsv(X_df, file.path(d_mat, "GSE57338_ML_matrix_sample_by_gene.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(data.frame(sample = meta$sample, group = meta$group,
                     group_binary = meta$group_binary, y = meta$y,
                     stringsAsFactors = FALSE),
          file.path(d_mat, "GSE57338_ML_labels.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 3.2 基因层面的描述统计（供报告与复现核对）
gene_stat <- data.frame(
  GeneSymbol = colnames(X_raw),
  mean = colMeans(X_raw),
  sd = apply(X_raw, 2, sd),
  min = apply(X_raw, 2, min),
  max = apply(X_raw, 2, max),
  mean_Case = colMeans(X_raw[meta$y == CASE_CODE, , drop = FALSE]),
  mean_Control = colMeans(X_raw[meta$y == CTRL_CODE, , drop = FALSE]),
  stringsAsFactors = FALSE
)
gene_stat$diff <- gene_stat$mean_Case - gene_stat$mean_Control
wcsv(gene_stat, file.path(d_mat, "GSE57338_candidate_gene_stats.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 3.3 零方差检查（零方差基因对任何模型都无信息，但必须报告而非静默删除）
zero_var <- gene_stat$GeneSymbol[gene_stat$sd < 1e-8]
if (length(zero_var) > 0) {
  add_issue("建模矩阵", "WARN", "零方差基因",
            sprintf("%d 个基因方差≈0：%s", length(zero_var), paste(zero_var, collapse = ", ")))
  keep <- setdiff(colnames(X_raw), zero_var)
  X_raw <- X_raw[, keep, drop = FALSE]
  msg(sprintf("剔除零方差基因后：%d 样本 × %d 基因", nrow(X_raw), ncol(X_raw)))
}

## 3.4 标准化矩阵（z-score，供 SVM 与可视化使用）
## 说明：此处使用全样本均值/标准差，属"特征筛选"语境下的常规做法；
##       本分析不构成独立验证集上的预测性能评估，局限已在报告中明确说明。
X_z <- scale(X_raw)
X_z[is.na(X_z)] <- 0
wcsv(cbind(sample = rownames(X_z), as.data.frame(X_z, check.names = FALSE)),
          file.path(d_mat, "GSE57338_ML_matrix_zscore.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(data.frame(GeneSymbol = colnames(X_raw),
                     center = attr(X_z, "scaled:center"),
                     scale = attr(X_z, "scaled:scale"),
                     stringsAsFactors = FALSE),
          file.path(d_mat, "GSE57338_ML_zscore_parameters.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

y_num <- meta$y
y_fac <- factor(ifelse(y_num == CASE_CODE, CASE_LABEL, CTRL_LABEL),
                levels = c(CTRL_LABEL, CASE_LABEL))
names(y_num) <- rownames(X_raw)

## 3.5 统一的分层交叉验证折（三种方法共用，保证比较公平）
reset_seed()
folds <- caret::createFolds(y_fac, k = N_FOLD, list = TRUE, returnTrain = FALSE)
fold_id <- integer(nrow(X_raw))
for (i in seq_along(folds)) fold_id[folds[[i]]] <- i
msg("分层 5 折样本分布：")
print_tab <- table(fold = fold_id, class = y_fac)
for (i in seq_len(nrow(print_tab))) {
  msg(sprintf("  Fold %d : Control %3d / Case %3d", i, print_tab[i, 1], print_tab[i, 2]))
}
wcsv(data.frame(sample = rownames(X_raw), fold = fold_id, y = y_num,
                     group = meta$group),
          file.path(d_mat, "GSE57338_CV_folds.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

###############################################################################
## 4. LASSO 逻辑回归（二分类 + 5 折 CV 选 lambda）
###############################################################################
msg("\n---- 4. LASSO ----")
reset_seed()

## glmnet 自带标准化（standardize=TRUE），系数会还原到原始尺度
lasso_fit <- glmnet(X_raw, y_num, family = "binomial", alpha = 1,
                    standardize = TRUE, type.measure = LASSO_MEASURE)
reset_seed()
lasso_cv <- cv.glmnet(X_raw, y_num, family = "binomial", alpha = 1,
                      standardize = TRUE, type.measure = LASSO_MEASURE,
                      nfolds = N_FOLD, foldid = fold_id)

lam_min <- lasso_cv$lambda.min
lam_1se <- lasso_cv$lambda.1se
msg(sprintf("lambda.min = %.6f（对应 CV %s = %.4f）", lam_min,
            LASSO_MEASURE, lasso_cv$cvm[lasso_cv$lambda == lam_min]))
msg(sprintf("lambda.1se = %.6f（对应 CV %s = %.4f）", lam_1se,
            LASSO_MEASURE, lasso_cv$cvm[lasso_cv$lambda == lam_1se]))
msg(sprintf("lambda 路径长度：%d（%.4f ~ %.4f）",
            length(lasso_cv$lambda), min(lasso_cv$lambda), max(lasso_cv$lambda)))

## 4.1 交叉验证曲线（glmnet >= 4.0 的 plot.cv.glmnet 已去掉 xvar 参数，直接用默认 lambda 轴）
nz_min <- sum(as.matrix(coef(lasso_cv, s = "lambda.min"))[-1, 1] != 0)
nz_1se <- sum(as.matrix(coef(lasso_cv, s = "lambda.1se"))[-1, 1] != 0)

draw_lasso_cv <- function() {
  par(mar = c(4.6, 4.8, 5.4, 1.4))
  plot(lasso_cv)
  title(main = sprintf("LASSO %d-fold CV (%s), seed = %d", N_FOLD, LASSO_MEASURE, SEED),
        line = 4.1, cex.main = 1.05)
  mtext(sprintf("lambda.min = %.4f (%d genes) | lambda.1se = %.4f (%d genes)",
                lam_min, nz_min, lam_1se, nz_1se),
        side = 3, line = 3.0, cex = 0.78, col = "grey30")
}
pdf(file.path(d_lasso, "GSE57338_LASSO_cv_curve.pdf"), width = 7, height = 5.8)
draw_lasso_cv(); dev.off()
png(file.path(d_lasso, "GSE57338_LASSO_cv_curve.png"), width = 7 * 300, height = 5.8 * 300, res = 300)
draw_lasso_cv(); dev.off()

## 4.2 系数路径图（自绘：glmnet 5.x 的 plot.glmnet 内部 x 轴语义已变更，
##     直接用 abline(v = log(lambda)) 会把参考线画到坐标范围之外，故手绘保证坐标准确）
draw_coef_path <- function() {
  beta <- as.matrix(lasso_fit$beta)      # p × nlambda
  lx <- log(lasso_fit$lambda)            # x 轴：log(lambda)
  nz <- colSums(beta != 0)               # 每个 lambda 下的非零系数个数
  pal <- colorRampPalette(c("#2C6E9B", "#7FA8C4", "#C0392B", "#E0A07A"))(nrow(beta))
  par(mar = c(4.6, 4.8, 4.6, 1.4))
  matplot(lx, t(beta), type = "l", lty = 1, lwd = 0.9,
          col = adjustcolor(pal, 0.5), xaxt = "n",
          xlab = expression(log(lambda)), ylab = "Coefficient",
          xlim = rev(range(lx)), main = "")
  axis(1, at = pretty(range(lx)))
  abline(h = 0, col = "grey55", lty = 3)
  abline(v = log(lam_min), lty = 2, col = "#C0392B", lwd = 1.3)
  abline(v = log(lam_1se), lty = 3, col = "#2C6E9B", lwd = 1.3)
  ## 顶部轴：非零系数个数（左侧 lambda 大、非零少）
  ats <- pretty(range(lx))
  nz_at <- approx(rev(lx), rev(nz), xout = ats, rule = 2)$y
  axis(3, at = ats, labels = round(nz_at), tick = FALSE, line = -0.7, cex.axis = 0.78)
  mtext("Number of non-zero coefficients", side = 3, line = 1.3, cex = 0.72)
  title(main = "LASSO coefficient path", line = 3.1, cex.main = 1.05)
  legend("topleft",
         c(sprintf("lambda.min = %.4f (%d genes)", lam_min, nz_min),
           sprintf("lambda.1se = %.4f (%d genes)", lam_1se, nz_1se)),
         lty = c(2, 3), col = c("#C0392B", "#2C6E9B"), bty = "n", cex = 0.76)
}
pdf(file.path(d_lasso, "GSE57338_LASSO_coef_path.pdf"), width = 7, height = 5.6)
draw_coef_path(); dev.off()
png(file.path(d_lasso, "GSE57338_LASSO_coef_path.png"), width = 7 * 300, height = 5.6 * 300, res = 300)
draw_coef_path(); dev.off()

## 4.3 系数表与非零基因
coef_min <- as.matrix(coef(lasso_cv, s = "lambda.min"))
coef_1se <- as.matrix(coef(lasso_cv, s = "lambda.1se"))
coef_tab <- data.frame(
  GeneSymbol = rownames(coef_min),
  Coef_lambda_min = as.numeric(coef_min[, 1]),
  Coef_lambda_1se = as.numeric(coef_1se[rownames(coef_min), 1]),
  stringsAsFactors = FALSE
)
coef_tab <- coef_tab[coef_tab$GeneSymbol != "(Intercept)", , drop = FALSE]
coef_tab$AbsCoef <- abs(coef_tab$Coef_lambda_min)
coef_tab$Selected_lambda_min <- coef_tab$Coef_lambda_min != 0
coef_tab$Selected_lambda_1se <- coef_tab$Coef_lambda_1se != 0
coef_tab$Direction <- ifelse(coef_tab$Coef_lambda_min > 0, "Risk(up in Case)",
                      ifelse(coef_tab$Coef_lambda_min < 0, "Protective(down in Case)", "None"))
coef_tab <- coef_tab[order(-coef_tab$AbsCoef), , drop = FALSE]

wcsv(coef_tab, file.path(d_lasso, "GSE57338_LASSO_coef_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

lasso_genes <- coef_tab$GeneSymbol[coef_tab$Selected_lambda_min]
lasso_genes_1se <- coef_tab$GeneSymbol[coef_tab$Selected_lambda_1se]
wcsv(coef_tab[coef_tab$Selected_lambda_min, , drop = FALSE],
          file.path(d_lasso, "GSE57338_LASSO_nonzero_genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(data.frame(GeneSymbol = lasso_genes_1se),
          file.path(d_lasso, "GSE57338_LASSO_nonzero_genes_lambda1se.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

msg(sprintf("lambda.min 非零系数基因：%d 个（正/风险 %d，负/保护 %d）",
            length(lasso_genes), sum(coef_tab$Coef_lambda_min > 0),
            sum(coef_tab$Coef_lambda_min < 0)))
msg(sprintf("lambda.1se 非零系数基因：%d 个", length(lasso_genes_1se)))
msg(sprintf("LASSO 全数据表观 AUC：%.3f",
            as.numeric(pROC::auc(pROC::roc(y_num, predict(lasso_cv, X_raw, s = "lambda.min", type = "response")[, 1],
                                           quiet = TRUE, direction = "<")))))

if (length(lasso_genes) == 0) {
  add_issue("LASSO", "ERROR", "lambda.min 无非零系数",
            "模型完全收缩，该阈值下无特征基因（不调整 lambda 强行取基因）")
}

## 4.4 系数条形图（top 20）
top_lasso <- head(coef_tab[coef_tab$Selected_lambda_min, , drop = FALSE], 20)
if (nrow(top_lasso) > 0) {
  p_lasso <- ggplot(top_lasso, aes(x = reorder(GeneSymbol, Coef_lambda_min), y = Coef_lambda_min,
                                   fill = Coef_lambda_min > 0)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C6E9B"),
                      labels = c(`TRUE` = "Up in Case", `FALSE` = "Down in Case"), name = NULL) +
    labs(x = NULL, y = "LASSO coefficient (lambda.min)",
         title = "LASSO selected genes (top 20 by |coef|)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top")
  ggsave(file.path(d_lasso, "GSE57338_LASSO_coef_barplot.pdf"), p_lasso,
         width = 6.2, height = 6.0, device = "pdf")
  ggsave(file.path(d_lasso, "GSE57338_LASSO_coef_barplot.png"), p_lasso,
         width = 6.2, height = 6.0, dpi = 300)
}

saveRDS(list(fit = lasso_fit, cv = lasso_cv, lambda.min = lam_min,
             lambda.1se = lam_1se, coef = coef_tab, genes = lasso_genes,
             seed = SEED, foldid = fold_id),
        file.path(d_lasso, "GSE57338_LASSO_model.rds"))

###############################################################################
## 5. SVM-RFE（线性核，5 折 CV 递归消除）
###############################################################################
msg("\n---- 5. SVM-RFE ----")

## 每一轮：用当前特征集做 5 折 CV 估计错误率；用全数据线性 SVM 的 |w| 排序，
## 剔除 |w| 最小的 10%。按"淘汰越晚越重要"给出 rank。
rfe_svm <- function(x, y, foldid, drop_prop, min_feat, cost) {
  feats <- colnames(x)
  k <- length(unique(foldid))
  hist_n <- integer(); hist_err <- numeric(); hist_se <- numeric()
  hist_feats <- list(); hist_w <- list(); elim_log <- list()
  step <- 0L
  repeat {
    if (length(feats) <= min_feat) break
    step <- step + 1L
    errs <- numeric(k)
    for (i in seq_len(k)) {
      te <- which(foldid == i); tr <- setdiff(seq_len(nrow(x)), te)
      m <- e1071::svm(x[tr, feats, drop = FALSE], y[tr], kernel = "linear",
                      cost = cost, scale = TRUE)
      p <- predict(m, x[te, feats, drop = FALSE])
      errs[i] <- mean(p != y[te])
    }
    hist_n   <- c(hist_n, length(feats))
    hist_err <- c(hist_err, mean(errs))
    hist_se  <- c(hist_se, sd(errs) / sqrt(k))
    hist_feats[[step]] <- feats

    m_all <- e1071::svm(x[, feats, drop = FALSE], y, kernel = "linear",
                        cost = cost, scale = TRUE)
    w <- drop(t(m_all$coefs) %*% m_all$SV)
    names(w) <- feats
    hist_w[[step]] <- w

    n_drop <- max(1L, floor(length(feats) * drop_prop))
    if (length(feats) - n_drop < min_feat) n_drop <- length(feats) - min_feat
    drop_idx <- order(abs(w))[seq_len(n_drop)]
    elim_log[[step]] <- names(w)[drop_idx]
    feats <- setdiff(feats, names(w)[drop_idx])
  }
  ## 最终剩余特征（|w| 从大到小 → rank 1 起）
  m_fin <- e1071::svm(x[, feats, drop = FALSE], y, kernel = "linear",
                      cost = cost, scale = TRUE)
  w_fin <- drop(t(m_fin$coefs) %*% m_fin$SV); names(w_fin) <- feats
  hist_n <- c(hist_n, length(feats)); hist_err <- c(hist_err, NA_real_); hist_se <- c(hist_se, NA_real_)
  hist_feats[[length(hist_feats) + 1]] <- feats
  hist_w[[length(hist_w) + 1]] <- w_fin
  list(n = hist_n, err = hist_err, se = hist_se, feats = hist_feats,
       w = hist_w, elim = elim_log, final_w = w_fin, final_feats = feats)
}

reset_seed()
rfe <- rfe_svm(X_z, y_fac, fold_id, SVM_DROP_PROP, SVM_MIN_FEAT, SVM_COST)

## 5.1 交叉验证错误率曲线
## 注意：idx 必须保留 rfe$feats 的原始位置。表格按 n_features 排序后，
##       若直接用排序后的行号去索引 feats，会取到错误轮次的特征集（曾经的 bug）。
rfe_df <- data.frame(idx = seq_along(rfe$n), n_features = rfe$n,
                     cv_error = rfe$err, cv_se = rfe$se)
rfe_df <- rfe_df[order(rfe_df$n_features), , drop = FALSE]
wcsv(rfe_df[, c("n_features", "cv_error", "cv_se")],
     file.path(d_svm, "GSE57338_SVM_RFE_cv_error_profile.csv"),
     row.names = FALSE, fileEncoding = "UTF-8")

## 最优特征数：最小 CV 错误率；若并列，取其中特征数最多者（保留信息最充分且最稳）
valid <- rfe_df[!is.na(rfe_df$cv_error), , drop = FALSE]
min_err <- min(valid$cv_error)
cand_n <- valid$n_features[valid$cv_error == min_err]
n_opt <- max(cand_n)
if (length(cand_n) > 1) {
  add_issue("SVM-RFE", "NOTE", "最优特征数并列",
            sprintf("CV 错误率同为 %.4f 的特征数：%s；按'取最大 n'规则选定 %d",
                    min_err, paste(cand_n, collapse = ", "), n_opt))
}
opt_pos <- which(!is.na(rfe_df$cv_error) & rfe_df$n_features == n_opt)[1]
svm_genes <- rfe$feats[[ rfe_df$idx[opt_pos] ]]
## 一致性断言：最优特征集的基因数必须等于 n_opt，否则说明索引仍错位
if (length(svm_genes) != n_opt) {
  add_issue("SVM-RFE", "ERROR", "最优特征集长度与 n 不一致",
            sprintf("n_opt = %d，实际取出 %d 个基因", n_opt, length(svm_genes)))
  stop("SVM-RFE 最优特征集索引错位，脚本终止（不输出错误结果）")
}
msg(sprintf("CV 错误率最低 = %.4f（±%.4f），对应最优特征数 = %d",
            min_err, rfe_df$cv_se[opt_pos], n_opt))
msg(sprintf("最优特征集基因数：%d", length(svm_genes)))

## 绘图时剔除最后一轮（终止轮无 CV 误差，为 NA）
rfe_plot_df <- rfe_df[!is.na(rfe_df$cv_error), , drop = FALSE]
p_rfe <- ggplot(rfe_plot_df, aes(x = n_features, y = cv_error)) +
  geom_line(color = "#2C6E9B", linewidth = 0.8) +
  geom_point(color = "#2C6E9B", size = 1.6) +
  geom_errorbar(aes(ymin = cv_error - cv_se, ymax = cv_error + cv_se),
                width = 1.2, color = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = n_opt, linetype = "dashed", color = "#C0392B") +
  annotate("point", x = n_opt, y = min_err, color = "#C0392B", size = 3) +
  ## 标注放在左下空白区，用箭头指向最优点，避免压住曲线
  annotate("text", x = 125, y = 0.0115,
           label = sprintf("optimal: n = %d features\nCV error = %.4f", n_opt, min_err),
           hjust = 0, vjust = 0.4, size = 3.1, color = "#C0392B") +
  annotate("segment", x = 118, y = 0.0135, xend = 55, yend = 0.0102,
           arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed"),
           color = "#C0392B", linewidth = 0.4) +
  scale_x_reverse() +
  labs(x = "Number of features", y = sprintf("%d-fold CV error rate", N_FOLD),
       title = "SVM-RFE recursive feature elimination") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(d_svm, "GSE57338_SVM_RFE_cv_error_curve.pdf"), p_rfe,
       width = 6.4, height = 5.0, device = "pdf")
ggsave(file.path(d_svm, "GSE57338_SVM_RFE_cv_error_curve.png"), p_rfe,
       width = 6.4, height = 5.0, dpi = 300)

## 5.2 完整排序表（rank 1 = 最重要）
all_feats <- colnames(X_z)
rank_tab <- data.frame(GeneSymbol = all_feats, Rank = NA_integer_,
                       stringsAsFactors = FALSE)
rownames(rank_tab) <- all_feats
fin_ord <- order(-abs(rfe$final_w))
rank_tab[rfe$final_feats[fin_ord], "Rank"] <- seq_along(rfe$final_feats)
next_rank <- length(rfe$final_feats) + 1L
if (length(rfe$elim) > 0) {
  for (s in rev(seq_along(rfe$elim))) {
    g <- rfe$elim[[s]]
    ws <- abs(rfe$w[[s]][g])
    g_ord <- g[order(-ws)]
    rank_tab[g_ord, "Rank"] <- next_rank:(next_rank + length(g) - 1L)
    next_rank <- next_rank + length(g)
  }
}
rank_tab <- rank_tab[order(rank_tab$Rank), , drop = FALSE]
rank_tab$Selected <- rank_tab$GeneSymbol %in% svm_genes
rank_tab$Elimination_round <- NA_integer_
for (s in seq_along(rfe$elim)) {
  rank_tab$Elimination_round[match(rfe$elim[[s]], rank_tab$GeneSymbol)] <- s
}
rank_tab$CV_error_at_optimal_n <- min_err   # 全局最优 n 下的 5 折平均错误率（非逐基因值）

wcsv(rank_tab, file.path(d_svm, "GSE57338_SVM_RFE_rank_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(rank_tab[rank_tab$Selected, , drop = FALSE],
          file.path(d_svm, "GSE57338_SVM_RFE_selected_genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

if (any(is.na(rank_tab$Rank))) {
  add_issue("SVM-RFE", "ERROR", "排序表存在未赋 rank 的基因",
            sprintf("%d 个基因 rank 缺失", sum(is.na(rank_tab$Rank))))
}

## 5.3 保存模型（最终最优特征集上的线性 SVM）
reset_seed()
svm_final <- e1071::svm(X_z[, svm_genes, drop = FALSE], y_fac,
                        kernel = "linear", cost = SVM_COST, scale = TRUE)
saveRDS(list(rfe = rfe, final_model = svm_final, n_opt = n_opt,
             min_err = min_err, genes = svm_genes, seed = SEED, foldid = fold_id),
        file.path(d_svm, "GSE57338_SVM_RFE_model.rds"))

## 5.4 权重条形图（最优特征集，top 20）
## 【注意】权重必须取自"最优特征集上重新训练的模型"(svm_final)；
##         rfe$final_w 是迭代终止时只剩 2 个特征的那一轮权重，直接用来画图会得到错误的 2 根柱子。
w_opt <- drop(t(svm_final$coefs) %*% svm_final$SV)
if (length(w_opt) != length(svm_genes)) {
  add_issue("SVM-RFE", "ERROR", "权重维度与最优特征集不一致",
            sprintf("模型权重 %d 个，最优特征 %d 个", length(w_opt), length(svm_genes)))
} else {
  names(w_opt) <- svm_genes
  w_df <- data.frame(GeneSymbol = names(w_opt), Weight = as.numeric(w_opt),
                     stringsAsFactors = FALSE)
  w_df <- w_df[order(-abs(w_df$Weight)), , drop = FALSE]
  wcsv(w_df, file.path(d_svm, "GSE57338_SVM_RFE_weights_table.csv"),
       row.names = FALSE, fileEncoding = "UTF-8")
  w_df_top <- head(w_df, 20)
  p_w <- ggplot(w_df_top, aes(x = reorder(GeneSymbol, Weight), y = Weight,
                              fill = Weight > 0)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C6E9B"),
                      labels = c(`TRUE` = "Up in Case", `FALSE` = "Down in Case"), name = NULL) +
    labs(x = NULL, y = "Linear SVM weight (|w| ranking)",
         title = sprintf("SVM-RFE selected genes (n = %d, top 20 by |w|)", length(svm_genes))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  ggsave(file.path(d_svm, "GSE57338_SVM_RFE_weight_barplot.pdf"), p_w,
         width = 6.2, height = 6.4, device = "pdf")
  ggsave(file.path(d_svm, "GSE57338_SVM_RFE_weight_barplot.png"), p_w,
         width = 6.2, height = 6.4, dpi = 300)
}

###############################################################################
## 6. 随机森林
###############################################################################
msg("\n---- 6. Random Forest ----")

## 6.1 5 折交叉验证（与 LASSO / SVM-RFE 共用同一套 foldid）
reset_seed()
cv_err_rf <- numeric(N_FOLD); oob_err_rf <- numeric(N_FOLD)
for (i in seq_len(N_FOLD)) {
  te <- which(fold_id == i); tr <- setdiff(seq_len(nrow(X_raw)), te)
  reset_seed()
  rf_i <- randomForest(x = X_raw[tr, , drop = FALSE], y = y_fac[tr],
                       ntree = RF_NTREE, importance = TRUE)
  pred_i <- predict(rf_i, X_raw[te, , drop = FALSE])
  cv_err_rf[i] <- mean(pred_i != y_fac[te])
  oob_err_rf[i] <- rf_i$err.rate[RF_NTREE, "OOB"]
}
msg(sprintf("5 折 CV 错误率：均值 %.4f（各折 %s）",
            mean(cv_err_rf), paste(sprintf("%.3f", cv_err_rf), collapse = ", ")))
msg(sprintf("各折 OOB 错误率：均值 %.4f", mean(oob_err_rf)))

## 6.2 全数据训练最终模型（用于重要性）
reset_seed()
rf_fit <- randomForest(x = X_raw, y = y_fac, ntree = RF_NTREE,
                       importance = TRUE, proximity = FALSE)
msg(sprintf("最终模型 OOB 错误率：%.4f（Control %.4f / Case %.4f）",
            rf_fit$err.rate[RF_NTREE, "OOB"],
            rf_fit$err.rate[RF_NTREE, CTRL_LABEL],
            rf_fit$err.rate[RF_NTREE, CASE_LABEL]))

imp <- importance(rf_fit, type = 1)   # type=1: MeanDecreaseAccuracy
imp_gini <- importance(rf_fit, type = 2)
imp_tab <- data.frame(
  GeneSymbol = rownames(imp),
  MeanDecreaseAccuracy = as.numeric(imp[, 1]),
  MeanDecreaseGini = as.numeric(imp_gini[rownames(imp), 1]),
  stringsAsFactors = FALSE
)
## 也导出按类别拆分的重要性（若可用）
if (ncol(imp) >= 3) {
  imp_tab$MDA_Control <- as.numeric(imp[, CTRL_LABEL])
  imp_tab$MDA_Case <- as.numeric(imp[, CASE_LABEL])
}
imp_tab <- imp_tab[order(-imp_tab$MeanDecreaseAccuracy), , drop = FALSE]
imp_tab$Rank_MDA <- seq_len(nrow(imp_tab))
imp_tab$Rank_Gini <- rank(-imp_tab$MeanDecreaseGini, ties.method = "first")

## 6.3 筛选规则（预先设定，不事后调参）
##   主规则：MDA > mean(MDA) + 1 SD（"显著高于平均重要性"）
##   补充规则：Top-N（N = RF_TOPN）
thr_mda <- mean(imp_tab$MeanDecreaseAccuracy) + sd(imp_tab$MeanDecreaseAccuracy)
rf_genes_main <- imp_tab$GeneSymbol[imp_tab$MeanDecreaseAccuracy > thr_mda]
rf_genes_topn <- head(imp_tab$GeneSymbol, min(RF_TOPN, nrow(imp_tab)))
msg(sprintf("MDA 阈值(mean+1SD) = %.4f → 入选 %d 个基因", thr_mda, length(rf_genes_main)))
if (length(rf_genes_main) == 0) {
  add_issue("随机森林", "WARN", "mean+1SD 规则筛出 0 个基因",
            sprintf("阈值 %.4f 高于全部基因的 MDA（最大 %.4f）；主结果改用 Top-%d 规则并明确标注",
                    thr_mda, max(imp_tab$MeanDecreaseAccuracy), RF_TOPN))
  rf_genes <- rf_genes_topn
  rf_rule_used <- sprintf("Top-%d（fallback：mean+1SD 规则为空）", RF_TOPN)
} else {
  rf_genes <- rf_genes_main
  rf_rule_used <- "MDA > mean(MDA) + 1 SD"
}
imp_tab$Selected_main <- imp_tab$GeneSymbol %in% rf_genes
imp_tab$Selected_topN <- imp_tab$GeneSymbol %in% rf_genes_topn
msg(sprintf("随机森林最终入选基因：%d 个（规则：%s）", length(rf_genes), rf_rule_used))

wcsv(imp_tab, file.path(d_rf, "GSE57338_RF_importance_table.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(imp_tab[imp_tab$Selected_main, , drop = FALSE],
          file.path(d_rf, "GSE57338_RF_selected_genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(imp_tab[imp_tab$Selected_topN, , drop = FALSE],
          file.path(d_rf, "GSE57338_RF_selected_genes_topN.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
wcsv(data.frame(fold = seq_len(N_FOLD), cv_error = cv_err_rf, oob_error = oob_err_rf),
          file.path(d_rf, "GSE57338_RF_cv_error.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  sprintf("ntree = %d", RF_NTREE),
  sprintf("mtry  = %d (default, floor(sqrt(p)) with p = %d)", rf_fit$mtry, ncol(X_raw)),
  sprintf("nodesize = %d", rf_fit$nodesize),
  sprintf("OOB error = %.4f", rf_fit$err.rate[RF_NTREE, "OOB"]),
  sprintf("5-fold CV error = %.4f (sd %.4f)", mean(cv_err_rf), sd(cv_err_rf)),
  sprintf("selection rule = %s", rf_rule_used),
  sprintf("MDA threshold (mean+1SD) = %.4f", thr_mda),
  sprintf("seed = %d", SEED)
), file.path(d_rf, "GSE57338_RF_model_summary.txt"))

## 6.4 重要性排序图
plot_n <- min(30L, nrow(imp_tab))
imp_plot <- head(imp_tab, plot_n)
imp_plot$GeneSymbol <- factor(imp_plot$GeneSymbol, levels = rev(imp_plot$GeneSymbol))
p_imp <- ggplot(imp_plot, aes(x = GeneSymbol, y = MeanDecreaseAccuracy,
                              fill = Selected_main)) +
  geom_col(width = 0.72) + coord_flip() +
  geom_hline(yintercept = thr_mda, linetype = "dashed", color = "#C0392B", linewidth = 0.5) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey65"),
                    labels = c(`TRUE` = "Selected", `FALSE` = "Not selected"), name = NULL) +
  labs(x = NULL, y = "Mean Decrease Accuracy",
       title = sprintf("Random Forest variable importance (top %d)", plot_n),
       subtitle = sprintf("dashed line = mean + 1 SD (%.3f); selected n = %d",
                          thr_mda, length(rf_genes))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.subtitle = element_text(size = 9, color = "grey35"))
ggsave(file.path(d_rf, "GSE57338_RF_importance_plot.pdf"), p_imp,
       width = 6.4, height = 7.2, device = "pdf")
ggsave(file.path(d_rf, "GSE57338_RF_importance_plot.png"), p_imp,
       width = 6.4, height = 7.2, dpi = 300)

## 6.5 重要性分布图（规则可解释性）
p_dist <- ggplot(imp_tab, aes(x = MeanDecreaseAccuracy)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = thr_mda, color = "#C0392B", linetype = "dashed") +
  labs(x = "Mean Decrease Accuracy", y = "Gene count",
       title = "Distribution of RF importance",
       subtitle = sprintf("threshold mean+1SD = %.3f", thr_mda)) +
  theme_bw(base_size = 11)
ggsave(file.path(d_rf, "GSE57338_RF_importance_distribution.pdf"), p_dist,
       width = 6.0, height = 4.2, device = "pdf")
ggsave(file.path(d_rf, "GSE57338_RF_importance_distribution.png"), p_dist,
       width = 6.0, height = 4.2, dpi = 300)

saveRDS(list(model = rf_fit, importance = imp_tab, genes = rf_genes,
             genes_topN = rf_genes_topn, threshold = thr_mda,
             cv_error = cv_err_rf, oob_error = oob_err_rf,
             rule = rf_rule_used, seed = SEED),
        file.path(d_rf, "GSE57338_RF_model.rds"))

###############################################################################
## 7. 三方法比较：Venn / UpSet / 交集
###############################################################################
msg("\n---- 7. 三方法比较 ----")

gene_sets <- list(LASSO = sort(unique(lasso_genes)),
                  `SVM-RFE` = sort(unique(svm_genes)),
                  `RandomForest` = sort(unique(rf_genes)))
n_sets <- sapply(gene_sets, length)
msg(sprintf("各方法入选基因数：LASSO %d | SVM-RFE %d | RandomForest %d",
            n_sets[1], n_sets[2], n_sets[3]))

## 7.1 并集长表（基因 × 方法 0/1）
all_genes <- sort(unique(unlist(gene_sets)))
membership <- data.frame(
  GeneSymbol = all_genes,
  LASSO = as.integer(all_genes %in% gene_sets$LASSO),
  SVM_RFE = as.integer(all_genes %in% gene_sets$`SVM-RFE`),
  RandomForest = as.integer(all_genes %in% gene_sets$RandomForest),
  stringsAsFactors = FALSE
)
membership$n_methods <- rowSums(membership[, c("LASSO", "SVM_RFE", "RandomForest")])
membership <- membership[order(-membership$n_methods, membership$GeneSymbol), , drop = FALSE]
wcsv(membership, file.path(d_cmp, "GSE57338_ML_gene_membership.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 7.2 交集（三方法共同）
common_all <- Reduce(intersect, gene_sets)
msg(sprintf("三方法共同交集基因：%d 个", length(common_all)))
if (length(common_all) == 0) {
  add_issue("比较", "ERROR", "三方法交集为空",
            "LASSO / SVM-RFE / RandomForest 无共同基因；已在报告中明确说明，未放宽任何阈值")
}

common_tab <- membership[membership$n_methods == 3, , drop = FALSE]
if (nrow(common_tab) > 0) {
  ## 附上各方法内部的排名/系数，便于解读
  common_tab$LASSO_coef <- coef_tab$Coef_lambda_min[match(common_tab$GeneSymbol, coef_tab$GeneSymbol)]
  common_tab$SVM_RFE_rank <- rank_tab$Rank[match(common_tab$GeneSymbol, rank_tab$GeneSymbol)]
  common_tab$RF_Rank_MDA <- imp_tab$Rank_MDA[match(common_tab$GeneSymbol, imp_tab$GeneSymbol)]
  common_tab$RF_MDA <- imp_tab$MeanDecreaseAccuracy[match(common_tab$GeneSymbol, imp_tab$GeneSymbol)]
}
wcsv(common_tab, file.path(d_cmp, "GSE57338_ML_intersection_genes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 7.3 两两交集（交集为空时仍有信息量）
pair_list <- list(
  LASSO_vs_SVM_RFE = intersect(gene_sets$LASSO, gene_sets$`SVM-RFE`),
  LASSO_vs_RF      = intersect(gene_sets$LASSO, gene_sets$`RandomForest`),
  SVM_RFE_vs_RF    = intersect(gene_sets$`SVM-RFE`, gene_sets$`RandomForest`)
)
pair_tab <- data.frame(
  comparison = names(pair_list),
  n_genes = sapply(pair_list, length),
  genes = sapply(pair_list, function(x) paste(x, collapse = ";")),
  stringsAsFactors = FALSE
)
wcsv(pair_tab, file.path(d_cmp, "GSE57338_ML_pairwise_intersection.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 方法一致性总览
overview <- data.frame(
  method = c("LASSO", "SVM-RFE", "RandomForest"),
  n_selected = as.integer(n_sets),
  rule = c(sprintf("lambda.min (%.4f), non-zero coef", lam_min),
           sprintf("min CV error %.4f at n = %d", min_err, n_opt),
           rf_rule_used),
  stringsAsFactors = FALSE
)
overview$overlap_with_others <- c(
  length(intersect(gene_sets$LASSO, union(gene_sets$`SVM-RFE`, gene_sets$`RandomForest`))),
  length(intersect(gene_sets$`SVM-RFE`, union(gene_sets$LASSO, gene_sets$`RandomForest`))),
  length(intersect(gene_sets$`RandomForest`, union(gene_sets$LASSO, gene_sets$`SVM-RFE`)))
)
overview$in_all_three <- length(common_all)
wcsv(overview, file.path(d_cmp, "GSE57338_ML_method_overview.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 7.4 Venn 图（PDF + PNG）
if (has_venn) {
  ## VennDiagram 通过 futile.logger 输出大量 INFO 日志，压到 ERROR 级别保持输出干净
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(futile.logger::ERROR)
  }
  venn_cols <- c("#C0392B", "#2C6E9B", "#6B8E5A")
  venn_labels <- c(sprintf("LASSO\n(%d)", n_sets[1]),
                   sprintf("SVM-RFE\n(%d)", n_sets[2]),
                   sprintf("Random Forest\n(%d)", n_sets[3]))
  make_venn <- function(fname) {
    VennDiagram::venn.diagram(
      x = gene_sets, filename = fname,
      category.names = venn_labels,
      fill = venn_cols, alpha = 0.45, col = "grey30", lwd = 1,
      cex = 1.1, fontfamily = "sans", cat.cex = 0.95, cat.fontfamily = "sans",
      cat.pos = c(-20, 20, 180), cat.dist = c(0.06, 0.06, 0.05),
      margin = 0.10, main = "Overlap of feature genes across three methods",
      main.cex = 1.0, main.fontfamily = "sans"
    )
    invisible(NULL)
  }
  ok_venn <- tryCatch({
    make_venn(file.path(d_cmp, "GSE57338_ML_venn.pdf"))
    make_venn(file.path(d_cmp, "GSE57338_ML_venn.png"))
    TRUE
  }, error = function(e) {
    add_issue("比较", "WARN", "Venn 图生成失败", conditionMessage(e))
    FALSE
  })
  if (ok_venn) msg("Venn 图已保存（PDF + PNG）")
}

## 7.5 UpSet 图（PDF + PNG）
if (has_upset) {
  draw_upset <- function() {
    print(UpSetR::upset(UpSetR::fromList(gene_sets), nsets = 3,
                        order.by = "freq", keep.order = TRUE,
                        sets.bar.color = "#2C6E9B",
                        main.bar.color = "#4A4A4A",
                        matrix.color = "#C0392B",
                        mainbar.y.label = "Intersection size",
                        sets.x.label = "Feature genes per method",
                        text.scale = c(1.2, 1.1, 1.0, 1.0, 1.1, 1.0),
                        mb.ratio = c(0.62, 0.38)))
  }
  ok_upset <- tryCatch({
    pdf(file.path(d_cmp, "GSE57338_ML_upset.pdf"), width = 8, height = 5)
    draw_upset(); dev.off()
    png(file.path(d_cmp, "GSE57338_ML_upset.png"), width = 8 * 300, height = 5 * 300, res = 300)
    draw_upset(); dev.off()
    TRUE
  }, error = function(e) {
    while (dev.cur() > 1) dev.off()
    add_issue("比较", "WARN", "UpSet 图生成失败（UpSetR 与当前 ggplot2 版本兼容性问题）",
              conditionMessage(e))
    FALSE
  })
  if (ok_upset) msg("UpSet 图已保存（PDF + PNG）")
}

###############################################################################
## 8. 结果汇总与 sessionInfo
###############################################################################
msg("\n---- 8. 结果汇总 ----")

summary_tab <- data.frame(
  item = c("候选基因数（输入）", "候选基因在矩阵中匹配数", "候选基因缺失数",
           "样本数（Case）", "样本数（Control）", "特征维度",
           "LASSO 入选基因数（lambda.min）", "LASSO 入选基因数（lambda.1se）",
           "SVM-RFE 最优特征数", "SVM-RFE 最低 CV 错误率",
           "随机森林 OOB 错误率", "随机森林 5 折 CV 错误率",
           "随机森林入选基因数（主规则）", "随机森林 Top-N 基因数",
           "三方法共同交集基因数"),
  value = c(length(cand_genes), length(cand_in_expr), length(cand_missing),
            sum(y_num == CASE_CODE), sum(y_num == CTRL_CODE),
            sprintf("%d samples x %d genes", nrow(X_raw), ncol(X_raw)),
            length(lasso_genes), length(lasso_genes_1se),
            n_opt, sprintf("%.4f", min_err),
            sprintf("%.4f", rf_fit$err.rate[RF_NTREE, "OOB"]),
            sprintf("%.4f", mean(cv_err_rf)),
            length(rf_genes_main), length(rf_genes_topn), length(common_all)),
  stringsAsFactors = FALSE
)
wcsv(summary_tab, file.path(out_root, "GSE57338_ML_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
for (i in seq_len(nrow(summary_tab))) {
  msg(sprintf("  %-32s : %s", summary_tab$item[i], summary_tab$value[i]))
}

## 8.1 各方法基因清单汇总（长表，便于后续取用）
list_long <- rbind(
  data.frame(method = "LASSO", GeneSymbol = gene_sets$LASSO, stringsAsFactors = FALSE),
  data.frame(method = "SVM-RFE", GeneSymbol = gene_sets$`SVM-RFE`, stringsAsFactors = FALSE),
  data.frame(method = "RandomForest", GeneSymbol = gene_sets$`RandomForest`, stringsAsFactors = FALSE)
)
wcsv(list_long, file.path(out_root, "GSE57338_ML_selected_genes_long.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 8.2 异常登记表落盘
if (nrow(issues) > 0) {
  wcsv(issues, file.path(out_root, "GSE57338_ML_issues.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  msg(sprintf("异常登记：共 %d 条，已写入 GSE57338_ML_issues.csv", nrow(issues)))
} else {
  wcsv(data.frame(step = character(), level = character(), item = character(),
                       detail = character()),
            file.path(out_root, "GSE57338_ML_issues.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  msg("异常登记：无异常")
}

## 8.3 sessionInfo
si_file <- file.path(d_rep, "sessionInfo_09_ml.txt")
writeLines(c(
  "================================================================================",
  "GSE57338 特征基因筛选 (LASSO / SVM-RFE / Random Forest) —— 运行环境",
  sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("R 版本  ：%s", R.version.string),
  sprintf("平台    ：%s", R.version$platform),
  sprintf("随机种子：%d", SEED),
  "================================================================================",
  "",
  capture.output(sessionInfo())
), si_file, useBytes = TRUE)
saveRDS(sessionInfo(), file.path(out_root, "GSE57338_ML_sessionInfo.rds"))

msg("\n全部完成：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
close(log_con)
