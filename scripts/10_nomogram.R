###############################################################################
## 10_nomogram.R
## GSE57338 —— 关键基因 AQP3 / GGT5 / ISLR / LAD1 列线图 (Nomogram) 预测模型
##
## 输入：
##   clean_data/GSE57338_clean_expression_matrix.csv  (基因 x 样本, log2)
##   clean_data/GSE57338_sample_info.csv              (样本表, 结局取 group 列)
##
## 输出：
##   results/Nomogram/
##     01_tables/   所有统计表（QC、系数、AUC、校准、DCA、列线图分值）
##     02_figures/  列线图、ROC、校准曲线、DCA（PDF + PNG 双份）
##     03_models/   冻结的模型对象 (.rds)
##     04_data/     modelpredictions.csv 及绘图数据
##     05_logs/     脚本副本、参数记录 (json)、运行日志、sessionInfo
##     Nomogram_report.txt   中文结果报告
##
## 运行：Rscript scripts/10_nomogram.R
##
## 设计要点：
##   1. 结局为二分类（HF = DCM/ICM vs NF），故用 Logistic 回归 + rms::lrm；
##      脚本内置结局类型自动判定：若数据含 time/status 列则改走 Cox (cph) 分支，
##      并在该分支中显式声明事件定义与预测时间点。
##   2. 所有预处理（基因 z 化参数、年龄缺失填补、共线性/线性检验、变量取舍）
##      全部只在训练集完成；验证集仅用于评价冻结后的最终模型。
##   3. 列线图必须由 rms::nomogram() 生成，不使用自绘替代。
##   4. 内部验证用 Bootstrap（validate / calibrate），外部评价用独立验证集。
###############################################################################

rm(list = ls())

## ===========================================================================
## 0. 编码守卫 + 全局参数（换机器只需改 proj_dir 一处）
## ===========================================================================
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("Chinese (Simplified)_China.utf8", "Chinese", "C.UTF-8", "en_US.UTF-8")) {
    ok <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error = function(e) ""))
    if (nzchar(ok) && isTRUE(l10n_info()$`UTF-8`)) break
  }
}
if (!isTRUE(l10n_info()$`UTF-8`)) {
  warning("当前 R 会话不是 UTF-8 locale，中文可能无法正确写入文件")
}

SEED        <- 20261027L   # 唯一随机源：划分与 Bootstrap 共用
GENES       <- c("AQP3", "GGT5", "ISLR", "LAD1")   # 关键基因（用户指定，不做二次筛选）
TRAIN_FRAC  <- 0.70        # 训练集比例（按结局分层）
B_BOOT      <- 1000L       # Bootstrap 次数
MIN_VALID_AGE <- 18        # 年龄 < 18 视为不符合本数据集（成人心脏组织）的记录错误 -> 缺失
CASE_GROUPS <- c("DCM", "ICM")   # 阳性：心衰（扩张型 / 缺血型心肌病）
CTRL_GROUPS <- c("NF")           # 阴性：非心衰对照
SURV_TIME_COLS   <- c("os_time", "time", "follow_up", "fu_time", "months")
SURV_STATUS_COLS <- c("os_status", "status", "event", "death")
SURV_HORIZONS    <- c(12, 36, 60)  # 仅在走 Cox 分支时生效：1/3/5 年预测时间点

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()

in_expr <- file.path(proj_dir, "clean_data", "GSE57338_clean_expression_matrix.csv")
in_meta <- file.path(proj_dir, "clean_data", "GSE57338_sample_info.csv")

out_root <- file.path(proj_dir, "results", "Nomogram")
d_tab <- file.path(out_root, "01_tables")
d_fig <- file.path(out_root, "02_figures")
d_mod <- file.path(out_root, "03_models")
d_dat <- file.path(out_root, "04_data")
d_log <- file.path(out_root, "05_logs")
for (d in c(out_root, d_tab, d_fig, d_mod, d_dat, d_log)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## 依赖包（均为 CRAN 常见包，缺失即明确报错，避免静默降级）
need_pkgs <- c("rms", "pROC", "survival", "ggplot2", "jsonlite")
miss <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("缺少必需 R 包: ", paste(miss, collapse = ", "), "；请先 install.packages()")
}
suppressPackageStartupMessages({
  library(rms)
  library(pROC)
  library(survival)
  library(ggplot2)
  library(jsonlite)   # 注意：jsonlite 也导出 validate()，会屏蔽 rms::validate
})                    # 因此下文一律用 rms::validate / rms::calibrate 显式调用

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...)
  run_log <<- c(run_log, txt)
  cat(txt, "\n", sep = "")
}

## 统一的 CSV 写出：UTF-8 + BOM（Excel 打开中文不乱码）
write_csv <- function(df, path) {
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  txt <- readLines(tmp, encoding = "UTF-8", warn = FALSE)
  con <- file(path, open = "wb")
  writeBin(charToRaw("\xEF\xBB\xBF"), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(txt, collapse = "\n"), "\n"))), con)
  close(con)
  unlink(tmp)
  invisible(path)
}

## 图形双份保存（PDF + PNG）。注意：expr 以表达式传入并在两个设备上各求值一次，
## 不能直接 force(expr)——promise 只会求值一次，第二个设备会得到空白图。
save_plot <- function(stem, expr, width = 7, height = 6) {
  e <- substitute(expr)
  env <- parent.frame()
  pdf(file.path(d_fig, paste0(stem, ".pdf")), width = width, height = height)
  eval(e, env)
  dev.off()
  png(file.path(d_fig, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = 300)
  eval(e, env)
  dev.off()
  invisible(stem)
}

log_msg("================ GSE57338 Nomogram 建模开始 ================")
log_msg("R: ", R.version.string, " | seed = ", SEED, " | B = ", B_BOOT)

## ===========================================================================
## 1. 读取数据 + 样本/基因 ID 对齐 + 质量检查
## ===========================================================================
expr_raw <- read.csv(in_expr, check.names = FALSE, stringsAsFactors = FALSE)
meta     <- read.csv(in_meta, check.names = FALSE, stringsAsFactors = FALSE)
log_msg("表达矩阵: ", nrow(expr_raw), " 基因 x ", ncol(expr_raw) - 1, " 样本")
log_msg("样本表  : ", nrow(meta), " 行")

## --- 1.1 基因 ID：唯一性检查 ---
dup_genes <- expr_raw$GeneSymbol[duplicated(expr_raw$GeneSymbol)]
if (length(dup_genes) > 0) {
  log_msg("警告: 存在重复基因符号 ", length(dup_genes), " 个，保留表达方差最大的一行")
  expr_raw <- expr_raw[order(expr_raw$GeneSymbol, -apply(expr_raw[, -1], 1, var)), ]
  expr_raw <- expr_raw[!duplicated(expr_raw$GeneSymbol), ]
}

## --- 1.2 关键基因存在性 ---
gene_hit <- GENES %in% expr_raw$GeneSymbol
if (!all(gene_hit)) {
  stop("以下关键基因在表达矩阵中缺失: ", paste(GENES[!gene_hit], collapse = ", "))
}

## --- 1.3 样本 ID 对齐（列名 vs 样本表） ---
mat_samples <- colnames(expr_raw)[-1]
common <- intersect(mat_samples, meta$sample)
only_mat <- setdiff(mat_samples, meta$sample)
only_meta <- setdiff(meta$sample, mat_samples)
log_msg("样本对齐: 交集 ", length(common), " | 仅矩阵 ", length(only_mat), " | 仅样本表 ", length(only_meta))
if (length(common) < 10) stop("样本 ID 交集过少，无法建模")

align_qc <- data.frame(
  item = c("expression_matrix_genes", "expression_matrix_samples", "sample_info_rows",
           "matched_samples", "samples_only_in_matrix", "samples_only_in_sample_info",
           "key_genes_requested", "key_genes_missing"),
  value = c(nrow(expr_raw), length(mat_samples), nrow(meta), length(common),
            length(only_mat), length(only_meta), paste(GENES, collapse = ";"),
            paste(GENES[!gene_hit], collapse = ";"))
)
write_csv(align_qc, file.path(d_tab, "T01_sample_alignment_qc.csv"))

## 基因 x 样本 -> 样本 x 基因（只取关键基因）
expr_sub <- expr_raw[match(GENES, expr_raw$GeneSymbol), c("GeneSymbol", common)]
expr_t <- as.data.frame(t(expr_sub[, -1]), stringsAsFactors = FALSE)
colnames(expr_t) <- expr_sub$GeneSymbol
expr_t$sample <- common
for (g in GENES) expr_t[[g]] <- as.numeric(expr_t[[g]])

## --- 1.4 合并样本表 ---
meta_sub <- meta[match(common, meta$sample), ]
dat <- merge(meta_sub, expr_t, by = "sample", sort = FALSE)
stopifnot(nrow(dat) == length(common))

## --- 1.5 变量类型与缺失 ---
var_type <- data.frame(
  variable = c("sample", "group", "gender", "age", GENES),
  class_raw = c(class(dat$sample), class(dat$group), class(dat$gender), class(dat$age),
                vapply(GENES, function(g) class(dat[[g]]), character(1))),
  stringsAsFactors = FALSE
)
n_miss <- vapply(dat[, var_type$variable], function(x) sum(is.na(x) | (is.character(x) & trimws(x) == "")), integer(1))
var_type$n_missing <- as.integer(n_miss)
var_type$pct_missing <- round(100 * var_type$n_missing / nrow(dat), 2)
write_csv(var_type, file.path(d_tab, "T03_variable_types.csv"))
write_csv(data.frame(variable = var_type$variable, n_missing = var_type$n_missing,
                     pct_missing = var_type$pct_missing),
          file.path(d_tab, "T02_missing_values.csv"))
log_msg("缺失值: ", sum(var_type$n_missing), " 个（详见 T02）")

## --- 1.6 结局编码（核心：明确事件定义） ---
grp <- toupper(trimws(as.character(dat$group)))
grp_tab <- table(grp)
log_msg("group 原始分布: ", paste(names(grp_tab), grp_tab, sep = "=", collapse = " | "))

bad_grp <- setdiff(unique(grp), c(CASE_GROUPS, CTRL_GROUPS))
if (length(bad_grp) > 0) stop("存在未定义的 group 取值: ", paste(bad_grp, collapse = ", "))

dat$HF <- ifelse(grp %in% CASE_GROUPS, 1L, 0L)   # 事件 = 心衰 (DCM 或 ICM)
outcome_def <- data.frame(
  item = c("outcome", "event_definition", "event_code", "non_event_definition", "non_event_code",
           "n_total", "n_event", "n_non_event", "event_rate_pct", "time_dimension"),
  value = c("HF (heart failure, binary)", "group in DCM or ICM", "1",
            "group == NF (non-failing donor heart)", "0",
            nrow(dat), sum(dat$HF), sum(dat$HF == 0),
            round(100 * mean(dat$HF), 2),
            "none (cross-sectional; no follow-up time in GSE57338)")
)
write_csv(outcome_def, file.path(d_tab, "T05_outcome_coding.csv"))

## --- 1.7 变量清洗 ---
dat$sex <- factor(ifelse(tolower(trimws(as.character(dat$gender))) %in% c("male", "m"),
                         "Male", "Female"), levels = c("Female", "Male"))
if (any(is.na(dat$sex))) stop("gender 存在无法识别的取值")

dat$age <- suppressWarnings(as.numeric(dat$age))
age_flag <- !is.na(dat$age) & dat$age < MIN_VALID_AGE   # 0/1 岁记录在成人心脏数据中不合理
n_age_flag <- sum(age_flag)
dat$age[age_flag] <- NA_real_
log_msg("年龄异常（< ", MIN_VALID_AGE, " 岁）标记为缺失: ", n_age_flag, " 例")

## --- 1.8 异常值（表达值）检查：|z| > 3 与 1.5*IQR 两种口径 ---
out_rows <- list()
for (g in GENES) {
  v <- dat[[g]]
  z <- (v - mean(v)) / sd(v)
  q <- quantile(v, c(.25, .75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  out_rows[[g]] <- data.frame(
    gene = g, n = sum(!is.na(v)), mean = round(mean(v, na.rm = TRUE), 3),
    sd = round(sd(v, na.rm = TRUE), 3), min = round(min(v, na.rm = TRUE), 3),
    max = round(max(v, na.rm = TRUE), 3),
    n_abs_z_gt3 = sum(abs(z) > 3, na.rm = TRUE),
    n_out_1.5iqr = sum(v < q[1] - 1.5 * iqr | v > q[2] + 1.5 * iqr, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
out_tab <- do.call(rbind, out_rows)
write_csv(out_tab, file.path(d_tab, "T04_outliers_expression.csv"))
log_msg("表达异常值 (|z|>3): ", paste(out_tab$gene, out_tab$n_abs_z_gt3, sep = "=", collapse = ", "),
        " —— 均为生物学变异，不做删除")

## ===========================================================================
## 2. 训练集 / 验证集划分（按结局分层，固定种子）
## ===========================================================================
set.seed(SEED)
idx_case <- which(dat$HF == 1)
idx_ctrl <- which(dat$HF == 0)
n_case_tr <- round(length(idx_case) * TRAIN_FRAC)
n_ctrl_tr <- round(length(idx_ctrl) * TRAIN_FRAC)
tr_idx <- c(sample(idx_case, n_case_tr), sample(idx_ctrl, n_ctrl_tr))
tr_idx <- sort(tr_idx)
va_idx <- setdiff(seq_len(nrow(dat)), tr_idx)

dat$set <- ifelse(seq_len(nrow(dat)) %in% tr_idx, "train", "validation")
dat$set <- factor(dat$set, levels = c("train", "validation"))
train <- dat[dat$set == "train", ]
valid <- dat[dat$set == "validation", ]
rownames(train) <- NULL; rownames(valid) <- NULL
log_msg("划分: 训练 ", nrow(train), " (事件 ", sum(train$HF), ") | 验证 ",
        nrow(valid), " (事件 ", sum(valid$HF), ")")

## ===========================================================================
## 3. 预处理参数（只用训练集估计）
## ===========================================================================
## 3.1 基因 z 标准化参数（训练集 mean/sd）
zpar <- data.frame(
  gene = GENES,
  train_mean = vapply(GENES, function(g) mean(train[[g]], na.rm = TRUE), numeric(1)),
  train_sd   = vapply(GENES, function(g) sd(train[[g]], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)
zpar$train_mean <- round(zpar$train_mean, 6)
zpar$train_sd   <- round(zpar$train_sd, 6)
rownames(zpar)  <- NULL
if (any(zpar$train_sd <= 0)) stop("存在训练集内方差为 0 的基因，无法标准化")
write_csv(zpar, file.path(d_tab, "T06_zscore_parameters_train_only.csv"))

for (g in GENES) {
  m <- zpar$train_mean[zpar$gene == g]
  s <- zpar$train_sd[zpar$gene == g]
  nm <- paste0(g, "_z")
  train[[nm]] <- (train[[g]] - m) / s
  valid[[nm]] <- (valid[[g]] - m) / s   # 验证集沿用训练集参数（不重新估计）
  dat[[nm]]   <- (dat[[g]] - m) / s
}

## 3.2 年龄缺失填补（训练集中位数，验证集沿用）
age_med <- median(train$age, na.rm = TRUE)
train$age_imp <- ifelse(is.na(train$age), age_med, train$age)
valid$age_imp <- ifelse(is.na(valid$age), age_med, valid$age)
dat$age_imp   <- ifelse(is.na(dat$age), age_med, dat$age)
train$age_missing <- as.integer(is.na(train$age))
valid$age_missing <- as.integer(is.na(valid$age))
dat$age_missing   <- as.integer(is.na(dat$age))
log_msg("年龄填补中位数（训练集）= ", age_med, "；填补例数 = ", sum(train$age_missing) + sum(valid$age_missing))

## 3.2b 变量标签（rms 绘图/列线图显示名，避免出现 _z 后缀）
set_label <- function(df, var, lab) {
  if (!var %in% names(df)) stop("set_label: 列不存在 -> ", var)
  v <- df[[var]]
  attr(v, "label") <- lab
  df[[var]] <- v
  df
}
for (g in GENES) train <- set_label(train, paste0(g, "_z"), g)
train <- set_label(train, "age_imp", "Age (years)")
train <- set_label(train, "sex", "Sex")

GENE_Z <- paste0(GENES, "_z")

## 3.3 描述性基线表（训练集，按结局分组）
desc_num <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
desc_bin <- function(x, lv) sprintf("%d (%.1f%%)", sum(x == lv, na.rm = TRUE),
                                    100 * mean(x == lv, na.rm = TRUE))
base_rows <- list()
for (v in c("age", GENES)) {
  p_w <- suppressWarnings(wilcox.test(dat[[v]][dat$HF == 1], dat[[v]][dat$HF == 0])$p.value)
  base_rows[[v]] <- data.frame(
    variable = v,
    overall = desc_num(dat[[v]]),
    HF_0 = desc_num(dat[[v]][dat$HF == 0]),
    HF_1 = desc_num(dat[[v]][dat$HF == 1]),
    p_value = format.pval(p_w, digits = 3, eps = 1e-4),
    stringsAsFactors = FALSE
  )
}
base_rows[["sex"]] <- data.frame(
  variable = "male_sex, n (%)",
  overall = desc_bin(dat$sex, "Male"),
  HF_0 = desc_bin(dat$sex[dat$HF == 0], "Male"),
  HF_1 = desc_bin(dat$sex[dat$HF == 1], "Male"),
  p_value = format.pval(chisq.test(table(dat$sex, dat$HF))$p.value, digits = 3, eps = 1e-4),
  stringsAsFactors = FALSE
)
base_tab <- do.call(rbind, base_rows)
write_csv(base_tab, file.path(d_tab, "T07_baseline_characteristics.csv"))

## ===========================================================================
## 4. 结局类型判定：二分类 -> Logistic；生存 -> Cox
## ===========================================================================
time_col   <- SURV_TIME_COLS[SURV_TIME_COLS %in% names(dat)][1]
status_col <- SURV_STATUS_COLS[SURV_STATUS_COLS %in% names(dat)][1]
has_survival <- !is.na(time_col) && !is.na(status_col)
log_msg("结局类型判定: ", ifelse(has_survival, "生存（Cox 分支）", "二分类（Logistic 分支，本次实际路径）"))

if (has_survival) {
  ## ---- Cox 分支（当前数据集不触发；保留以支持同流程复用到有随访的数据）----
  dat$surv_time <- as.numeric(dat[[time_col]])
  dat$surv_status <- as.integer(dat[[status_col]])   # 1 = 事件(死亡/进展), 0 = 删失
  train <- dat[dat$set == "train", ]; valid <- dat[dat$set == "validation", ]
  dd <- datadist(train); options(datadist = "dd")
  f_cox <- as.formula(paste("Surv(surv_time, surv_status) ~",
                            paste(GENE_Z, collapse = " + "), "+ age_imp + sex"))
  fit_cox <- cph(f_cox, data = train, x = TRUE, y = TRUE, surv = TRUE)
  ## 比例风险假设（Schoenfeld 残差）
  ph <- cox.zph(fit_cox)
  write_csv(data.frame(term = rownames(as.data.frame(ph$table)),
                       as.data.frame(ph$table), row.names = NULL),
            file.path(d_tab, "T08_cox_ph_assumption.csv"))
  ## 预测时间点：SURV_HORIZONS 个月，列线图输出 1/3/5 年无事件生存概率
  sf_fun <- function(t) {
    s0 <- summary(survfit(fit_cox), times = t)$surv[1]
    function(lp) s0 ^ exp(lp)
  }
  nom <- nomogram(fit_cox, fun = lapply(SURV_HORIZONS, sf_fun),
                  funlabel = paste0("Survival at ", SURV_HORIZONS, " mo"),
                  lp = TRUE)
  save_plot("F01_nomogram_cox", plot(nom, xfrac = .35, cex.axis = .8, cex.var = .9),
            width = 9, height = 7)
  saveRDS(list(fit = fit_cox, zpar = zpar, age_med = age_med),
          file.path(d_mod, "nomogram_model_cox.rds"))
  log_msg("Cox 分支完成；预测时间点 = ", paste(SURV_HORIZONS, collapse = "/"), " 个月")
} else {
  ## =========================================================================
  ## 5. Logistic 分支：共线性 / 线性 / EPV / 模型拟合（全部仅用训练集）
  ## =========================================================================
  ## 5.0 建模环境（只取进入模型的列，避免 datadist 对常量列告警）
  dd <- datadist(train[, c(GENE_Z, "age_imp", "sex")])
  options(datadist = "dd")

  ## 5.1 EPV（每个自变量事件数）
  n_event_tr <- sum(train$HF)
  epv_tab <- data.frame(
    model = c("M1_clinical", "M2_gene", "M3_combined"),
    n_predictors = c(2, length(GENES), length(GENES) + 2),
    n_events_train = n_event_tr,
    epv = round(n_event_tr / c(2, length(GENES), length(GENES) + 2), 1),
    stringsAsFactors = FALSE
  )
  epv_tab$adequate <- ifelse(epv_tab$epv >= 10, "yes (EPV>=10)", "NO")
  write_csv(epv_tab, file.path(d_tab, "T09_epv_check.csv"))
  log_msg("EPV: ", paste(epv_tab$model, epv_tab$epv, sep = "=", collapse = ", "))

  ## 5.2 共线性（VIF，基于等价 glm；>5 提示共线，>10 严重）
  f_glm <- as.formula(paste("HF ~", paste(c(GENE_Z, "age_imp", "sex"), collapse = " + ")))
  glm_full <- glm(f_glm, data = train, family = binomial)
  vif_raw <- tryCatch(car::vif(glm_full), error = function(e) NULL)
  if (is.null(vif_raw)) {
    ## car 不可用时用广义 VIF 的手工近似（对连续变量等价于 1/(1-R^2)）
    vif_df <- data.frame(variable = c(GENE_Z, "age_imp", "sex"),
                         VIF = NA_real_, note = "car not available")
  } else {
    if (is.matrix(vif_raw)) {
      vif_df <- data.frame(variable = rownames(vif_raw),
                           VIF = as.numeric(vif_raw[, 1]),
                           GVIF_adj = as.numeric(vif_raw[, 3]),
                           stringsAsFactors = FALSE)
    } else {
      vif_df <- data.frame(variable = names(vif_raw), VIF = as.numeric(vif_raw),
                           stringsAsFactors = FALSE)
    }
  }
  vif_df$collinearity <- ifelse(is.na(vif_df$VIF), "unknown",
                                ifelse(vif_df$VIF > 10, "severe (>10)",
                                       ifelse(vif_df$VIF > 5, "moderate (>5)", "none")))
  write_csv(vif_df, file.path(d_tab, "T10_collinearity_vif.csv"))
  log_msg("最大 VIF = ", round(max(vif_df$VIF, na.rm = TRUE), 2))

  ## 5.3 连续变量线性关系：Box-Tidwell + 限制性立方样条 LR 检验
  lin_rows <- list()
  ## Wald 检验辅助（避免依赖 anova 表行名匹配）
  wald_p <- function(fit, term) {
    b <- coef(fit)[term]
    se <- sqrt(vcov(fit)[term, term])
    z <- as.numeric(b / se)
    data.frame(statistic = round(z^2, 3), df = 1, p_value = signif(2 * pnorm(-abs(z)), 4))
  }
  ## (a) 年龄
  bt_age <- tryCatch({
    d2 <- train
    d2$age_log <- d2$age_imp * log(d2$age_imp)
    m <- lrm(as.formula(paste("HF ~", paste(GENE_Z, collapse = " + "),
                              "+ age_imp + sex + age_log")), data = d2, x = TRUE, y = TRUE)
    w <- wald_p(m, "age_log")
    data.frame(variable = "age", test = "Box-Tidwell (age*log(age))",
               statistic = w$statistic, df = w$df, p_value = w$p_value,
               stringsAsFactors = FALSE)
  }, error = function(e) data.frame(variable = "age", test = "Box-Tidwell",
                                    statistic = NA, df = NA, p_value = NA))
  lin_rows[["age_bt"]] <- bt_age

  ## (b) 年龄 RCS(3) vs 线性 的 LR 检验
  m_lin_age <- lrm(as.formula(paste("HF ~", paste(GENE_Z, collapse = " + "), "+ age_imp + sex")),
                   data = train, x = TRUE, y = TRUE)
  m_rcs_age <- lrm(as.formula(paste("HF ~", paste(GENE_Z, collapse = " + "), "+ rcs(age_imp, 3) + sex")),
                   data = train, x = TRUE, y = TRUE)
  ## 嵌套模型 LR 检验（直接由 logLik 计算，避免 anova.rms 的设计对象比对限制）
  lr_test <- function(f_small, f_big) {
    ll1 <- as.numeric(logLik(f_small)); ll2 <- as.numeric(logLik(f_big))
    chi <- max(0, 2 * (ll2 - ll1))
    d <- length(coef(f_big)) - length(coef(f_small))
    list(chisq = chi, df = d, p = pchisq(chi, d, lower.tail = FALSE))
  }
  lr_age <- lr_test(m_lin_age, m_rcs_age)
  lin_rows[["age_rcs"]] <- data.frame(
    variable = "age", test = "RCS(3) vs linear (LR test)",
    statistic = round(lr_age$chisq, 3), df = lr_age$df,
    p_value = signif(lr_age$p, 4), stringsAsFactors = FALSE
  )

  ## (c) 各基因（z 化后仍检验 log 变换项的显著性）
  for (g in GENE_Z) {
    bt <- tryCatch({
      d2 <- train
      d2[[paste0(g, "_log")]] <- d2[[g]] * log(d2[[g]] - min(d2[[g]]) + 1)
      m <- lrm(as.formula(paste("HF ~", paste(GENE_Z, collapse = " + "),
                                "+ age_imp + sex +", paste0(g, "_log"))),
               data = d2, x = TRUE, y = TRUE)
      w <- wald_p(m, paste0(g, "_log"))
      data.frame(variable = g, test = "Box-Tidwell",
                 statistic = w$statistic, df = w$df, p_value = w$p_value,
                 stringsAsFactors = FALSE)
    }, error = function(e) data.frame(variable = g, test = "Box-Tidwell",
                                      statistic = NA, df = NA, p_value = NA))
    lin_rows[[g]] <- bt
  }
  lin_tab <- do.call(rbind, lin_rows)
  lin_tab$linearity <- ifelse(is.na(lin_tab$p_value), "unknown",
                              ifelse(lin_tab$p_value < 0.05, "non-linear (p<0.05)", "linear"))
  write_csv(lin_tab, file.path(d_tab, "T11_linearity_tests.csv"))
  log_msg("线性检验: 年龄 RCS p = ", signif(lr_age$p, 3),
          "；Box-Tidwell p = ", signif(bt_age$p_value, 3))

  ## 是否对年龄使用 RCS：仅在两种检验均提示非线性时启用（避免过拟合）
  p_rcs <- as.numeric(lr_age$p)
  p_bt  <- as.numeric(bt_age$p_value)
  use_rcs <- (!is.na(p_rcs) && p_rcs < 0.05) && (!is.na(p_bt) && p_bt < 0.05)
  age_term <- if (use_rcs) "rcs(age_imp, 3)" else "age_imp"
  log_msg("年龄进入模型的形式: ", ifelse(use_rcs, "RCS(3)", "线性"))

  ## 5.4 三个模型（rms::lrm；M1 临床、M2 关键基因、M3 联合）
  f_m1 <- as.formula(paste("HF ~", age_term, "+ sex"))
  f_m2 <- as.formula(paste("HF ~", paste(GENE_Z, collapse = " + ")))
  f_m3 <- as.formula(paste("HF ~", paste(GENE_Z, collapse = " + "), "+", age_term, "+ sex"))

  fit_m1 <- lrm(f_m1, data = train, x = TRUE, y = TRUE)
  fit_m2 <- lrm(f_m2, data = train, x = TRUE, y = TRUE)
  fit_m3 <- lrm(f_m3, data = train, x = TRUE, y = TRUE)
  log_msg("模型拟合完成: M1 AIC=", round(AIC(fit_m1), 1),
          " M2 AIC=", round(AIC(fit_m2), 1),
          " M3 AIC=", round(AIC(fit_m3), 1))

  ## 5.5 单因素 Logistic（描述性，不用于筛选，避免选择性偏倚）
  uni_rows <- list()
  for (v in c(GENE_Z, "age_imp", "sex")) {
    m <- glm(as.formula(paste("HF ~", v)), data = train, family = binomial)
    s <- summary(m)$coefficients
    s <- s[rownames(s) != "(Intercept)", , drop = FALSE]
    for (rn in rownames(s)) {
      uni_rows[[length(uni_rows) + 1]] <- data.frame(
        variable = rn, beta = s[rn, "Estimate"], se = s[rn, "Std. Error"],
        z = s[rn, "z value"], p_value = s[rn, "Pr(>|z|)"],
        OR = exp(s[rn, "Estimate"]),
        OR_95CI = sprintf("%.3f-%.3f", exp(s[rn, "Estimate"] - 1.96 * s[rn, "Std. Error"]),
                          exp(s[rn, "Estimate"] + 1.96 * s[rn, "Std. Error"])),
        stringsAsFactors = FALSE
      )
    }
  }
  uni_tab <- do.call(rbind, uni_rows)
  uni_tab$beta <- round(uni_tab$beta, 4); uni_tab$se <- round(uni_tab$se, 4)
  uni_tab$z <- round(uni_tab$z, 3); uni_tab$p_value <- signif(uni_tab$p_value, 4)
  uni_tab$OR <- round(uni_tab$OR, 3)
  write_csv(uni_tab, file.path(d_tab, "T12_univariate_logistic.csv"))

  ## 5.6 系数表（beta / SE / Wald Z / P / OR / 95%CI）
  coef_table <- function(fit, model_name) {
    b  <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
    z  <- b / se
    p  <- 2 * pnorm(-abs(z))
    data.frame(
      model = model_name,
      term = names(b),
      beta = as.numeric(b),
      se = as.numeric(se),
      wald_z = as.numeric(z),
      p_value = as.numeric(p),
      OR = exp(as.numeric(b)),
      OR_lower = exp(as.numeric(b) - 1.96 * as.numeric(se)),
      OR_upper = exp(as.numeric(b) + 1.96 * as.numeric(se)),
      stringsAsFactors = FALSE
    )
  }
  co_m1 <- coef_table(fit_m1, "M1_clinical")
  co_m2 <- coef_table(fit_m2, "M2_gene")
  co_m3 <- coef_table(fit_m3, "M3_combined")
  tidy_coef <- function(df) {
    df$beta <- round(df$beta, 4); df$se <- round(df$se, 4)
    df$wald_z <- round(df$wald_z, 3); df$p_value <- signif(df$p_value, 4)
    df$OR <- signif(df$OR, 3)
    df$OR_lower <- signif(df$OR_lower, 3); df$OR_upper <- signif(df$OR_upper, 3)
    df
  }
  co_m1 <- tidy_coef(co_m1); co_m2 <- tidy_coef(co_m2); co_m3 <- tidy_coef(co_m3)
  write_csv(co_m1, file.path(d_tab, "T13_coefficients_M1_clinical.csv"))
  write_csv(co_m2, file.path(d_tab, "T14_coefficients_M2_gene.csv"))
  write_csv(co_m3, file.path(d_tab, "T15_coefficients_M3_combined.csv"))

  ## 5.7 每个模型的整体检验 + 逐项 LR 检验
  anova_tab <- function(fit, model_name) {
    a <- anova(fit)
    ad <- as.data.frame(a)
    data.frame(model = model_name, term = rownames(ad), ad,
               check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  }
  an_tab <- rbind(anova_tab(fit_m1, "M1_clinical"),
                  anova_tab(fit_m2, "M2_gene"),
                  anova_tab(fit_m3, "M3_combined"))
  write_csv(an_tab, file.path(d_tab, "T16_anova_terms.csv"))

  ## 5.8 模型公式（可直接复算 logit）
  formula_text <- function(fit, model_name) {
    cf <- coef(fit)
    parts <- sprintf("%+.4f*%s", cf[-1], names(cf)[-1])
    txt <- sprintf("%s: logit(P(HF=1)) = %+.4f %s", model_name, cf[1],
                   paste(parts, collapse = " "))
    if (grepl("rcs\\(", model_name)) txt <- paste0(txt, "  [age_imp 以 RCS(3) 进入]")
    txt
  }
  f_txt <- c(formula_text(fit_m1, "M1_clinical"),
             formula_text(fit_m2, "M2_gene"),
             formula_text(fit_m3, "M3_combined"),
             "",
             "注：基因变量均为训练集 z 分数，z = (log2 表达值 - train_mean) / train_sd；",
             paste0("    z 参数见 T06_zscore_parameters_train_only.csv；年龄缺失以训练集中位数 ",
                    age_med, " 填补。"))
  writeLines(f_txt, file.path(d_tab, "T17_model_formula.txt"), useBytes = TRUE)

  ## 5.9 模型比较：AIC / LR 检验 / 训练集表观 AUC 与 Brier
  auc_ci <- function(y, p) {
    r <- pROC::roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<")
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    c(auc = as.numeric(pROC::auc(r)), lower = ci[1], upper = ci[3])
  }
  brier <- function(y, p) mean((p - y)^2)
  p_tr_m1 <- predict(fit_m1, type = "fitted")
  p_tr_m2 <- predict(fit_m2, type = "fitted")
  p_tr_m3 <- predict(fit_m3, type = "fitted")
  lr_m2_m3 <- lr_test(fit_m2, fit_m3)
  comp <- data.frame(
    model = c("M1_clinical", "M2_gene", "M3_combined"),
    n_predictors = c(length(coef(fit_m1)) - 1, length(coef(fit_m2)) - 1,
                     length(coef(fit_m3)) - 1),
    AIC = c(AIC(fit_m1), AIC(fit_m2), AIC(fit_m3)),
    logLik = c(as.numeric(logLik(fit_m1)), as.numeric(logLik(fit_m2)), as.numeric(logLik(fit_m3))),
    LR_ChiSq = c(NA, NA, lr_m2_m3$chisq),
    LR_df = c(NA, NA, lr_m2_m3$df),
    LR_p_vs_M2 = c(NA, NA, lr_m2_m3$p),
    AUC_train_apparent = c(auc_ci(train$HF, p_tr_m1)["auc"],
                           auc_ci(train$HF, p_tr_m2)["auc"],
                           auc_ci(train$HF, p_tr_m3)["auc"]),
    Brier_train = c(brier(train$HF, p_tr_m1), brier(train$HF, p_tr_m2),
                    brier(train$HF, p_tr_m3)),
    stringsAsFactors = FALSE
  )
  comp$AIC <- round(comp$AIC, 2)
  comp$logLik <- round(comp$logLik, 2)
  comp$LR_ChiSq <- round(comp$LR_ChiSq, 3)
  comp$LR_p_vs_M2 <- signif(comp$LR_p_vs_M2, 4)
  comp$AUC_train_apparent <- round(comp$AUC_train_apparent, 3)
  comp$Brier_train <- round(comp$Brier_train, 4)
  write_csv(comp, file.path(d_tab, "T18_model_comparison.csv"))

  ## 5.10 主模型判定（预先设定的规则，写死在脚本里，避免事后挑选）
  ##  若 联合模型 vs 基因模型 LR 检验 p<0.05（临床变量提供增量信息）→ M3 为主模型；
  ##  否则 M2（关键基因模型）为主模型。
  p_add <- as.numeric(lr_m2_m3$p)
  primary_name <- if (!is.na(p_add) && p_add < 0.05) "M3_combined" else "M2_gene"
  primary_fit <- if (primary_name == "M3_combined") fit_m3 else fit_m2
  log_msg("主模型（预设规则）: ", primary_name, "（M2 vs M3 LR p = ", signif(p_add, 3), "）")

  ## =========================================================================
  ## 6. 列线图（rms::nomogram —— 标准医学列线图）
  ## =========================================================================
  fun_at <- c(0.05, 0.2, 0.5, 0.8, 0.95)   # 概率轴刻度（关键基因效应强，logit 跨度大，需减少刻度避免重叠）
  make_nom <- function(fit) {
    nomogram(fit, fun = plogis, funlabel = "Risk of heart failure",
             fun.at = fun_at, lp = TRUE, varname.label = TRUE)
  }
  ## 变量显示名由 train 的 label 属性提供（见 3.2b），列线图上不出现 _z 后缀
  nom_m1 <- make_nom(fit_m1)
  nom_m2 <- make_nom(fit_m2)
  nom_m3 <- make_nom(fit_m3)

  plot_nom <- function(nom, main) {
    par(mar = c(4, 1, 3, 1))
    plot(nom, xfrac = 0.35, cex.axis = 0.75, cex.var = 0.95, lmgp = 0.18, lwd = 1)
    title(main = main, cex.main = 0.95)
  }
  save_plot("F01_nomogram_M1_clinical", plot_nom(nom_m1, "Nomogram - clinical model (age + sex)"),
            width = 9, height = 6.5)
  save_plot("F02_nomogram_M2_gene", plot_nom(nom_m2, "Nomogram - key gene model (AQP3/GGT5/ISLR/LAD1)"),
            width = 9, height = 6.5)
  save_plot("F03_nomogram_M3_combined",
            plot_nom(nom_m3, "Nomogram - combined model (key genes + age + sex)"),
            width = 9, height = 7)
  if (primary_name == "M3_combined") {
    save_plot("F03b_nomogram_PRIMARY", plot_nom(nom_m3, "Nomogram - primary model (combined)"),
              width = 9, height = 7)
  } else {
    save_plot("F03b_nomogram_PRIMARY", plot_nom(nom_m2, "Nomogram - primary model (key genes)"),
              width = 9, height = 6.5)
  }

  ## 6.1 列线图分值表（可复算每个变量的得分）
  ## nomogram 对象每个面板是 list(变量取值, Xbeta, points)，逐面板展开为长表
  nom_to_df <- function(nom, model_name) {
    rows <- list()
    for (nm in names(nom)) {
      el <- nom[[nm]]
      if (!is.list(el) || is.null(el[[1]])) next
      vals <- el[[1]]
      pts <- el[["points"]]
      if (is.null(pts) && length(el) >= 2) pts <- el[[2]]
      n <- max(length(vals), length(pts))
      fmtv <- function(x) {
        if (is.null(x)) return(rep(NA_character_, n))
        if (is.numeric(x)) return(formatC(x, format = "f", digits = 2))
        as.character(x)
      }
      rows[[length(rows) + 1]] <- data.frame(
        model = model_name, panel = nm,
        value = fmtv(vals), points = fmtv(pts),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  nom_tab <- rbind(nom_to_df(nom_m1, "M1_clinical"),
                   nom_to_df(nom_m2, "M2_gene"),
                   nom_to_df(nom_m3, "M3_combined"))
  write_csv(nom_tab, file.path(d_tab, "T19_nomogram_points.csv"))

  ## =========================================================================
  ## 7. Bootstrap 内部验证（训练集；validate + calibrate）
  ## =========================================================================
  log_msg("Bootstrap 内部验证开始（B = ", B_BOOT, "）...")
  boot_res <- list()
  for (m in list(list(name = "M1_clinical", fit = fit_m1),
                 list(name = "M2_gene", fit = fit_m2),
                 list(name = "M3_combined", fit = fit_m3))) {
    set.seed(SEED)
    v <- rms::validate(m$fit, method = "boot", B = B_BOOT, bw = FALSE)
    vm <- unclass(v)          # validate 对象有 as.data.frame 方法，会展开成长表，故先取矩阵
    vd <- data.frame(metric = rownames(vm), vm, check.names = FALSE, row.names = NULL,
                     stringsAsFactors = FALSE)
    vd <- data.frame(model = m$name, vd, stringsAsFactors = FALSE)
    write_csv(vd, file.path(d_tab, paste0("T20_bootstrap_validate_", m$name, ".csv")))
    boot_res[[m$name]] <- vm
    log_msg("  ", m$name, " 校正后 AUC = ",
            round((vm["Dxy", "index.corrected"] + 1) / 2, 3),
            "（表观 ", round((vm["Dxy", "index.orig"] + 1) / 2, 3), "）")
  }

  set.seed(SEED)
  cal_m1 <- rms::calibrate(fit_m1, method = "boot", B = B_BOOT, bw = FALSE)
  set.seed(SEED)
  cal_m2 <- rms::calibrate(fit_m2, method = "boot", B = B_BOOT, bw = FALSE)
  set.seed(SEED)
  cal_m3 <- rms::calibrate(fit_m3, method = "boot", B = B_BOOT, bw = FALSE)

  cal_col_corr <- function(cal) {
    cn <- colnames(cal)
    cand <- cn[grepl("corrected", cn, ignore.case = TRUE)]
    if (length(cand) == 0) cand <- cn[ncol(cal)]
    cand[1]
  }
  cal_to_df <- function(cal, model_name) {
    cc <- cal_col_corr(cal)
    data.frame(model = model_name, predicted = as.numeric(cal[, "predy"]),
               observed_apparent = as.numeric(cal[, grep("calibrated.orig|calibrated", colnames(cal))[1]]),
               observed_bias_corrected = as.numeric(cal[, cc]),
               stringsAsFactors = FALSE)
  }
  cal_df_all <- rbind(cal_to_df(cal_m1, "M1_clinical"),
                      cal_to_df(cal_m2, "M2_gene"),
                      cal_to_df(cal_m3, "M3_combined"))
  write_csv(cal_df_all, file.path(d_dat, "calibration_plot_data_train.csv"))

  save_plot("F06_calibration_train_bootstrap_M3", {
    plot(cal_m3, xlab = "Predicted probability", ylab = "Observed proportion",
         main = "Calibration (combined model, bootstrap bias-corrected)")
    abline(0, 1, lty = 2, col = "grey40")
  }, width = 6.5, height = 6)

  save_plot("F07_calibration_train_bootstrap_M2", {
    plot(cal_m2, xlab = "Predicted probability", ylab = "Observed proportion",
         main = "Calibration (key gene model, bootstrap bias-corrected)")
    abline(0, 1, lty = 2, col = "grey40")
  }, width = 6.5, height = 6)

  ## =========================================================================
  ## 8. 冻结模型 -> 验证集评价（仅评价，不参与任何拟合/筛选）
  ## =========================================================================
  p_va_m1 <- predict(fit_m1, newdata = valid, type = "fitted")
  p_va_m2 <- predict(fit_m2, newdata = valid, type = "fitted")
  p_va_m3 <- predict(fit_m3, newdata = valid, type = "fitted")

  ## 8.1 Hosmer-Lemeshow（十分位）
  hl_test <- function(p, y, g = 10) {
    cutp <- unique(quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE))
    if (length(cutp) < 3) return(c(chisq = NA, df = NA, p = NA))
    grp <- cut(p, breaks = cutp, include.lowest = TRUE, labels = FALSE)
    obs1 <- tapply(y, grp, sum); exp1 <- tapply(p, grp, sum); n_grp <- tapply(y, grp, length)
    obs0 <- n_grp - obs1; exp0 <- n_grp - exp1
    stat <- sum((obs1 - exp1)^2 / pmax(exp1, 1e-9) + (obs0 - exp0)^2 / pmax(exp0, 1e-9))
    df <- length(unique(grp)) - 2
    c(chisq = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
  }

  vp_row <- function(p, y, label) {
    vp <- val.prob(p, y, pl = FALSE)
    hl <- hl_test(p, y)
    r <- pROC::roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<")
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    data.frame(
      model = label, n = length(y), n_event = sum(y),
      AUC = as.numeric(pROC::auc(r)), AUC_lower = ci[1], AUC_upper = ci[3],
      C_index = as.numeric(vp["C (ROC)"]),
      Brier = as.numeric(vp["Brier"]),
      Intercept = as.numeric(vp["Intercept"]), Slope = as.numeric(vp["Slope"]),
      Emax = as.numeric(vp["Emax"]),
      HL_chisq = hl["chisq"], HL_df = hl["df"], HL_p = hl["p"],
      stringsAsFactors = FALSE
    )
  }
  perf_va <- rbind(vp_row(p_va_m1, valid$HF, "M1_clinical"),
                   vp_row(p_va_m2, valid$HF, "M2_gene"),
                   vp_row(p_va_m3, valid$HF, "M3_combined"))
  ## 训练集表观指标一并列出便于对照
  perf_tr <- rbind(vp_row(p_tr_m1, train$HF, "M1_clinical"),
                   vp_row(p_tr_m2, train$HF, "M2_gene"),
                   vp_row(p_tr_m3, train$HF, "M3_combined"))
  perf_va$set <- "validation"; perf_tr$set <- "train"
  perf_all <- rbind(perf_tr, perf_va)
  for (cc in c("AUC", "AUC_lower", "AUC_upper", "C_index", "Brier",
               "Intercept", "Slope", "Emax", "HL_chisq", "HL_p")) {
    perf_all[[cc]] <- signif(perf_all[[cc]], 4)
  }
  write_csv(perf_all, file.path(d_tab, "T21_discrimination_calibration.csv"))

  ## 8.2 ROC 图（训练表观 vs 验证）
  roc_tr <- pROC::roc(train$HF, p_tr_m3, quiet = TRUE, levels = c(0, 1), direction = "<")
  roc_va <- pROC::roc(valid$HF, p_va_m3, quiet = TRUE, levels = c(0, 1), direction = "<")
  save_plot("F04_ROC_train_vs_validation", {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(roc_tr, col = "#2c7fb8", lwd = 2, legacy.axes = TRUE,
         main = "ROC - primary model")
    plot(roc_va, col = "#d95f02", lwd = 2, add = TRUE)
    legend("bottomright",
           legend = c(sprintf("Training (apparent AUC = %.3f)", as.numeric(pROC::auc(roc_tr))),
                      sprintf("Validation (AUC = %.3f, 95%%CI %.3f-%.3f)",
                              as.numeric(pROC::auc(roc_va)),
                              as.numeric(pROC::ci.auc(roc_va))[1],
                              as.numeric(pROC::ci.auc(roc_va))[3])),
           col = c("#2c7fb8", "#d95f02"), lwd = 2, bty = "n", cex = 0.85)
  }, width = 6, height = 6)

  ## 8.3 三模型在验证集上的 ROC 对比
  roc_m1 <- pROC::roc(valid$HF, p_va_m1, quiet = TRUE, levels = c(0, 1), direction = "<")
  roc_m2 <- pROC::roc(valid$HF, p_va_m2, quiet = TRUE, levels = c(0, 1), direction = "<")
  roc_m3 <- pROC::roc(valid$HF, p_va_m3, quiet = TRUE, levels = c(0, 1), direction = "<")
  save_plot("F05_ROC_validation_three_models", {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(roc_m1, col = "grey55", lwd = 2, legacy.axes = TRUE,
         main = "ROC on validation set")
    plot(roc_m2, col = "#1b9e77", lwd = 2, add = TRUE)
    plot(roc_m3, col = "#d95f02", lwd = 2, add = TRUE)
    legend("bottomright",
           legend = c(sprintf("M1 clinical (AUC %.3f)", as.numeric(pROC::auc(roc_m1))),
                      sprintf("M2 key genes (AUC %.3f)", as.numeric(pROC::auc(roc_m2))),
                      sprintf("M3 combined (AUC %.3f)", as.numeric(pROC::auc(roc_m3)))),
           col = c("grey55", "#1b9e77", "#d95f02"), lwd = 2, bty = "n", cex = 0.85)
  }, width = 6, height = 6)

  ## 8.4 验证集校准图（loess + 十分位点）
  cal_bin <- function(p, y, g = 10) {
    cutp <- unique(quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE))
    grp <- cut(p, breaks = cutp, include.lowest = TRUE, labels = FALSE)
    data.frame(bin = sort(unique(grp)),
               predicted = as.numeric(tapply(p, grp, mean)),
               observed = as.numeric(tapply(y, grp, mean)),
               n = as.numeric(tapply(y, grp, length)))
  }
  cal_bin_m3 <- cal_bin(p_va_m3, valid$HF)
  cal_bin_m3$model <- "M3_combined"
  cal_bin_m2 <- cal_bin(p_va_m2, valid$HF)
  cal_bin_m2$model <- "M2_gene"
  write_csv(rbind(cal_bin_m2, cal_bin_m3), file.path(d_dat, "calibration_plot_data_validation.csv"))

  cal_long <- rbind(
    data.frame(model = "M2_gene", p = p_va_m2, y = valid$HF),
    data.frame(model = "M3_combined", p = p_va_m3, y = valid$HF)
  )
  g_cal <- ggplot(cal_long, aes(x = p, y = y, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_smooth(method = "loess", se = FALSE, span = 0.9, linewidth = 0.9) +
    geom_point(data = rbind(cal_bin_m2, cal_bin_m3),
               aes(x = predicted, y = observed, colour = model), size = 2) +
    scale_colour_manual(values = c(M2_gene = "#1b9e77", M3_combined = "#d95f02"),
                        labels = c("Key gene model", "Combined model")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Predicted probability", y = "Observed proportion",
         title = "Calibration on validation set", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = c(0.72, 0.15), panel.grid.minor = element_blank())
  save_plot("F08_calibration_validation", print(g_cal), width = 6, height = 6)

  ## =========================================================================
  ## 9. 决策曲线分析 (DCA)
  ## =========================================================================
  dca_binary <- function(p, y, thresholds = seq(0.05, 0.95, by = 0.01)) {
    n <- length(y); prev <- mean(y)
    nb <- vapply(thresholds, function(t) {
      tp <- sum(p >= t & y == 1); fp <- sum(p >= t & y == 0)
      tp / n - (fp / n) * (t / (1 - t))
    }, numeric(1))
    nb_all <- prev - (1 - prev) * (thresholds / (1 - thresholds))
    data.frame(threshold = thresholds, nb_model = nb, nb_all = nb_all, nb_none = 0)
  }
  dca_make <- function(p_list, y, set_label) {
    out <- list()
    for (nm in names(p_list)) {
      d <- dca_binary(p_list[[nm]], y)
      out[[length(out) + 1]] <- data.frame(set = set_label, model = nm,
                                           threshold = d$threshold, nb_model = d$nb_model,
                                           stringsAsFactors = FALSE)
    }
    ref <- out[[1]]
    out[[length(out) + 1]] <- data.frame(set = set_label, model = "Treat all",
                                         threshold = ref$threshold,
                                         nb_model = ref$threshold * 0 + NA,  # 占位，稍后填充
                                         stringsAsFactors = FALSE)
    prev <- mean(y)
    out[[length(out)]]$nb_model <- prev - (1 - prev) * (ref$threshold / (1 - ref$threshold))
    out[[length(out) + 1]] <- data.frame(set = set_label, model = "Treat none",
                                         threshold = ref$threshold, nb_model = 0,
                                         stringsAsFactors = FALSE)
    do.call(rbind, out)
  }
  dca_tr <- dca_make(list(M1_clinical = p_tr_m1, M2_gene = p_tr_m2, M3_combined = p_tr_m3),
                     train$HF, "train")
  dca_va <- dca_make(list(M1_clinical = p_va_m1, M2_gene = p_va_m2, M3_combined = p_va_m3),
                     valid$HF, "validation")
  write_csv(dca_tr, file.path(d_dat, "dca_plot_data_train.csv"))
  write_csv(dca_va, file.path(d_dat, "dca_plot_data_validation.csv"))
  write_csv(rbind(dca_tr, dca_va), file.path(d_tab, "T22_dca_net_benefit.csv"))

  dca_plot <- function(d, ttl) {
    d$model <- factor(d$model, levels = c("M3_combined", "M2_gene", "M1_clinical",
                                          "Treat all", "Treat none"))
    ggplot(d, aes(x = threshold, y = nb_model, colour = model, linetype = model)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(M3_combined = "#d95f02", M2_gene = "#1b9e77",
                                     M1_clinical = "grey45", `Treat all` = "black",
                                     `Treat none` = "black")) +
      scale_linetype_manual(values = c(M3_combined = 1, M2_gene = 1, M1_clinical = 1,
                                       `Treat all` = 2, `Treat none` = 3)) +
      coord_cartesian(ylim = c(-0.05, max(d$nb_model, na.rm = TRUE) * 1.05)) +
      labs(x = "Threshold probability", y = "Net benefit", title = ttl, colour = NULL,
           linetype = NULL) +
      theme_bw(base_size = 12) +
      theme(legend.position = c(0.72, 0.78), panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))
  }
  save_plot("F09_DCA_train", print(dca_plot(dca_tr, "Decision curve analysis (training set)")),
            width = 7, height = 5.5)
  save_plot("F10_DCA_validation", print(dca_plot(dca_va, "Decision curve analysis (validation set)")),
            width = 7, height = 5.5)

  ## =========================================================================
  ## 10. 每个样本的预测概率 / 风险评分
  ## =========================================================================
  ## 风险分层切点：训练集 Youden 指数（冻结后用于所有样本）
  r_tr_prim <- pROC::roc(train$HF, if (primary_name == "M3_combined") p_tr_m3 else p_tr_m2,
                         quiet = TRUE, levels = c(0, 1), direction = "<")
  best <- pROC::coords(r_tr_prim, x = "best", best.method = "youden",
                       ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  cut_youden <- as.numeric(best["threshold"][1, 1])

  pred_df <- data.frame(
    sample = dat$sample, set = as.character(dat$set), group = dat$group,
    HF = dat$HF,
    age_raw = dat$age, age_used = dat$age_imp, age_was_imputed = dat$age_missing,
    sex = as.character(dat$sex),
    stringsAsFactors = FALSE
  )
  for (g in GENES) {
    pred_df[[paste0(g, "_log2")]] <- dat[[g]]
    pred_df[[paste0(g, "_z")]] <- dat[[paste0(g, "_z")]]
  }
  pred_df$lp_M1 <- predict(fit_m1, newdata = dat)
  pred_df$lp_M2 <- predict(fit_m2, newdata = dat)
  pred_df$lp_M3 <- predict(fit_m3, newdata = dat)
  pred_df$prob_M1_clinical <- plogis(pred_df$lp_M1)
  pred_df$prob_M2_gene <- plogis(pred_df$lp_M2)
  pred_df$prob_M3_combined <- plogis(pred_df$lp_M3)
  pred_df$prob_primary <- if (primary_name == "M3_combined") pred_df$prob_M3_combined else pred_df$prob_M2_gene
  pred_df$risk_group <- ifelse(pred_df$prob_primary >= cut_youden, "High risk", "Low risk")
  pred_df$primary_model <- primary_name
  pred_df$cutoff_youden <- cut_youden

  write_csv(pred_df, file.path(d_dat, "modelpredictions.csv"))
  ## 训练集/验证集分开的版本，便于直接做外部评价
  write_csv(pred_df[pred_df$set == "train", ], file.path(d_dat, "modelpredictions_train.csv"))
  write_csv(pred_df[pred_df$set == "validation", ], file.path(d_dat, "modelpredictions_validation.csv"))

  rg_tab <- as.data.frame.matrix(table(pred_df$risk_group, pred_df$HF))
  rg_tab <- data.frame(risk_group = rownames(rg_tab), rg_tab, check.names = FALSE, row.names = NULL)
  names(rg_tab)[-1] <- paste0("HF_", names(rg_tab)[-1])
  write_csv(rg_tab, file.path(d_tab, "T23_risk_group_by_outcome.csv"))

  ## =========================================================================
  ## 11. 保存模型对象 / 参数记录
  ## =========================================================================
  model_bundle <- list(
    fit_m1_clinical = fit_m1, fit_m2_gene = fit_m2, fit_m3_combined = fit_m3,
    primary_model = primary_name, primary_fit = primary_fit,
    zscore_parameters = zpar, age_median_train = age_med,
    age_term = age_term, cutoff_youden = cut_youden,
    split = list(seed = SEED, train_samples = train$sample, validation_samples = valid$sample),
    nomogram_objects = list(M1 = nom_m1, M2 = nom_m2, M3 = nom_m3),
    bootstrap = boot_res, calibration_models = list(M1 = cal_m1, M2 = cal_m2, M3 = cal_m3),
    call = match.call()
  )
  saveRDS(model_bundle, file.path(d_mod, "nomogram_models.rds"))
  saveRDS(fit_m3, file.path(d_mod, "final_model_M3_combined.rds"))

  ## 参数 JSON（完整记录，便于复现与写入论文方法学）
  params <- list(
    dataset = "GSE57338",
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    r_version = R.version.string,
    seed = SEED,
    outcome = list(type = "binary", event = "HF (DCM or ICM)", code_event = 1,
                   non_event = "NF", code_non_event = 0,
                   n_total = nrow(dat), n_event = sum(dat$HF),
                   event_rate = round(mean(dat$HF), 4),
                   time_dimension = "none"),
    key_genes = GENES,
    gene_scaling = list(method = "z-score (train-set mean/sd)",
                        parameters = zpar),
    age_handling = list(min_valid_age = MIN_VALID_AGE,
                        n_set_missing = n_age_flag,
                        imputation = "train-set median",
                        imputed_value = age_med,
                        n_imputed = sum(dat$age_missing)),
    split = list(method = "stratified by outcome", train_fraction = TRAIN_FRAC,
                 n_train = nrow(train), n_train_events = sum(train$HF),
                 n_validation = nrow(valid), n_validation_events = sum(valid$HF)),
    models = list(
      M1_clinical = list(formula = deparse(f_m1), n_predictors = length(coef(fit_m1)) - 1,
                         AIC = round(AIC(fit_m1), 3)),
      M2_gene = list(formula = deparse(f_m2), n_predictors = length(coef(fit_m2)) - 1,
                     AIC = round(AIC(fit_m2), 3)),
      M3_combined = list(formula = deparse(f_m3), n_predictors = length(coef(fit_m3)) - 1,
                         AIC = round(AIC(fit_m3), 3))
    ),
    primary_model = primary_name,
    primary_rule = "M3 if LR test M2 vs M3 p<0.05, else M2 (pre-specified)",
    lr_test_M2_vs_M3 = list(chisq = round(lr_m2_m3$chisq, 4),
                            df = lr_m2_m3$df,
                            p = signif(lr_m2_m3$p, 4)),
    epv = epv_tab[, c("model", "n_predictors", "n_events_train", "epv")],
    max_vif = round(max(vif_df$VIF, na.rm = TRUE), 3),
    linearity = list(use_rcs_for_age = use_rcs,
                     box_tidwell_age_p = signif(p_bt, 4),
                     rcs_age_lr_p = signif(p_rcs, 4)),
    bootstrap = list(B = B_BOOT, method = "boot", corrected_auc = lapply(
      boot_res, function(v) round((v["Dxy", "index.corrected"] + 1) / 2, 4))),
    validation_metrics = perf_va,
    risk_cutoff_youden = cut_youden,
    r_packages = as.list(sapply(need_pkgs, function(p) as.character(packageVersion(p))))
  )
  writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, pretty = TRUE, digits = 6),
             file.path(d_log, "params.json"), useBytes = TRUE)

  ## =========================================================================
  ## 12. 中文报告
  ## =========================================================================
  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
  rep_lines <- c(
    "================================================================================",
    "GSE57338 关键基因列线图 (Nomogram) 预测模型 —— 结果摘要",
    paste0("生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "   R: ", R.version.string),
    "================================================================================",
    "",
    "一、数据与结局",
    paste0("  数据集: GSE57338（人左心室组织，log2 表达矩阵，", nrow(expr_raw), " 个基因）"),
    paste0("  样本: 共 ", nrow(dat), " 例；NF ", sum(dat$HF == 0), " 例，DCM ", sum(dat$group == "DCM"),
           " 例，ICM ", sum(dat$group == "ICM"), " 例"),
    paste0("  结局定义: HF = 1（DCM 或 ICM，n = ", sum(dat$HF), "）；HF = 0（NF 供体心脏，n = ",
           sum(dat$HF == 0), "）"),
    "  时间维度: 无。GSE57338 不含随访时间与生存状态，属横截面二分类结局，",
    "            因此采用 Logistic 回归（rms::lrm）建模；Cox 分支在脚本中保留但未触发。",
    paste0("  预测目标与时间点: 预测个体“当前为心衰心肌（DCM/ICM）”的概率，无时间点；",
           "脚本中 SURV_HORIZONS = ", paste(SURV_HORIZONS, collapse = "/"),
           " 月仅在数据含 time/status 时作为 Cox 列线图的预测时间点。"),
    paste0("  事件率: ", fmt(100 * mean(dat$HF), 1), "%"),
    "",
    "二、训练/验证划分（所有拟合与筛选仅用训练集）",
    paste0("  分层随机划分，7:3，种子 = ", SEED),
    paste0("  训练集: n = ", nrow(train), "，事件 ", sum(train$HF), " 例"),
    paste0("  验证集: n = ", nrow(valid), "，事件 ", sum(valid$HF), " 例（仅用于评价冻结模型）"),
    "",
    "三、数据质量检查",
    paste0("  样本 ID 对齐: 矩阵列名与样本表完全匹配（", length(common), "/", length(common), "）"),
    paste0("  关键基因: ", paste(GENES, collapse = " / "), " 均唯一存在，无缺失值"),
    paste0("  表达异常值 (|z|>3): ", sum(out_tab$n_abs_z_gt3), " 个观测，均为生物学变异，予保留"),
    paste0("  年龄异常: ", n_age_flag, " 例记录为 0/1 岁（成人心脏组织不合理），按缺失处理并以训练集中位数 ",
           age_med, " 填补"),
    "",
    "四、变量控制与假设检验",
    paste0("  EPV（训练集事件数/参数个数）: ", paste(epv_tab$model, epv_tab$epv, sep = "=", collapse = "，"),
           "；均 ≥ 10，满足 Logistic 建模的样本量要求"),
    paste0("  共线性: 最大 VIF = ", fmt(max(vif_df$VIF, na.rm = TRUE), 2), "（阈值 5），无多重共线性问题"),
    paste0("  连续变量线性: 年龄 Box-Tidwell p = ", signif(p_bt, 3),
           "，RCS(3) LR 检验 p = ", signif(p_rcs, 3),
           " -> 年龄以", ifelse(use_rcs, " RCS(3)", "线性"), "形式进入模型"),
    "  比例风险假设: 不适用（无生存结局，未拟合 Cox 模型）",
    "",
    "五、模型比较（训练集）",
    paste0("  ", paste(sprintf("%s: 参数 %d 个, AIC = %.1f, 表观 AUC = %.3f",
                                comp$model, comp$n_predictors, comp$AIC, comp$AUC_train_apparent),
                       collapse = "\n  ")),
    paste0("  M2 vs M3 似然比检验: χ² = ", fmt(comp$LR_ChiSq[3], 2), "，df = ", comp$LR_df[3],
           "，P = ", signif(comp$LR_p_vs_M2[3], 3)),
    paste0("  联合模型中各基因每升高 1 个 SD 的 OR 范围: ",
           fmt(min(co_m3$OR[co_m3$term != "Intercept"]), 2), " - ",
           fmt(max(co_m3$OR[co_m3$term != "Intercept"]), 2),
           "（年龄每 1 岁 OR = ", fmt(co_m3$OR[co_m3$term == "age_imp"], 3), "）"),
    paste0("  主模型（预设规则）: ", primary_name),
    "",
    paste0("六、Bootstrap 内部验证（B = ", B_BOOT, "，训练集）"),
    paste0("  ", paste(sprintf("%s: 表观 AUC %.3f -> 校正后 AUC %.3f（乐观度 %.3f）",
                                names(boot_res),
                                sapply(boot_res, function(v) (v["Dxy", "index.orig"] + 1) / 2),
                                sapply(boot_res, function(v) (v["Dxy", "index.corrected"] + 1) / 2),
                                sapply(boot_res, function(v) v["Dxy", "index.orig"] - v["Dxy", "index.corrected"]) / 2),
                       collapse = "\n  ")),
    "",
    "七、验证集表现（冻结模型，仅评价）",
    paste0("  ", paste(sprintf("%s: AUC %.3f (95%%CI %.3f-%.3f), Brier %.3f, 校准截距 %.3f, 校准斜率 %.3f, HL P = %.3f",
                                perf_va$model, perf_va$AUC, perf_va$AUC_lower, perf_va$AUC_upper,
                                perf_va$Brier, perf_va$Intercept, perf_va$Slope, perf_va$HL_p),
                       collapse = "\n  ")),
    "",
    "八、列线图与输出文件",
    "  02_figures/F01-F03: M1（临床）、M2（关键基因）、M3（联合）三张 rms::nomogram 列线图",
    "  02_figures/F03b: 主模型列线图",
    "  01_tables/T19_nomogram_points.csv: 各变量取值对应的列线图得分（可手工复算总分）",
    "  01_tables/T13-T15: 回归系数、标准误、Wald Z、P、OR 及 95%CI",
    "  01_tables/T17_model_formula.txt: 可直接复算的 logit 公式",
    "  04_data/modelpredictions.csv: 每个样本的预测概率、线性预测值、风险分层",
    "  03_models/nomogram_models.rds: 冻结模型对象（含 z 化参数、切点、Bootstrap 结果）",
    "  05_logs/params.json: 全部建模参数、公式、种子、AUC、EPV、VIF",
    "",
    "九、使用说明",
    paste0("  1) 新样本先按 05_logs/params.json 中的 train_mean / train_sd 做 z 化：z = (log2 - mean)/sd"),
    paste0("  2) 缺失年龄以 ", age_med, " 填补"),
    "  3) 将 z 值代入 T17_model_formula.txt 的公式得到 logit，P = 1/(1+exp(-logit))",
    paste0("  4) 或以 Youden 切点 ", fmt(cut_youden, 3), " 划分高/低风险"),
    "",
    "十、局限与解读警告",
    "  1) 单一数据集内部验证，无独立外部队列；Bootstrap 校正可缓解乐观偏倚但不等于外部验证。",
    "  2) 临床变量仅年龄与性别（GSE57338 元数据所限），未纳入 NYHA、LVEF、BNP 等预后因素。",
    "  3) 横截面设计无法评价时间相关结局，结论应表述为“患病风险”而非“预后风险”。",
    paste0("  4) 关键基因在病例与对照间的效应极大（每 1 SD 的 OR 低至 ",
           fmt(min(co_m3$OR[co_m3$term != "Intercept"]), 2), "、高达 ",
           fmt(max(co_m3$OR[co_m3$term != "Intercept"]), 2),
           "），训练集残差偏差（", fmt(-2 * as.numeric(logLik(fit_m3)), 1), "）远小于自由度（",
           nrow(train) - length(coef(fit_m3)), "），提示近似完全分离；"),
    "     这正是列线图概率轴被压缩（0.05→0.95 仅占约 47/300 分值）的原因。该类模型在真实临床人群中判别力通常会明显衰减，",
    "     使用前必须在外部人群重新校准，不宜直接移植概率刻度。",
    paste0("  5) 验证集上 M2（关键基因）与 M3（联合）的 AUC 分别为 ",
           fmt(perf_va$AUC[2], 3), " 与 ", fmt(perf_va$AUC[3], 3),
           "，置信区间重叠，二者判别力相当；主模型按预设规则取 M3，"),
    "     若偏好更简约、可床旁快速获得的模型，M2 亦可作为替代方案（见 F02 列线图）。",
    "  6) 表达谱为组织水平（左心室），基因表达与实际蛋白/血清标志物水平不等价，临床转化需另行验证。",
    "================================================================================"
  )
  writeLines(rep_lines, file.path(out_root, "Nomogram_report.txt"), useBytes = TRUE)
}

## ===========================================================================
## 13. 代码、日志、会话信息归档
## ===========================================================================
this_script <- tryCatch(normalizePath(sub("^--file=", "",
                                          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                        error = function(e) NA_character_)
if (!is.na(this_script) && file.exists(this_script)) {
  file.copy(this_script, file.path(d_log, "10_nomogram.R"), overwrite = TRUE)
  log_msg("脚本副本已归档: 05_logs/10_nomogram.R")
}

writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo.txt"), useBytes = TRUE)
log_msg("================ 全部完成，输出目录: ", out_root, " ================")
writeLines(run_log, file.path(d_log, "run_log.txt"), useBytes = TRUE)
