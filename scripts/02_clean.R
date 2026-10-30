# =============================================================================
# 02_clean.R  ——  GSE57338 表达矩阵清洗（探针注释 + 去重 + 低表达过滤）
#
# 输入（raw_data/，只读）：
#   GSE57338_series_matrix.txt.gz   表达矩阵（33,297 探针 × 313 样本，已 RMA+log2）
#                                   + 样本信息（!Sample_* 元数据行，44 列）
#   GPL11532.annot.gz               平台注释（探针 → Gene symbol / Entrez ID）
#   download_log.txt                下载日志（非数据文件，仅记录）
#
# 输出（clean_data/）：
#   GSE57338_clean_expression_matrix.csv   基因 × 样本（行名=Gene symbol）
#   GSE57338_sample_info.csv               样本信息表（含分组）
#   GSE57338_gene_annotation.csv           探针→基因 映射明细（含被平均的探针数）
#   GSE57338_qc_summary.txt                清洗质控汇总（机器可读）
# 输出（report/）：
#   sessionInfo_02_clean.txt
#
# 用法：cd 到项目根目录后
#   Rscript scripts/02_clean.R
#
# 设计说明：
#   1) 表达值已是 RMA 归一化 + log2 转换（range 0.9~13.5），因此本脚本**不再做**
#      log 转换或分位数归一化，避免二次标准化引入假象；只做注释、合并、过滤。
#   2) 过滤采用「阈值 + 样本比例」双条件：log2 表达 > MIN_EXPR 的样本占比
#      >= MIN_FRAC 才保留，比单纯 rowMeans 阈值更稳健（不受少数极端样本拉高均值）。
#   3) 多探针映射同一基因时按均值合并（用户指定），而非取最大表达/最大方差探针。
# =============================================================================

rm(list = ls())

# ---------- 0. 路径与参数 ----------
proj_dir  <- getwd()
raw_dir   <- file.path(proj_dir, "raw_data")
clean_dir <- file.path(proj_dir, "clean_data")
rep_dir   <- file.path(proj_dir, "report")

if (!dir.exists(raw_dir)) stop("找不到 raw_data/，请先 cd 到项目根目录再运行。")
for (d in c(clean_dir, rep_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

accession <- "GSE57338"
platform  <- "GPL11532"

MIN_EXPR     <- 3.0    # log2 表达阈值（约等于原始信号 8，接近背景水平之上）
MIN_FRAC     <- 0.25   # 至少在 25% 的样本中超过阈值才保留
NA_FRAC_MAX  <- 0.20   # 单个探针缺失比例上限，超过即丢弃
MULTI_SYMBOL <- "first" # 一个探针注释到多个 symbol（"A///B"）时："first"取第一个，"drop"丢弃
ROUND_DIGITS <- 4      # 输出矩阵保留小数位

t_start <- Sys.time()

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 68), "\n", sep = "")

# 统一用 UTF-8 写文本，避免 Windows 中文乱码
write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con)
}

# ---------- 1. 识别 raw_data/ 中每个文件的用途 ----------
msg(">>> [1/7] 识别 raw_data/ 中的文件 ...")

all_files <- list.files(raw_dir, full.names = TRUE)
is_matrix <- grepl("_series_matrix\\.txt(\\.gz)?$", basename(all_files), ignore.case = TRUE)
is_annot  <- grepl("^GPL\\d+\\.annot(\\.gz)?$",      basename(all_files), ignore.case = TRUE)
is_log    <- grepl("^(download|readme|.*log).*\\.(txt|md)$", basename(all_files), ignore.case = TRUE)

file_role <- ifelse(is_matrix, "表达矩阵 + 样本信息（series matrix）",
             ifelse(is_annot,  "平台注释（probe → gene）",
             ifelse(is_log,    "下载/说明日志（非数据）", "未识别")))

file_tab <- data.frame(
  文件      = basename(all_files),
  大小_MB   = round(file.info(all_files)$size / 1024^2, 2),
  用途      = file_role,
  stringsAsFactors = FALSE
)
print(file_tab, row.names = FALSE)

matrix_file <- all_files[is_matrix][1]
annot_file  <- all_files[is_annot][1]
if (is.na(matrix_file)) stop("未找到 series matrix 文件（*_series_matrix.txt.gz）。")
msg("    表达矩阵文件 : ", basename(matrix_file))
msg("    平台注释文件 : ", ifelse(is.na(annot_file), "未找到（将退回 ID 原样保留）", basename(annot_file)))

# ---------- 2. 读取表达矩阵与样本信息 ----------
msg(">>> [2/7] 读取表达矩阵与样本信息 ...")

suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })

eset <- getGEO(filename = matrix_file, getGPL = FALSE)
if (is.list(eset) && !is(eset, "ExpressionSet")) eset <- eset[[1]]

expr_raw <- exprs(eset)
pd_raw   <- pData(eset)

n_probe_raw <- nrow(expr_raw)
n_sample    <- ncol(expr_raw)
msg("    原始矩阵：", n_probe_raw, " 探针 x ", n_sample, " 样本")
msg("    表达值范围：", paste(round(range(expr_raw, na.rm = TRUE), 3), collapse = " ~ "),
    "（已 log2，无需再转换）")
msg("    缺失值总数：", sum(is.na(expr_raw)))

# ---------- 3. 判断基因 ID 类型 ----------
msg(">>> [3/7] 判断表达矩阵中的 ID 类型 ...")

ids <- rownames(expr_raw)
id_type <- if (all(grepl("^ENS[A-Z]*G\\d+", ids))) {
  "Ensembl ID"
} else if (all(grepl("^\\d+$", ids))) {
  "芯片探针 ID（Affymetrix transcript cluster，纯数字）"
} else if (mean(grepl("^[A-Za-z][A-Za-z0-9\\-\\.]*$", ids)) > 0.9) {
  "Gene symbol"
} else {
  "芯片探针 ID（含字母数字混合）"
}
msg("    判定结果：", id_type)
msg("    ID 示例：", paste(head(ids, 3), collapse = ", "))

# ---------- 4. 探针 ID → Gene symbol ----------
msg(">>> [4/7] 探针 ID → Gene symbol 注释 ...")

read_gpl <- function(path) {
  # .annot.gz 头部有 ^Annotation / !Annotation_* / #ID 说明行，数据表头以 "ID\t" 开始
  con <- gzfile(path, "r")
  head_lines <- readLines(con, n = 200, warn = FALSE)
  close(con)
  skip_n <- max(grep("^ID\t", head_lines)) - 1L
  read.delim(gzfile(path), header = TRUE, sep = "\t", quote = "\"",
             comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
             skip = skip_n)
}

ann <- NULL
if (!is.na(annot_file)) {
  ann <- read_gpl(annot_file)
  msg("    注释表：", nrow(ann), " 行 x ", ncol(ann), " 列")
  msg("    注释列示例：", paste(head(colnames(ann), 6), collapse = " | "))
}

# 统一取列名
pick_col <- function(df, cands) {
  hit <- cands[cands %in% colnames(df)][1]
  if (is.na(hit)) NULL else hit
}
col_id  <- if (!is.null(ann)) pick_col(ann, c("ID", "ID_ref", "Probe_Set_ID", "SPOT_ID")) else NULL
col_sym <- if (!is.null(ann)) pick_col(ann, c("Gene symbol", "GENE_SYMBOL", "Symbol", "gene_assignment")) else NULL
col_ens <- if (!is.null(ann)) pick_col(ann, c("Ensembl", "ENSEMBL", "Ensembl ID")) else NULL

sym_map <- NULL   # 最终 probe -> symbol 映射

if (id_type == "Gene symbol") {
  # 本身就是 symbol：直接使用，仅做格式清洗
  sym_map <- data.frame(probe = ids, symbol = ids, n_candidates = 1L, stringsAsFactors = FALSE)
  msg("    ID 已是 Gene symbol，无需转换")

} else if (grepl("^Ensembl", id_type)) {
  # Ensembl ID：去版本号 → 转 symbol
  probe_nover <- sub("\\.[0-9]+$", "", ids)
  msg("    去除 Ensembl 版本号：", sum(probe_nover != ids), " 个 ID 带版本号后缀")
  map_tab <- NULL
  if (!is.null(ann) && !is.null(col_ens) && !is.null(col_sym)) {
    map_tab <- data.frame(
      probe  = sub("\\.[0-9]+$", "", as.character(ann[[col_ens]])),
      symbol = as.character(ann[[col_sym]]), stringsAsFactors = FALSE)
  } else if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(probe_nover),
                               keytype = "ENSEMBL", columns = "SYMBOL")
    map_tab <- data.frame(probe = m$ENSEMBL, symbol = m$SYMBOL, stringsAsFactors = FALSE)
  }
  if (is.null(map_tab)) stop("Ensembl ID 需要注释文件或 org.Hs.eg.db 才能转换，二者均不可用。")
  map_tab <- map_tab[!duplicated(map_tab$probe), ]
  sym_map <- data.frame(probe = ids,
                        symbol = map_tab$symbol[match(probe_nover, map_tab$probe)],
                        n_candidates = 1L, stringsAsFactors = FALSE)

} else {
  # 芯片探针 ID：必须依赖平台注释
  if (is.null(ann) || is.null(col_sym)) stop("检测到探针 ID，但缺少可用的平台注释文件。")
  m <- data.frame(probe  = as.character(ann[[col_id]]),
                  symbol = as.character(ann[[col_sym]]),
                  stringsAsFactors = FALSE)
  m <- m[!duplicated(m$probe), ]
  sym_map <- data.frame(probe = ids,
                        symbol = m$symbol[match(ids, m$probe)],
                        stringsAsFactors = FALSE)
  # 记录一个探针对应几个候选 symbol（"A///B" 形式）
  n_cand <- ifelse(is.na(sym_map$symbol), 0L,
                   lengths(strsplit(sym_map$symbol, "///", fixed = TRUE)))
  sym_map$n_candidates <- n_cand
}

# 处理 "A///B" 多基因注释
multi_probe <- sum(sym_map$n_candidates > 1, na.rm = TRUE)
if (multi_probe > 0) {
  if (MULTI_SYMBOL == "first") {
    sym_map$symbol <- trimws(sub("///.*$", "", sym_map$symbol))
  } else {
    sym_map$symbol[sym_map$n_candidates > 1] <- NA
  }
}
sym_map$symbol <- trimws(sym_map$symbol)

# 识别无效注释
bad_sym <- is.na(sym_map$symbol) | sym_map$symbol == "" | sym_map$symbol %in% c("---", "NA", "null")
n_bad <- sum(bad_sym)
sym_map$annotated <- !bad_sym

annot_rate <- 100 * sum(!bad_sym) / nrow(sym_map)
msg("    可注释探针：", sum(!bad_sym), " / ", nrow(sym_map),
    sprintf("（%.2f%%）", annot_rate))
msg("    无法注释（空/---）：", n_bad)
msg("    多基因注释探针（A///B）：", multi_probe,
    ifelse(MULTI_SYMBOL == "first", " → 取第一个 symbol", " → 丢弃"))

# 仅保留可注释探针
keep_annot <- sym_map$annotated
expr_ann   <- expr_raw[keep_annot, , drop = FALSE]
sym_ann    <- sym_map$symbol[keep_annot]

# ---------- 5. 多探针 → 单基因（均值） ----------
msg(">>> [5/7] 同一基因多探针合并（取均值）...")

sums   <- rowsum(expr_ann, group = sym_ann, reorder = TRUE)
counts <- rowsum(rep(1, nrow(expr_ann)), group = sym_ann, reorder = TRUE)[, 1]
gene_mat <- sums / counts
n_gene_before_filter <- nrow(gene_mat)

probe_per_gene <- table(sym_ann)
msg("    合并后基因数：", n_gene_before_filter,
    "（", nrow(expr_ann), " 个探针，平均 ", round(nrow(expr_ann) / n_gene_before_filter, 2),
    " 探针/基因，最多 ", max(probe_per_gene), " 个）")

# ---------- 6. 数据清洗 ----------
msg(">>> [6/7] 数据清洗（缺失 / 低表达 / 零方差）...")

# 6.1 缺失值：探针层面缺失比例过高的已在上游丢弃；这里检查合并后
na_frac_gene <- rowMeans(is.na(gene_mat))
drop_na <- na_frac_gene > NA_FRAC_MAX
# 剩余少量缺失以该基因行中位数填补（仅影响极少数位点，报告中记录条数）
n_na_cell <- sum(is.na(gene_mat))
if (n_na_cell > 0) {
  for (i in which(rowSums(is.na(gene_mat)) > 0)) {
    v <- gene_mat[i, ]
    v[is.na(v)] <- median(v, na.rm = TRUE)
    gene_mat[i, ] <- v
  }
}

# 6.2 低表达过滤
pass_expr <- rowSums(gene_mat > MIN_EXPR, na.rm = TRUE) >= MIN_FRAC * ncol(gene_mat)

# 6.3 零方差（完全无变化，无法用于差异分析）
gene_var <- apply(gene_mat, 1, var)
pass_var <- gene_var > 0

keep_gene <- !drop_na & pass_expr & pass_var
gene_clean <- gene_mat[keep_gene, , drop = FALSE]

msg("    缺失比例过高剔除：", sum(drop_na), " 个基因（合并后缺失位点 ", n_na_cell, " 个）")
msg("    低表达剔除（log2 > ", MIN_EXPR, " 的样本 < ", MIN_FRAC * 100, "%）：",
    sum(!pass_expr), " 个基因")
msg("    零方差剔除：", sum(!pass_var), " 个基因")

# ---------- 7. 整理样本信息 ----------
msg(">>> [7/7] 整理样本信息 ...")

pd <- pd_raw
get_col <- function(df, cands) {
  hit <- cands[cands %in% colnames(df)][1]
  if (is.na(hit)) rep(NA_character_, nrow(df)) else as.character(df[[hit]])
}

sample_info <- data.frame(
  sample        = get_col(pd, c("geo_accession")),
  title         = get_col(pd, c("title")),
  disease       = get_col(pd, c("disease status:ch1", "disease state:ch1", "source_name_ch1")),
  heart_failure = get_col(pd, c("heart failure:ch1")),
  gender        = get_col(pd, c("gender:ch1", "Sex:ch1")),
  age           = suppressWarnings(as.numeric(get_col(pd, c("age:ch1", "Age:ch1")))),
  tissue        = get_col(pd, c("tissue:ch1")),
  platform      = get_col(pd, c("platform_id")),
  stringsAsFactors = FALSE
)
sample_info$disease <- trimws(sample_info$disease)

# 疾病分组缩写： idiopathic dilated CMP -> DCM；ischemic -> ICM；non-failing -> NF
sample_info$group <- ifelse(grepl("non-failing", sample_info$disease, ignore.case = TRUE), "NF",
                     ifelse(grepl("ischemic",    sample_info$disease, ignore.case = TRUE), "ICM",
                     ifelse(grepl("dilated",     sample_info$disease, ignore.case = TRUE), "DCM", NA)))
sample_info$group <- factor(sample_info$group, levels = c("NF", "DCM", "ICM"))
# 保证样本顺序与表达矩阵列一致
sample_info <- sample_info[match(colnames(gene_clean), sample_info$sample), ]
stopifnot(!any(is.na(sample_info$sample)))

# 与 02b_sample_info.R 保持同一 schema：sample / group 固定在前两列。
# 完整的协变量表（含 group_source、patient_id、一致性标记等）由 02b 生成，
# 若已跑过 02b，请在其后重跑 02b 以恢复完整版本。
sample_info <- sample_info[, c("sample", "group",
                               setdiff(names(sample_info), c("sample", "group")))]

msg("    分组统计：")
print(table(sample_info$group, useNA = "ifany"))
msg("    性别统计：")
print(table(sample_info$gender, useNA = "ifany"))

# ---------- 8. 输出 ----------
msg(">>> 写出结果 ...")

gene_clean_out <- round(gene_clean, ROUND_DIGITS)

mat_path <- file.path(clean_dir, paste0(accession, "_clean_expression_matrix.csv"))
mat_df <- data.frame(GeneSymbol = rownames(gene_clean_out), gene_clean_out,
                     check.names = FALSE, stringsAsFactors = FALSE)
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fwrite(mat_df, mat_path, quote = FALSE, na = "NA")
} else {
  write.csv(mat_df, mat_path, row.names = FALSE, quote = FALSE)
}
msg("    表达矩阵：", basename(mat_path), "（", nrow(mat_df), " 基因 x ", ncol(mat_df) - 1, " 样本，",
    round(file.info(mat_path)$size / 1024^2, 1), " MB）")

info_path <- file.path(clean_dir, paste0(accession, "_sample_info.csv"))
write.csv(sample_info, info_path, row.names = FALSE, fileEncoding = "UTF-8")
msg("    样本信息：", basename(info_path))

# 探针→基因 映射明细（仅可注释探针）
map_out <- sym_map[sym_map$annotated, c("probe", "symbol", "n_candidates")]
map_out$n_probes_in_gene <- as.integer(probe_per_gene[map_out$symbol])
map_out <- map_out[order(map_out$symbol), ]
map_path <- file.path(clean_dir, paste0(accession, "_gene_annotation.csv"))
write.csv(map_out, map_path, row.names = FALSE, fileEncoding = "UTF-8")
msg("    注释明细：", basename(map_path))

# ---------- 9. 质控汇总 ----------
n_probe_kept_raw <- sum(keep_gene[sym_map$annotated])  # 保留基因对应的原始探针数（近似）
qc <- c(
  "================================================================",
  paste0(accession, " 表达矩阵清洗 质控汇总"),
  "================================================================",
  paste0("生成时间     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " (", format(Sys.time(), "%Z"), ")"),
  paste0("平台         : ", platform, "  [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array"),
  paste0("原始数据类型 : 芯片（RMA 归一化 + log2，GEO 已处理）"),
  "",
  "-------------------------- 维度 --------------------------",
  paste0("原始矩阵     : ", n_probe_raw, " 探针 x ", n_sample, " 样本"),
  paste0("ID 类型      : ", id_type),
  paste0("注释后基因   : ", n_gene_before_filter, " 基因 x ", n_sample, " 样本（多探针取均值）"),
  paste0("清洗后矩阵   : ", nrow(gene_clean), " 基因 x ", ncol(gene_clean), " 样本"),
  "",
  "-------------------------- 注释 --------------------------",
  paste0("可注释探针   : ", sum(!bad_sym), " / ", nrow(sym_map),
         sprintf("（%.2f%%）", annot_rate)),
  paste0("无法注释探针 : ", n_bad, sprintf("（%.2f%%）", 100 * n_bad / nrow(sym_map))),
  paste0("多基因探针   : ", multi_probe, "（处理方式：", MULTI_SYMBOL, "）"),
  paste0("唯一基因符号 : ", n_gene_before_filter),
  paste0("平均探针/基因: ", round(nrow(expr_ann) / n_gene_before_filter, 2)),
  paste0("单探针基因数 : ", sum(probe_per_gene == 1),
         sprintf("（%.1f%%）", 100 * mean(probe_per_gene == 1))),
  "",
  "-------------------------- 过滤 --------------------------",
  paste0("剔除-无法注释: ", n_bad, " 探针"),
  paste0("剔除-缺失过多: ", sum(drop_na), " 基因（阈值 > ", NA_FRAC_MAX * 100, "% 缺失）"),
  paste0("剔除-低表达  : ", sum(!pass_expr), " 基因（log2 > ", MIN_EXPR,
         " 的样本 < ", MIN_FRAC * 100, "%）"),
  paste0("剔除-零方差  : ", sum(!pass_var), " 基因"),
  paste0("保留基因     : ", nrow(gene_clean)),
  "",
  "-------------------------- 比例 --------------------------",
  paste0("探针保留率   : ", sprintf("%.2f%%", 100 * nrow(expr_ann) / n_probe_raw),
         "（", nrow(expr_ann), " / ", n_probe_raw, "）"),
  paste0("基因保留率   : ", sprintf("%.2f%%", 100 * nrow(gene_clean) / n_gene_before_filter),
         "（相对注释后基因数 ", n_gene_before_filter, "）"),
  paste0("基因保留率   : ", sprintf("%.2f%%", 100 * nrow(gene_clean) / n_probe_raw),
         "（相对原始探针数 ", n_probe_raw, "，含去冗余效果）"),
  paste0("整体清洗比例 : 原始 ", n_probe_raw, " 探针 → 最终 ", nrow(gene_clean), " 基因，",
         sprintf("减少 %.2f%%", 100 * (1 - nrow(gene_clean) / n_probe_raw))),
  "",
  "-------------------------- 样本 --------------------------",
  paste0("样本数       : ", ncol(gene_clean)),
  paste0("缺失值       : 原始 ", sum(is.na(expr_raw)), " 个；合并后 ", n_na_cell, " 个（已按行中位数填补）"),
  paste0("分组分布     : ", paste(names(table(sample_info$group)),
                              as.integer(table(sample_info$group)), sep = "=", collapse = ", ")),
  paste0("性别分布     : ", paste(names(table(sample_info$gender)),
                              as.integer(table(sample_info$gender)), sep = "=", collapse = ", ")),
  paste0("年龄范围     : ", paste(range(sample_info$age, na.rm = TRUE), collapse = " ~ ")),
  paste0("表达值范围   : ", paste(round(range(gene_clean), 3), collapse = " ~ "),
         "（log2，未做二次标准化）"),
  "",
  "-------------------------- 参数 --------------------------",
  paste0("MIN_EXPR     : ", MIN_EXPR),
  paste0("MIN_FRAC     : ", MIN_FRAC),
  paste0("NA_FRAC_MAX  : ", NA_FRAC_MAX),
  paste0("MULTI_SYMBOL : ", MULTI_SYMBOL),
  paste0("ROUND_DIGITS : ", ROUND_DIGITS),
  "",
  "================================================================"
)
qc_path <- file.path(clean_dir, paste0(accession, "_qc_summary.txt"))
write_utf8(qc, qc_path)
msg("    质控汇总：", basename(qc_path))

# sessionInfo
si_path <- file.path(rep_dir, "sessionInfo_02_clean.txt")
write_utf8(capture.output(sessionInfo()), si_path)

# ---------- 10. 控制台摘要 ----------
t_end <- Sys.time()
cat("\n", strrep("=", 68), "\n", sep = "")
cat("清洗完成\n")
cat(strrep("=", 68), "\n", sep = "")
cat(sprintf("%-16s %s\n", "原始矩阵", paste0(n_probe_raw, " 探针 x ", n_sample, " 样本")))
cat(sprintf("%-16s %s\n", "注释后",   paste0(n_gene_before_filter, " 基因（注释率 ", sprintf("%.2f%%", annot_rate), "）")))
cat(sprintf("%-16s %s\n", "清洗后",   paste0(nrow(gene_clean), " 基因 x ", ncol(gene_clean), " 样本")))
cat(sprintf("%-16s %s\n", "耗时",     paste0(round(as.numeric(difftime(t_end, t_start, units = "secs")), 1), " 秒")))
cat("\n输出目录：", normalizePath(clean_dir, winslash = "/"), "\n", sep = "")
cat(paste0("  - ", basename(mat_path), "\n"))
cat(paste0("  - ", basename(info_path), "\n"))
cat(paste0("  - ", basename(map_path), "\n"))
cat(paste0("  - ", basename(qc_path), "\n"))
