# =============================================================================
# 02b_sample_info.R  ——  GSE57338 样本分组信息标准化
#
# 目标：从 series matrix 的样本元数据中提取样本编号与实验分组，生成标准样本
#       信息表 <accession>_sample_info.csv（至少含 sample / group 两列），
#       并与 clean expression matrix 的列名做一致性校验、自动修正。
#
# 输入（raw_data/，只读）：
#   GSE57338_series_matrix.txt.gz    表达矩阵 + 样本元数据（313 样本 × 44 列 pData）
# 输入（clean_data/）：
#   GSE57338_clean_expression_matrix.csv   清洗后矩阵（仅读取表头，用于列名校验）
#
# 输出（clean_data/）：
#   GSE57338_sample_info.csv            标准样本信息表（sample, group, + 协变量）
#   GSE57338_group_counts.csv           每组样本数量统计
#   GSE57338_sample_name_check.csv      列名一致性校验明细（逐样本修正规则）
#   GSE57338_sample_info_qc.txt         分组质控汇总（机器可读）
# 输出（report/）：
#   sessionInfo_02b_sample_info.txt
#
# 用法：cd 到项目根目录后
#   Rscript scripts/02b_sample_info.R
#
# 关键设计：
#   1) 分组不是单字段硬编码，而是「四级规则链 + 多字段交叉验证」：
#      优先级 disease status:ch1 > source_name_ch1 > title > heart failure:ch1，
#      并对每个样本记录最终生效字段（group_source）与字段间是否矛盾（group_conflict）。
#   2) 样本名与矩阵列名比对采用分阶段规则（R0 精确 → R1 去空白 → R2 去引号
#      → R3 去 R 的 X 前缀 → R4 统一大小写 → R5 去非字母数字），命中哪一级就
#      记哪条规则，未命中即报错停止，不做静默丢弃。
#   3) 最终输出顺序严格等于表达矩阵列顺序，避免下游按位置取样本时错位。
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

# 分组缩写（顺序即因子水平：对照在前，便于下游 DCM-NF / ICM-NF 对比）
GROUP_LEVELS <- c("NF", "DCM", "ICM")
GROUP_FULL <- c(
  NF  = "non-failing（非衰竭对照）",
  DCM = "idiopathic dilated cardiomyopathy（特发性扩张型心肌病）",
  ICM = "ischemic cardiomyopathy（缺血性心肌病）"
)

t_start <- Sys.time()
msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 70), "\n", sep = "")

write_utf8 <- function(lines, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con)
}

LOG <- character(0)          # 全程日志，最后落盘
log_add <- function(...) LOG <<- c(LOG, paste0(...))

# ---------- 1. 读取 series matrix（离线，不下载 GPL） ----------
msg(">>> [1/6] 读取 series matrix ...")
suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })

matrix_file <- list.files(raw_dir, pattern = "_series_matrix\\.txt(\\.gz)?$",
                          full.names = TRUE, ignore.case = TRUE)[1]
if (is.na(matrix_file)) stop("未找到 *_series_matrix.txt.gz")

# getGPL = FALSE：样本信息与表达值都在 series matrix 内，无需联网下载平台注释
eset <- getGEO(filename = matrix_file, getGPL = FALSE)
if (is.list(eset) && !is(eset, "ExpressionSet")) eset <- eset[[1]]

pd        <- pData(eset)
expr_raw  <- exprs(eset)
ids_expr  <- colnames(expr_raw)          # 来源 A：原始矩阵列名
ids_pd    <- rownames(pd)                # 来源 B：pData 行名
ids_geo   <- as.character(pd$geo_accession)  # 来源 C：geo_accession 字段

msg("    pData：", nrow(pd), " 样本 x ", ncol(pd), " 列")
msg("    原始表达矩阵：", nrow(expr_raw), " 探针 x ", ncol(expr_raw), " 样本")

# 三源一致性：A vs B vs C 必须完全相同，否则说明 GEOquery 解析异常
stopifnot(identical(ids_expr, ids_pd), identical(ids_expr, ids_geo))
msg("    三处样本编号完全一致（colnames / rownames / geo_accession）：",
    length(ids_expr), " 个")
log_add("三源样本编号一致：colnames(exprs) = rownames(pData) = geo_accession = ",
        length(ids_expr), " 个，identical() 校验 TRUE")
log_add("样本编号唯一性：重复 ", sum(duplicated(ids_expr)), " 个；缺失 ",
        sum(is.na(ids_expr) | ids_expr == ""), " 个")

# ---------- 2. 字段抽取 ----------
msg(">>> [2/6] 抽取分组相关字段 ...")

get_col <- function(df, cands) {
  hit <- cands[cands %in% colnames(df)][1]
  if (is.na(hit)) rep(NA_character_, nrow(df)) else as.character(df[[hit]])
}

# characteristics_ch1 在该数据集中被拆成 5 个独立列（tissue / heart failure /
# disease status / gender / age），同时也保留了被 :ch1 后缀命名的同义列。
# 优先用 :ch1 显名列，回退到从 characteristics_ch1* 系列按 "key: value" 解析。
parse_characteristic <- function(df, key) {
  ch_cols <- grep("^characteristics_ch1", colnames(df), value = TRUE)
  out <- rep(NA_character_, nrow(df))
  for (cn in ch_cols) {
    v <- as.character(df[[cn]])
    hit <- grepl(paste0("^\\s*", key, "\\s*:"), v, ignore.case = TRUE)
    out[hit & is.na(out)] <- sub(paste0("^\\s*", key, "\\s*:\\s*"), "", v[hit & is.na(out)],
                                 ignore.case = TRUE)
  }
  out
}

pick_with_fallback <- function(df, cands, key = NULL) {
  v <- get_col(df, cands)
  if (all(is.na(v)) && !is.null(key)) v <- parse_characteristic(df, key)
  trimws(v)
}

disease_raw   <- pick_with_fallback(pd, c("disease status:ch1", "disease state:ch1"), "disease status")
source_name   <- pick_with_fallback(pd, c("source_name_ch1"))
title_raw     <- pick_with_fallback(pd, c("title"))
heart_failure <- pick_with_fallback(pd, c("heart failure:ch1"), "heart failure")
gender        <- pick_with_fallback(pd, c("gender:ch1", "Sex:ch1"), "gender")
age_raw       <- pick_with_fallback(pd, c("age:ch1", "Age:ch1"), "age")
tissue        <- pick_with_fallback(pd, c("tissue:ch1"), "tissue")
organism      <- pick_with_fallback(pd, c("organism_ch1"))
molecule      <- pick_with_fallback(pd, c("molecule_ch1"))
platform_col  <- pick_with_fallback(pd, c("platform_id"))
patient_code  <- pick_with_fallback(pd, c("description"))
age           <- suppressWarnings(as.numeric(age_raw))

msg("    disease status 取值：", paste(sort(unique(disease_raw)), collapse = " | "))
msg("    source_name    取值：", paste(sort(unique(source_name)), collapse = " | "))
msg("    heart failure 取值：", paste(sort(unique(heart_failure)), collapse = " | "))

# ---------- 3. 分组规则链 ----------
msg(">>> [3/6] 按规则链判定分组 ...")

# 单值分类器：把一个自由文本字段映射到 NF / DCM / ICM，无法识别返回 NA
# 注意顺序：必须先判 non-failing，否则 "non-failing" 不含 dilated/ischemic 不会误伤，
# 但 "non-ischemic" 这类表述只能落到 NA（由规则链下一级接手），不做猜测。
classify_text <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[grepl("non[-_ ]?failing|nonfailing|\\bnormal\\b|healthy|donor|control", s)] <- "NF"
  out[is.na(out) & grepl("dilated|\\bdcm\\b|idiopathic", s)]                     <- "DCM"
  out[is.na(out) & grepl("ischa?emic|\\bicm\\b", s)]                             <- "ICM"
  out
}

# heart failure 字段只能区分「是否衰竭」，无法区分 DCM/ICM，故单独处理
classify_hf <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[grepl("^no$|^n$|false|^0$", s)] <- "NF"   # 非衰竭
  out[grepl("^yes$|^y$|true|^1$", s)] <- "HF"   # 衰竭但亚型未知
  out
}

# 规则链：按优先级依次尝试，先命中者生效
rules <- list(
  list(field = "disease status:ch1", values = disease_raw,    fn = classify_text, note = "GEO 表型字段，最权威"),
  list(field = "source_name_ch1",    values = source_name,    fn = classify_text, note = "样本来源描述"),
  list(field = "title",              values = title_raw,      fn = classify_text, note = "样本标题"),
  list(field = "heart failure:ch1",  values = heart_failure,  fn = classify_hf,   note = "仅能区分衰竭/非衰竭")
)

group        <- rep(NA_character_, nrow(pd))
group_source <- rep(NA_character_, nrow(pd))

for (r in rules) {
  cand <- r$fn(r$values)
  need <- is.na(group) & !is.na(cand)
  group[need]        <- cand[need]
  group_source[need] <- r$field
}
# 第四级给出的是 "HF"（亚型未知），不能当作最终分组：降级为 NA 并标记待人工处理
hf_only <- !is.na(group) & group == "HF"
group[hf_only]        <- NA_character_
group_source[hf_only] <- "heart failure:ch1（仅判为衰竭，亚型未定）"

# 每级规则的独立判定结果，用于交叉验证
per_field <- data.frame(
  from_disease = classify_text(disease_raw),
  from_source  = classify_text(source_name),
  from_title   = classify_text(title_raw),
  from_hf      = classify_hf(heart_failure),
  stringsAsFactors = FALSE
)
per_field$from_hf[per_field$from_hf == "HF"] <- "DCM/ICM(unknown)"

# 冲突检测：把所有能判出亚型的字段拉平比较，出现 >=2 种不同 NF/DCM/ICM 即冲突
subtype_mat <- as.matrix(per_field[, c("from_disease", "from_source", "from_title")])
n_agree <- apply(subtype_mat, 1, function(v) {
  v <- v[!is.na(v)]; if (!length(v)) 0L else max(tabulate(match(v, unique(v))))
})
n_conflict <- apply(subtype_mat, 1, function(v) {
  v <- v[!is.na(v)]; length(unique(v)) > 1
})
# heart failure 与亚型的一致性：NF 应为 "no"，DCM/ICM 应为 "yes"
hf_consistent <- ifelse(is.na(group), NA,
                 ifelse(group == "NF", heart_failure == "no", heart_failure == "yes"))

msg("    最终生效字段分布：")
print(table(group_source, useNA = "ifany"))
msg("    字段间冲突样本数：", sum(n_conflict, na.rm = TRUE))

if (any(is.na(group))) {
  bad <- ids_expr[is.na(group)]
  stop("以下样本无法判定分组，请人工核查：", paste(bad, collapse = ", "))
}

group <- factor(group, levels = GROUP_LEVELS)
msg("    分组结果：", paste(names(table(group)), as.integer(table(group)), sep = "=", collapse = ", "))

# ---------- 4. 样本名与表达矩阵列名一致性校验 ----------
msg(">>> [4/6] 校验 sample 与 clean matrix 列名 ...")

mat_path <- file.path(clean_dir, paste0(accession, "_clean_expression_matrix.csv"))
if (!file.exists(mat_path)) {
  msg("    警告：未找到 clean matrix，跳过列名校验")
  mat_cols <- NULL
} else {
  # 只读取表头一行，避免把 38 MB 矩阵整体读入
  con <- file(mat_path, open = "r")
  hdr <- readLines(con, n = 1, warn = FALSE)
  close(con)
  hdr <- sub("^\ufeff", "", hdr)                       # 去 BOM
  mat_cols <- strsplit(hdr, ",", fixed = TRUE)[[1]]
  mat_cols <- gsub('^"|"$', "", mat_cols)              # 去引号
  mat_cols <- mat_cols[-1]                             # 首列为 GeneSymbol
  msg("    clean matrix 列数：", length(mat_cols), "（首列 GeneSymbol 已排除）")
}

# 分阶段匹配规则：每一级都是「对矩阵列名做变换后与 pData 样本名精确比对」
norm_rules <- list(
  list(id = "R0", desc = "无需修正（精确匹配）",          fn = function(x) x),
  list(id = "R1", desc = "去除首尾空白",                  fn = function(x) trimws(x)),
  list(id = "R2", desc = "去除包裹引号",                  fn = function(x) gsub('^["\']|["\']$', "", x)),
  list(id = "R3", desc = "去除 R read.csv 自动添加的 X 前缀", fn = function(x) sub("^X(?=GSM)", "", x)),
  list(id = "R4", desc = "统一为大写",                    fn = function(x) toupper(x)),
  list(id = "R5", desc = "去除非字母数字字符后匹配",       fn = function(x) gsub("[^A-Za-z0-9]", "", toupper(x)))
)

match_one <- function(tok, ref, ref_norm5) {
  for (r in norm_rules) {
    v <- r$fn(tok)
    if (r$id == "R5") {
      hit <- match(gsub("[^A-Za-z0-9]", "", toupper(v)), ref_norm5)
    } else {
      hit <- match(v, ref)
    }
    if (!is.na(hit)) return(list(idx = hit, rule = r$id, desc = r$desc))
  }
  list(idx = NA_integer_, rule = NA_character_, desc = "未匹配")
}

if (!is.null(mat_cols)) {
  ref <- ids_expr
  ref_norm5 <- gsub("[^A-Za-z0-9]", "", toupper(ref))
  res <- lapply(mat_cols, match_one, ref = ref, ref_norm5 = ref_norm5)
  map_idx  <- vapply(res, function(z) z$idx,  integer(1))
  map_rule <- vapply(res, function(z) z$rule, character(1))
  map_desc <- vapply(res, function(z) z$desc, character(1))

  n_unmatched <- sum(is.na(map_idx))
  if (n_unmatched > 0) {
    stop("有 ", n_unmatched, " 个矩阵列名无法映射到 series matrix 样本：",
         paste(mat_cols[is.na(map_idx)], collapse = ", "))
  }
  if (anyDuplicated(map_idx)) {
    stop("映射结果出现重复样本，请检查列名规范化规则：",
         paste(mat_cols[duplicated(map_idx)], collapse = ", "))
  }

  name_check <- data.frame(
    matrix_col   = mat_cols,
    matched_sample = ref[map_idx],
    rule         = map_rule,
    rule_desc    = map_desc,
    changed      = mat_cols != ref[map_idx],
    stringsAsFactors = FALSE
  )
  msg("    矩阵列 → 样本 映射：", nrow(name_check), " 列全部匹配，未匹配 ",
      n_unmatched, " 个")
  msg("    实际需要改写的列名：", sum(name_check$changed), " 个")
  print(table(name_check$rule, name_check$changed))

  # 集合层面双向校验
  only_mat <- setdiff(mat_cols, ref)
  only_ref <- setdiff(ref, mat_cols)
  msg("    仅存在于矩阵：", length(only_mat), " 个；仅存在于样本表：", length(only_ref), " 个")
  log_add("列名校验：矩阵 ", length(mat_cols), " 列 vs 样本 ", length(ref),
          " 个；集合差异 矩阵独有=", length(only_mat), "，样本独有=", length(only_ref))
  log_add("需要改写列名的数量：", sum(name_check$changed),
          "；按规则分布：", paste(names(table(name_check$rule[name_check$changed])),
                                as.integer(table(name_check$rule[name_check$changed])),
                                sep = "=", collapse = ", ", recycle0 = TRUE))
  # 顺序一致性：矩阵列顺序是否等于 series matrix 样本顺序
  order_same <- identical(mat_cols[map_idx], mat_cols) && identical(ref[map_idx], mat_cols)
  msg("    矩阵列顺序与 series matrix 顺序一致：", order_same)
  log_add("矩阵列顺序与 series matrix 样本顺序一致：", order_same)
} else {
  map_idx  <- seq_along(ids_expr)
  name_check <- data.frame(
    matrix_col = ids_expr, matched_sample = ids_expr, rule = "NO_MATRIX",
    rule_desc = "未找到矩阵，未校验", changed = FALSE, stringsAsFactors = FALSE)
}

# ---------- 5. 组装标准样本信息表 ----------
msg(">>> [5/6] 组装标准样本信息表 ...")

info <- data.frame(
  sample        = ids_expr[map_idx],                  # 与矩阵列名严格一致
  group         = as.character(group)[map_idx],       # NF / DCM / ICM
  group_full    = unname(GROUP_FULL[as.character(group)[map_idx]]),
  title         = title_raw[map_idx],
  patient_id    = patient_code[map_idx],
  disease_raw   = disease_raw[map_idx],
  heart_failure = heart_failure[map_idx],
  gender        = gender[map_idx],
  age           = age[map_idx],
  tissue        = tissue[map_idx],
  organism      = organism[map_idx],
  platform      = platform_col[map_idx],
  molecule      = molecule[map_idx],
  group_source  = group_source[map_idx],
  group_conflict = n_conflict[map_idx],
  hf_consistent = hf_consistent[map_idx],
  stringsAsFactors = FALSE
)
info$group <- factor(info$group, levels = GROUP_LEVELS)

# 硬校验：sample 必须与矩阵列名逐一相同、且顺序一致
if (!is.null(mat_cols)) {
  stopifnot(identical(as.character(info$sample), mat_cols))
}
stopifnot(!any(is.na(info$sample)), !anyDuplicated(info$sample), !any(is.na(info$group)))
msg("    sample 与矩阵列名逐一 identical：TRUE")

# ---------- 6. 输出 ----------
msg(">>> [6/6] 写出结果 ...")

fwrite_safe <- function(df, path, quote = FALSE) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, path, quote = quote, na = "NA")
  } else {
    write.csv(df, path, row.names = FALSE, quote = quote, na = "NA", fileEncoding = "UTF-8")
  }
}

# 6.1 标准样本信息表
info_path <- file.path(clean_dir, paste0(accession, "_sample_info.csv"))
fwrite_safe(info, info_path)
msg("    样本信息表：", basename(info_path), "（", nrow(info), " 行 x ", ncol(info), " 列）")

# 6.2 分组统计
tab  <- table(info$group)
cnt  <- data.frame(
  group      = names(tab),
  group_full = unname(GROUP_FULL[names(tab)]),
  n          = as.integer(tab),
  percent    = round(100 * as.integer(tab) / nrow(info), 2),
  stringsAsFactors = FALSE
)
cnt_path <- file.path(clean_dir, paste0(accession, "_group_counts.csv"))
fwrite_safe(cnt, cnt_path)
print(cnt, row.names = FALSE)

# 6.3 列名校验明细
chk_path <- file.path(clean_dir, paste0(accession, "_sample_name_check.csv"))
fwrite_safe(name_check, chk_path)
msg("    列名校验明细：", basename(chk_path))

# 6.4 质控汇总
xt_gender <- table(info$group, info$gender)
xt_hf     <- table(info$group, info$heart_failure)
age_stat  <- tapply(info$age, info$group, function(v)
  sprintf("%.0f (%.0f-%.0f)", median(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))

qc <- c(
  "================================================================",
  paste0(accession, " 样本分组信息 质控汇总"),
  "================================================================",
  paste0("生成时间   : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " (", format(Sys.time(), "%Z"), ")"),
  paste0("平台       : ", platform),
  paste0("样本数     : ", nrow(info)),
  paste0("分组字段   : disease status:ch1 / source_name_ch1 / title / heart failure:ch1"),
  "",
  "-------------------------- 分组 ----------------------------",
  paste0(sprintf("  %-4s %4d  %6.2f%%  %s", cnt$group, cnt$n, cnt$percent, cnt$group_full)),
  paste0("  合计 ", sum(cnt$n)),
  "",
  "-------------------------- 判定来源 ------------------------",
  paste0(sprintf("  %-45s %4d", names(table(info$group_source)),
                 as.integer(table(info$group_source)))),
  "",
  "-------------------------- 一致性 --------------------------",
  paste0("字段间亚型冲突样本 : ", sum(info$group_conflict)),
  paste0("heart failure 与亚型不一致 : ", sum(!info$hf_consistent, na.rm = TRUE)),
  paste0("group 缺失 : ", sum(is.na(info$group))),
  paste0("sample 重复 : ", sum(duplicated(info$sample))),
  paste0("sample 与矩阵列名 identical : ", !is.null(mat_cols)),
  paste0("需改写列名数 : ", sum(name_check$changed)),
  "",
  "-------------------------- 协变量 --------------------------",
  "性别（group x gender）：",
  paste0("  ", capture.output(print(xt_gender))),
  "心衰状态（group x heart_failure）：",
  paste0("  ", capture.output(print(xt_hf))),
  "年龄 中位(范围)：",
  paste0("  ", names(age_stat), ": ", as.character(age_stat)),
  "",
  "-------------------------- 规则 ----------------------------",
  "group 判定优先级 : disease status:ch1 > source_name_ch1 > title > heart failure:ch1",
  "样本名匹配规则   : R0 精确 > R1 去空白 > R2 去引号 > R3 去X前缀 > R4 统一大写 > R5 去非字母数字",
  "================================================================"
)
qc_path <- file.path(clean_dir, paste0(accession, "_sample_info_qc.txt"))
write_utf8(qc, qc_path)
msg("    质控汇总：", basename(qc_path))

si_path <- file.path(rep_dir, "sessionInfo_02b_sample_info.txt")
write_utf8(capture.output(sessionInfo()), si_path)

# 6.5 详细统计落盘，供报告引用
age_tab <- do.call(rbind, lapply(split(info$age, info$group), function(v) data.frame(
  n = length(v), mean = round(mean(v, na.rm = TRUE), 1), median = median(v, na.rm = TRUE),
  sd = round(sd(v, na.rm = TRUE), 1), min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))))
age_tab <- data.frame(group = rownames(age_tab), age_tab, row.names = NULL)
age_lines <- c(paste(colnames(age_tab), collapse = "\t"),
               apply(age_tab, 1, function(r) paste(r, collapse = "\t")))

detail <- c(
  "=== GROUP_COUNTS ===",
  "group\tn\tpercent\tgroup_full",
  sprintf("%s\t%d\t%.2f\t%s", cnt$group, cnt$n, cnt$percent, cnt$group_full),
  "",
  "=== GROUP_SOURCE ===",
  paste0(names(table(info$group_source)), "\t", as.integer(table(info$group_source))),
  "",
  "=== CONFLICT ===",
  paste0("conflict=", sum(info$group_conflict),
         "; hf_inconsistent=", sum(!info$hf_consistent, na.rm = TRUE)),
  "",
  "=== GENDER ===",
  capture.output(print(xt_gender)),
  "",
  "=== HEART_FAILURE ===",
  capture.output(print(xt_hf)),
  "",
  "=== AGE ===",
  paste0(names(age_stat), "\t", as.character(age_stat)),
  "",
  "=== AGE_STAT ===",
  age_lines,
  "",
  "=== NAME_RULE ===",
  capture.output(print(table(name_check$rule, name_check$changed))),
  "",
  "=== LOG ===",
  LOG,
  "",
  "=== HEAD ===",
  capture.output(print(utils::head(info[, c("sample", "group", "title", "patient_id",
                                            "heart_failure", "gender", "age")], 5)))
)
write_utf8(detail, file.path(clean_dir, paste0(accession, "_sample_info_detail.txt")))

# ---------- 控制台摘要 ----------
t_end <- Sys.time()
cat("\n", strrep("=", 70), "\n", sep = "")
cat("样本信息整理完成\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("%-18s %s\n", "样本数", nrow(info)))
cat(sprintf("%-18s %s\n", "分组", paste(cnt$group, cnt$n, sep = "=", collapse = ", ")))
cat(sprintf("%-18s %s\n", "字段冲突", sum(info$group_conflict)))
cat(sprintf("%-18s %s\n", "列名需改写", sum(name_check$changed)))
cat(sprintf("%-18s %s\n", "耗时", paste0(round(as.numeric(difftime(t_end, t_start, units = "secs")), 1), " 秒")))
cat("\n输出：")
cat(paste0("\n  - ", basename(info_path)), "\n")
cat(paste0("  - ", basename(cnt_path)), "\n")
cat(paste0("  - ", basename(chk_path)), "\n")
cat(paste0("  - ", basename(qc_path)), "\n")
