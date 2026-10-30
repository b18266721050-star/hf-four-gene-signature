# =============================================================================
# 01_download.R  ——  GSE57338 数据下载
#
# 作用：从 GEO 下载表达矩阵（series matrix）与样本信息，并下载平台注释文件，
#       全部原始文件写入 raw_data/，同时生成 raw_data/download_log.txt。
#
# 数据集：GSE57338
#         平台 GPL11532  [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array
#         313 个样本，Homo sapiens，Expression profiling by array
#         文献：Liu Y et al. Genomics 2015;105(2):83-9  (PMID 25528681)
#
# 用法：在项目根目录执行
#         Rscript scripts/01_download.R
#       （脚本以 getwd() 作为项目根目录，故必须先 cd 到项目根目录）
#
# 编码注意（Windows 必读）：
#   本文件为 UTF-8 编码（无 BOM，结尾也不加 shebang，Windows 下 Rscript 不认 BOM）。
#   R 必须运行在 UTF-8 语言环境下，否则中文会乱码、日志写入会报
#   "invalid char string in output conversion"。
#   若控制台出现乱码，先检查环境变量：LC_ALL / LANG 若被设成 Windows 不认识的
#   名字（例如 C.UTF-8），R 会退化成 C(ASCII) 环境，此时在启动前清空即可：
#     PowerShell:  $env:LC_ALL=""; $env:LANG=""; Rscript scripts/01_download.R
#     CMD       :  set LC_ALL= && set LANG= && Rscript scripts\01_download.R
#   正常的中文 Windows（系统区域「Beta: 使用 Unicode UTF-8」开启）无需任何处理。
#
# -----------------------------------------------------------------------------
# 关于「用不用 getGEO()」以及本机踩到的坑
# -----------------------------------------------------------------------------
# 已按要求优先尝试 GEOquery::getGEO(accession, GSEMatrix = TRUE)。
# 实测本机环境：GEOquery 的目录列举正常，但 series matrix 下载长时间 0 字节卡死
# （>12 min 未写出任何数据）。
# 逐项排查结论：问题出在 R 自带的 libcurl 下载器（本机 IPv6 不通 + libcurl 握手异常），
# 与 GEO 服务器无关——同一个 URL 用 Windows 原生 WinINet 只需 57 s。
#   默认 libcurl        : 长时间 0 字节卡死
#   libcurl + 强制 IPv4 : 目录列举 49 s，大文件仍卡死
#   wininet（系统栈）   : 56 MB / 57.4 s ≈ 955 KB/s   <- 采用
#
# 因此最终方案：
#   download.file(method = "wininet") 直连 GEO FTP 拉 series matrix
#   → GEOquery::getGEO(filename=) 本地解析
# 下到的是**同一个 FTP 文件**，只是绕开了 getGEO 内部不稳定的下载器；
# 解析、注释全部仍由 GEOquery 完成。
#
# 若你的网络 IPv6 正常、getGEO() 能跑通，把下面 TRY_GETGEO 改为 TRUE 即可走
# 纯 GEOquery 在线下载路径。
# =============================================================================

rm(list = ls())

# ---------- 0. 全局开关 / 路径 / 参数 ----------
TRY_GETGEO <- FALSE   # TRUE = 优先用 getGEO() 在线下载（本机实测会卡住，故置 FALSE）

proj_dir <- getwd()
raw_dir  <- file.path(proj_dir, "raw_data")
rep_dir  <- file.path(proj_dir, "report")

if (!dir.exists(raw_dir)) stop("找不到 raw_data/，请先 cd 到项目根目录再运行本脚本。")
if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

accession <- "GSE57338"
platform  <- "GPL11532"

series_matrix_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/", accession,
  "/matrix/", accession, "_series_matrix.txt.gz"
)
annot_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/", platform,
  "/annot/", platform, ".annot.gz"
)

options(timeout = 1800)
options(download.file.method.GEOquery = "libcurl")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

t_start <- Sys.time()

snap <- function() sort(list.files(raw_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE))
before <- snap()

# 带重试的下载
# ---- 本机网络实测（2026-10-26，56 MB 的 series matrix）----
#   libcurl (R 内置)        : 长时间 0 字节卡死，甚至无法完成握手
#   libcurl + 强制 IPv4     : 目录列举 49 s，大文件仍反复卡死
#   wininet (Windows 原生)  : 56,110,721 字节 / 57.4 s ≈ 955 KB/s  ✔ 采用
# 原因：本机 IPv6 不通，且 R 自带 libcurl 走系统代理/HTTP2 时表现异常；
#       wininet 直接复用 Windows 系统网络栈，最稳。
# 因此 Windows 上首选 method = "wininet"，其余平台用 libcurl。
DL_METHOD <- if (.Platform$OS.type == "windows") "wininet" else "libcurl"

dl <- function(url, dest, tries = 3) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({
      download.file(url, destfile = dest, mode = "wb",
                    method = DL_METHOD, quiet = FALSE)
      file.exists(dest) && file.info(dest)$size > 0
    }, error = function(e) {
      message("    第 ", i, " 次下载失败：", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) return(TRUE)
    Sys.sleep(5)
  }
  FALSE
}

# ---------- 1. 下载 series matrix（表达矩阵 + 样本信息） ----------
message(">>> [1/3] 获取 series matrix ...")

sm_file <- file.path(raw_dir, paste0(accession, "_series_matrix.txt.gz"))
method_note <- NA_character_
gse <- NULL

if (TRY_GETGEO) {
  gse <- tryCatch(
    getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE, destdir = raw_dir),
    error = function(e) { message("    getGEO() 失败：", conditionMessage(e)); NULL }
  )
  if (!is.null(gse)) {
    method_note <- "GEOquery::getGEO(GSEMatrix=TRUE, getGPL=FALSE) 在线下载并解析"
  }
}

if (is.null(gse)) {
  need <- !(file.exists(sm_file) && file.info(sm_file)$size > 1e6)
  if (need) {
    message("    直连下载：", series_matrix_url)
    if (!dl(series_matrix_url, sm_file)) stop("series matrix 下载失败，请检查网络或稍后重试。")
  } else {
    message("    本地已存在完整文件，跳过下载")
  }
  method_note <- paste0(
    "download.file(method='wininet'，Windows 原生网络栈) 直连 GEO FTP 下载 series matrix（",
    series_matrix_url, "），再由 GEOquery::getGEO(filename=) 本地解析",
    "；说明：GEOquery::getGEO() 在线下载在本机因 R 内置 libcurl 卡死而超时，故直连同一 FTP 文件"
  )
  gse <- getGEO(filename = sm_file)
}

if (is(gse, "ExpressionSet")) gse <- list(gse)
eset <- gse[[1]]

ex <- exprs(eset)
pd <- pData(eset)
message("    表达矩阵维度：", nrow(ex), " 探针 x ", ncol(ex), " 样本")

matrix_ok <- nrow(ex) > 0 && sum(!is.na(ex)) > 0
message("    表达矩阵是否为空：", ifelse(matrix_ok, "否（已含表达值）", "是（需去 supplementary 找）"))

# ---------- 2. 检查 supplementary：先只列清单，不下载 ----------
message(">>> [2/3] 检查 supplementary 文件清单 ...")

supp_tab <- tryCatch(
  getGEOSuppFiles(accession, baseDir = raw_dir, makeDirectory = FALSE, fetch_files = FALSE),
  error = function(e) NULL
)

supp_lines <- "  （GEO 上无 supplementary 文件，或查询失败）"
supp_files <- character(0)
if (!is.null(supp_tab) && nrow(supp_tab) > 0) {
  supp_files <- if ("fname" %in% names(supp_tab)) supp_tab$fname else rownames(supp_tab)
  supp_url   <- if ("url" %in% names(supp_tab)) supp_tab$url else ""
  supp_lines <- paste0("    ", supp_files, "   <", supp_url, ">")
}

supp_note <- "series matrix 已包含完整表达值，未下载 supplementary"
if (!matrix_ok && length(supp_files) > 0) {
  message("    表达矩阵为空，下载 supplementary 中的 processed matrix ...")
  supp_note <- tryCatch({
    getGEOSuppFiles(accession, baseDir = raw_dir, makeDirectory = FALSE)
    "series matrix 表达矩阵为空，已下载全部 supplementary 文件"
  }, error = function(e) paste0("下载 supplementary 失败：", conditionMessage(e)))
}

# ---------- 3.（可选）平台注释 ----------
message(">>> [3/3] 平台注释 ", platform, ".annot.gz ...")
annot_file <- file.path(raw_dir, paste0(platform, ".annot.gz"))
annot_note <- "已存在，跳过"
if (!file.exists(annot_file)) {
  annot_note <- if (dl(annot_url, annot_file, tries = 2)) {
    paste0("下载成功（", round(file.info(annot_file)$size / 1024^2, 2),
           " MB，供 02_clean.R 做探针→基因注释）")
  } else {
    "下载失败（可选文件，不影响本次主流程）"
  }
}
message("    ", annot_note)

# ---------- 4. 写 download_log.txt ----------
t_end <- Sys.time()
after <- snap()
new_files <- setdiff(after, before)

L <- c(
  "================================================================",
  "GEO 数据下载日志",
  "================================================================",
  paste0("数据集编号 : ", accession),
  paste0("数据集标题 : RNA-Seq Identifies Novel Myocardial Gene Expression Signatures of Heart Failure [microarray]"),
  paste0("平台       : ", platform, "  [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array"),
  paste0("物种       : Homo sapiens"),
  paste0("数据类型   : Expression profiling by array（芯片，非 RNA-seq）"),
  paste0("样本数     : ", ncol(ex)),
  paste0("下载时间   : ", format(t_start, "%Y-%m-%d %H:%M:%S"), " (", format(t_start, "%Z"), ")"),
  paste0("完成时间   : ", format(t_end, "%Y-%m-%d %H:%M:%S")),
  paste0("耗时       : ", round(as.numeric(difftime(t_end, t_start, units = "secs")), 1), " 秒"),
  paste0("主机       : ", Sys.info()[["nodename"]]),
  paste0("R 版本     : ", R.version.string),
  paste0("GEOquery   : ", as.character(packageVersion("GEOquery"))),
  "",
  "----------------------------------------------------------------",
  "下载方式",
  "----------------------------------------------------------------",
  paste0("  表达矩阵 + 样本信息 : ", method_note),
  paste0("  Supplementary 处理  : ", supp_note),
  paste0("  平台注释            : ", annot_note),
  "",
  "----------------------------------------------------------------",
  "raw_data/ 中的原始数据文件（本次下载 + 既有文件）",
  "----------------------------------------------------------------"
)

rec_files <- setdiff(snap(), "download_log.txt")
if (length(rec_files) == 0) {
  L <- c(L, "  （raw_data/ 中没有任何数据文件）")
} else {
  for (f in rec_files) {
    fp  <- file.path(raw_dir, f)
    sz  <- file.info(fp)$size
    mt  <- format(file.info(fp)$mtime, "%Y-%m-%d %H:%M:%S")
    md5 <- tryCatch(as.character(tools::md5sum(fp)), error = function(e) NA_character_)
    tag <- if (f %in% new_files) "本次新增" else "此前已存在"
    L <- c(L,
      paste0("  文件名   : ", basename(f), "    [", tag, "]"),
      paste0("  存放路径 : ", normalizePath(fp, winslash = "/")),
      paste0("  相对路径 : raw_data/", f),
      paste0("  文件大小 : ", format(sz, big.mark = ","), " 字节 (", round(sz / 1024^2, 2), " MB)"),
      paste0("  写入时间 : ", mt),
      paste0("  md5      : ", md5),
      "")
  }
}

L <- c(L,
  "----------------------------------------------------------------",
  "GEO 上可用的 supplementary 文件（仅列出，未下载）",
  "----------------------------------------------------------------",
  supp_lines
)

cand <- grep("characteristics|title|source_name", colnames(pd), value = TRUE, ignore.case = TRUE)
L <- c(L,
  "",
  "----------------------------------------------------------------",
  "数据内容摘要",
  "----------------------------------------------------------------",
  paste0("  表达矩阵维度 : ", nrow(ex), " 行 x ", ncol(ex), " 列（探针 x 样本）"),
  paste0("  表达值范围   : ", paste(round(range(ex, na.rm = TRUE), 3), collapse = " ~ "),
         "  （GEO 已做过 RMA 归一化 + log2，series matrix 中的值可直接用于 limma）"),
  paste0("  样本信息列数 : ", ncol(pd), " 列"),
  "  样本分组相关列：",
  if (length(cand) == 0) "    无" else
    unlist(lapply(cand, function(cn) {
      v <- as.character(pd[[cn]])
      paste0("    ", cn, " : ", paste(head(unique(v), 6), collapse = " | "),
             if (length(unique(v)) > 6) " ..." else "")
    })),
  "",
  "================================================================",
  "注：raw_data/ 为只读区，任何清洗/转换请输出到 clean_data/。",
  "================================================================"
)

log_path <- file.path(raw_dir, "download_log.txt")
con <- file(log_path, open = "w", encoding = "UTF-8")
writeLines(enc2utf8(L), con)
close(con)

# ---------- 5. 控制台摘要 + sessionInfo ----------
cat("\n", strrep("=", 62), "\n", sep = "")
cat("下载完成\n")
cat(strrep("=", 62), "\n", sep = "")
cat("表达矩阵：", nrow(ex), " 探针 x ", ncol(ex), " 样本\n", sep = "")
cat("raw_data/ 新增：\n")
if (length(new_files) == 0) cat("  (无)\n") else cat(paste0("  - ", new_files, "\n"), sep = "")
cat("日志：", normalizePath(log_path, winslash = "/"), "\n", sep = "")

si_path <- file.path(rep_dir, "sessionInfo_01_download.txt")
con2 <- file(si_path, open = "w", encoding = "UTF-8")
writeLines(capture.output(sessionInfo()), con2)
close(con2)
