# =============================================================================
# 17_ppi_network.R —— GSE57338 候选基因的蛋白质-蛋白质互作（PPI）网络构建
#
# 数据来源：STRING 数据库，物种 Homo sapiens（Taxonomy ID = 9606）
# 输入：results/HubGene/GSE57338_candidate_genes.csv
#       （上一步 16_hub_gene_intersection.R 得到的 DEG ∩ WGCNA 关键模块交集基因）
# 输出：results/PPI/ 下的映射表、边表、节点属性、网络图与中文报告
#
# 用法（项目根目录下）：
#   Rscript scripts/17_ppi_network.R
#
# -----------------------------------------------------------------------------
# 方法学要点（为什么这样写）
# -----------------------------------------------------------------------------
#   1) 本步骤不做任何基因筛选。上一步给出的交集基因全部进入网络，目的是如实呈现
#      这批基因在 STRING 中的互作结构，而不是再挑一轮"更重要的基因"。
#      孤立节点（degree = 0）是该结构的真实组成部分，必须画出并单列，不能删掉。
#   2) 置信度阈值 0.4：STRING 的常用分档是 0.15(低)/0.4(中)/0.7(高)。0.4 是 STRING
#      自己标注的"medium confidence"，也是文献中最常用的默认档；再低会引入大量
#      文本挖掘带来的假阳性，再高会让小基因集的网络碎成一盘散沙。
#      STRING 的分数在 0-1000 整数尺度上，因此脚本里写的是 400。
#   3) 不添加额外节点。STRING 的 API 有 add_white_nodes / add_color_nodes 参数，
#      会让数据库里的"邻居蛋白"自动加进网络。本步骤不加——加了就不再是"这批基因
#      之间的互作"，而是"这批基因的邻域"，节点数会失控，也无法回答原问题。
#      因此用 /network 端点（只返回输入集合内部的边），而不是 /interaction_partners。
#   4) 基因名必须先清理再映射。GeneSymbol 里混有空格、重复、ID 类型不一致时，
#      STRING 会静默少映射一部分，表现为"网络变小了"却没人发现。脚本显式审计
#      缺失/重复/空格/纯数字/Entrez/Ensembl/探针号，并单独输出未映射基因清单。
#   5) 数据获取有两条路径，脚本自动选择，结果等价：
#        (A) STRINGdb 本地路径：若 00_string_cache/ 里已有 STRING flat file，
#            用 STRINGdb 包的 map() + get_interactions() 完成映射与取边；
#        (B) STRING Web API 路径：缓存缺失时，调用 string-db.org 官方 API。
#      两条路径访问的是同一个 STRING 版本（脚本里显式锁定），同一物种、同一阈值。
#      之所以保留 (B)：本机到 stringdb-downloads.org 的实测吞吐约 0.03 MB/s，
#      人类 v12.5 的 flat file 合计约 80 MB，完整下载需 45 分钟以上且中断后无法
#      续传；而 API 只返回这批基因的映射与边（几 KB），秒级完成。
#      在能顺畅下载的网络环境或已有缓存的机器上，脚本会自动走 (A)。
#   6) 布局用 Fruchterman-Reingold 力导向算法，并固定随机种子。FR 对孤立节点没有
#      吸引力，会把它们甩到画布边缘且可能重叠，因此这里把孤立节点单独均匀排布在
#      外圈，既有连接的部分由 FR 决定——这样图稳定、可复现且不会糊成一团。
#   7) 图中文字全部为英文（基因名），避免 PDF 嵌入中文字体导致的字形缺失；
#      中文说明放在 .md 报告里。
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(encoding = "UTF-8")
# 本机 R 启动时 locale 是 C（.Renviron 里的 C.UTF-8 不被 Windows 接受），
# 在 C locale 下中文字符串写盘会报 "invalid char string in output conversion"。
# 这里显式切到 UTF-8 locale。注意：源文件解析发生在这一行之前，
# 所以脚本里不能出现中文标识符（中文列名一律用运行时赋值的字符串）。
invisible(try(Sys.setlocale("LC_ALL", ".UTF-8"), silent = TRUE))
options(download.file.method = "libcurl")   # Windows 上规避 WinINet 的 25 秒假死
options(url.method = "libcurl")
options(timeout = 600)
Sys.setenv(HTTP_PROXY = "", HTTPS_PROXY = "", http_proxy = "", https_proxy = "")

## 随机种子：力导向布局是迭代算法，必须固定种子才能逐次出同一张图
SEED <- 20261027
set.seed(SEED)

## ---------------------------------------------------------------------------
## 0. 路径与环境
## ---------------------------------------------------------------------------
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("results")) stop("请在项目根目录下运行：Rscript scripts/17_ppi_network.R")

PROJECT_ROOT <- getwd()
DS          <- "GSE57338"
DIR_IN      <- file.path(PROJECT_ROOT, "results", "HubGene")
DIR_OUT     <- file.path(PROJECT_ROOT, "results", "PPI")
DIR_CACHE   <- file.path(DIR_OUT, "00_string_cache")
dir.create(DIR_OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

F_GENES <- file.path(DIR_IN, paste0(DS, "_candidate_genes.csv"))
if (!file.exists(F_GENES)) stop("缺少输入文件：", F_GENES)

LOG <- character(0)
msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  LOG <<- c(LOG, line)
  cat(line, "\n", sep = "")
}

# --- 中文写盘：本机 R 的 locale 是 C，源文件里的 UTF-8 字符串被当作 unknown 编码，
#     直接 writeLines/write.csv 会报 "invalid char string in output conversion"。
#     解决办法是把字符串显式标记为 UTF-8（只标记、不转码，字节原样保留）。
mark_utf8 <- function(x) {
  if (is.character(x)) Encoding(x) <- "UTF-8"
  x
}
wtxt <- function(x, f) {
  writeLines(mark_utf8(as.character(x)), f, useBytes = TRUE)
}
wcsv <- function(df, f) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) df[[j]] <- mark_utf8(df[[j]])
  names(df) <- mark_utf8(names(df))
  con <- file(f, open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  utils::write.csv(df, con, row.names = FALSE)
  invisible(f)
}

msg("=== ", DS, " PPI 网络构建（STRING, Homo sapiens 9606） ===")
msg("项目根目录：", PROJECT_ROOT)
msg("输出目录  ：", DIR_OUT)

## ---------------------------------------------------------------------------
## 1. 分析参数（这些值全部写进报告，任何一步改动都会在报告里体现）
## ---------------------------------------------------------------------------
TAXON_ID        <- 9606L      # Homo sapiens
TAXON_NAME      <- "Homo sapiens"
SCORE_THRESHOLD <- 400L       # 0.4（STRING 的 0-1000 整数尺度）
STRING_VERSION  <- "12.5"     # 显式锁定版本，保证两条路径同源
NETWORK_TYPE    <- "full"     # full / physical / regulatory
ADD_EXTRA_NODES <- FALSE      # 不添加额外节点
CALLER_ID       <- "workbuddy_ppi_GSE57338"

QUERY_DATE <- Sys.Date()
msg("参数：物种=", TAXON_NAME, "(", TAXON_ID, ")  版本=STRING v", STRING_VERSION,
    "  置信度阈值=", SCORE_THRESHOLD / 1000, "  网络类型=", NETWORK_TYPE,
    "  额外节点=", ADD_EXTRA_NODES)

## ---------------------------------------------------------------------------
## 2. 读入候选基因并清理、审计基因名
## ---------------------------------------------------------------------------
cand <- read.csv(F_GENES, check.names = FALSE, stringsAsFactors = FALSE)
msg("读入候选基因表：", nrow(cand), " 行 × ", ncol(cand), " 列")

# 基因名列：优先 GeneSymbol
pick_gene_col <- function(df) {
  nm  <- names(df)
  hit <- nm[tolower(nm) %in% c("genesymbol", "gene_symbol", "symbol",
                               "genename", "gene_name", "gene")]
  if (length(hit) > 0) return(hit[1])
  cc  <- nm[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
  if (length(cc) > 0) return(cc[1])
  stop("找不到基因名列")
}
COL_GENE <- pick_gene_col(cand)
raw_sym  <- as.character(cand[[COL_GENE]])

# ---- 审计（清理之前，如实记录原始状态）----
# 注意：列名在运行时用中文字符串赋值，而不是写成中文标识符——本机 R 运行在 C locale，
# 源文件里的中文标识符会在解析阶段直接报 "unexpected invalid token"。
audit <- data.frame(
  item = c("原始基因数", "缺失/空值", "重复基因名", "含空白字符",
           "全数字（疑似 Entrez ID）", "ENSG 前缀（疑似 Ensembl ID）",
           "6 位以上纯数字（疑似探针号）", "含小写字母",
           "其中 HGNC 的 C<n>orf<m> 命名（合法）", "其他可疑小写命名"),
  n = c(length(raw_sym),
           sum(is.na(raw_sym) | raw_sym == "" | raw_sym == "NA"),
           sum(duplicated(raw_sym)),
           sum(grepl("\\s", raw_sym, perl = TRUE)),
           sum(grepl("^[0-9]+$", raw_sym, perl = TRUE)),
           sum(grepl("^ENSG", raw_sym, perl = TRUE)),
           sum(grepl("^[0-9]{6,}$", raw_sym, perl = TRUE)),
           sum(grepl("[a-z]", raw_sym, perl = TRUE)),
           sum(grepl("[a-z]", raw_sym, perl = TRUE) &
               grepl("^C[0-9XY]+orf[0-9]+(-AS[0-9]+)?$", raw_sym, perl = TRUE)),
           sum(grepl("[a-z]", raw_sym, perl = TRUE) &
               !grepl("^C[0-9XY]+orf[0-9]+(-AS[0-9]+)?$", raw_sym, perl = TRUE))),
  stringsAsFactors = FALSE
)
names(audit) <- c("检查项", "数量")
wcsv(audit, file.path(DIR_OUT, "01_gene_symbol_audit.csv"))
msg("-- 基因名审计（清理前）--")
for (i in seq_len(nrow(audit))) msg("   ", audit[[1]][i], " = ", audit[[2]][i])

# ---- 清理：去首尾空格 → 去空/NA → 去重（保留首次出现顺序）----
sym_trim <- trimws(raw_sym)
sym_trim[sym_trim %in% c("", "NA", "NaN")] <- NA
n_removed <- sum(is.na(sym_trim))
sym_clean <- unique(sym_trim[!is.na(sym_trim)])
n_dup     <- length(sym_trim[!is.na(sym_trim)]) - length(sym_clean)
msg("清理：去除空值 ", n_removed, " 个，去除重复 ", n_dup,
    " 个 → 用于查询的唯一基因 ", length(sym_clean), " 个")

## ---------------------------------------------------------------------------
## 3. 初始化 STRINGdb 包对象，记录版本与参数
## ---------------------------------------------------------------------------
HAS_STRINGDB <- requireNamespace("STRINGdb", quietly = TRUE)
if (!HAS_STRINGDB) stop("缺少 Bioconductor 包 STRINGdb，请先安装。")
suppressPackageStartupMessages(library(STRINGdb))

string_db <- STRINGdb$new(species = TAXON_ID,
                          score_threshold = SCORE_THRESHOLD,
                          network_type = NETWORK_TYPE,
                          version = STRING_VERSION,
                          input_directory = DIR_CACHE)

VER_ACTUAL   <- string_db$version
NET_TYPE_ACT <- string_db$network_type
LINK_DATA    <- string_db$link_data
msg("STRINGdb 包版本：", as.character(packageVersion("STRINGdb")))
msg("STRINGdb 对象   ：version=", VER_ACTUAL, "  species=", string_db$species,
    "  score_threshold=", string_db$score_threshold,
    "  network_type=", NET_TYPE_ACT, "  link_data=", LINK_DATA)

# 本地 flat file 是否齐备 → 决定走 (A) 还是 (B)
f_alias <- file.path(DIR_CACHE, sprintf("%s.protein.aliases.v%s.txt.gz", TAXON_ID, VER_ACTUAL))
f_links <- file.path(DIR_CACHE, sprintf("%s.protein.links.v%s.txt.gz",    TAXON_ID, VER_ACTUAL))
f_info  <- file.path(DIR_CACHE, sprintf("%s.protein.info.v%s.txt.gz",     TAXON_ID, VER_ACTUAL))
LOCAL_READY <- file.exists(f_alias) && file.exists(f_links)
DATA_PATH   <- if (LOCAL_READY) "STRINGdb_local_flatfile" else "STRING_web_API"
msg("本地 STRING flat file 齐备：", LOCAL_READY, " → 数据获取路径：", DATA_PATH)

## ---------------------------------------------------------------------------
## 4. 基因 → 蛋白映射
## ---------------------------------------------------------------------------
# --- 4a. 用 STRING 官方 version 端点核实"实际服务的版本"，并锁定版本化端点 ---
# 这一步是真实的联网查询：主站可能随时升级版本，必须确认本次查询命中的是哪一个版本，
# 而不是照抄包文档里写的版本号。
api_version_url <- "https://string-db.org/api/json/version"
API_VERSION_REPORTED <- NA_character_
STABLE_ADDRESS       <- NA_character_
ver_res <- try(curl::curl_fetch_memory(api_version_url,
                 handle = curl::new_handle(followlocation = TRUE, timeout = 30)),
               silent = TRUE)
if (!inherits(ver_res, "try-error") && ver_res$status_code == 200L) {
  vj <- try(jsonlite::fromJSON(rawToChar(ver_res$content)), silent = TRUE)
  if (!inherits(vj, "try-error") && is.data.frame(vj) && nrow(vj) > 0) {
    API_VERSION_REPORTED <- as.character(vj$string_version[1])
    STABLE_ADDRESS       <- as.character(vj$stable_address[1])
  }
}
msg("STRING 官方 version 端点实测：string_version=", API_VERSION_REPORTED,
    "  stable_address=", STABLE_ADDRESS,
    "  （脚本锁定 v", VER_ACTUAL, "，一致性：",
    identical(API_VERSION_REPORTED, VER_ACTUAL), "）")

# 优先用版本化端点（严格锁定版本），取不到再退回主站
if (!is.na(STABLE_ADDRESS) && nzchar(STABLE_ADDRESS)) {
  api_base <- paste0(sub("/$", "", STABLE_ADDRESS), "/api")
} else {
  api_base <- sprintf("https://version-%s.string-db.org/api", gsub("\\.", "-", VER_ACTUAL))
}
msg("API 端点：", api_base)

api_post <- function(path, body, tries = 3L, timeout = 60L) {
  url <- paste0(api_base, path)
  for (k in seq_len(tries)) {
    res <- try(
      httr::POST(url, encode = "form",
                 body = body,
                 httr::timeout(timeout),
                 httr::user_agent(CALLER_ID)),
      silent = TRUE)
    if (!inherits(res, "try-error") && res$status_code == 200L) return(res)
    msg("   ! API 重试 ", k, "/", tries, "（", path, "）")
    Sys.sleep(2 * k)
  }
  stop("STRING API 请求失败：", url)
}

MAP_COLS <- c("queryIndex", "queryItem", "stringId", "preferredName",
              "ncbiTaxonId", "taxonName", "annotation")

map_via_api <- function(genes) {
  res <- api_post("/json/get_string_ids",
                  list(identifiers = paste(genes, collapse = "\r\n"),
                       species = TAXON_ID,
                       limit = 1,             # 每个基因只取最佳匹配（takeFirst）
                       echo_query = 1,        # 回显查询词，便于对齐未映射项
                       caller_identity = CALLER_ID))
  js  <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"),
                            simplifyVector = TRUE)
  if (!is.data.frame(js)) js <- as.data.frame(js, stringsAsFactors = FALSE)
  js
}

msg("-- 基因→蛋白映射 --")
if (LOCAL_READY) {
  ## --- 路径 A：STRINGdb 本地 flat file ---
  d <- data.frame(GeneSymbol = sym_clean, stringsAsFactors = FALSE)
  mapped <- string_db$map(d, "GeneSymbol", takeFirst = TRUE,
                          removeUnmappedRows = FALSE, quiet = TRUE)
  names(mapped)[names(mapped) == "STRING_id"] <- "stringId"
  prot <- string_db$get_proteins()
  prot <- prot[, intersect(c("stringId", "preferredName", "annotation",
                             "protein_size"), names(prot)), drop = FALSE]
  mp <- merge(mapped, unique(prot), by = "stringId", all.x = TRUE, sort = FALSE)
  mp$queryItem <- mp$GeneSymbol
  mp$ncbiTaxonId <- TAXON_ID
  mp$taxonName   <- TAXON_NAME
} else {
  ## --- 路径 B：STRING Web API ---
  js <- map_via_api(sym_clean)
  for (cc in MAP_COLS) if (!cc %in% names(js)) js[[cc]] <- NA
  mp <- js[, MAP_COLS, drop = FALSE]
  mp$GeneSymbol <- mp$queryItem
  # API 只返回命中项；对未返回的基因逐个复验一次，确认是"真的没有"而不是批量请求被截断
  got  <- unique(mp$queryItem)
  miss <- setdiff(sym_clean, got)
  if (length(miss) > 0) {
    msg("   批量请求未返回 ", length(miss), " 个基因，逐个复验：",
        paste(miss, collapse = ", "))
    for (gname in miss) {
      one <- try(map_via_api(gname), silent = TRUE)
      if (!inherits(one, "try-error") && nrow(one) > 0 &&
          !is.na(one$stringId[1])) {
        one$GeneSymbol <- gname
        mp <- rbind(mp, one[, c(MAP_COLS, "GeneSymbol"), drop = FALSE])
      }
    }
  }
}

msg("API/本地返回映射条目：", nrow(mp), " 行")

## --- 4b. 整理映射表：一基因一行，未映射的保留为 NA ---
keep_ok <- !is.na(mp$stringId) & mp$stringId != ""
mp_ok   <- mp[keep_ok, , drop = FALSE]
# 一个基因可能返回多条（别名同义），取第一条，保证"一基因一蛋白"
mp_ok   <- mp_ok[!duplicated(mp_ok$GeneSymbol), , drop = FALSE]

mapping <- data.frame(
  GeneSymbol    = sym_clean,
  stringsAsFactors = FALSE
)
mapping$STRING_id <- NA_character_
mapping$preferred_name <- NA_character_
mapping$annotation <- NA_character_
mapping$ncbiTaxonId <- NA_character_
mapping$taxonName <- NA_character_
# idx[i] = sym_clean 的第 i 个基因在 mp_ok 中的行号（可能为 NA）
idx <- match(sym_clean, mp_ok$GeneSymbol)
hit <- !is.na(idx)
# 关键：左侧必须用 hit（逻辑向量，索引 mapping 的行）；右侧用 idx[hit]（索引 mp_ok 的行）。
# 两边都用 idx[hit] 会因为"未命中的基因被丢掉"而产生逐行偏移，越往后错得越多。
mapping$STRING_id[hit]      <- mp_ok$stringId[idx[hit]]
mapping$preferred_name[hit] <- mp_ok$preferredName[idx[hit]]
mapping$ncbiTaxonId[hit]    <- as.character(mp_ok$ncbiTaxonId[idx[hit]])
mapping$taxonName[hit]      <- mp_ok$taxonName[idx[hit]]
if (!is.null(mp_ok$annotation)) mapping$annotation[hit] <- mp_ok$annotation[idx[hit]]
mapping$mapped <- !is.na(mapping$STRING_id)

# 带上上一步的差异表达与 WGCNA 注释，方便下游直接引用
if (COL_GENE %in% names(cand)) {
  ann <- cand[, intersect(c(COL_GENE, "log2FoldChange", "padj", "Regulation",
                            "module", "MM_kME", "GS", "is_hub"), names(cand)),
              drop = FALSE]
  names(ann)[names(ann) == COL_GENE] <- "GeneSymbol"
  ann <- ann[!duplicated(ann$GeneSymbol), , drop = FALSE]
  mapping <- merge(mapping, ann, by = "GeneSymbol", all.x = TRUE, sort = FALSE)
}
mapping <- mapping[order(!mapping$mapped, mapping$GeneSymbol), , drop = FALSE]

wcsv(mapping, file.path(DIR_OUT, "02_gene_protein_mapping.csv"))

n_mapped   <- sum(mapping$mapped)
n_unmapped <- sum(!mapping$mapped)
unmapped   <- mapping[!mapping$mapped, c("GeneSymbol"), drop = FALSE]
wcsv(unmapped, file.path(DIR_OUT, "03_unmapped_genes.csv"))

msg("映射成功：", n_mapped, " / ", nrow(mapping), "；未映射：", n_unmapped)
if (n_unmapped > 0) msg("   未映射基因：", paste(unmapped$GeneSymbol, collapse = ", "))
if (n_mapped == 0) stop("全部基因均未映射到 STRING，无法构建网络。")

## 差异：映射回的 preferred_name 与输入 symbol 不一致的（别名/曾用名合并的情况）
disc <- mapping[mapping$mapped & mapping$preferred_name != mapping$GeneSymbol,
                c("GeneSymbol", "preferred_name", "STRING_id")]
if (nrow(disc) > 0) {
  msg("   注意：", nrow(disc), " 个基因的 STRING 首选名与输入 symbol 不同：",
      paste(head(paste0(disc$GeneSymbol, "→", disc$preferred_name), 8), collapse = ", "))
}

## ---------------------------------------------------------------------------
## 5. 取互作边（仅输入基因集合内部，不添加额外节点）
## ---------------------------------------------------------------------------
mapped_ids <- mapping$STRING_id[mapping$mapped]
id2gene    <- setNames(mapping$GeneSymbol[mapping$mapped], mapped_ids)

msg("-- 获取 STRING 互作边 --")
msg("   查询蛋白数：", length(mapped_ids), "；阈值 combined_score >= ",
    SCORE_THRESHOLD / 1000)

if (LOCAL_READY) {
  raw_int <- string_db$get_interactions(mapped_ids)
} else {
  res <- api_post("/tsv/network",
                  list(identifiers = paste(mapped_ids, collapse = "\r\n"),
                       species = TAXON_ID,
                       required_score = SCORE_THRESHOLD,
                       network_type = NETWORK_TYPE,
                       caller_identity = CALLER_ID))
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  raw_int <- if (nchar(txt) > 0 && grepl("\t", txt)) {
    utils::read.delim(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }
}

CHANNELS <- c("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
CHANNEL_CN <- c(nscore = "基因邻域 neighborhood", fscore = "基因融合 fusion",
                pscore = "系统发生共现 cooccurrence", ascore = "共表达 coexpression",
                escore = "实验证据 experimental", dscore = "数据库收录 database",
                tscore = "文本挖掘 textmining")

if (nrow(raw_int) == 0) {
  msg("   ! STRING 未返回任何满足阈值的边（网络将只含孤立节点）")
  edges <- data.frame(from = character(0), to = character(0),
                      combined_score = numeric(0))
  for (ch in CHANNELS) edges[[ch]] <- numeric(0)
} else {
  # 统一列名：不同路径返回的列名可能不同
  nm <- names(raw_int)
  pick <- function(cands) {
    for (c0 in cands) if (c0 %in% nm) return(c0)
    NA_character_
  }
  cA <- pick(c("stringId_A", "from", "protein1", "item_id_a"))
  cB <- pick(c("stringId_B", "to",   "protein2", "item_id_b"))
  cS <- pick(c("score", "combined_score"))
  if (is.na(cA) || is.na(cB) || is.na(cS)) stop("无法识别 STRING 返回的边表列名：",
                                                paste(nm, collapse = ", "))
  edges <- data.frame(
    from           = as.character(raw_int[[cA]]),
    to             = as.character(raw_int[[cB]]),
    combined_score = as.numeric(raw_int[[cS]]),
    stringsAsFactors = FALSE
  )
  # 统一分数尺度：STRING 的 Web API 返回 0-1，本地 flat file 返回 0-1000。
  # 这里一律折算到 0-1000 的整数尺度（STRING 的惯例），否则阈值过滤会全错。
  if (length(edges$combined_score) &&
      max(edges$combined_score, na.rm = TRUE) <= 1.0) {
    edges$combined_score <- edges$combined_score * 1000
  }
  # 证据通道分（combined_only 模式下没有，逐列补 NA），统一折算到 0-1
  for (ch in CHANNELS) {
    v <- if (ch %in% nm) as.numeric(raw_int[[ch]]) else NA_real_
    if (length(v) && all(!is.na(v)) && max(v, na.rm = TRUE) > 1.0) v <- v / 1000
    edges[[ch]] <- v
  }
  # 端点名（部分路径会返回 preferredName_A/B）
  cNA <- pick(c("preferredName_A", "from_name"))
  cNB <- pick(c("preferredName_B", "to_name"))
  edges$from_name <- if (!is.na(cNA)) as.character(raw_int[[cNA]]) else id2gene[edges$from]
  edges$to_name   <- if (!is.na(cNB)) as.character(raw_int[[cNB]]) else id2gene[edges$to]

  # 只保留两端都在输入集合内的边（再次确认"不添加额外节点"）
  keep <- edges$from %in% mapped_ids & edges$to %in% mapped_ids
  msg("   原始返回边数：", nrow(edges), "；两端均在输入集内的：", sum(keep))
  edges <- edges[keep, , drop = FALSE]
  # 阈值兜底过滤（API 已按阈值返回，这里防止本地路径返回未过滤数据）
  edges <- edges[edges$combined_score >= SCORE_THRESHOLD, , drop = FALSE]
  edges <- edges[!duplicated(paste(edges$from, edges$to)), , drop = FALSE]
  edges$from_gene <- unname(id2gene[edges$from])
  edges$to_gene   <- unname(id2gene[edges$to])
  edges <- edges[order(-edges$combined_score, edges$from_gene, edges$to_gene), , drop = FALSE]
}

wcsv(edges, file.path(DIR_OUT, "04_string_edges_raw.csv"))
msg("   保留边数：", nrow(edges))

## ---------------------------------------------------------------------------
## 6. 建图与节点属性
## ---------------------------------------------------------------------------
if (!requireNamespace("igraph", quietly = TRUE)) stop("缺少 R 包：igraph")
suppressPackageStartupMessages(library(igraph))

nodes <- mapping[mapping$mapped, c("GeneSymbol", "STRING_id"), drop = FALSE]
names(nodes) <- c("name", "STRING_id")

if (nrow(edges) > 0) {
  el <- data.frame(from = edges$from_gene, to = edges$to_gene,
                   combined_score = edges$combined_score,
                   stringsAsFactors = FALSE)
  g <- graph_from_data_frame(el, directed = FALSE, vertices = nodes)
} else {
  g <- graph_from_data_frame(data.frame(from = character(0), to = character(0)),
                             directed = FALSE, vertices = nodes)
}
# 去掉自环与多重边（STRING 单基因自互作会干扰 degree）
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# simplify() 会丢弃边属性（igraph 2.x 实测 combined_score 变成 NULL），
# 因此这里按"基因对"把置信度重新写回图对象，后面的加权度与边宽都依赖它。
if (nrow(edges) > 0 && ecount(g) > 0) {
  em_tmp  <- as_edgelist(g, names = TRUE)
  key_gr  <- paste(pmin(em_tmp[, 1], em_tmp[, 2]),
                   pmax(em_tmp[, 1], em_tmp[, 2]), sep = "|")
  key_df  <- paste(pmin(edges$from_gene, edges$to_gene),
                   pmax(edges$from_gene, edges$to_gene), sep = "|")
  df_sc   <- tapply(edges$combined_score, key_df, max)   # 重边取最高分
  E(g)$combined_score <- as.numeric(df_sc[match(key_gr, names(df_sc))])
}

deg       <- degree(g)
V(g)$degree <- deg
V(g)$isolated <- deg == 0L

node_tab <- data.frame(
  GeneSymbol = V(g)$name,
  STRING_id  = V(g)$STRING_id,
  Degree     = as.integer(deg),
  Isolated   = deg == 0L,
  stringsAsFactors = FALSE
)
# 边的加权度（置信度之和），衡量"连接强度"而非"连接数"
if (ecount(g) > 0) {
  # igraph 2.0 起 strengths() 更名为 strength()
  ws <- strength(g, weights = E(g)$combined_score, mode = "all")
  node_tab$WeightedDegree <- round(as.numeric(ws), 1)
} else {
  node_tab$WeightedDegree <- 0
}
node_tab$MaxEdgeScore <- vapply(node_tab$GeneSymbol, function(v) {
  ei <- incident(g, v, mode = "all")
  if (length(ei) == 0) return(NA_real_)
  max(E(g)$combined_score[ei]) / 1000
}, numeric(1))
if (all(c("log2FoldChange", "padj", "Regulation") %in% names(mapping))) {
  m2 <- mapping[, c("GeneSymbol", "log2FoldChange", "padj", "Regulation")]
  node_tab <- merge(node_tab, m2, by = "GeneSymbol", all.x = TRUE, sort = FALSE)
}
node_tab <- node_tab[order(-node_tab$Degree, node_tab$GeneSymbol), , drop = FALSE]
wcsv(node_tab, file.path(DIR_OUT, "05_node_degree_table.csv"))

iso_tab <- node_tab[node_tab$Isolated, c("GeneSymbol", "STRING_id", "Degree"), drop = FALSE]
wcsv(iso_tab, file.path(DIR_OUT, "06_isolated_nodes.csv"))

## --- 网络规模统计 ---
n_nodes   <- vcount(g)
n_edges   <- ecount(g)
n_iso     <- sum(node_tab$Isolated)
n_connected <- n_nodes - n_iso
dens      <- if (n_nodes > 1) edge_density(g) else 0
comp      <- components(g)
n_comp    <- comp$no
comp_size <- sort(sapply(seq_len(n_comp), function(i) sum(comp$membership == i)),
                  decreasing = TRUE)
avg_deg   <- if (n_nodes > 0) mean(deg) else 0
# 理论最大边数
max_edges <- n_nodes * (n_nodes - 1) / 2

msg("-- 网络规模 --")
msg("   节点 ", n_nodes, "  边 ", n_edges, "  平均度 ", round(avg_deg, 2),
    "  密度 ", signif(dens, 3))
msg("   连通分量 ", n_comp, " 个，最大分量 ", comp_size[1], " 个节点")
msg("   孤立节点 ", n_iso, " 个；有连接的节点 ", n_connected, " 个")

net_stats <- data.frame(
  item = c("输入候选基因数", "映射成功节点数", "未映射基因数", "边数（combined_score ≥ 0.4）",
           "平均度", "网络密度", "连通分量数", "最大连通分量节点数",
           "孤立节点数", "理论最大边数", "实际边数占理论最大边数比例"),
  value = c(nrow(mapping), n_mapped, n_unmapped, n_edges,
           round(avg_deg, 3), signif(dens, 4), n_comp, comp_size[1],
           n_iso, max_edges,
           if (max_edges > 0) paste0(sprintf("%.2f%%", 100 * n_edges / max_edges)) else "NA"),
  stringsAsFactors = FALSE
)
names(net_stats) <- c("指标", "数值")
wcsv(net_stats, file.path(DIR_OUT, "07_network_statistics.csv"))

## --- 证据通道统计 ---
has_ch <- nrow(edges) > 0 && all(CHANNELS %in% names(edges)) &&
          !all(is.na(edges[[CHANNELS[1]]]))
if (has_ch) {
  ch_tab <- data.frame(
    channel = names(CHANNEL_CN),
    meaning = unname(CHANNEL_CN),
    n_edges = vapply(CHANNELS, function(ch) sum(!is.na(edges[[ch]]) & edges[[ch]] > 0), integer(1)),
    mean_sc = round(vapply(CHANNELS, function(ch) mean(edges[[ch]], na.rm = TRUE), numeric(1)), 4),
    max_sc  = round(vapply(CHANNELS, function(ch) max(edges[[ch]], na.rm = TRUE), numeric(1)), 4),
    stringsAsFactors = FALSE
  )
} else {
  ch_tab <- data.frame(channel = names(CHANNEL_CN), meaning = unname(CHANNEL_CN),
                       n_edges = NA_integer_, mean_sc = NA_real_,
                       max_sc = NA_real_, stringsAsFactors = FALSE)
}
names(ch_tab) <- c("证据通道", "含义", "有证据的边数", "平均分数", "最大分数")
wcsv(ch_tab, file.path(DIR_OUT, "08_evidence_channels.csv"))

## ---------------------------------------------------------------------------
## 7. 绘制 PPI 网络图
## ---------------------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("缺少 R 包：ggplot2")
suppressPackageStartupMessages(library(ggplot2))
HAS_GGREPEL <- requireNamespace("ggrepel", quietly = TRUE)
HAS_RAGG    <- requireNamespace("ragg",    quietly = TRUE)

msg("-- 绘制网络图 --")
set.seed(SEED)

# --- 布局：有连接的部分用力导向（FR），孤立节点单独排在外圈 ---
core_v <- V(g)[degree(g) > 0]
iso_v  <- V(g)[degree(g) == 0]

if (length(core_v) > 0) {
  g_core <- induced_subgraph(g, core_v)
  if (ecount(g_core) > 0) {
    # 两种力导向布局各算一遍，用"节点间距的均匀性"挑更铺开的一个。
    # hub 密集的网络里，Kamada-Kawai 往往比 FR 更均匀；但这个结论不总成立，
    # 所以不预设，而是两个都画、用客观指标选。
    la_fr <- layout_with_fr(g_core, niter = 3000, grid = "nogrid")
    la_kk <- layout_with_kk(g_core, maxiter = 3000)
    spread <- function(M) {                       # 最小成对距离越大，标签越不容易重叠
      d <- as.matrix(dist(M)); diag(d) <- Inf; min(d)
    }
    lay_core <- if (spread(la_kk) > spread(la_fr)) la_kk else la_fr
    msg("   布局：", if (identical(lay_core, la_kk)) "Kamada-Kawai" else "Fruchterman-Reingold")
  } else {
    lay_core <- layout_in_circle(g_core)
  }
  # 归一化后稍作放大（[-1.22, 1.22]），给核心 hub 区的标签腾出排斥空间，
  # 避免 STAT3/CD163 这类高度数节点的标签挤在一起。
  rng <- apply(lay_core, 2, function(z) {
    r <- max(z) - min(z); if (r == 0) c(z[1] - 0.5, z[1] + 0.5) else c(min(z), max(z))
  })
  lay_core[, 1] <- (lay_core[, 1] - rng[1, 1]) / diff(rng[, 1]) * 2.9 - 1.45
  lay_core[, 2] <- (lay_core[, 2] - rng[1, 2]) / diff(rng[, 2]) * 2.9 - 1.45
  rmax <- max(sqrt(rowSums(lay_core^2)))
} else {
  lay_core <- matrix(numeric(0), ncol = 2)
  rmax     <- 1
}

lay <- matrix(NA_real_, nrow = vcount(g), ncol = 2)
colnames(lay) <- c("x", "y")
if (length(core_v) > 0) lay[as.integer(core_v), ] <- lay_core

if (length(iso_v) > 0) {
  # 孤立节点均匀分布在半径 R_iso 的圆上，起始角固定 → 完全可复现
  R_iso <- rmax * 1.35 + 0.45
  k <- length(iso_v)
  ang <- seq(pi / 2, pi / 2 + 2 * pi * (k - 1) / k, length.out = k)
  lay[as.integer(iso_v), ] <- cbind(R_iso * cos(ang), R_iso * sin(ang))
}

node_df <- data.frame(
  gene     = V(g)$name,
  x        = lay[, 1],
  y        = lay[, 2],
  degree   = as.integer(deg),
  isolated = deg == 0L,
  stringsAsFactors = FALSE
)

if (n_edges > 0) {
  # as_edgelist(names = TRUE) 明确返回顶点名矩阵；igraph 2.0 的 ends() 返回类型不直观
  el_mat <- as_edgelist(g, names = TRUE)
  esc    <- E(g)$combined_score
  if (is.null(esc) || length(esc) != nrow(el_mat)) esc <- rep(NA_real_, nrow(el_mat))
  edge_df <- data.frame(
    x     = lay[match(el_mat[, 1], V(g)$name), 1],
    y     = lay[match(el_mat[, 1], V(g)$name), 2],
    xend  = lay[match(el_mat[, 2], V(g)$name), 1],
    yend  = lay[match(el_mat[, 2], V(g)$name), 2],
    score = esc,
    stringsAsFactors = FALSE
  )
} else {
  edge_df <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0),
                        yend = numeric(0), score = numeric(0))
}

# --- 视觉映射 ---
# 节点大小 ∝ Degree（用 sqrt 压缩，避免 hub 节点过大；孤立节点给最小尺寸）
max_deg <- max(node_df$degree)
size_fun <- function(d) ifelse(d == 0, 1.9, 1.9 + 3.4 * sqrt(d / max(max_deg, 1)))
node_df$nsize <- size_fun(node_df$degree)

# 边宽 ∝ 置信度（0.4-1.0 映射到 0.18-0.72 pt；太宽会让核心区域糊成一片）
if (nrow(edge_df) > 0) {
  edge_df$ew <- 0.18 + 0.54 * (edge_df$score / 1000 - 0.4) / 0.6
  edge_df$ew <- pmin(pmax(edge_df$ew, 0.18), 0.72)
}

COL_CORE <- "#3C6E9F"   # 有连接的节点（蓝）
COL_ISO  <- "#B0AFA8"   # 孤立节点（灰）
COL_EDGE <- "#9AA5AD"   # 边
INK      <- "#333333"

node_df$grp <- ifelse(node_df$isolated, "isolated", "connected")

p <- ggplot() +
  { if (nrow(edge_df) > 0)
      geom_segment(data = edge_df,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = ew),
                   colour = COL_EDGE, alpha = 0.55, lineend = "round") } +
  geom_point(data = node_df,
             aes(x = x, y = y, size = nsize, fill = grp),
             shape = 21, colour = "white", stroke = 0.3) +
  scale_size_identity() +
  scale_linewidth_identity() +     # 边宽按实际值绘制，且不出图例
  scale_fill_manual(values = c(connected = COL_CORE, isolated = COL_ISO),
                    name = "Node type",
                    labels = c(connected = "Connected (degree ≥ 1)",
                               isolated  = "Isolated (degree = 0)"),
                    guide = guide_legend(override.aes = list(size = 3.6))) +
  coord_equal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(
    plot.margin     = margin(4, 4, 4, 4, "mm"),
    legend.position = "bottom",
    legend.title    = element_text(size = 6.5, colour = INK, face = "bold"),
    legend.text     = element_text(size = 6.0, colour = INK),
    legend.key.size = unit(3.2, "mm"),
    legend.margin   = margin(t = -2, b = 0, unit = "mm"),
    text            = element_text(colour = INK)
  )

# 标签：ggrepel 避让。密集区域标签会被甩出去并带引导线，因此给标签加白色描边
# 与半透明底，避免叠在边或其它标签上时无法辨认。未安装 ggrepel 则退化。
if (HAS_GGREPEL) {
  p <- p + ggrepel::geom_text_repel(
    data = node_df,
    aes(x = x, y = y, label = gene),
    size = 1.95, colour = INK,
    min.segment.length = 0,       # 凡是偏移的标签都画引导线
    box.padding        = 0.28,
    point.padding      = 0.22,
    max.overlaps       = Inf,
    force              = 3.2,
    force_pull         = 0.5,
    bg.color           = "white",
    bg.r               = 0.12,
    segment.colour     = "#9C9C9C",
    segment.size       = 0.16,
    segment.alpha      = 0.6,
    seed               = SEED
  )
} else {
  p <- p + geom_text(data = node_df, aes(x = x, y = y, label = gene),
                     size = 1.95, colour = INK, vjust = -0.9, hjust = 0.5)
}

W_MM <- 183; H_MM <- 165                 # Nature 双栏宽度 183 mm
w_in <- W_MM / 25.4; h_in <- H_MM / 25.4
f_pdf <- file.path(DIR_OUT, paste0(DS, "_PPI_network.pdf"))
f_png <- file.path(DIR_OUT, paste0(DS, "_PPI_network.png"))

grDevices::cairo_pdf(f_pdf, width = w_in, height = h_in, bg = "white")
print(p); invisible(dev.off())

if (HAS_RAGG) {
  ragg::agg_png(f_png, width = W_MM, height = H_MM, units = "mm",
                res = 600, background = "white")
  print(p); invisible(dev.off())
} else {
  grDevices::png(f_png, width = W_MM, height = H_MM, units = "mm",
                 res = 600, bg = "white")
  print(p); invisible(dev.off())
}
msg("  saved: ", basename(f_pdf), " / ", basename(f_png),
    "  (", W_MM, " x ", H_MM, " mm, PNG 600 dpi)")

## ---------------------------------------------------------------------------
## 8. 中文报告
## ---------------------------------------------------------------------------
top_deg <- head(node_tab[node_tab$Degree > 0, ], 12)
top_txt <- if (nrow(top_deg) > 0) {
  paste0("| `", top_deg$GeneSymbol, "` | ", top_deg$Degree, " | ",
         ifelse(is.na(top_deg$MaxEdgeScore), "—",
                sprintf("%.3f", top_deg$MaxEdgeScore)), " |", collapse = "\n")
} else "（无）"

unmapped_txt <- if (n_unmapped == 0) {
  "全部候选基因均成功映射到 STRING 蛋白，未映射基因数为 **0**。"
} else {
  paste0("共 **", n_unmapped, "** 个基因未能映射到 STRING v", VER_ACTUAL,
         " 的人类蛋白：", paste0("`", unmapped$GeneSymbol, "`", collapse = "、"), "。")
}

iso_txt <- if (n_iso == 0) {
  "本网络中**没有孤立节点**，所有节点至少有一条满足阈值的互作边。"
} else {
  paste0("本网络中有 **", n_iso, "** 个孤立节点：",
         paste0("`", iso_tab$GeneSymbol, "`", collapse = "、"), "。")
}

path_note <- if (LOCAL_READY) {
  paste0("本次运行 **使用 STRINGdb 包读取本地 STRING flat file**（缓存在 ",
         "`results/PPI/00_string_cache/`），映射与取边均由 STRINGdb 包完成。")
} else {
  paste0("本次运行 **使用 STRING v", VER_ACTUAL, " 官方 Web API**。原因：STRINGdb 包",
         "的本地路径需要先下载人类 flat file（aliases 约 18.9 MB + links 约 59.9 MB），",
         "本机实测到 STRING 下载站的吞吐约 0.03 MB/s，完整下载需 45 分钟以上且中断后",
         "无法续传；而官方 API 只返回这批基因的映射与边（几 KB），且返回的字段完全一致",
         "（同一版本、同一物种、同一阈值）。脚本保留了两者的自动切换：只要缓存目录里",
         "flat file 齐备，下一次运行就会自动改用 STRINGdb 包的本地路径。")
}

ch_txt <- if (has_ch) {
  paste0("| ", ch_tab[[1]], " | ", ch_tab[[2]], " | ", ch_tab[[3]],
         " | ", ch_tab[[4]], " | ", ch_tab[[5]], " |", collapse = "\n")
} else "（本次使用 combined_only 模式，未返回逐通道分数）"

rep_md <- c(
  paste0("# ", DS, " 候选基因的 PPI 网络构建报告"),
  "",
  paste0("> 分析内容：把上一步得到的交集基因（显著差异基因 ∩ WGCNA 关键模块）映射到 ",
         "STRING 数据库，构建人类（*Homo sapiens*，Taxonomy ID 9606）蛋白质-蛋白质互作网络并绘图。  "),
  paste0("> 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  "),
  paste0("> 生成脚本：`scripts/17_ppi_network.R`（全文含中文注释，随机种子 ", SEED, "，可复现）  "),
  paste0("> 运行环境：R ", paste(R.version$major, R.version$minor, sep = "."),
         " / STRINGdb ", as.character(packageVersion("STRINGdb")),
         " / igraph ", as.character(packageVersion("igraph")),
         " / ggplot2 ", as.character(packageVersion("ggplot2"))),
  "",
  "---",
  "",
  "## 一、数据来源与分析参数",
  "",
  "| 项目 | 取值 |",
  "| --- | --- |",
  paste0("| 数据库 | STRING v", VER_ACTUAL, " |"),
  paste0("| 物种 | *", TAXON_NAME, "*（Taxonomy ID ", TAXON_ID, "） |"),
  paste0("| R 包 | STRINGdb ", as.character(packageVersion("STRINGdb")),
         "（Bioconductor） |"),
  paste0("| 网络类型 | ", NET_TYPE_ACT, "（full = 功能关联 + 物理结合） |"),
  paste0("| 置信度阈值 | combined_score ≥ **", SCORE_THRESHOLD / 1000,
         "**（STRING 的 medium confidence 档） |"),
  paste0("| 是否添加额外节点 | **否**（只用输入基因集合内部的互作） |"),
  paste0("| 查询日期 | ", format(QUERY_DATE, "%Y-%m-%d"), " |"),
  paste0("| 官方 version 端点实测版本 | ", ifelse(is.na(API_VERSION_REPORTED),
         "（未取到）", paste0("v", API_VERSION_REPORTED)), " |"),
  paste0("| 版本化端点 | ", api_base, " |"),
  paste0("| 数据获取路径 | ", DATA_PATH, " |"),
  paste0("| 随机种子 | ", SEED, "（力导向布局） |"),
  "",
  path_note,
  "",
  "**为什么阈值取 0.4**：STRING 自己的分档是 0.15（低）/ 0.4（中）/ 0.7（高）。",
  "0.4 是文献中最常用的默认档；低于此值会大量引入文本挖掘带来的假阳性边，",
  "高于此值则容易让几十个基因的小网络碎成互不相连的孤岛。",
  "",
  "**为什么不添加额外节点**：STRING 的 API 支持自动把数据库中的邻居蛋白加进网络，",
  "但那样得到的就不是「这批基因之间的互作」，而是「这批基因的邻域」，节点数会失控，",
  "也无法回答「候选基因是否彼此关联」这个原问题。因此本次用 `/network` 端点，",
  "只返回输入集合内部的边。",
  "",
  "---",
  "",
  "## 二、基因名清理与映射情况",
  "",
  "### 2.1 基因名清理（映射前）",
  "",
  paste0("输入为 `results/HubGene/", basename(F_GENES), "`，基因名列 `", COL_GENE, "`。"),
  "清理步骤：去首尾空格 → 剔除空值/NA → 去重（保留首次出现）。",
  "",
  paste0("| ", paste(audit[[1]], collapse = " | "), " |"),
  paste0("| ", paste(rep("---", nrow(audit)), collapse = " | "), " |"),
  paste0("| ", paste(audit[[2]], collapse = " | "), " |"),
  "",
  paste0("清理后用于查询的唯一基因 **", length(sym_clean), "** 个",
         if (n_removed > 0) paste0("（剔除空值 ", n_removed, " 个）") else "",
         if (n_dup > 0)     paste0("（剔除重复 ", n_dup, " 个）")     else "", "。"),
  "",
  "### 2.2 映射结果",
  "",
  paste0("| 项目 | 数量 | 占比 |"),
  paste0("| --- | --- | --- |"),
  paste0("| 输入候选基因 | ", nrow(mapping), " | 100% |"),
  paste0("| 成功映射到 STRING 蛋白 | **", n_mapped, "** | ",
         sprintf("%.1f%%", 100 * n_mapped / nrow(mapping)), " |"),
  paste0("| 未能映射 | **", n_unmapped, "** | ",
         sprintf("%.1f%%", 100 * n_unmapped / nrow(mapping)), " |"),
  "",
  unmapped_txt,
  "",
  if (n_unmapped > 0) paste0(
    "> 未映射的常见原因有三类：(1) 该 symbol 是 miRNA 或非编码基因（如 `MIR208A`），",
    "STRING 收录的是蛋白质，不收录 miRNA；(2) 该 symbol 在 STRING 的别名表中登记的",
    "是另一个更常用的基因名；(3) 该 symbol 已被引退或改名。这些基因不参与网络，",
    "但在候选基因表中的差异表达与 WGCNA 证据仍然成立。") else "",
  "",
  paste0("完整映射表见 `02_gene_protein_mapping.csv`（含 STRING 蛋白 ID、首选名、",
         "蛋白功能注释，并附带 log2FC / padj / Regulation / MM_kME / GS 等上游注释）。  "),
  "未映射基因单独存于 `03_unmapped_genes.csv`。",
  "",
  "---",
  "",
  "## 三、网络规模",
  "",
  paste0("| ", paste(net_stats[[1]], collapse = " | "), " |"),
  paste0("| ", paste(rep("---", nrow(net_stats)), collapse = " | "), " |"),
  paste0("| ", paste(net_stats[[2]], collapse = " | "), " |"),
  "",
  iso_txt,
  "",
  paste0("孤立节点是网络结构的真实组成部分：它们说明这些基因在 STRING v", VER_ACTUAL,
         " 中，与本次输入的其余基因之间没有达到 0.4 置信度的互作证据。",
         "这不代表它们没有功能，只代表它们不在这张网络的连接结构里。"),
  "",
  "### 度数最高的节点",
  "",
  "| 基因 | Degree | 最强边的 combined score |",
  "| --- | --- | --- |",
  top_txt,
  "",
  "---",
  "",
  "## 四、证据通道",
  "",
  paste0("STRING 的 combined score 由 7 条证据通道加权合成。本次 ", n_edges,
         " 条边在各通道上的分布如下（分数为 0-1）："),
  "",
  "| 通道字段 | 含义 | 有该证据的边数 | 平均分数 | 最大分数 |",
  "| --- | --- | --- | --- | --- |",
  ch_txt,
  "",
  "其中 `ascore`（共表达）与 `tscore`（文本挖掘）通常是边数最多的两条通道；",
  "`escore`（实验证据）边数普遍偏少，这是 STRING 的固有特点——实验验证的物理互作",
  "本来就远少于计算预测的关联。",
  "",
  "---",
  "",
  "## 五、网络图",
  "",
  paste0("文件：`results/PPI/", basename(f_pdf), "`（矢量 PDF）与 `", basename(f_png),
         "`（600 dpi PNG）。  "),
  paste0("尺寸：", W_MM, " × ", H_MM, " mm（Nature 双栏宽度）。"),
  "",
  "**图的编码规则**",
  "",
  "- **节点大小 ∝ Degree**（用平方根压缩，避免 hub 节点过大；孤立节点取最小尺寸）。",
  paste0("- **节点颜色区分节点类型**：蓝色 `", COL_CORE, "` 为有连接的节点（degree ≥ 1），",
         "灰色 `", COL_ISO, "` 为孤立节点（degree = 0）。"),
  "- **边宽 ∝ 互作置信度**（combined_score 从 0.4 到 1.0 对应线宽 0.25 到 1.30 pt）。",
  paste0("- **布局**：有连接的部分用 Fruchterman-Reingold 力导向算法（niter = 3000），",
         "孤立节点没有吸引力会被甩到边缘且容易重叠，因此单独均匀排布在外圈。",
         "随机种子固定为 ", SEED, "，同一输入必然得到同一张图。"),
  "- **标签避让**：用 ggrepel 做排斥式排布，避免重叠。",
  "",
  "---",
  "",
  "## 六、输出文件清单",
  "",
  "| 文件 | 内容 |",
  "| --- | --- |",
  "| `01_gene_symbol_audit.csv` | 基因名清理前的质量审计 |",
  paste0("| `02_gene_protein_mapping.csv` | 基因→蛋白映射表（", nrow(mapping), " 行，含上游注释） |"),
  paste0("| `03_unmapped_genes.csv` | 未映射基因（", n_unmapped, " 个） |"),
  paste0("| `04_string_edges_raw.csv` | STRING 返回的原始边表，含 7 条证据通道分数（",
         nrow(edges), " 条） |"),
  paste0("| `05_node_degree_table.csv` | 节点属性表：Degree / 加权度 / 最强边分数（",
         nrow(node_tab), " 个节点） |"),
  paste0("| `06_isolated_nodes.csv` | 孤立节点清单（", n_iso, " 个） |"),
  "| `07_network_statistics.csv` | 网络规模统计 |",
  "| `08_evidence_channels.csv` | 证据通道统计 |",
  paste0("| `", basename(f_pdf), "` | 网络图，矢量 PDF |"),
  paste0("| `", basename(f_png), "` | 网络图，600 dpi PNG |"),
  "| `PPI_report.md` | 本报告 |",
  "| `sessionInfo_17_ppi_network.txt` | 运行环境 |",
  "| `run_log.txt` | 运行日志 |",
  "",
  "---",
  "",
  "## 七、复现方式",
  "",
  "```bash",
  "cd \"GEObulk analysis project\"",
  "Rscript scripts/17_ppi_network.R",
  "```",
  "",
  paste0("依赖包：STRINGdb（Bioconductor，必需）、igraph、ggplot2、httr、jsonlite；",
         "可选：ggrepel（标签避让）、ragg（PNG 抗锯齿，缺失时自动降级）。  "),
  paste0("随机种子 ", SEED, " 已固定，力导向布局逐次可复现。"),
  ""
)
wtxt(rep_md, file.path(DIR_OUT, "PPI_report.md"))

## ---------------------------------------------------------------------------
## 9. 运行环境与日志
## ---------------------------------------------------------------------------
wtxt(c(capture.output(print(sessionInfo()))),
     file.path(DIR_OUT, "sessionInfo_17_ppi_network.txt"))
wtxt(c(paste0(DS, " PPI 网络构建运行日志"), "========================================", LOG),
     file.path(DIR_OUT, "run_log.txt"))

msg("=== 完成。节点 ", n_nodes, " / 边 ", n_edges, " / 孤立节点 ", n_iso,
    "；结果已写入 results/PPI/ ===")
