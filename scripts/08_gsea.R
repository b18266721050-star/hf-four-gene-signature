# =============================================================================
# 08_gsea.R  ——  GSE57338 全基因排序的 GSEA（MSigDB Hallmark, H）
#
# 输入（results/DEG/）：
#   GSE57338_DEGall.csv                     全部检验基因（18,726 行），Disease vs Control
#   GSE57338_DEGall_DCM_vs_Control.csv      亚组：DCM vs Control（补充分析）
#   GSE57338_DEGall_ICM_vs_Control.csv      亚组：ICM vs Control（补充分析）
#
# 输出（results/GSEA/，关键文件同步一份到 results/ 根目录）：
#   GSE57338_GSEA_result.csv                显著富集通路（padj<0.05）—— 主交付表
#   GSE57338_GSEA_allSets.csv               全部 50 条 Hallmark 集的结果（不过滤）
#   GSE57338_DCM_vs_Control_GSEA_result.csv 亚组补充表
#   GSE57338_ICM_vs_Control_GSEA_result.csv 亚组补充表
#   GSEA_NES_comparison.csv                 三套分析的 NES 对照（共享 / 亚型特异）
#   GSE57338_GSEA_ranking_sensitivity.csv   排序指标敏感性（log2FC vs limma t）
#   GSEA_gene_ranking.csv                   实际送入 GSEA 的排序向量（可复核）
#   GSEA_id_mapping.csv                     SYMBOL -> ENTREZID 映射明细
#   GSEA_bubble_top5.pdf/.png               总览气泡图（最显著 5 条）
#   GSEA_bubble_overview.pdf/.png           总览气泡图（激活 top5 + 抑制 top5）
#   GSEA_runningscore_top5.pdf/.png         前 5 条通路的 running score 图（拼版）
#   GSEA_runningscore_<i>_<pathway>.pdf/png 前 5 条通路各自的 running score 图
#   GSEA_summary.txt                        参数、条目数、top5 激活/抑制汇总
#   GSEA_figure_legends.txt                 图注（英文，可直接引用）
# 输出（report/）：
#   sessionInfo_08_gsea.txt
#
# 用法（项目根目录下）：
#   Rscript scripts/08_gsea.R
#
# -----------------------------------------------------------------------------
# 方法学要点（为什么这样设参数）
# -----------------------------------------------------------------------------
#   1) 物种：GSE57338 平台 GPL11532（Affymetrix HuGene-1.1-ST），组织为人左心室
#      心肌 -> Homo sapiens。基因集用 msigdbr(species="Homo sapiens",
#      collection="H")，ID 取 ncbi_gene（Entrez），与 org.Hs.eg.db 的 ENTREZID
#      同一命名空间；脚本内会显式断言 db_target_species 全为 Homo sapiens，
#      并报告基因集覆盖度，防止"人基因集套小鼠 ID"这类物种错配。
#   2) 输入是**全部检验基因**（18,726），不是 269 个显著基因。GSEA 的前提就是
#      无阈值全基因组排序；用显著基因做 GSEA 会把排序截断在阈值处，结果无意义。
#   3) 排序依据：任务指定 log2FoldChange（从大到小）。排序指标的选择会决定结果
#      ——这是 GSEA 最容易被忽略的一点。log2FC 只带方向与幅度、不含方差校准，
#      低表达/低计数基因的极端倍数变化可能占据 leading edge。因此本脚本额外用
#      limma moderated t 作排序做一次敏感性分析（见第 9 节），报告两者 NES 的
#      Spearman 相关与显著集重叠率，让读者自行判断排序指标的影响。
#      同一 log2FC 值的并列基因用 -log10(P) 做 1e-7 量级的确定性微扰打破
#      （远小于 log2FC 的 1e-4 分辨率，不改变真实排序，只是让并列可复现）。
#   4) 基因 ID 转换：SYMBOL -> ENTREZID（bitr）。一对多时保留第一个；同一
#      ENTREZID 被多个 symbol 命中时，保留 |log2FC| 最大的那条（避免重复基因在
#      running sum 中被重复计数）。GSEA 要求命名向量严格降序且基因唯一。
#   5) 参数：exponent=1（加权 KS，按 |统计量| 加权，Subramanian 2005 默认）、
#      minGSSize=10 / maxGSSize=500（过小集 p 不稳，过大集几乎必然"显著"）、
#      eps=0（替代已废弃的 nPerm，解析到极小的精确 p）、BH 校正、
#      set.seed(123)+seed=TRUE（fgseaMultilevel 是蒙特卡洛，必须固定种子）。
#   6) 显著性：pvalueCutoff 设为 1 先取全部条目，再在 R 里显式按 p.adjust<0.05
#      过滤——阈值透明、可复核（pvalueCutoff 过滤的是 p.adjust 而非原始 p）。
#   7) 方向：NES>0 = 该通路基因整体偏向排序顶端 = 在疾病组被激活；NES<0 = 偏向
#      底端 = 被抑制。注意这是**基因集水平**的协同变化，与"某个基因上调"不是
#      一回事，报告里必须区分。
#   8) 已知局限（写进报告）：preranked GSEA 只能用基因置换（gene permutation），
#      它假设基因独立，会低估共调控基因集的方差，对相关基因集偏反保守；若需
#      严格的 I 型错误控制，应在有表达矩阵与设计矩阵时用 limma::camera。
#      本脚本严格按任务要求使用 clusterProfiler/fgsea 的 preranked 实现，并在
#      报告中声明置换类型。
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 0. 路径
## ---------------------------------------------------------------------------
get_this_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  return(getwd())
}

wd <- getwd()
proj_dir <- if (dir.exists(file.path(wd, "results"))) normalizePath(wd) else dirname(get_this_dir())
deg_dir  <- file.path(proj_dir, "results", "DEG")
res_root <- file.path(proj_dir, "results")
out_dir  <- file.path(res_root, "GSEA")
rep_dir  <- file.path(proj_dir, "report")

for (d in c(out_dir, rep_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function(ch = "-", n = 76) cat(strrep(ch, n), "\n", sep = "")

## ---------------------------------------------------------------------------
## 1. 参数
## ---------------------------------------------------------------------------
ACCESSION  <- "GSE57338"
ORGANISM   <- "Homo sapiens"
ORGDB      <- "org.Hs.eg.db"
COLLECTION <- "H"                 # MSigDB Hallmark
RANK_COL   <- "log2FoldChange"    # 任务指定排序依据
PADJ_CUT   <- 0.05
ADJ_METHOD <- "BH"
MIN_GS     <- 10
MAX_GS     <- 500
EXPONENT   <- 1
EPS        <- 0
SEED       <- 123
TOP_N      <- 5
FONT_FAMILY <- "Arial"

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
})

## ---------------------------------------------------------------------------
## 2. 待分析的数据集（主分析 + 亚组补充）
## ---------------------------------------------------------------------------
ds_def <- data.frame(
  tag   = c("GSE57338",
            "GSE57338_DCM_vs_Control",
            "GSE57338_ICM_vs_Control"),
  file  = c(file.path(deg_dir, "GSE57338_DEGall.csv"),
            file.path(deg_dir, "GSE57338_DEGall_DCM_vs_Control.csv"),
            file.path(deg_dir, "GSE57338_DEGall_ICM_vs_Control.csv")),
  label = c("Disease (DCM+ICM) vs Control",
            "DCM vs Control",
            "ICM vs Control"),
  main  = c(TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE)

miss <- ds_def$file[!file.exists(ds_def$file)]
if (length(miss) > 0) stop("missing DEG table: ", paste(miss, collapse = ", "))

# 兜底：把 results/ 下其余 *_DEGall*.csv 列出来，说明为何未纳入
extra <- setdiff(list.files(deg_dir, pattern = "_DEGall.*\\.csv$", full.names = TRUE), ds_def$file)

msg("== Input ==")
for (i in seq_len(nrow(ds_def))) {
  d <- read.csv(ds_def$file[i], check.names = FALSE)
  msg(sprintf("  %-28s %6d genes  (%s)", ds_def$tag[i], nrow(d), ds_def$label[i]))
}
if (length(extra) > 0) {
  msg("\n  [note] other *_DEGall* tables found but not analysed (same contrast, sensitivity variants):")
  for (e in extra) msg("         - ", basename(e))
}

## ---------------------------------------------------------------------------
## 3. Hallmark 基因集（物种核验）
## ---------------------------------------------------------------------------
msg("\n== MSigDB Hallmark (H), species check ==")
hm <- msigdbr(species = ORGANISM, collection = COLLECTION)
hm <- as.data.frame(hm)

sp_code <- unique(as.character(hm$db_target_species))
# msigdbr 26.x 用两字母物种码（HS = Homo sapiens），旧版用全名；两种都接受
if (!all(sp_code %in% c("HS", "Homo sapiens", "Homo_sapiens")))
  stop("Hallmark gene sets are NOT human: ", paste(sp_code, collapse = ", "))
if (!all(grepl("^[0-9]+$", as.character(hm$ncbi_gene)))) stop("ncbi_gene is not Entrez-like")

t2g <- unique(data.frame(term = as.character(hm$gs_name),
                         gene = as.character(hm$ncbi_gene),
                         stringsAsFactors = FALSE))
t2n <- unique(data.frame(term = as.character(hm$gs_name),
                         name = as.character(hm$gs_name),
                         stringsAsFactors = FALSE))

msg("  msigdbr version      : ", as.character(packageVersion("msigdbr")))
msg("  db_target_species    : ", paste(sp_code, collapse = ", "), "  (HS = Homo sapiens, verified)")
msg("  collection           : ", paste(unique(hm$gs_collection), collapse = ", "),
    " (", unique(hm$gs_collection_name), ")")
msg("  gene sets            : ", length(unique(t2g$term)))
msg("  term2gene rows       : ", nrow(t2g), " (Entrez, unique term-gene pairs)")
msg("  MSigDB db_version    : ", paste(unique(hm$db_version), collapse = ", "))
msg("  Entrez namespace     : ncbi_gene -> matches org.Hs.eg.db ENTREZID")

## ---------------------------------------------------------------------------
## 4. 构建排序向量（SYMBOL -> ENTREZID，去重，降序）
## ---------------------------------------------------------------------------
id_map_all <- list()

build_rank <- function(deg, rank_col = RANK_COL, tag = "") {
  d <- deg
  d$GeneSymbol <- as.character(d$GeneSymbol)
  d <- d[!is.na(d$GeneSymbol) & nzchar(d$GeneSymbol), , drop = FALSE]
  d <- d[!is.na(d[[rank_col]]), , drop = FALSE]

  map <- suppressWarnings(bitr(unique(d$GeneSymbol), fromType = "SYMBOL",
                               toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  map <- map[!duplicated(map$SYMBOL), , drop = FALSE]

  idx  <- match(d$GeneSymbol, map$SYMBOL)
  keep <- !is.na(idx)
  d    <- d[keep, , drop = FALSE]
  d$ENTREZID <- map$ENTREZID[idx[keep]]

  # 同一 ENTREZID 多 symbol 命中：保留 |rank| 最大的一条（避免重复计数）
  d <- d[order(-abs(d[[rank_col]])), , drop = FALSE]
  n_dup_entrez <- sum(duplicated(d$ENTREZID))
  d <- d[!duplicated(d$ENTREZID), , drop = FALSE]

  v <- setNames(as.numeric(d[[rank_col]]), d$ENTREZID)

  # ---- 并列处理 -------------------------------------------------------------
  # log2FC 只保留 4 位小数（分辨率 1e-4），1.7 万个基因落在约 4400 个取值上，
  # 必然产生大量并列（本例 >1 万）。fgsea 遇到并列会随机决定顺序，结果不可复现。
  # 这里按"主键 = log2FC 降序，次键 = P 值升序，末键 = 基因符号"排定唯一顺序，
  # 再把组内次序折算成 <5e-5 的微小偏移（小于 1e-4 的分辨率，不改变真实排序，
  # 只保证严格单调且可复现）。
  pv <- if ("P.Value" %in% names(d)) d$P.Value else rep(NA_real_, nrow(d))
  pv[is.na(pv)] <- 1
  ord <- order(-v, pv, names(v), method = "radix")
  n_tie <- sum(duplicated(v))
  v_ord <- v[ord]
  grp <- cumsum(c(TRUE, diff(v_ord) != 0))
  rank_in_tie <- stats::ave(seq_along(ord), grp, FUN = seq_along)
  grp_size    <- stats::ave(seq_along(ord), grp, FUN = length)
  offset <- ((grp_size - rank_in_tie) / (grp_size + 1)) * 5e-5
  v2 <- v_ord + offset
  names(v2) <- names(v)[ord]
  v2 <- sort(v2, decreasing = TRUE)
  n_residual_tie <- sum(duplicated(v2))

  list(vec = v2, n_input = nrow(deg), n_named = length(idx),
       n_mapped = sum(keep), n_dup_entrez = n_dup_entrez,
       n_tie = n_tie, n_residual_tie = n_residual_tie, map = map)
}

rank_info <- list()

## ---------------------------------------------------------------------------
## 5. GSEA 主函数
## ---------------------------------------------------------------------------
run_gsea <- function(deg, tag, label, rank_col = RANK_COL) {
  msg("\n", strrep("=", 76))
  msg("== GSEA: ", tag, "  (", label, ") ==")
  hr()

  br <- build_rank(deg, rank_col = rank_col, tag = tag)
  gl <- br$vec
  rank_info[[tag]] <<- br

  ov <- sum(names(gl) %in% t2g$gene) / length(gl)
  msg(sprintf("  genes in DEG table        : %d", br$n_input))
  msg(sprintf("  SYMBOL -> ENTREZID mapped : %d (%.1f%%)", br$n_mapped, 100 * br$n_mapped / br$n_input))
  msg(sprintf("  duplicate ENTREZID removed: %d", br$n_dup_entrez))
  msg(sprintf("  ties on %-16s : %d (broken by P value, residual ties: %d)",
              rank_col, br$n_tie, br$n_residual_tie))
  msg(sprintf("  ranked vector length      : %d", length(gl)))
  msg(sprintf("  range                     : %.4f .. %.4f (sorted decreasing: %s)",
              min(gl), max(gl),
              identical(unname(gl), unname(sort(gl, decreasing = TRUE)))))
  msg(sprintf("  covered by Hallmark sets  : %d (%.1f%%)", sum(names(gl) %in% t2g$gene), 100 * ov))

  set.seed(SEED)
  t0 <- Sys.time()
  gse <- GSEA(geneList     = gl,
              TERM2GENE    = t2g,
              TERM2NAME    = t2n,
              exponent     = EXPONENT,
              minGSSize    = MIN_GS,
              maxGSSize    = MAX_GS,
              eps          = EPS,
              pvalueCutoff = 1,          # 先全取，随后显式按 p.adjust 过滤
              pAdjustMethod = ADJ_METHOD,
              seed         = TRUE,
              verbose      = FALSE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  df <- as.data.frame(gse)
  msg(sprintf("  sets tested               : %d  (%.1fs, fgseaMultilevel, eps=%s)",
              nrow(df), dt, EPS))

  # ENTREZ -> SYMBOL（core_enrichment 可读化）
  gse_r <- tryCatch(setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                    error = function(e) gse)
  dfr <- as.data.frame(gse_r)
  core_entrez <- df$core_enrichment[match(dfr$ID, df$ID)]

  out <- data.frame(
    Pathway          = gsub("_", " ", sub("^HALLMARK_", "", dfr$ID)),
    Hallmark_ID      = dfr$ID,
    Direction        = ifelse(dfr$NES > 0, "Activated", "Suppressed"),
    SetSize          = dfr$setSize,
    ES               = round(dfr$enrichmentScore, 4),
    NES              = round(dfr$NES, 4),
    P.Value          = signif(dfr$pvalue, 4),
    padj             = signif(dfr$p.adjust, 4),
    qvalue           = signif(dfr$qvalue, 4),
    LeadingEdge_n    = vapply(strsplit(as.character(dfr$core_enrichment), "/", fixed = TRUE), length, integer(1)),
    Rank_at_max_ES   = dfr$rank,
    PeakES_pct       = round(100 * dfr$rank / length(gl), 2),
    Leading_edge     = dfr$leading_edge,
    Core_Enrichment  = dfr$core_enrichment,
    Core_Enrichment_Entrez = core_entrez,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$padj, out$P.Value), , drop = FALSE]
  rownames(out) <- NULL

  sig <- out[out$padj < PADJ_CUT, , drop = FALSE]
  up  <- sig[sig$NES > 0, , drop = FALSE]; up <- up[order(up$padj), , drop = FALSE]
  dn  <- sig[sig$NES < 0, , drop = FALSE]; dn <- dn[order(dn$padj), , drop = FALSE]

  msg(sprintf("  significant (padj<%.2f)   : %d  (activated %d / suppressed %d)",
              PADJ_CUT, nrow(sig), nrow(up), nrow(dn)))

  list(gse = gse_r, all = out, sig = sig, up = up, dn = dn, gl = gl)
}

## ---------------------------------------------------------------------------
## 6. 逐个数据集运行
## ---------------------------------------------------------------------------
res_list <- list()
for (i in seq_len(nrow(ds_def))) {
  deg <- read.csv(ds_def$file[i], check.names = FALSE)
  res_list[[ds_def$tag[i]]] <- run_gsea(deg, ds_def$tag[i], ds_def$label[i])
}

# ID 映射明细（主分析）
write.csv(rank_info[["GSE57338"]]$map, file.path(out_dir, "GSEA_id_mapping.csv"), row.names = FALSE)

# 实际送入 GSEA 的排序向量（可复核）
rk <- rank_info[["GSE57338"]]
write.csv(data.frame(ENTREZID = names(rk$vec), log2FoldChange = as.numeric(rk$vec)),
          file.path(out_dir, "GSEA_gene_ranking.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 7. 写结果表
## ---------------------------------------------------------------------------
msg("\n== Output tables ==")
for (tag in names(res_list)) {
  r  <- res_list[[tag]]
  f1 <- file.path(out_dir, paste0(tag, "_GSEA_result.csv"))
  f2 <- file.path(out_dir, paste0(tag, "_GSEA_allSets.csv"))
  write.csv(r$sig, f1, row.names = FALSE)
  write.csv(r$all, f2, row.names = FALSE)
  if (ds_def$main[ds_def$tag == tag]) {
    write.csv(r$sig, file.path(res_root, paste0(tag, "_GSEA_result.csv")), row.names = FALSE)
  }
  msg(sprintf("  %-58s sig=%2d / all=%2d", basename(f1), nrow(r$sig), nrow(r$all)))
}

# 三套分析的 NES 对照
cmp <- Reduce(function(a, b) merge(a, b, by = c("Hallmark_ID", "Pathway"), all = TRUE),
              list(
                data.frame(Hallmark_ID = res_list[[1]]$all$Hallmark_ID,
                           Pathway = res_list[[1]]$all$Pathway,
                           NES_Disease = res_list[[1]]$all$NES,
                           padj_Disease = res_list[[1]]$all$padj),
                data.frame(Hallmark_ID = res_list[[2]]$all$Hallmark_ID,
                           Pathway = res_list[[2]]$all$Pathway,
                           NES_DCM = res_list[[2]]$all$NES,
                           padj_DCM = res_list[[2]]$all$padj),
                data.frame(Hallmark_ID = res_list[[3]]$all$Hallmark_ID,
                           Pathway = res_list[[3]]$all$Pathway,
                           NES_ICM = res_list[[3]]$all$NES,
                           padj_ICM = res_list[[3]]$all$padj)))
cmp$Direction <- with(cmp, ifelse(NES_Disease > 0, "Activated", "Suppressed"))
cmp$Sig_in <- with(cmp, paste0(ifelse(padj_Disease < PADJ_CUT, "D", "-"),
                               ifelse(padj_DCM < PADJ_CUT, "M", "-"),
                               ifelse(padj_ICM < PADJ_CUT, "I", "-")))
cmp <- cmp[order(-abs(cmp$NES_Disease)), ]
write.csv(cmp, file.path(out_dir, "GSEA_NES_comparison.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 8. 显著通路汇总文本
## ---------------------------------------------------------------------------
sf <- file.path(out_dir, "GSEA_summary.txt")
con <- file(sf, open = "wt", encoding = "UTF-8")
w <- function(...) writeLines(paste0(...), con)

w("GSEA (MSigDB Hallmark, H) summary - ", ACCESSION)
w("run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("organism: ", ORGANISM, " | msigdbr ", as.character(packageVersion("msigdbr")),
  " | org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")),
  " | clusterProfiler ", as.character(packageVersion("clusterProfiler")),
  " | fgsea ", as.character(packageVersion("fgsea")))
w("ranking metric: ", RANK_COL, " (decreasing) | exponent: ", EXPONENT,
  " | gene permutation (preranked) | eps=", EPS, " | seed=", SEED, " | BH")
w("cutoff: p.adjust < ", PADJ_CUT, " | minGSSize=", MIN_GS, " | maxGSSize=", MAX_GS)
w("")
for (tag in names(res_list)) {
  r <- res_list[[tag]]
  w(sprintf("[%s] ranked genes=%d  sets tested=%d  significant=%d  (activated=%d, suppressed=%d)",
            tag, length(r$gl), nrow(r$all), nrow(r$sig), nrow(r$up), nrow(r$dn)))
  w("  -- top 5 ACTIVATED (NES > 0) --")
  if (nrow(r$up) == 0) w("     none") else
    for (j in seq_len(min(TOP_N, nrow(r$up))))
      w(sprintf("  %d. %-44s NES=%6.2f  p=%.2e  padj=%.2e  n=%d",
                j, r$up$Pathway[j], r$up$NES[j], r$up$P.Value[j], r$up$padj[j], r$up$SetSize[j]))
  w("  -- top 5 SUPPRESSED (NES < 0) --")
  if (nrow(r$dn) == 0) w("     none") else
    for (j in seq_len(min(TOP_N, nrow(r$dn))))
      w(sprintf("  %d. %-44s NES=%6.2f  p=%.2e  padj=%.2e  n=%d",
                j, r$dn$Pathway[j], r$dn$NES[j], r$dn$P.Value[j], r$dn$padj[j], r$dn$SetSize[j]))
  w("")
}
w("Note: NES > 0 = pathway genes shifted to the top of the log2FC ranking = activated in disease;")
w("      NES < 0 = shifted to the bottom = suppressed. Preranked GSEA uses GENE permutation,")
w("      which is anti-conservative for co-regulated gene sets; see the report for caveats.")
close(con)

## ---------------------------------------------------------------------------
## 9. 排序指标敏感性：log2FC vs limma moderated t（仅主分析）
## ---------------------------------------------------------------------------
sens_file <- file.path(out_dir, paste0(ACCESSION, "_GSEA_ranking_sensitivity.csv"))
deg_main <- read.csv(ds_def$file[1], check.names = FALSE)
if ("t" %in% names(deg_main)) {
  msg("\n== Ranking sensitivity: log2FC vs limma t ==")
  r_t <- run_gsea(deg_main, paste0(ACCESSION, "_rankT"),
                  "Disease vs Control (ranked by limma moderated t)", rank_col = "t")
  a <- res_list[[1]]$all; b <- r_t$all
  m <- merge(a[, c("Hallmark_ID", "NES", "padj")], b[, c("Hallmark_ID", "NES", "padj")],
             by = "Hallmark_ID", suffixes = c("_log2FC", "_t"))
  m <- m[order(-abs(m$NES_log2FC)), ]
  rho <- suppressWarnings(cor(m$NES_log2FC, m$NES_t, method = "spearman"))
  s1 <- m$Hallmark_ID[m$padj_log2FC < PADJ_CUT]; s2 <- m$Hallmark_ID[m$padj_t < PADJ_CUT]
  ov <- length(intersect(s1, s2))
  m$Sign_agree <- sign(m$NES_log2FC) == sign(m$NES_t)
  write.csv(m, sens_file, row.names = FALSE)
  msg(sprintf("  Spearman rho(NES) = %.3f | sign agreement = %.1f%% | significant overlap = %d/%d (log2FC) vs %d (t)",
              rho, 100 * mean(m$Sign_agree), ov, length(s1), length(s2)))
  sens_txt <- c(
    sprintf("Ranking sensitivity check (%s)", ACCESSION),
    sprintf("  metric A: log2FoldChange  -> %d significant sets", length(s1)),
    sprintf("  metric B: limma moderated t -> %d significant sets", length(s2)),
    sprintf("  Spearman rho of NES = %.3f", rho),
    sprintf("  sign agreement = %.1f%%", 100 * mean(m$Sign_agree)),
    sprintf("  significant sets shared = %d", ov),
    "  detail: GSE57338_GSEA_ranking_sensitivity.csv")
  writeLines(sens_txt, file.path(out_dir, "GSEA_ranking_sensitivity.txt"))
} else {
  msg("\n  [skip] no 't' column; ranking sensitivity not run")
}

## ---------------------------------------------------------------------------
## 10. 可视化
## ---------------------------------------------------------------------------
MORANDI_DIV <- c(low = "#3D6E8F", mid = "#F4F1EA", high = "#A6402D")   # 蓝=抑制 红=激活

nature_theme <- function(base_size = 7) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", linewidth = 0.35, fill = NA),
      axis.line          = element_blank(),
      axis.ticks         = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length  = unit(1.1, "mm"),
      axis.text          = element_text(colour = "black", size = base_size),
      axis.title         = element_text(colour = "black", size = base_size + 0.5),
      legend.title       = element_text(size = base_size, face = "plain"),
      legend.text        = element_text(size = base_size - 0.5),
      legend.key.height  = unit(3.0, "mm"),
      legend.key.width   = unit(2.4, "mm"),
      legend.background  = element_blank(),
      plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin        = margin(3, 3, 3, 3)
    )
}

wrap_lab <- function(x, width = 34) vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))

# 统一配色：激活=红，抑制=蓝
dir_cols <- c(Activated = "#A6402D", Suppressed = "#3D6E8F")

make_bubble <- function(d, title, order_by = c("padj", "NES"), x_pad = 0.08, nes_lim = 2.5) {
  order_by <- match.arg(order_by)
  d <- d[!is.na(d$NES), , drop = FALSE]
  d$lab <- wrap_lab(d$Pathway)
  d$NLP <- -log10(pmax(d$padj, 1e-300))
  if (order_by == "padj") {
    d <- d[order(d$padj), , drop = FALSE]
    lev <- d$lab[order(d$padj, decreasing = TRUE)]
  } else {
    d <- d[order(d$NES), , drop = FALSE]
    lev <- d$lab[order(d$NES)]
  }
  d$lab <- factor(d$lab, levels = unique(lev))

  ggplot(d, aes(x = NES, y = lab)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey40") +
    geom_segment(aes(x = 0, xend = NES, y = lab, yend = lab, colour = NES),
                 linewidth = 0.45) +
    geom_point(aes(size = NLP, colour = NES), stroke = 0) +
    geom_text(data = d[d$NES > 0, , drop = FALSE], aes(label = sprintf("%.2f", NES)),
              hjust = -0.45, size = 2.0, family = FONT_FAMILY) +
    geom_text(data = d[d$NES < 0, , drop = FALSE], aes(label = sprintf("%.2f", NES)),
              hjust = 1.45, size = 2.0, family = FONT_FAMILY) +
    scale_colour_gradient2(low = MORANDI_DIV["low"], mid = MORANDI_DIV["mid"],
                           high = MORANDI_DIV["high"], midpoint = 0,
                           limits = c(-nes_lim, nes_lim), oob = scales::squish,
                           name = "NES", guide = guide_colourbar(order = 1)) +
    scale_size_continuous(range = c(1.0, 4.6), breaks = scales::pretty_breaks(n = 3),
                          name = expression(-log[10]~italic(P)[adj])) +
    scale_x_continuous(expand = expansion(mult = c(x_pad, x_pad))) +
    labs(title = title, x = "Normalized enrichment score (NES)", y = NULL) +
    nature_theme(7) + theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25))
}

save_fig <- function(p, stem, width_mm, height_mm, dir = out_dir) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(dir, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  msg(sprintf("  saved: %s.pdf/.png  (%.0f x %.0f mm, 600 dpi)", stem, width_mm, height_mm))
}

msg("\n== Figures ==")
r1 <- res_list[[1]]
sig1 <- r1$sig

# 两张气泡图共用同一 NES 色标范围，颜色才可跨图比较（对称、含 0）
NES_LIM <- ceiling(max(abs(sig1$NES), na.rm = TRUE) * 2) / 2

# 10.1 总览气泡图 —— 最显著 5 条（任务要求）
top5 <- head(sig1[order(sig1$padj), , drop = FALSE], TOP_N)
p_top5 <- make_bubble(top5, sprintf("Top 5 Hallmark pathways (%s, padj < %.2f)", ACCESSION, PADJ_CUT),
                      order_by = "padj", nes_lim = NES_LIM)
save_fig(p_top5, "GSEA_bubble_top5", 168, 62)
ggsave(file.path(res_root, "GSEA_bubble_top5.pdf"), p_top5, width = 168/25.4, height = 62/25.4, device = cairo_pdf, bg = "white")
ggsave(file.path(res_root, "GSEA_bubble_top5.png"), p_top5, width = 168/25.4, height = 62/25.4, dpi = 600, device = ragg::agg_png, bg = "white")

# 10.2 总览气泡图 —— 激活 top5 + 抑制 top5（方向对照，信息量更大）
ov <- rbind(head(r1$up[, ], TOP_N), head(r1$dn[, ], TOP_N))
p_ov <- make_bubble(ov, sprintf("Hallmark GSEA overview: top 5 activated / top 5 suppressed (%s)", ACCESSION),
                    order_by = "NES", nes_lim = NES_LIM)
save_fig(p_ov, "GSEA_bubble_overview", 168, 104)
ggsave(file.path(res_root, "GSEA_bubble_overview.pdf"), p_ov, width = 168/25.4, height = 104/25.4, device = cairo_pdf, bg = "white")
ggsave(file.path(res_root, "GSEA_bubble_overview.png"), p_ov, width = 168/25.4, height = 104/25.4, dpi = 600, device = ragg::agg_png, bg = "white")

# 10.3 前 5 条通路的 running score 图（经典三段式：RS 曲线 / 命中位置 / 排序指标）
# 说明：这里按 Subramanian 2005 的定义自行复算加权 running sum，而不是调用
# enrichplot::gseaplot2。原因有二：(1) gseaplot2 返回 aplot 复合对象，多图拼版
# 时无法被 patchwork/cowplot 再组合；(2) 自绘可统一配色与 Nature 风格主题。
# 复算的 ES 会与 fgsea 报告的 enrichmentScore 逐条核对（见日志 ES_check）。
running_sum <- function(gl, genes, exponent = 1) {
  N <- length(gl)
  hit <- names(gl) %in% genes
  nH <- sum(hit)
  if (nH == 0) return(NULL)
  w <- abs(gl)^exponent
  ph <- cumsum(ifelse(hit, w, 0)) / sum(w[hit])
  pm <- cumsum(ifelse(!hit, 1 / (N - nH), 0))
  ph - pm
}

plot_running_score <- function(gl, row, set_ids, tag = NULL, base_size = 7) {
  # set_ids 必须是**完整基因集**（不是 leading edge），否则 running sum 的三段结构不成立
  rs <- running_sum(gl, set_ids)
  if (is.null(rs)) return(NULL)
  N <- length(gl)
  i_max <- which.max(abs(rs))
  es_calc <- rs[i_max]
  col <- if (row$NES > 0) unname(dir_cols["Activated"]) else unname(dir_cols["Suppressed"])

  d <- data.frame(rank = seq_len(N), rs = rs)
  hit_pos <- which(names(gl) %in% set_ids)
  # leading edge：正向富集取 ES 峰之前的命中，负向取峰之后的命中
  le_pos  <- if (row$NES > 0) hit_pos[hit_pos <= i_max] else hit_pos[hit_pos >= i_max]

  # 面板 a：running enrichment score
  pa <- ggplot(d, aes(x = rank, y = rs)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey45") +
    geom_area(data = d[d$rs >= 0, , drop = FALSE], aes(y = rs), fill = col, alpha = 0.18) +
    geom_area(data = d[d$rs <= 0, , drop = FALSE], aes(y = rs), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 0.55) +
    geom_vline(xintercept = i_max, linetype = "22", linewidth = 0.3, colour = "grey30") +
    geom_point(data = d[i_max, , drop = FALSE], aes(x = rank, y = rs), size = 1.1, colour = col) +
    annotate("text", x = i_max,
             y = if (row$NES > 0) min(rs) else max(rs),
             hjust = if (i_max > 0.5 * N) 1.15 else -0.15,
             vjust = if (row$NES > 0) 1.3 else -0.7,
             label = sprintf("ES = %.3f", row$ES), size = 2.1, family = FONT_FAMILY) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, N + 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.18)),
                       labels = function(x) sprintf("%.2f", x)) +
    labs(title = tag, x = NULL, y = "Running ES") +
    nature_theme(base_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = base_size, face = "bold", lineheight = 1.05))

  # 面板 b：基因集成员在排序中的位置（leading edge 加深）
  pb <- ggplot() +
    geom_segment(data = data.frame(x = setdiff(hit_pos, le_pos)),
                 aes(x = x, xend = x, y = 0, yend = 1), linewidth = 0.28, colour = "grey62") +
    geom_segment(data = data.frame(x = le_pos),
                 aes(x = x, xend = x, y = 0, yend = 1), linewidth = 0.4, colour = col) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, N + 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size) +
    theme(plot.margin = margin(1, 3, 0, 3))

  # 面板 c：底层排序指标（log2FC）
  dc <- data.frame(rank = seq_len(N), v = as.numeric(gl),
                   sgn = factor(ifelse(gl >= 0, "up", "down"), levels = c("up", "down")))
  pc <- ggplot(dc, aes(x = rank, y = v, fill = sgn)) +
    geom_col(width = 1) +
    scale_fill_manual(values = c(up = "#A6402D", down = "#3D6E8F"), guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey45") +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, N + 0.5)) +
    labs(x = sprintf("Gene rank (log2FC, decreasing; 1 - %s)", format(N, big.mark = ",")),
         y = "log2FC") +
    nature_theme(base_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  list(plot = pa / pb / pc + patchwork::plot_layout(heights = c(3.0, 0.55, 1.25)),
       es_calc = es_calc, es_reported = row$ES)
}

gse1 <- r1$gse
top5_ids <- top5$Hallmark_ID
plist <- list()
es_chk <- data.frame()
for (j in seq_along(top5_ids)) {
  tg <- sprintf("(%s) %s\nNES = %.2f, padj = %.1e",
                letters[j], top5$Pathway[j], top5$NES[j], top5$padj[j])
  pj <- tryCatch(plot_running_score(r1$gl, top5[j, , drop = FALSE],
                                    set_ids = t2g$gene[t2g$term == top5_ids[j]], tag = tg),
                 error = function(e) { msg("  [warn] running-score plot failed for ",
                                           top5_ids[j], ": ", conditionMessage(e)); NULL })
  if (is.null(pj)) next
  es_chk <- rbind(es_chk, data.frame(Hallmark_ID = top5_ids[j],
                                     ES_reported = pj$es_reported,
                                     ES_recomputed = round(pj$es_calc, 4)))
  nm <- sprintf("GSEA_runningscore_%d_%s", j,
                gsub("[^A-Za-z0-9]+", "_", sub("^HALLMARK_", "", top5_ids[j])))
  save_fig(pj$plot, nm, 132, 112)
  plist[[length(plist) + 1]] <- pj$plot
}
if (nrow(es_chk) > 0) {
  es_chk$Diff <- with(es_chk, round(ES_recomputed - ES_reported, 4))
  write.csv(es_chk, file.path(out_dir, "GSEA_ES_recompute_check.csv"), row.names = FALSE)
  msg("  ES check (recomputed vs fgsea): max |diff| = ", max(abs(es_chk$Diff)))
}

if (length(plist) >= 2) {
  p_comb <- patchwork::wrap_plots(plist, ncol = 2)
  save_fig(p_comb, "GSEA_runningscore_top5", 180, 300)
  ggsave(file.path(res_root, "GSEA_runningscore_top5.pdf"), p_comb, width = 180/25.4,
         height = 300/25.4, device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(file.path(res_root, "GSEA_runningscore_top5.png"), p_comb, width = 180/25.4,
         height = 300/25.4, dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

## ---------------------------------------------------------------------------
## 11. 图注
## ---------------------------------------------------------------------------
leg <- file.path(out_dir, "GSEA_figure_legends.txt")
con2 <- file(leg, open = "wt", encoding = "UTF-8")
w2 <- function(...) writeLines(paste0(...), con2)

w2("Figure legends - GSEA against MSigDB Hallmark gene sets (", ACCESSION, ")")
w2("")
w2("Common methods. Gene set enrichment analysis (Subramanian et al., PNAS 2005) was run with")
w2("clusterProfiler ", as.character(packageVersion("clusterProfiler")), " / fgsea ",
   as.character(packageVersion("fgsea")), " on the preranked list of all ",
   length(r1$gl), " tested genes, ordered by log2 fold change (Disease vs Control).")
w2("Gene identifiers were converted from HGNC symbols to Entrez Gene IDs before ranking;")
w2("duplicated identifiers were collapsed to the entry with the largest |log2FC|.")
w2("Gene sets were the 50 MSigDB Hallmark (H) sets (msigdbr ",
   as.character(packageVersion("msigdbr")), "; Homo sapiens), restricted to ",
   MIN_GS, "-", MAX_GS, " genes. The weighted enrichment score used exponent = 1;")
w2("significance was assessed by gene permutation with the multilevel algorithm (eps = 0,")
w2("seed = 123) and adjusted by Benjamini-Hochberg. NES > 0 indicates enrichment at the top of")
w2("the ranking (pathway activated in disease); NES < 0 indicates enrichment at the bottom")
w2("(pathway suppressed). ", nrow(sig1), " of ", nrow(r1$all), " sets reached padj < ", PADJ_CUT, ".")
w2("")
w2("Figure 1 | Hallmark pathways most strongly enriched in failing myocardium.")
w2("Bubble plot of the ", nrow(top5), " most significant Hallmark sets. Horizontal position is the NES;")
w2("colour encodes NES on a diverging scale centred at zero (red, activated; blue, suppressed);")
w2("bubble area is -log10(adjusted P). Numbers next to bubbles are NES values. The dashed line")
w2("marks NES = 0.")
w2("")
w2("Figure 2 | Direction-resolved overview of Hallmark enrichment.")
w2("Top five activated (NES > 0) and top five suppressed (NES < 0) sets among the ", nrow(sig1))
w2("sets with padj < ", PADJ_CUT, "; all other encodings as in Figure 1. Pathways are ordered by NES.")
w2("")
w2("Figure 3 | Running enrichment score of the five most significant Hallmark pathways.")
w2("For each set: (top) the weighted running-sum enrichment score across the ranked gene list,")
w2("with the position of the maximum ES marked; (middle) the rank positions of the set members")
w2("(leading-edge genes in the darker shade); (bottom) the ranked log2FC metric. Positive running")
w2("scores indicate accumulation of set members at the top of the ranking. Each set is also")
w2("supplied as an individual file (GSEA_runningscore_<n>_<pathway>.pdf/png).")
close(con2)

## ---------------------------------------------------------------------------
## 12. sessionInfo
## ---------------------------------------------------------------------------
si <- file.path(rep_dir, "sessionInfo_08_gsea.txt")
sink(si); cat("Script: 08_gsea.R\nRun at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n"); print(sessionInfo()); sink()

msg("\n== Done ==  output: ", out_dir)
