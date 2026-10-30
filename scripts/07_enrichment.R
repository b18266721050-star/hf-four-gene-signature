# =============================================================================
# 07_enrichment.R  ——  GSE57338 显著差异基因的 GO / KEGG 富集分析（ORA）
#
# 输入（results/DEG/）：
#   GSE57338_DEG_significant.csv   显著差异基因（269 个，|log2FC|>0.585 且 padj<0.05）
#   GSE57338_DEGall.csv            全部检验基因（18,726 行）—— 用作背景 universe
#
# 输出（results/Enrich/）：
#   GO_BP.csv / GO_CC.csv / GO_MF.csv / KEGG.csv   富集结果（padj<0.05）
#   另存 results/ 根目录下的 GOBP.csv / GOCC.csv / GO_MF.csv / KEGG.csv 副本
#   GO_*_simplified.csv          GO 语义去冗余后的补充表（simplify, Wang 相似度）
#   GO_dotplot.pdf/.png          GO 三层面气泡图（BP/CC/MF，各 top 10）
#   KEGG_dotplot.pdf/.png        KEGG 气泡图（top 10）
#   enrichment_summary.txt       参数、基因数、条目数、top10 汇总
#   kegg_hsa_<date>.gson         KEGG 通路集快照（固定版本，保证可复现）
#   id_mapping.csv               基因 ID 转换明细
# 输出（report/）：
#   sessionInfo_07_enrich.txt
#
# 用法（在项目根目录下）：
#   Rscript scripts/07_enrichment.R
#
# -----------------------------------------------------------------------------
# 方法学要点（为什么这样设参数）
# -----------------------------------------------------------------------------
#   1) 物种：GSE57338 平台为 GPL11532（Affymetrix HuGene-1.1-ST），样本为人
#      左心室心肌 -> Homo sapiens，注释库用 org.Hs.eg.db（Entrez ID）。
#   2) 富集方法：ORA（过表征分析），clusterProfiler::enrichGO / enrichKEGG，
#      单侧超几何检验（等价于单侧 Fisher 精确检验）。
#   3) 背景 universe：**必须**设为"本研究中真正被检验过的基因"（18,726 个），
#      而不是全基因组（~20,000+ 注释基因）。用全基因组作背景会把组织不表达、
#      低表达基因算进分母，使组织特异 / 低表达基因相关的条目假性显著。
#      前景与背景使用同一套 SYMBOL→ENTREZID 映射，避免 ID 层面不一致。
#   4) 富集函数的 pvalueCutoff 过滤的是 **p.adjust**（不是原始 p），且 qvalueCutoff
#      是第二道闸。为严格遵守"padj<0.05"这一唯一判据，这里把两个 cutoff 都设为 1，
#      取出全部条目后在 R 里显式按 p.adjust < 0.05 过滤 —— 阈值透明、可复核。
#   5) 多重检验：BH（Benjamini-Hochberg）控制 FDR。GO 条目间因 DAG 真路径传播
#      而正相关，BH 在此条件下仍然有效（Bonferroni 过于保守）。
#   6) GO 冗余：GO 的父子层级会让同一信号点亮一整条祖先谱系，因此额外输出
#      simplify()（Wang 语义相似度，cutoff=0.7）去冗余后的补充表。
#   7) KEGG 是**在线数据库**，enrichKEGG 每次调用都查 KEGG REST API，不同时间
#      跑结果会变。脚本会把当次 KEGG 通路集写成 gson 快照并记录访问日期，
#      保证结果可追溯；这是"结果可复现"的必要条件。
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
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  return(getwd())
}

wd <- getwd()
proj_dir <- if (dir.exists(file.path(wd, "clean_data"))) normalizePath(wd) else dirname(get_this_dir())
deg_dir  <- file.path(proj_dir, "results", "DEG")
res_root <- file.path(proj_dir, "results")
out_dir  <- file.path(res_root, "Enrich")
rep_dir  <- file.path(proj_dir, "report")

for (d in c(out_dir, rep_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function(ch = "-", n = 74) cat(strrep(ch, n), "\n", sep = "")

## ---------------------------------------------------------------------------
## 1. 参数
## ---------------------------------------------------------------------------
ACCESSION <- "GSE57338"
ORGANISM  <- "Homo sapiens"          # 人；平台 GPL11532，组织 = 左心室心肌
ORGDB     <- "org.Hs.eg.db"
KEGG_CODE <- "hsa"                   # KEGG 物种代码：hsa = Homo sapiens

PADJ_CUT  <- 0.05                    # 唯一显著性判据
ADJ_METHOD <- "BH"
MIN_GS <- 10                         # 过小的基因集 p 值不稳定
MAX_GS <- 500                        # 过大的基因集几乎必然"显著"，信息量低
TOP_N  <- 10                         # 气泡图每类展示条目数
SIMPLIFY_CUTOFF <- 0.7               # GO 语义去冗余阈值
SIMP_MEASURE <- "Wang"               # 图拓扑相似度，跨注释版本稳定
SIMPLIFY_MAX <- 600                  # 去冗余最多处理的条目数（O(n^2) 计算，需限流）

FONT_FAMILY <- "Arial"               # Nature 惯用无衬线字体

## ---------------------------------------------------------------------------
## 2. 读入差异基因表
## ---------------------------------------------------------------------------
sig_file <- file.path(deg_dir, paste0(ACCESSION, "_DEG_significant.csv"))
all_file <- file.path(deg_dir, paste0(ACCESSION, "_DEGall.csv"))
if (!file.exists(sig_file)) stop("DEG significant table not found: ", sig_file)
if (!file.exists(all_file)) stop("DEG all table not found: ", all_file)

deg_sig <- read.csv(sig_file, check.names = FALSE)
deg_all <- read.csv(all_file, check.names = FALSE)

msg("== Input ==")
msg("significant DEGs : ", nrow(deg_sig),
    "  (Up ", sum(deg_sig$Regulation == "Up"), " / Down ", sum(deg_sig$Regulation == "Down"), ")")
msg("tested genes     : ", nrow(deg_all))

## ---------------------------------------------------------------------------
## 3. ID 转换（SYMBOL -> ENTREZID）
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

fg_symbols <- unique(deg_sig$GeneSymbol)
# universe = 真正被检验过的基因（padj 非 NA 的行），不是全基因组
bg_symbols <- unique(deg_all$GeneSymbol[!is.na(deg_all$P.Value)])

fg_map <- suppressWarnings(bitr(fg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
bg_map <- suppressWarnings(bitr(bg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))

fg_map$set <- "foreground"
bg_map$set <- "universe"
write.csv(rbind(fg_map, bg_map), file.path(out_dir, "id_mapping.csv"), row.names = FALSE)

gene_fg <- unique(fg_map$ENTREZID)      # 去重：一对多映射会虚增 Count
gene_bg <- unique(bg_map$ENTREZID)

conv_rate <- length(gene_fg) / length(fg_symbols)
msg("\n== ID mapping ==")
msg("foreground: ", length(fg_symbols), " symbols -> ", length(gene_fg),
    " Entrez IDs (conversion ", sprintf("%.1f%%", conv_rate * 100), ")")
msg("universe  : ", length(bg_symbols), " symbols -> ", length(gene_bg), " Entrez IDs")
if (conv_rate < 0.85) warning("ID conversion loss > 15%; foreground may be too eroded to trust.")

## ---------------------------------------------------------------------------
## 4. GO 富集（BP / CC / MF 分别做）
## ---------------------------------------------------------------------------
run_go <- function(ont) {
  ego <- enrichGO(gene          = gene_fg,
                  universe      = gene_bg,          # 关键：实验可检测到的基因集
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont,              # 必须显式指定（默认是 MF）
                  pAdjustMethod = ADJ_METHOD,
                  pvalueCutoff  = 1,                # 先全取，随后按 p.adjust 显式过滤
                  qvalueCutoff  = 1,
                  minGSSize     = MIN_GS,
                  maxGSSize     = MAX_GS,
                  readable      = TRUE)             # ENTREZ -> SYMBOL，输出可读
  ego
}

ratio2num <- function(x) {
  s <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(s, function(y) as.numeric(y[1]) / as.numeric(y[2]), numeric(1))
}

tidy_result <- function(obj, ontology) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(df)
  df$FoldEnrichment <- round(ratio2num(df$GeneRatio) / ratio2num(df$BgRatio), 3)
  df$GeneRatioNum   <- round(ratio2num(df$GeneRatio), 3)
  df$Ontology       <- ontology
  df <- df[, c("ID", "Description", "Ontology", "GeneRatio", "GeneRatioNum",
               "BgRatio", "FoldEnrichment", "pvalue", "p.adjust", "qvalue",
               "geneID", "Count")]
  df <- add_direction(df)
  df <- df[order(df$p.adjust, df$pvalue), ]
  rownames(df) <- NULL
  df
}

# ORA 本身是"方向不可知"的：它只问这些基因是否在条目里过表征，不问上调还是下调。
# 把前景基因里的 up/down 构成一并算出来，读者才能判断某条目是被"上调驱动"还是
# 被"下调驱动"——否则很容易把"炎症相关条目显著"误读成"炎症通路被激活"。
add_direction <- function(df) {
  if (nrow(df) == 0) return(df)
  fc <- stats::setNames(deg_sig$log2FoldChange, deg_sig$GeneSymbol)
  gl <- strsplit(as.character(df$geneID), "/", fixed = TRUE)
  df$N_Up <- vapply(gl, function(g) sum(fc[g] > 0, na.rm = TRUE), numeric(1))
  df$N_Down <- vapply(gl, function(g) sum(fc[g] < 0, na.rm = TRUE), numeric(1))
  df$Median_log2FC <- round(vapply(gl, function(g) stats::median(fc[g], na.rm = TRUE),
                                   numeric(1)), 3)
  df$Dominant <- ifelse(df$N_Up >= df$N_Down, "Up", "Down")
  df
}

go_res <- list()
msg("\n== GO enrichment (padj < ", PADJ_CUT, ") ==")
for (ont in c("BP", "CC", "MF")) {
  t0 <- Sys.time()
  ego_full <- run_go(ont)
  df_all   <- tidy_result(ego_full, ont)
  df_sig   <- df_all[df_all$p.adjust < PADJ_CUT, ]
  msg(sprintf("  [%s] enrich done: %d tested / %d sig  (%.1fs)",
              ont, nrow(df_all), nrow(df_sig),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  go_res[[ont]] <- list(obj = ego_full, all = df_all, sig = df_sig)

  # 语义去冗余（补充表）：simplify 只在单一本体内部有意义
  # 注意：simplify 的复杂度是 O(n^2) 的语义相似度矩阵，若对全部条目（BP 可达
  #       数千条）做会跑上数小时。这里只对**显著条目**做，并限制在 SIMPLIFY_MAX 内。
  df_simp <- tryCatch({
    keep <- df_all$p.adjust < PADJ_CUT
    if (sum(keep) >= 2) {
      ego_sub <- ego_full
      res_sub <- ego_full@result[keep, , drop = FALSE]
      res_sub <- res_sub[order(res_sub$p.adjust), , drop = FALSE]
      if (nrow(res_sub) > SIMPLIFY_MAX) res_sub <- res_sub[seq_len(SIMPLIFY_MAX), , drop = FALSE]
      res_sub$ID <- as.character(res_sub$ID)
      rownames(res_sub) <- res_sub$ID
      ego_sub@result <- res_sub

      es <- simplify(ego_sub, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust",
                     select_fun = min, measure = SIMP_MEASURE)
      d <- tidy_result(es, ont)
      d[d$p.adjust < PADJ_CUT, , drop = FALSE]
    } else {
      NULL
    }
  }, error = function(e) {
    msg("  [warn] simplify failed for ", ont, ": ", conditionMessage(e))
    NULL
  })
  go_res[[ont]]$simp <- df_simp

  msg(sprintf("  %s: %4d terms tested, %3d significant (padj<%.2f)%s",
              ont, nrow(df_all), nrow(df_sig), PADJ_CUT,
              if (is.null(df_simp)) "" else sprintf(", %d after simplify", nrow(df_simp))))
}

## ---------------------------------------------------------------------------
## 5. KEGG 富集 + 快照固定版本
## ---------------------------------------------------------------------------
msg("\n== KEGG enrichment (padj < ", PADJ_CUT, ") ==")
kegg_full <- tryCatch(
  enrichKEGG(gene         = gene_fg,
             organism     = KEGG_CODE,
             keyType      = "ncbi-geneid",
             universe     = gene_bg,
             pvalueCutoff = 1,
             qvalueCutoff = 1,
             pAdjustMethod = ADJ_METHOD,
             minGSSize    = MIN_GS,
             maxGSSize    = MAX_GS),
  error = function(e) stop("enrichKEGG failed (network?): ", conditionMessage(e))
)
kegg_full <- setReadable(kegg_full, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_all  <- tidy_result(kegg_full, "KEGG")
kegg_sig  <- kegg_all[kegg_all$p.adjust < PADJ_CUT, ]
msg(sprintf("  KEGG: %d pathways tested, %d significant (padj<%.2f)",
            nrow(kegg_all), nrow(kegg_sig), PADJ_CUT))

# KEGG 在线库会随时间变化 —— 存一份快照并记录访问日期
kegg_date <- format(Sys.Date(), "%Y%m%d")
kegg_pin  <- tryCatch({
  g <- gson_KEGG(species = KEGG_CODE, keyType = "kegg")
  g@accessed_date <- as.character(Sys.Date())
  gson::write.gson(g, file.path(out_dir, paste0("kegg_", KEGG_CODE, "_", kegg_date, ".gson")))
  as.character(Sys.Date())
}, error = function(e) {
  msg("  [warn] KEGG snapshot failed: ", conditionMessage(e)); NA_character_
})

## ---------------------------------------------------------------------------
## 6. 输出 CSV
## ---------------------------------------------------------------------------
out_files <- c(BP = "GOBP.csv", CC = "GOCC.csv", MF = "GO_MF.csv")
for (ont in names(out_files)) {
  d <- go_res[[ont]]$sig
  write.csv(d, file.path(out_dir, out_files[[ont]]), row.names = FALSE)
  # 需求书另要求放在 results 根目录，这里同步一份副本
  write.csv(d, file.path(res_root, out_files[[ont]]), row.names = FALSE)

  s <- go_res[[ont]]$simp
  if (!is.null(s) && nrow(s) > 0) {
    write.csv(s, file.path(out_dir, paste0("GO_", ont, "_simplified.csv")), row.names = FALSE)
  }
}
write.csv(kegg_sig, file.path(out_dir, "KEGG.csv"), row.names = FALSE)
write.csv(kegg_sig, file.path(res_root, "KEGG.csv"), row.names = FALSE)

msg("\n== CSV written ==")
msg("  ", paste(c(out_files, "KEGG.csv"), collapse = ", "), "  ->  results/Enrich/ (+ results/)")

## ---------------------------------------------------------------------------
## 7. 可视化：Nature 风格气泡图（莫兰迪配色）
## ---------------------------------------------------------------------------
# 莫兰迪色系：低饱和、灰调；由浅到深对应 -log10(padj) 由小到大
MORANDI <- c("#E4DFD3", "#C3CCBF", "#9AAEA4", "#6F877D", "#42554C")

# 单一图内共用一个色标；跨图尺度差异极大（GO-CC 最小 padj 达 1e-26，而 BP 只有
# 1e-4 量级），若三张图硬套同一绝对尺度，BP 的点会全部挤在色标最浅端而看不清。
# 处理方式：每张图按自身取值范围定色标，并对极显著值做上限截断（squish），
# 截断处会在图注中明确写出（见 enrichment_figure_legends.txt）。
COLOR_CAP_GO   <- 10
COLOR_CAP_KEGG <- 4

make_color_scale <- function(v, cap, min_span = 2.2) {
  lo <- min(v)
  hi <- min(max(v), cap)
  if (hi - lo < min_span) hi <- lo + min_span
  scale_colour_gradientn(
    colours = MORANDI,
    limits  = c(lo, hi),
    oob     = scales::squish,
    name    = expression(bold(-log[10]~italic(p)[adj])),
    guide   = guide_colourbar(order = 1, ticks.colour = "white", frame.colour = NA)
  )
}

wrap_lab <- function(x, width = 46) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

nature_theme <- function(base_size = 6.5) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", linewidth = 0.35, fill = NA),
      axis.line          = element_blank(),
      axis.ticks         = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length  = unit(1.2, "mm"),
      axis.text          = element_text(colour = "black", size = base_size),
      axis.title         = element_text(colour = "black", size = base_size + 0.5),
      axis.title.x       = element_text(margin = margin(t = 2)),
      axis.title.y       = element_text(margin = margin(r = 3)),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", hjust = 0,
                                        size = base_size + 1,
                                        margin = margin(b = 1.5, t = 1)),
      legend.title       = element_text(size = base_size, face = "plain"),
      legend.text        = element_text(size = base_size - 0.5),
      legend.key.height  = unit(3.2, "mm"),
      legend.key.width   = unit(2.6, "mm"),
      legend.background  = element_blank(),
      legend.margin      = margin(0, 0, 0, 0),
      plot.margin        = margin(2, 2, 2, 2)
    )
}

make_dotdata <- function(df, ont_label) {
  d <- head(df, TOP_N)
  # 最显著的排在最上：y 因子顺序 = 由下往上，故按 p.adjust 降序排列
  d <- d[order(d$p.adjust, d$Count, decreasing = TRUE), ]
  d$mlog10p <- -log10(d$p.adjust)
  d$label   <- paste0(ont_label, "___", wrap_lab(d$Description))
  d
}

## ---- 7.1 GO 三层面气泡图（BP / CC / MF 竖向排列） ---------------------------
go_plot_df <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(ont) {
  d <- go_res[[ont]]$sig
  if (nrow(d) == 0) return(NULL)
  make_dotdata(d, ont)
}))
go_plot_df$panel <- factor(sub("___.*$", "", go_plot_df$label),
                           levels = c("BP", "CC", "MF"))
go_plot_df$label <- factor(go_plot_df$label, levels = go_plot_df$label)

panel_labs <- c(BP = "a  Biological process (GO-BP)",
                CC = "b  Cellular component (GO-CC)",
                MF = "c  Molecular function (GO-MF)")

p_go <- ggplot(go_plot_df, aes(x = GeneRatioNum, y = label)) +
  geom_point(aes(size = Count, colour = mlog10p), stroke = 0) +
  make_color_scale(go_plot_df$mlog10p, COLOR_CAP_GO) +
  scale_size_continuous(range = c(0.9, 4.2), name = "Gene count",
                        breaks = scales::pretty_breaks(n = 3)) +
  scale_y_discrete(labels = function(x) sub("^[A-Z]+___", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.06))) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y",
             labeller = as_labeller(panel_labs)) +
  labs(x = "Gene ratio", y = NULL) +
  nature_theme(6.5) +
  theme(panel.spacing.y = unit(2.2, "mm"))

## ---- 7.2 KEGG 气泡图 -------------------------------------------------------
p_kegg <- NULL
if (nrow(kegg_sig) > 0) {
  kd <- make_dotdata(kegg_sig, "KEGG")
  kd$label <- factor(kd$label, levels = kd$label)
  p_kegg <- ggplot(kd, aes(x = GeneRatioNum, y = label)) +
    geom_point(aes(size = Count, colour = mlog10p), stroke = 0) +
    make_color_scale(kd$mlog10p, COLOR_CAP_KEGG) +
    scale_size_continuous(range = c(0.9, 4.2), name = "Gene count",
                          breaks = scales::pretty_breaks(n = 3)) +
    scale_y_discrete(labels = function(x) sub("^[A-Z]+___", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    labs(x = "Gene ratio", y = NULL) +
    nature_theme(6.5)
}

save_fig <- function(p, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  msg("  saved: ", stem, ".pdf / .png  (", width_mm, " x ", height_mm, " mm, 600 dpi)")
}

msg("\n== Figures ==")
save_fig(p_go, "GO_dotplot", 183, 178)
if (!is.null(p_kegg)) save_fig(p_kegg, "KEGG_dotplot", 183, 82)

## ---- 7.3 图注（英文，Nature 投稿格式，可直接引用） --------------------------
leg <- file.path(out_dir, "enrichment_figure_legends.txt")
con2 <- file(leg, open = "wt", encoding = "UTF-8")
w2 <- function(...) writeLines(paste0(...), con2)

w2("Figure legends - GO / KEGG enrichment (GSE57338, ", nrow(deg_sig), " DEGs)")
w2("")
w2("Figure 1 | Gene Ontology over-representation analysis of the differentially expressed genes.")
w2("Dot plot of the ten most significant terms in each of the three ontologies: ",
   "(a) biological process, (b) cellular component and (c) molecular function. ",
   "The analysis used the ", nrow(deg_sig), " DEGs (Disease vs Control, |log2FC| > 0.585, ",
   "BH-adjusted P < 0.05; ", length(gene_fg), " of them mapped to Entrez Gene IDs) as the ",
   "foreground and the ", nrow(deg_all), " genes tested in the differential analysis (",
   length(gene_bg), " with Entrez IDs) as the background universe. ",
   "Dot size indicates the number of DEGs annotated to the term; ",
   "dot colour indicates the BH-adjusted P value on a -log10 scale ",
   "(darker = more significant), truncated at -log10(Padj) = ", COLOR_CAP_GO,
   " so that the far more significant extracellular-matrix terms of the cellular-component ",
   "panel do not compress the colour range of the other panels. ",
   "Terms are ordered by adjusted P value, most significant at the top. ",
   "Horizontal position is the gene ratio (annotated DEGs / all annotated DEGs). ",
   "Numbers of significant terms: BP ", nrow(go_res$BP$sig),
   ", CC ", nrow(go_res$CC$sig), ", MF ", nrow(go_res$MF$sig),
   ". GO terms are DAG-nested, so a single biological theme can occupy several rows; ",
   "semantically de-redundified lists are provided in GO_BP/CC/MF_simplified.csv.")
w2("")
if (!is.null(p_kegg)) {
  w2("Figure 2 | KEGG pathway over-representation analysis of the differentially expressed genes.")
  w2("Dot plot of the ten most significant KEGG pathways (of ", nrow(kegg_sig),
     " pathways with BH-adjusted P < 0.05). Foreground, background, dot size and ordering ",
     "are as in Figure 1. Dot colour indicates -log10(adjusted P), truncated at ",
     COLOR_CAP_KEGG, ". KEGG pathway definitions were retrieved from the KEGG REST API on ",
     ifelse(is.na(kegg_pin), "an unspecified date", kegg_pin),
     "; a frozen copy of the gene sets is supplied as kegg_", KEGG_CODE, "_", kegg_date,
     ".gson to make this panel exactly reproducible. Several hits are KEGG disease maps ",
     "(e.g. malaria, Chagas disease, systemic lupus erythematosus, rheumatoid arthritis) whose ",
     "shared content is the complement cascade, cytokine signalling and antigen presentation; ",
     "they are reported as inflammation modules, not as evidence of infection.")
}
close(con2)

## ---------------------------------------------------------------------------
## 8. 汇总文本
## ---------------------------------------------------------------------------
sum_txt <- file.path(out_dir, "enrichment_summary.txt")
con <- file(sum_txt, open = "wt", encoding = "UTF-8")
w <- function(...) writeLines(paste0(...), con)

w("GSE57338 GO / KEGG enrichment summary")
w("date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("organism: ", ORGANISM, " | OrgDb: ", ORGDB, " ", as.character(packageVersion(ORGDB)))
w("KEGG organism code: ", KEGG_CODE, " | KEGG accessed: ", ifelse(is.na(kegg_pin), "unknown", kegg_pin))
w("method: over-representation analysis (one-sided hypergeometric), BH adjusted p")
w("foreground: ", length(gene_fg), " Entrez IDs (", length(fg_symbols), " symbols -> ",
  sprintf("%.1f%%", conv_rate * 100), " mapped)")
w("universe  : ", length(gene_bg), " Entrez IDs (tested genes, NOT whole genome)")
w("cutoff    : p.adjust < ", PADJ_CUT, " | minGSSize = ", MIN_GS, " | maxGSSize = ", MAX_GS)
w("DEG direction: ", sum(deg_sig$Regulation == "Up"), " up / ",
  sum(deg_sig$Regulation == "Down"), " down",
  " (ORA is direction-agnostic; see N_Up / N_Down / Median_log2FC columns)")
w("")
for (ont in c("BP", "CC", "MF")) {
  d <- go_res[[ont]]
  w(sprintf("[GO-%s] tested=%d  significant=%d  after_simplify=%s",
            ont, nrow(d$all), nrow(d$sig),
            if (is.null(d$simp)) "NA" else nrow(d$simp)))
  if (nrow(d$sig) > 0) {
    tt <- head(d$sig, TOP_N)
    for (i in seq_len(nrow(tt))) {
      w(sprintf("  %2d. %-58s p.adj=%.2e  FE=%5.2f  n=%2d (up %2d / down %2d)",
                i, substr(tt$Description[i], 1, 58), tt$p.adjust[i],
                tt$FoldEnrichment[i], tt$Count[i], tt$N_Up[i], tt$N_Down[i]))
    }
  }
  w("")
}
w(sprintf("[KEGG] tested=%d  significant=%d", nrow(kegg_all), nrow(kegg_sig)))
if (nrow(kegg_sig) > 0) {
  tt <- head(kegg_sig, TOP_N)
  for (i in seq_len(nrow(tt))) {
    w(sprintf("  %2d. %-58s p.adj=%.2e  FE=%5.2f  n=%2d (up %2d / down %2d)",
              i, substr(tt$Description[i], 1, 58), tt$p.adjust[i],
              tt$FoldEnrichment[i], tt$Count[i], tt$N_Up[i], tt$N_Down[i]))
  }
}
w("")
w("Note: GO dotplot shows the 10 most significant terms per ontology as requested;")
w("      GO terms are DAG-nested, so the top list may contain one lineage several times.")
w("      See GO_*_simplified.csv for the semantic de-redundified short list.")
close(con)

## ---------------------------------------------------------------------------
## 9. 会话信息
## ---------------------------------------------------------------------------
si <- file.path(rep_dir, "sessionInfo_07_enrich.txt")
sink(si)
cat("Script: 07_enrichment.R\n")
cat("Run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

msg("\n== Done ==")
msg("output dir: ", out_dir)
