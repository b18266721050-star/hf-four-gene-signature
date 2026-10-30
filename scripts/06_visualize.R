# =============================================================================
# 06_visualize.R  ——  GSE57338 差异分析可视化（Nature 风格火山图 + 热图）
#
# 输入（results/DEG/）：
#   GSE57338_DEGall.csv           全部基因结果（含 log2FoldChange / padj / Regulation）
# 输入（clean_data/）：
#   GSE57338_clean_expression_matrix.csv   log2 表达矩阵（热图、PCA 用）
#   GSE57338_sample_info.csv               样本信息
#
# 输出（results/DEG/）：
#   GSE57338_volcano.pdf / .png            火山图（主图，单栏 89 mm）
#   GSE57338_heatmap.pdf / .png            Top 50 差异基因热图（双栏 183 mm）
#   GSE57338_QC_diagnostics.pdf / .png     诊断图（PCA + P 值分布，附图）
#
# 用法（在项目根目录下，先跑 04_deg.R）：
#   Rscript scripts/06_visualize.R
#
# -----------------------------------------------------------------------------
# 作图约定（Nature 投稿风格）
# -----------------------------------------------------------------------------
#   1) 尺寸按 Nature 栏宽：单栏 89 mm，双栏 183 mm；字号 5–8 pt，
#      PDF 为矢量（cairo_pdf，字体嵌入），PNG 为 600 dpi 位图预览。
#   2) 无背景网格、无多余边框，只保留左侧/底部轴线（theme_classic）；
#      轴线与刻度 0.35 pt，接近 Nature 要求的细线。
#   3) 配色不用彩虹色：上调砖红 #B2182B，下调钢蓝 #2166AC，不显著浅灰 #C9C9C9
#      （红蓝发散方案对红绿色觉缺陷也友好；热图用同一套蓝-白-红）。
#   4) 火山图 y 轴用 -log10(padj)（与显著性判据同一尺度，阈值线位置自洽），
#      横轴 log2FC；标签用 ggrepel，max.overlaps = Inf 防止标签被静默丢弃。
#   5) 热图按基因做 z-score（scale = "row"）以突出变化模式；样本**不做**列聚类
#      全混，而是"组内聚类 + 组间分块"，避免 313 个样本的树状图糊成一团。
# =============================================================================

rm(list = ls())

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
clean_dir <- file.path(proj_dir, "clean_data")
res_dir   <- file.path(proj_dir, "results", "DEG")
rep_dir   <- file.path(proj_dir, "report")
if (!dir.exists(res_dir)) stop("results/DEG not found. Run 04_deg.R first.")

## ---------------------------------------------------------------------------
## 1. 参数
## ---------------------------------------------------------------------------
ACCESSION <- "GSE57338"
LFC_CUT   <- 0.585        # 与 04_deg.R 保持一致
PADJ_CUT  <- 0.05
N_LABEL   <- 5            # 火山图每侧标注的基因数
N_HEAT    <- 25           # 热图每侧取的基因数（共 50）
TOP_VAR   <- 2000         # PCA 用方差最大的前 N 个基因
DPI       <- 600

# 配色（全脚本统一，图之间视觉一致）
COL_UP   <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_NS   <- "#C9C9C9"
COL_GRP  <- c(Control = "#9E9E9E", Disease = "#B2182B")
# 注意：subtype 的因子水平是 Control / DCM / ICM（对照在亚型里同样叫 Control），
# 配色向量的名字必须与因子水平逐一对应，否则 pheatmap 会报错
COL_SUB  <- c(Control = "#D9D9D9", DCM = "#F4A582", ICM = "#92C5DE")

## ---------------------------------------------------------------------------
## 2. 包与主题
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
  library(patchwork)
  library(grid)
  library(scales)
})

# Nature 风格主题：细轴线、无网格、黑色文字、小字号
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(1.2, "pt"),
      axis.text        = element_text(colour = "black", size = base_size),
      axis.title       = element_text(colour = "black", size = base_size + 0.5),
      legend.title     = element_text(size = base_size - 0.5),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(6, "pt"),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey35", hjust = 0),
      plot.margin      = margin(3, 3, 3, 3, "pt")
    )
}
theme_set(theme_nature())

# 同时输出 PDF（矢量）与 PNG（600 dpi）
save_both <- function(plot, name, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggsave(file.path(res_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(res_dir, paste0(name, ".png")), plot,
         width = w, height = h, dpi = DPI, device = ragg::agg_png, bg = "white")
  invisible(TRUE)
}

# grid 类图形（ComplexHeatmap 对象）在下面 5. 节用专用函数保存

## ---------------------------------------------------------------------------
## 3. 读入
## ---------------------------------------------------------------------------
deg  <- fread(file.path(res_dir, paste0(ACCESSION, "_DEGall.csv")), data.table = FALSE)
meta <- read.csv(file.path(clean_dir, paste0(ACCESSION, "_sample_info.csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)
expr_dt <- fread(file.path(clean_dir, paste0(ACCESSION, "_clean_expression_matrix.csv")),
                 data.table = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

# 与 04_deg.R 相同的样本顺序校验
if (!identical(colnames(expr), meta$sample)) {
  meta <- meta[match(colnames(expr), meta$sample), , drop = FALSE]
}
meta$group2 <- factor(ifelse(meta$group == "NF", "Control", "Disease"),
                      levels = c("Control", "Disease"))
meta$subtype <- factor(ifelse(meta$group == "NF", "Control", meta$group),
                       levels = c("Control", "DCM", "ICM"))

# 统一 Regulation 因子顺序（图例顺序：Up / NS / Down）
deg$Regulation <- factor(deg$Regulation, levels = c("Up", "NS", "Down"))
n_up   <- sum(deg$Regulation == "Up")
n_down <- sum(deg$Regulation == "Down")

## ---------------------------------------------------------------------------
## 4. 火山图
## ---------------------------------------------------------------------------
# y = -log10(padj)；padj 理论最小值为 0，取最小正数截断避免 Inf
min_padj <- min(deg$padj[deg$padj > 0], na.rm = TRUE)
deg$neglog10_padj <- -log10(pmax(deg$padj, min_padj))

# 标注基因：上/下调各取 padj 最小的 N_LABEL 个
up_df   <- deg[deg$Regulation == "Up", , drop = FALSE]
down_df <- deg[deg$Regulation == "Down", , drop = FALSE]
lab_df  <- rbind(head(up_df[order(up_df$padj), ], N_LABEL),
                 head(down_df[order(down_df$padj), ], N_LABEL))

y_max  <- max(deg$neglog10_padj, na.rm = TRUE)
x_lim  <- max(abs(deg$log2FoldChange), na.rm = TRUE) * 1.05
y_line <- -log10(PADJ_CUT)

p_vol <- ggplot(deg, aes(x = log2FoldChange, y = neglog10_padj)) +
  # 阈值参考线
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed",
             linewidth = 0.25, colour = "grey55") +
  geom_hline(yintercept = y_line, linetype = "dashed",
             linewidth = 0.25, colour = "grey55") +
  # 散点：先画不显著，再画显著，保证显著点在上层
  geom_point(data = deg[deg$Regulation == "NS", ], colour = COL_NS,
             size = 0.55, stroke = 0, alpha = 0.55) +
  geom_point(data = deg[deg$Regulation != "NS", ],
             aes(colour = Regulation), size = 0.9, stroke = 0, alpha = 0.95) +
  scale_colour_manual(
    values = c(Up = COL_UP, NS = COL_NS, Down = COL_DOWN),
    breaks = c("Up", "NS", "Down"),
    labels = c(paste0("Up (", n_up, ")"), paste0("NS (", nrow(deg) - n_up - n_down, ")"),
               paste0("Down (", n_down, ")")),
    name = NULL, drop = FALSE
  ) +
  # 基因标签
  geom_text_repel(data = lab_df, aes(label = GeneSymbol, colour = Regulation),
                  size = 2.1, fontface = "italic", show.legend = FALSE,
                  max.overlaps = Inf, box.padding = 0.18, point.padding = 0.10,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey40", seed = 42) +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(0, y_max * 1.06),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(
    x = expression(log[2]~"fold change (Disease / Control)"),
    y = expression(-log[10]~italic(P)[adj]),
    title = "GSE57338: Disease vs Control",
    subtitle = paste0("limma | |log2FC| > ", LFC_CUT, " & padj < ", PADJ_CUT,
                      " | n = 313 (177 disease / 136 control)")
  ) +
  theme(legend.position = c(0.015, 0.02),
        legend.justification = c(0, 0),
        legend.direction = "vertical",
        legend.spacing.y = unit(1, "pt")) +
  guides(colour = guide_legend(override.aes = list(size = 1.6)))

# 只在右下角标注校正后 P 的阈值线（阈值数值在副标题里已给出，避免重复堆文字）
p_vol <- p_vol +
  annotate("text", x = x_lim * 0.97, y = y_line * 1.6,
           label = paste0("padj = ", PADJ_CUT), hjust = 1, vjust = 0,
           size = 1.9, colour = "grey40")

save_both(p_vol, paste0(ACCESSION, "_volcano"), 89, 78)

## ---------------------------------------------------------------------------
## 5. 热图（Top 50 差异基因，ComplexHeatmap）
## ---------------------------------------------------------------------------
sig <- deg[deg$Regulation != "NS", , drop = FALSE]
sig_up   <- sig[sig$Regulation == "Up", , drop = FALSE]
sig_down <- sig[sig$Regulation == "Down", , drop = FALSE]
top_up   <- head(sig_up[order(sig_up$padj), ], N_HEAT)
top_down <- head(sig_down[order(sig_down$padj), ], N_HEAT)

# 下调在前、上调在后，颜色块在图上就是"左蓝右红 / 左红右蓝"两段，一眼可分
heat_df <- rbind(top_down, top_up)
heat_df <- heat_df[heat_df$GeneSymbol %in% rownames(expr), , drop = FALSE]

genes <- heat_df$GeneSymbol
hm    <- expr[genes, , drop = FALSE]
ok_row <- apply(hm, 1, sd) > 0            # 去掉零方差行，避免 z-score 出现 Inf
hm <- hm[ok_row, , drop = FALSE]
heat_df <- heat_df[ok_row, , drop = FALSE]

# 行 z-score：突出每个基因在各样本间的相对高低（不保留绝对表达水平）
hm_z <- t(scale(t(hm)))
hm_z[hm_z > 2]  <- 2      # 截断极端值，避免个别样本把色标拉满
hm_z[hm_z < -2] <- -2

# 列顺序：对照组内部聚类 → 疾病组内部聚类，两组之间由 column_split 留出空隙
ord_within <- function(idx) {
  sub <- hm_z[, idx, drop = FALSE]
  if (ncol(sub) < 3) return(idx)
  # 313 个样本列用欧氏距离聚类计算量大且易被离群样本带偏，改用 1-correlation
  hc <- hclust(as.dist(1 - cor(sub, method = "pearson")), method = "average")
  idx[hc$order]
}
idx_ctrl <- which(meta$group2 == "Control")
idx_dis  <- which(meta$group2 == "Disease")
col_ord  <- c(ord_within(idx_ctrl), ord_within(idx_dis))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

# 顶部注释条：分组 + 亚型（细条，2 mm 高，符合 Nature 的紧凑版式）
ha <- HeatmapAnnotation(
  Group   = meta$group2[col_ord],
  Subtype = meta$subtype[col_ord],
  col = list(Group = COL_GRP, Subtype = COL_SUB),
  annotation_name_gp = gpar(fontsize = 5.2),
  annotation_name_side = "left",
  simple_anno_size = unit(2.0, "mm"),
  gap = unit(0.5, "mm"),
  annotation_legend_param = list(
    Group   = list(title_gp = gpar(fontsize = 5.2, fontface = "plain"),
                   labels_gp = gpar(fontsize = 5),
                   grid_height = unit(2.2, "mm"), grid_width = unit(2.2, "mm")),
    Subtype = list(title_gp = gpar(fontsize = 5.2, fontface = "plain"),
                   labels_gp = gpar(fontsize = 5),
                   grid_height = unit(2.2, "mm"), grid_width = unit(2.2, "mm"))
  )
)

ht <- Heatmap(
  hm_z[, col_ord, drop = FALSE],
  name = "Row z-score",
  col  = colorRamp2(c(-2, 0, 2), c(COL_DOWN, "white", COL_UP)),
  cluster_rows = FALSE,                # 行已按 log2FC 排序，保持可解释的顺序
  cluster_columns = FALSE,             # 列顺序见上面的组内聚类
  row_split = factor(heat_df$Regulation, levels = c("Down", "Up")),
  column_split = factor(meta$group2[col_ord], levels = c("Control", "Disease")),
  row_title = NULL,
  column_title_gp = gpar(fontsize = 5.8),
  row_gap = unit(1.2, "mm"),
  column_gap = unit(1.2, "mm"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 4.8, fontface = "italic"),
  top_annotation = ha,
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5.2, fontface = "plain"),
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(14, "mm"),
    grid_width = unit(2.2, "mm"),
    at = c(-2, -1, 0, 1, 2),
    border = NA
  ),
  width  = unit(142, "mm"),
  height = unit(88, "mm")
)

# 保存：ComplexHeatmap 用 grid 绘制，同样输出 PDF + 600 dpi PNG
save_heatmap <- function(ht_obj, name, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  cairo_pdf(file.path(res_dir, paste0(name, ".pdf")), width = w, height = h, bg = "white")
  draw(ht_obj, merge_legend = TRUE, ht_gap = unit(1.5, "mm"),
       column_title = paste0("Top ", nrow(hm_z), " DE genes (",
                             sum(heat_df$Regulation == "Down"), " down / ",
                             sum(heat_df$Regulation == "Up"), " up)"),
       column_title_gp = gpar(fontsize = 6.5, fontface = "bold"))
  dev.off()
  ragg::agg_png(file.path(res_dir, paste0(name, ".png")),
                width = w, height = h, units = "in", res = DPI, background = "white")
  draw(ht_obj, merge_legend = TRUE, ht_gap = unit(1.5, "mm"),
       column_title = paste0("Top ", nrow(hm_z), " DE genes (",
                             sum(heat_df$Regulation == "Down"), " down / ",
                             sum(heat_df$Regulation == "Up"), " up)"),
       column_title_gp = gpar(fontsize = 6.5, fontface = "bold"))
  dev.off()
  invisible(TRUE)
}

save_heatmap(ht, paste0(ACCESSION, "_heatmap"), 183, 104)

## ---------------------------------------------------------------------------
## 6. 诊断图：PCA + P 值分布
## ---------------------------------------------------------------------------
# 6.1 PCA：全部基因中取方差最大的 TOP_VAR 个（log2 表达已是归一化值，无需再变换）
v  <- rowVars(expr)
top_var_genes <- order(v, decreasing = TRUE)[seq_len(min(TOP_VAR, length(v)))]
pca <- prcomp(t(expr[top_var_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
ve  <- 100 * (pca$sdev^2) / sum(pca$sdev^2)

pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  Group = meta$group2, Subtype = meta$subtype
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, fill = Group)) +
  geom_point(shape = 21, size = 1.15, stroke = 0.15, colour = "white", alpha = 0.9) +
  stat_ellipse(aes(colour = Group), type = "norm", linetype = "dashed",
               linewidth = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = COL_GRP) +
  scale_colour_manual(values = COL_GRP) +
  labs(x = sprintf("PC1 (%.1f%% variance)", ve[1]),
       y = sprintf("PC2 (%.1f%% variance)", ve[2]),
       title = "PCA (top 2000 variable genes)") +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 1.8)))

# 6.2 P 值分布：模型是否校准（理想形态 = 均匀分布 + 接近 0 处的峰）
y_null <- nrow(deg) / 50
p_hist <- ggplot(deg, aes(x = P.Value)) +
  geom_histogram(aes(y = after_stat(count)), bins = 50,
                 fill = "grey35", colour = "white", linewidth = 0.15) +
  geom_hline(yintercept = y_null, linetype = "dashed",
             linewidth = 0.3, colour = COL_UP) +
  annotate("text", x = 0.55, y = y_null * 2.2,
           label = "uniform null", size = 1.9, colour = "grey40", hjust = 0) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "Raw P value", y = "Gene count",
       title = "P value distribution (all 18,726 genes)") +
  theme(legend.position = "none")

p_diag <- (p_pca | p_hist) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "",
                  theme = theme(plot.tag = element_text(size = 8, face = "bold"))) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_both(p_diag, paste0(ACCESSION, "_QC_diagnostics"), 183, 80)

## ---------------------------------------------------------------------------
## 7. 图注素材（写入 txt，方便直接粘到论文/报告）
## ---------------------------------------------------------------------------
caption <- c(
  "Figure legends",
  "==============",
  "",
  "Figure 1 | Volcano plot of differential expression between failing (DCM + ICM)",
  "and non-failing (NF) human left ventricular myocardium (GSE57338; 177 vs 136).",
  "Linear models were fitted with limma and moderated t-statistics were obtained",
  "with empirical Bayes shrinkage (trend = TRUE, robust = TRUE), adjusting for sex",
  "and age. Dashed lines mark |log2FC| = 0.585 (fold change 1.5) and padj = 0.05",
  "(Benjamini-Hochberg). Genes are coloured red (up-regulated in disease) or blue",
  "(down-regulated); grey denotes non-significant genes. The top five genes per",
  "direction (by adjusted P) are labelled.",
  "",
  "Figure 2 | Heatmap of the 50 most significant differentially expressed genes",
  "(25 up, 25 down). Values are row-wise z-scores of log2 expression. Samples are",
  "clustered within group and separated by a vertical gap (control left, disease",
  "right); the horizontal gap separates down- from up-regulated genes.",
  "",
  "Supplementary Figure | (a) Principal component analysis of the 2,000 most",
  "variable genes; (b) distribution of raw P values across all tested genes,",
  "expected to be uniform apart from an enrichment near zero under a",
  "well-calibrated null model.",
  "",
  "Colours: up #B2182B, down #2166AC, non-significant #C9C9C9.",
  "Output: vector PDF + 600 dpi PNG."
)
con <- file(file.path(res_dir, paste0(ACCESSION, "_figure_legends.txt")), "w", encoding = "UTF-8")
writeLines(caption, con)
close(con)

cat("Done. Figures written to:", res_dir, "\n")
cat("  volcano  : 89 x 78 mm\n")
cat("  heatmap  : 183 x 120 mm\n")
cat("  QC diag  : 183 x 80 mm\n")
