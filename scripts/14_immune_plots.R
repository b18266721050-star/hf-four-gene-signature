###############################################################################
# 04_immune_plots.R
# Figures + plot data for immune infiltration analysis
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(scales)
})

ROOT   <- "E:/another data file/GEObulk analysis project"
OUT    <- file.path(ROOT, "results", "Immune")
DIR_S  <- file.path(OUT, "02_ssgsea")
DIR_T  <- file.path(OUT, "04_stats")
DIR_P  <- file.path(OUT, "03_plots")
DIR_PD <- file.path(OUT, "06_plotdata")
DIR_LOG<- file.path(OUT, "00_logs")
for (d in c(DIR_P, DIR_PD)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sink(file.path(DIR_LOG, "04_plots_console.log"), split = TRUE)
cat("=== 04_immune_plots.R started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n\n")

st   <- readRDS(file.path(DIR_T, "immune_stats.rds"))
obj  <- readRDS(file.path(OUT, "01_qc", "qc_input.rds"))
es   <- readRDS(file.path(DIR_S, "GSE57338_immune_ssGSEA.rds"))$es
si   <- obj$si
expr <- obj$expr

cells <- rownames(es)
grp   <- factor(si$group, levels = c("NF", "DCM", "ICM"))
n_grp <- table(grp)
key_genes <- c("AQP3", "GGT5", "ISLR", "LAD1")

grp_cols <- c(NF = "#4C72B0", DCM = "#DD8452", ICM = "#55A868")

## helper: save both formats
save_fig <- function(p, name, width, height) {
  ggsave(file.path(DIR_P, paste0(name, ".pdf")), p, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(DIR_P, paste0(name, ".png")), p, width = width, height = height,
         dpi = 300, limitsize = FALSE)
  cat("[fig]", name, "\n")
}

## ---------------------------------------------------------------------------
## FIG 1  relative composition stacked bar
## ---------------------------------------------------------------------------
long <- melt(as.data.table(es, keep.rownames = "cell_type"),
             id.vars = "cell_type", variable.name = "sample", value.name = "score")
long[, group := si$group[match(sample, si$sample)]]

## per-sample relative proportion of ssGSEA scores (documented transformation)
long[, score_pos := score - min(score) + 1e-3, by = cell_type]
long[, prop := score_pos / sum(score_pos), by = sample]

ord <- si[order(factor(group, levels = c("NF", "DCM", "ICM")))]
long[, sample := factor(sample, levels = ord$sample)]
long[, group  := factor(group,  levels = c("NF", "DCM", "ICM"))]

## order cell types by mean relative abundance so the legend matches the stacking
cell_ord <- long[, .(m = mean(prop)), by = cell_type][order(-m)]$cell_type
long[, cell_type := factor(cell_type, levels = cell_ord)]

## 28 distinguishable fills: even hue rotation + alternating lightness
pal28 <- grDevices::hcl(h = seq(12, 372, length.out = 28),
                        c = 78,
                        l = rep(c(72, 42, 58, 82, 50, 66, 36), length.out = 28))

p1 <- ggplot(long, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(width = 1) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal28, breaks = cell_ord) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Samples (ordered by group)", y = "Relative proportion of ssGSEA score (%)",
       fill = "Immune cell type",
       title = "Immune cell infiltration composition across 313 left-ventricle samples",
       subtitle = "ssGSEA enrichment scores rescaled per sample; NOT absolute cell fractions") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 7))
save_fig(p1, "Fig1_immune_composition_stacked", 15, 7.2)
fwrite(long[, .(sample, group, cell_type, ssgsea_score = score, rel_prop = prop)],
       file.path(DIR_PD, "Fig1_composition_plotdata.csv"))

## ---------------------------------------------------------------------------
## FIG 2  infiltration heatmap (row z-score)
## ---------------------------------------------------------------------------
z <- t(scale(t(es)))
z[z >  2.5] <-  2.5
z[z < -2.5] <- -2.5

ha <- HeatmapAnnotation(
  Group = grp,
  col = list(Group = grp_cols),
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 9),
                                 labels_gp = gpar(fontsize = 8))
)

ht2 <- Heatmap(
  z, name = "z-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#3B4CC0", "white", "#B40426")),
  top_annotation = ha,
  column_split = grp, cluster_column_slices = FALSE,
  column_title_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 4),
  show_column_names = FALSE,
  cluster_rows = TRUE, cluster_columns = TRUE,
  row_title = NULL,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)),
  use_raster = TRUE, raster_quality = 3
)
pdf(file.path(DIR_P, "Fig2_immune_infiltration_heatmap.pdf"), width = 9, height = 8)
draw(ht2, merge_legend = TRUE)
dev.off()
png(file.path(DIR_P, "Fig2_immune_infiltration_heatmap.png"), width = 9, height = 8,
    units = "in", res = 300)
draw(ht2, merge_legend = TRUE)
dev.off()
cat("[fig] Fig2_immune_infiltration_heatmap\n")
write.csv(data.frame(cell_type = rownames(z), z, check.names = FALSE),
          file.path(DIR_PD, "Fig2_heatmap_zscore_data.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## FIG 3  immune-cell correlation heatmap
## ---------------------------------------------------------------------------
cm  <- st$cm[cells, cells]
cmp <- st$cmp[cells, cells]
star <- matrix("", nrow(cmp), ncol(cmp), dimnames = dimnames(cmp))
padj <- p.adjust(cmp[upper.tri(cmp)], "BH")
pu <- cmp[upper.tri(cmp)]
star[upper.tri(star)] <- ifelse(pu < 0.001, "***", ifelse(pu < 0.01, "**", ifelse(pu < 0.05, "*", "")))
star[lower.tri(star)] <- t(star)[lower.tri(star)]

ht3 <- Heatmap(
  cm, name = "Spearman rho",
  col = colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                   c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
  row_names_gp = gpar(fontsize = 7.5),
  column_names_gp = gpar(fontsize = 7.5),
  row_names_side = "left",
  cluster_rows = TRUE, cluster_columns = TRUE,
  row_title = NULL, column_title = "Immune cell type",
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)),
  width  = unit(13, "cm"), height = unit(13, "cm")
)
pdf(file.path(DIR_P, "Fig3_immune_cell_correlation_heatmap.pdf"), width = 11, height = 10)
draw(ht3, merge_legend = TRUE)
dev.off()
png(file.path(DIR_P, "Fig3_immune_cell_correlation_heatmap.png"), width = 11, height = 10,
    units = "in", res = 300)
draw(ht3, merge_legend = TRUE)
dev.off()
cat("[fig] Fig3_immune_cell_correlation_heatmap\n")
write.csv(data.frame(cell_type = rownames(cm), cm, check.names = FALSE),
          file.path(DIR_PD, "Fig3_correlation_plotdata.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## FIG 4  group comparison violin + box
## ---------------------------------------------------------------------------
kw <- st$kw; pw <- st$pw
sig_cells <- kw[FDR < 0.05][order(p_value)]$cell_type
if (length(sig_cells) < 6) sig_cells <- kw[order(p_value)][1:min(12, nrow(kw))]$cell_type
cat("[fig4] cell types plotted:", length(sig_cells), "\n")

df4 <- long[cell_type %in% sig_cells, .(sample, group, cell_type, score)]
df4[, cell_type := factor(cell_type, levels = sig_cells)]

ann <- rbindlist(lapply(sig_cells, function(cl) {
  k  <- kw[cell_type == cl]
  d1 <- pw[cell_type == cl & contrast == "NF_vs_DCM"]
  d2 <- pw[cell_type == cl & contrast == "NF_vs_ICM"]
  ## Cliff's delta < 0 means the second group (disease) is higher
  dir1 <- if (d1$cliffs_delta < 0) "DCM higher" else "DCM lower"
  dir2 <- if (d2$cliffs_delta < 0) "ICM higher" else "ICM lower"
  data.table(cell_type = cl,
             label = paste0("KW: P=", formatC(k$p_value, format = "e", digits = 1),
                            ", FDR=", formatC(k$FDR, format = "e", digits = 1),
                            "\nNF-DCM: P=", formatC(d1$p_value, format = "e", digits = 1),
                            ", FDR=", formatC(d1$FDR, format = "e", digits = 1),
                            " (", dir1, ")",
                            "\nNF-ICM: P=", formatC(d2$p_value, format = "e", digits = 1),
                            ", FDR=", formatC(d2$FDR, format = "e", digits = 1),
                            " (", dir2, ")"))
}))
ann[, cell_type := factor(cell_type, levels = sig_cells)]

ymax <- df4[, .(ymax = max(score) + 0.26 * diff(range(score))), by = cell_type]
ann  <- merge(ann, ymax, by = "cell_type")

p4 <- ggplot(df4, aes(x = group, y = score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.95, linewidth = 0.3) +
  geom_jitter(width = 0.14, size = 0.35, alpha = 0.35, colour = "grey25") +
  geom_text(data = ann, aes(x = 2, y = ymax, label = label), inherit.aes = FALSE,
            size = 1.95, lineheight = 1.02, vjust = 1) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = grp_cols) +
  scale_x_discrete(labels = c(paste0("NF\n(n=", n_grp["NF"], ")"),
                              paste0("DCM\n(n=", n_grp["DCM"], ")"),
                              paste0("ICM\n(n=", n_grp["ICM"], ")"))) +
  labs(x = NULL, y = "ssGSEA enrichment score",
       title = "Differentially infiltrated immune cell types across NF / DCM / ICM",
       subtitle = paste("Kruskal-Wallis (global) + pairwise Wilcoxon rank-sum, BH-FDR;",
                        "effect direction is stated relative to NF;\n",
                        "violin = density, box = median/IQR, dots = individual samples")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 7.5, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7.5),
        plot.subtitle = element_text(size = 7.2))
save_fig(p4, "Fig4_immune_group_comparison_violin", 12, 2.6 * ceiling(length(sig_cells) / 4) + 1.4)

fwrite(df4, file.path(DIR_PD, "Fig4_group_comparison_plotdata.csv"))
fwrite(ann[, .(cell_type, annotation = label)], file.path(DIR_PD, "Fig4_annotation_text.csv"))

## ---------------------------------------------------------------------------
## FIG 5  key gene x immune cell correlation heatmap
## ---------------------------------------------------------------------------
cor_tab <- st$cor_tab
ct5 <- copy(cor_tab)
ct5[, sig := ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**",
              ifelse(FDR < 0.05, "*", "")))]
mat_r <- dcast(ct5, cell_type ~ gene, value.var = "rho")
mat_s <- dcast(ct5, cell_type ~ gene, value.var = "sig")
r_mat <- as.matrix(mat_r[, -1]); rownames(r_mat) <- mat_r$cell_type
s_mat <- as.matrix(mat_s[, -1]); rownames(s_mat) <- mat_s$cell_type
r_mat <- r_mat[, key_genes, drop = FALSE]; s_mat <- s_mat[, key_genes, drop = FALSE]

## order rows by hierarchical clustering
hc <- hclust(as.dist(1 - cor(t(r_mat), use = "pairwise.complete.obs")))
r_mat <- r_mat[hc$order, , drop = FALSE]
s_mat <- s_mat[hc$order, , drop = FALSE]

ht5 <- Heatmap(
  r_mat, name = "Spearman rho",
  col = colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                   c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(s_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 0,
  row_names_side = "left",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title = "Immune cell type", column_title = "Key gene",
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)),
  width = unit(6.2, "cm"), height = unit(11.5, "cm")
)
pdf(file.path(DIR_P, "Fig5_keygene_immune_correlation_heatmap.pdf"), width = 7.6, height = 9.2)
draw(ht5, merge_legend = TRUE)
dev.off()
png(file.path(DIR_P, "Fig5_keygene_immune_correlation_heatmap.png"), width = 7.6, height = 9.2,
    units = "in", res = 300)
draw(ht5, merge_legend = TRUE)
dev.off()
cat("[fig] Fig5_keygene_immune_correlation_heatmap\n")

write.csv(data.frame(cell_type = rownames(r_mat), r_mat, check.names = FALSE),
          file.path(DIR_PD, "Fig5_keygene_correlation_matrix.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## FIG 6  scatter: strongest key gene - immune cell pairs
## ---------------------------------------------------------------------------
top <- cor_tab[FDR < 0.05][order(-abs(rho))][1:min(9, .N)]
df6 <- rbindlist(lapply(seq_len(nrow(top)), function(i) {
  g  <- top$gene[i]; cl <- top$cell_type[i]
  data.table(gene = g, cell_type = cl,
             expr = as.numeric(expr[g, ]),
             score = as.numeric(es[cl, ]),
             group = si$group,
             rho = top$rho[i], FDR = top$FDR[i])
}))
df6[, panel := factor(paste0(gene, " vs\n", cell_type),
                      levels = unique(paste0(gene, " vs\n", cell_type)))]

lab6 <- df6[, .(label = paste0("rho=", sprintf("%.2f", rho[1]),
                               "\nFDR=", formatC(FDR[1], format = "e", digits = 1))),
            by = panel]

p6 <- ggplot(df6, aes(expr, score, colour = group)) +
  geom_point(size = 0.9, alpha = 0.65) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
              linewidth = 0.6, fill = "grey80") +
  geom_text(data = lab6, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.15, size = 2.3, lineheight = 0.95) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_colour_manual(values = grp_cols) +
  labs(x = "Gene expression (log2)", y = "ssGSEA enrichment score",
       title = "Key gene expression versus immune cell infiltration (Spearman)",
       colour = "Group") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(size = 7.2),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
save_fig(p6, "Fig6_keygene_immune_scatter", 10, 3.2 * ceiling(nrow(top) / 3) + 1.0)
fwrite(df6, file.path(DIR_PD, "Fig6_scatter_plotdata.csv"))

## ---------------------------------------------------------------------------
## FIG 7  heatmap of group medians (overview of all 28 cells)
## ---------------------------------------------------------------------------
med_mat <- sapply(levels(grp), function(g) matrixStats::rowMedians(es[, grp == g, drop = FALSE]))
rownames(med_mat) <- rownames(es)
colnames(med_mat) <- levels(grp)
med <- data.table(cell_type = rownames(med_mat), as.data.table(med_mat))
zf <- t(scale(t(med_mat)))
zf[is.na(zf)] <- 0

ht7 <- Heatmap(
  zf, name = "z of median",
  col = colorRamp2(c(-1.5, 0, 1.5), c("#4575B4", "white", "#D73027")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 0,
  row_names_side = "left",
  cluster_rows = TRUE, cluster_columns = FALSE,
  row_title = NULL, column_title = "Group",
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey90", lwd = 0.4))
    grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6.2))
  },
  width = unit(6, "cm"), height = unit(12.5, "cm"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8))
)
pdf(file.path(DIR_P, "Fig7_group_median_immune_heatmap.pdf"), width = 7.6, height = 9.5)
draw(ht7, merge_legend = TRUE)
dev.off()
png(file.path(DIR_P, "Fig7_group_median_immune_heatmap.png"), width = 7.6, height = 9.5,
    units = "in", res = 300)
draw(ht7, merge_legend = TRUE)
dev.off()
cat("[fig] Fig7_group_median_immune_heatmap\n")
write.csv(data.frame(cell_type = rownames(med_mat), med_mat, check.names = FALSE),
          file.path(DIR_PD, "Fig7_group_median_matrix.csv"), row.names = FALSE)

cat("\n=== 04_immune_plots.R finished ===\n")
sink()
