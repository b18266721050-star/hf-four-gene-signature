## re-draw marker plots using the score-based ranking
## (score = avg_log2FC * (pct.1 - pct.2)) which avoids lowly-expressed lncRNAs
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})
proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
fig_d <- file.path(sdir, "figures", "04_markers")
tab_d <- file.path(sdir, "tables")

so <- readRDS(file.path(sdir, "objects", "03_merged", "04_integrated.rds"))
so$cluster <- factor(as.character(so@meta.data[["harmony_snn_res.0.6"]]),
                     levels = sort(unique(as.integer(as.character(so@meta.data[["harmony_snn_res.0.6"]])))))
Idents(so) <- "cluster"

top <- read.csv(file.path(tab_d, "04_markers_top10_by_score.csv"), stringsAsFactors = FALSE)
top$cluster <- as.character(top$cluster)
sym <- setNames(ifelse(is.na(top$gene_symbol), top$gene, top$gene_symbol), top$gene)
sym <- sym[!duplicated(names(sym))]
lab_fun <- function(x) { y <- unname(sym[x]); ifelse(is.na(y), x, y) }

top3  <- top %>% group_by(cluster) %>% slice_head(n = 3) %>% ungroup()
top10 <- top %>% group_by(cluster) %>% slice_head(n = 10) %>% ungroup()

## DotPlot, 4 clusters per page
lev <- levels(so$cluster)
pages <- split(lev, ceiling(seq_along(lev) / 4))
for (k in seq_along(pages)) {
  cc <- pages[[k]]
  feats <- top3$gene[top3$cluster %in% cc]
  dp <- DotPlot(so, features = feats, group.by = "cluster", assay = "RNA") +
    RotatedAxis() +
    scale_x_discrete(labels = lab_fun) +
    labs(title = paste0("Top 3 markers (score-ranked) - clusters ", paste(cc, collapse = ",")),
         x = NULL, y = "cluster") +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(fig_d, sprintf("dotplot_score_page%d.png", k)), dp,
         width = 9.5, height = 4.6, dpi = 300)
}

## FeaturePlot, 6 genes per page
fp_genes <- top3$gene
fp_pages <- split(fp_genes, ceiling(seq_along(fp_genes) / 6))
for (k in seq_along(fp_pages)) {
  fp <- FeaturePlot(so, features = fp_pages[[k]], ncol = 3, pt.size = 0.25,
                    alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300),
                    order = TRUE, reduction = "umap.harmony")
  fp <- fp & scale_colour_gradientn(colours = c("grey90", "#FEE0B6", "#D7191C"))
  ggsave(file.path(fig_d, sprintf("featureplot_score_page%d.png", k)), fp,
         width = 11, height = 7, dpi = 250)
}

## heatmap of top10 by score on downsampled cells
set.seed(20261028)
sub <- subset(so, downsample = 100)
hm <- unique(top10$gene)
sub <- ScaleData(sub, features = hm, verbose = FALSE)
ph <- DoHeatmap(sub, features = hm, group.by = "cluster", size = 2.6,
                raster = FALSE, angle = 90) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C")
ph2 <- tryCatch(ph + scale_y_discrete(labels = lab_fun), error = function(e) ph)
ggsave(file.path(fig_d, "heatmap_top10_score.png"), ph2, width = 15, height = 13, dpi = 250)
ggsave(file.path(fig_d, "heatmap_top10_score.pdf"), ph2, width = 15, height = 13)

cat("done:", length(pages), "dotplot pages,", length(fp_pages), "featureplot pages\n")
