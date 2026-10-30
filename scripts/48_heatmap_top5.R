## cleaner heatmap: top-5 score-ranked markers per cluster (85 rows instead of 170)
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})
proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
fig_d <- file.path(sdir, "figures", "04_markers")
tab_d <- file.path(sdir, "tables")

so <- readRDS(file.path(sdir, "objects", "03_merged", "04_integrated.rds"))
cl <- as.character(so@meta.data[["harmony_snn_res.0.6"]])
so$cluster <- factor(cl, levels = sort(unique(as.integer(cl))))
Idents(so) <- "cluster"

top <- read.csv(file.path(tab_d, "04_markers_top10_by_score.csv"), stringsAsFactors = FALSE)
top$cluster <- as.character(top$cluster)
top5 <- top %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup()
sym <- setNames(ifelse(is.na(top5$gene_symbol), top5$gene, top5$gene_symbol), top5$gene)
sym <- sym[!duplicated(names(sym))]
lab_fun <- function(x) { y <- unname(sym[x]); ifelse(is.na(y), x, y) }

set.seed(20261028)
sub <- subset(so, downsample = 120)
sub <- ScaleData(sub, features = unique(top5$gene), verbose = FALSE)
ph <- DoHeatmap(sub, features = unique(top5$gene), group.by = "cluster",
                size = 4, raster = FALSE, angle = 0) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C")
ph2 <- tryCatch(ph + scale_y_discrete(labels = lab_fun), error = function(e) ph)
ggsave(file.path(fig_d, "heatmap_top5_score.png"), ph2, width = 12, height = 14, dpi = 300)
ggsave(file.path(fig_d, "heatmap_top5_score.pdf"), ph2, width = 12, height = 14)
cat("ok\n")
