suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
so <- readRDS("E:/another data file/GEObulk analysis project/results/Singlecel/objects/03_merged_validation/04_integrated_val2.rds")
out <- "E:/another data file/GEObulk analysis project/_tmp/plot_test.png"
cols <- c(GSM8848306 = "#E41A1C", GSM8848307 = "#377EB8")
mk <- function(ps, dpi, al, ttl) {
  DimPlot(so, reduction = "umap.harmony", group.by = "sample", cols = cols,
          pt.size = ps, alpha = al, raster = TRUE, raster.dpi = c(dpi, dpi)) +
    ggtitle(ttl) + theme(legend.position = "none")
}
p <- mk(0.1, 700, 1, "pt0.1 raster700") | mk(0.5, 300, 0.8, "pt0.5 raster300") |
     mk(1, 300, 0.5, "pt1 raster300 a0.5")
ggsave(out, p, width = 13, height = 4.4, dpi = 300)
cat("OK\n", file = "E:/another data file/GEObulk analysis project/_tmp/plot_test.txt")
