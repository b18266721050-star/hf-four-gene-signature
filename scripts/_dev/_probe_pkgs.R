pk <- c("Seurat", "SeuratObject", "harmony", "mclust", "presto", "clustree",
        "scDblFinder", "uwot", "Rcpp", "RcppArmadillo", "ggplot2", "patchwork",
        "pheatmap", "dplyr", "future", "BPCells")
res <- data.frame(pkg = pk,
                  installed = pk %in% rownames(installed.packages()),
                  stringsAsFactors = FALSE)
res$version <- NA
for (i in which(res$installed)) res$version[i] <- tryCatch(
  as.character(packageVersion(res$pkg[i])), error = function(e) NA)
write.csv(res, "E:/another data file/GEObulk analysis project/_tmp/probe_pkgs.csv", row.names = FALSE)
