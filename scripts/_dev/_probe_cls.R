logf <- "E:/another data file/GEObulk analysis project/_tmp/probe_cls.txt"
cat("", file = logf)
lg <- function(...) cat(paste0(..., "\n"), file = logf, append = TRUE)
suppressPackageStartupMessages(library(Seurat))
so <- readRDS("E:/another data file/GEObulk analysis project/results/Singlecel/objects/03_merged_validation/04_integrated_val2.rds")
for (v in c("preint_cluster", "sample", "donor", "condition", "integrated_cluster",
            "harmony_snn_res.0.2", "harmony_snn_res.0.6")) {
  x <- tryCatch(so[[v]], error = function(e) NULL)
  lg(sprintf("%-22s class=%s  len=%s", v, paste(class(x), collapse = "/"),
             if (is.null(x)) "NULL" else length(x)))
}
lg("--- via $ ---")
for (v in c("sample", "donor", "condition", "preint_cluster")) {
  x <- so[[v]]
  lg(sprintf("%-16s class=%s len=%s", v, paste(class(x), collapse = "/"), length(x)))
}
lg("meta colnames: ", paste(colnames(so@meta.data), collapse = " | "))
lg("END")
