logf <- "E:/another data file/GEObulk analysis project/_tmp/probe_meta.txt"
cat("", file = logf)
lg <- function(...) cat(paste0(..., "\n"), file = logf, append = TRUE)
suppressPackageStartupMessages({ library(Seurat) })
d <- "E:/another data file/GEObulk analysis project/results/Singlecel/objects/02_filtered"
o <- readRDS(file.path(d, "GSM8848306_filtered.rds"))
lg("dim: ", nrow(o), " x ", ncol(o))
lg("assays: ", paste(Assays(o), collapse = ","), " | layers: ", paste(Layers(o[["RNA"]]), collapse = ","))
lg("meta cols: ", paste(colnames(o@meta.data), collapse = " | "))
lg("--- head(meta,3) ---")
m <- o@meta.data[1:3, , drop = FALSE]
lg(paste(capture.output(print(m)), collapse = "\n"))
lg("--- unique values ---")
for (cc in c("sample", "donor", "condition", "batch", "tissue", "technology")) {
  if (cc %in% colnames(o@meta.data)) {
    v <- unique(as.character(o@meta.data[[cc]]))
    lg(cc, ": ", paste(head(v, 5), collapse = "/"), " (n=", length(v), ")", sep = "")
  } else lg(cc, ": MISSING")
}
lg("reductions: ", paste(Reductions(o), collapse = ","))
lg("misc: ", paste(names(o@misc), collapse = ","))
lg("END")
