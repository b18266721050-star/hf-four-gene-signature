logf <- "E:/another data file/GEObulk analysis project/_tmp/probe_mem.txt"
cat("", file = logf)
lg <- function(...) { cat(paste0(..., "\n"), file = logf, append = TRUE); flush.console() }
lg("start ", format(Sys.time()), " pid=", Sys.getpid())
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
lg("Seurat ", as.character(packageVersion("Seurat")))
proj <- "E:/another data file/GEObulk analysis project"
d <- file.path(proj, "results", "Singlecel", "objects", "02_filtered")
for (s in c("GSM8848306", "GSM8848314", "GSM8848324")) {
  lg("reading ", s, " ...")
  f <- file.path(d, paste0(s, "_filtered.rds"))
  lg("  fileMB=", round(file.size(f) / 1024^2, 1))
  t0 <- Sys.time()
  o <- tryCatch(readRDS(f), error = function(e) { lg("  ERROR ", conditionMessage(e)); NULL })
  if (is.null(o)) next
  lg("  loaded cells=", ncol(o), " genes=", nrow(o),
     " objMB=", round(as.numeric(object.size(o)) / 1024^2, 1),
     " s=", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  m <- LayerData(o[["RNA"]], layer = "counts")
  lg("  nnz=", length(m@x), " nnzMB=", round(as.numeric(object.size(m)) / 1024^2, 1),
     " layers=", paste(Layers(o[["RNA"]]), collapse = "|"))
  rm(o, m); invisible(gc())
}
lg("Rmem_MB=", round(sum(gc()[, 2])), "  freeGB=",
   round(as.numeric(system("wmic OS get FreePhysicalMemory /value", intern = TRUE)[2]) / 1e6, 2))
lg("END")
