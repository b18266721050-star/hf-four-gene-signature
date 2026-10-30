# =============================================================================
# 41_sc_stage2_merge_norm_pca.R
# Stage 2 - Step 2: merge the 25 filtered objects, LogNormalize, HVG, scale,
# PCA, and an UNINTEGRATED neighbour graph / UMAP / clustering.
#
# Design notes (memory-constrained machine, ~16 GB RAM):
#   * objects are merged sequentially and released immediately (peak ~2x final)
#   * scale.data is dropped right after RunPCA (frees ~2.4 GB) because the
#     downstream steps only need the PCA embedding
#   * every stage is checkpointed so a later failure does not repeat the merge
#
# Env vars:
#   SC_SUBSET=n   -> use only the first n samples (validation run, separate dirs)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2)
  library(patchwork); library(RColorBrewer)
})

SEED <- 20261028
set.seed(SEED)

proj   <- "E:/another data file/GEObulk analysis project"
sdir   <- file.path(proj, "results", "Singlecel")
obj_in <- file.path(sdir, "objects", "02_filtered")
tab_d  <- file.path(sdir, "tables")
fig_d  <- file.path(sdir, "figures")
log_d  <- file.path(sdir, "logs")

n_sub  <- Sys.getenv("SC_SUBSET", "")
tag    <- if (nzchar(n_sub)) paste0("_val", n_sub) else ""
if (nzchar(n_sub)) {
  obj_d <- file.path(sdir, "objects", "03_merged_validation")
  fig_d <- file.path(fig_d, "02_preint_validation")
  tab_d <- file.path(sdir, "tables", "validation")
} else {
  obj_d <- file.path(sdir, "objects", "03_merged")
  fig_d <- file.path(fig_d, "02_preint")
}
for (d in c(obj_d, fig_d, tab_d, log_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file(file.path(log_d, paste0("41_merge_norm_pca", tag, ".log")), "wt")
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf); flush(lgf)
}
options(error = function() {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ERROR: ", geterrmessage(), "\n"),
      file = lgf, append = TRUE); flush(lgf)
  q(status = 1)
})
mem <- function(tag2 = "") {
  g <- gc(verbose = FALSE)
  lg(sprintf("  mem%s: R=%d MB  peak=%d MB", tag2,
             round(sum(g[, 2])), round(sum(g[, 6]))))
}

f_merged <- file.path(obj_d, paste0("03a_merged_counts", tag, ".rds"))
f_norm   <- file.path(obj_d, paste0("03b_norm_pca", tag, ".rds"))
f_pre    <- file.path(obj_d, paste0("03_preintegration", tag, ".rds"))

lg("=== 41 merge / normalize / PCA  (seed=", SEED, ") ===")
lg("R ", R.version.string, " | Seurat ", as.character(packageVersion("Seurat")),
   " | SeuratObject ", as.character(packageVersion("SeuratObject")))

# ---------------------------------------------------------------- 1. merge ---
files <- sort(list.files(obj_in, pattern = "_filtered[.]rds$", full.names = TRUE))
if (nzchar(n_sub)) files <- files[seq_len(min(as.integer(n_sub), length(files)))]
lg("input objects: ", length(files))

if (file.exists(f_merged)) {
  lg("SKIP merge (checkpoint exists)")
  so <- readRDS(f_merged)
} else {
  stopifnot(length(files) >= 2)
  t0 <- Sys.time()
  so <- readRDS(files[1])
  lg(sprintf("[ 1/%d] %s  %d cells", length(files), basename(files[1]), ncol(so)))
  for (i in seq_along(files)[-1]) {
    tmp <- readRDS(files[i])
    so <- merge(so, y = tmp)
    lg(sprintf("[%2d/%d] %s  +%d cells -> %d cells", i, length(files),
               basename(files[i]), ncol(tmp), ncol(so)))
    rm(tmp); invisible(gc(verbose = FALSE))
  }
  lg("merged layers: ", length(Layers(so[["RNA"]])), " | cells: ", ncol(so))
  so <- JoinLayers(so)
  lg("after JoinLayers: ", paste(Layers(so[["RNA"]]), collapse = ","))
  Idents(so) <- "sample"
  lg(sprintf("merge done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  mem("after-merge")
  saveRDS(so, f_merged, compress = FALSE)
  lg("checkpoint written: ", basename(f_merged),
     sprintf(" (%.1f GB)", file.size(f_merged) / 1024^3))
}

so$sample    <- factor(so$sample)
so$condition <- factor(so$condition, levels = intersect(c("Control", "Disease"), unique(so$condition)))
so$donor     <- factor(so$donor)

# --------------------------------------------------- 2. normalization / HVG ---
if (file.exists(f_norm)) {
  lg("SKIP normalize/scale/PCA (checkpoint exists)")
  so <- readRDS(f_norm)
} else {
  t0 <- Sys.time()
  so <- NormalizeData(so, normalization.method = "LogNormalize",
                      scale.factor = 1e4, verbose = FALSE)
  lg(sprintf("NormalizeData(LogNormalize, 1e4) done in %.1f min",
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  mem("after-norm")

  so <- FindVariableFeatures(so, selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE)
  hv <- VariableFeatures(so)
  lg("HVG (vst, top 2000): ", length(hv),
     " | top10 = ", paste(head(hv, 10), collapse = ", "))
  hvi <- HVFInfo(so)[hv, , drop = FALSE]
  write.csv(data.frame(gene = hv, rank = seq_along(hv),
                       mean = hvi$mean, variance = hvi$variance,
                       variance.standardized = hvi$variance.standardized),
            file.path(tab_d, paste0("02_hvg_top2000", tag, ".csv")), row.names = FALSE)

  t0 <- Sys.time()
  so <- ScaleData(so, features = hv, verbose = FALSE)
  lg(sprintf("ScaleData on %d HVG done in %.1f min", length(hv),
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  mem("after-scale")

  t0 <- Sys.time()
  so <- RunPCA(so, features = hv, npcs = 30, seed.use = SEED, verbose = FALSE)
  lg(sprintf("RunPCA(npcs=30) done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ## variance table for the elbow decision
  sd <- so[["pca"]]@stdev
  vex <- sd^2 / sum(sd^2) * 100
  pv <- data.frame(PC = seq_along(sd), stdev = sd, var_pct = vex,
                   cum_var_pct = cumsum(vex))
  write.csv(pv, file.path(tab_d, paste0("02_pca_variance", tag, ".csv")), row.names = FALSE)
  lg("PC1-5 var%: ", paste(sprintf("%.2f", vex[1:5]), collapse = "/"),
     " | cum PC20 = ", sprintf("%.1f%%", pv$cum_var_pct[20]))

  ## free the dense scale.data before the graph steps
  if ("scale.data" %in% Layers(so[["RNA"]])) {
    ok <- tryCatch({ LayerData(so[["RNA"]], layer = "scale.data") <- NULL; TRUE },
                   error = function(e) FALSE)
    lg("scale.data dropped: ", ok, " | remaining layers: ",
       paste(Layers(so[["RNA"]]), collapse = ","))
    invisible(gc(verbose = FALSE))
  }
  mem("after-drop-scale")
  saveRDS(so, f_norm, compress = FALSE)
  lg("checkpoint written: ", basename(f_norm),
     sprintf(" (%.1f GB)", file.size(f_norm) / 1024^3))
}

# --------------------------------- 3. unintegrated graph / UMAP / clusters ---
DIMS <- 1:20
PT <- if (ncol(so) > 80000) 0.25 else 0.6   # point size tuned to cell number
lg("unintegrated graph on dims ", min(DIMS), "-", max(DIMS))
so <- FindNeighbors(so, reduction = "pca", dims = DIMS, k.param = 20, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = DIMS, n.neighbors = 30,
              min.dist = 0.3, seed.use = SEED, verbose = FALSE,
              reduction.name = "umap.unintegrated", reduction.key = "UMAPun_")
so <- FindClusters(so, resolution = 0.5, random.seed = SEED, verbose = FALSE,
                   graph.name = "RNA_snn")
so$preint_cluster <- so$seurat_clusters
lg("unintegrated clusters (res=0.5): ", length(unique(so$preint_cluster)))

# ------------------------------------------------------------ 4. figures ----
pal25 <- tryCatch(DiscretePalette(n = nlevels(so$sample), palette = "glasbey"),
                  error = function(e) DiscretePalette(n = nlevels(so$sample), palette = "polychrome"))
names(pal25) <- levels(so$sample)

p1 <- DimPlot(so, reduction = "umap.unintegrated", group.by = "sample",
              cols = pal25, pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Unintegrated UMAP - sample") + theme(legend.text = element_text(size = 6))
p2 <- DimPlot(so, reduction = "umap.unintegrated", group.by = "condition",
              cols = c(Control = "#4C72B0", Disease = "#C44E52"),
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Unintegrated UMAP - condition")
ggsave(file.path(fig_d, paste0("umap_pca_sample_condition", tag, ".png")),
       p1 | p2, width = 15, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_pca_sample_condition", tag, ".pdf")),
       p1 | p2, width = 15, height = 6.4)

p3 <- DimPlot(so, reduction = "umap.unintegrated", group.by = "preint_cluster",
              label = TRUE, repel = TRUE, pt.size = PT, alpha = 0.8,
              raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Unintegrated UMAP - cluster (res 0.5)") + NoLegend()
ggsave(file.path(fig_d, paste0("umap_pca_cluster", tag, ".png")), p3,
       width = 7.4, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_pca_cluster", tag, ".pdf")), p3,
       width = 7.4, height = 6.4)

p4 <- ElbowPlot(so, ndims = 30, reduction = "pca") +
  geom_vline(xintercept = max(DIMS), linetype = 2, colour = "red") +
  ggtitle("PCA elbow (red = dims used: 20)")
ggsave(file.path(fig_d, paste0("pca_elbow", tag, ".png")), p4,
       width = 6, height = 4.6, dpi = 300)

top10 <- head(VariableFeatures(so), 10)
p5 <- LabelPoints(plot = VariableFeaturePlot(so), points = top10, repel = TRUE)
ggsave(file.path(fig_d, paste0("hvg_top10", tag, ".png")), p5,
       width = 6.6, height = 4.6, dpi = 300)

# ------------------------------------------------------------- 5. tables ----
comp <- as.data.frame.matrix(table(so$preint_cluster, so$sample))
comp <- cbind(data.frame(cluster = rownames(comp),
                         n_cells = as.integer(rowSums(comp)),
                         stringsAsFactors = FALSE), comp)
comp_pct <- sweep(as.matrix(comp[, as.character(levels(so$sample)), drop = FALSE]),
                  1, rowSums(comp[, as.character(levels(so$sample)), drop = FALSE]), "/") * 100
out <- data.frame(cluster = comp$cluster, n_cells = comp$n_cells,
                  n_samples_present = rowSums(comp[, as.character(levels(so$sample))] > 0),
                  max_sample_pct = round(apply(comp_pct, 1, max), 1),
                  stringsAsFactors = FALSE)
write.csv(out, file.path(tab_d, paste0("02_preint_cluster_summary", tag, ".csv")), row.names = FALSE)
write.csv(comp, file.path(tab_d, paste0("02_preint_cluster_by_sample", tag, ".csv")), row.names = FALSE)

ct <- table(so$preint_cluster, so$condition)
write.csv(as.data.frame.matrix(ct),
          file.path(tab_d, paste0("02_preint_cluster_by_condition", tag, ".csv")))

so$preint_cluster <- factor(so$preint_cluster, levels = as.character(sort(as.integer(unique(so$preint_cluster)))))
saveRDS(so, f_pre, compress = FALSE)
lg("saved: ", basename(f_pre), sprintf(" (%.1f GB)", file.size(f_pre) / 1024^3))

writeLines(capture.output(sessionInfo()),
           file.path(log_d, paste0("sessionInfo_41", tag, ".txt")))
write.csv(data.frame(package = c("Seurat", "SeuratObject", "Matrix", "ggplot2", "uwot", "R"),
                     version = c(as.character(packageVersion("Seurat")),
                                 as.character(packageVersion("SeuratObject")),
                                 as.character(packageVersion("Matrix")),
                                 as.character(packageVersion("ggplot2")),
                                 as.character(packageVersion("uwot")),
                                 paste(R.version$major, R.version$minor, sep = "."))),
          file.path(tab_d, paste0("02_software_versions", tag, ".csv")), row.names = FALSE)
lg("=== 41 END ===")
close(lgf)
