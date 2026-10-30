# =============================================================================
# 42_sc_stage2_harmony_cluster.R
# Stage 2 - Step 3: Harmony integration on the PCA embedding, neighbour graph,
# UMAP and clustering in the integrated space, multi-resolution comparison,
# and before/after composition diagnostics (including an over-correction check).
#
# Env vars:
#   SC_SUBSET=n   -> validation mode (reads the _valN checkpoints)
#   SC_RES=0.6    -> primary resolution used for composition tables / markers
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork)
  library(dplyr); library(mclust); library(harmony)
})

SEED <- 20261028
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
tab_d <- file.path(sdir, "tables")
fig_d <- file.path(sdir, "figures")
log_d <- file.path(sdir, "logs")

n_sub <- Sys.getenv("SC_SUBSET", "")
tag   <- if (nzchar(n_sub)) paste0("_val", n_sub) else ""
if (nzchar(n_sub)) {
  obj_d <- file.path(sdir, "objects", "03_merged_validation")
  fig_d <- file.path(fig_d, "02_preint_validation")
  tab_d <- file.path(sdir, "tables", "validation")
} else {
  obj_d <- file.path(sdir, "objects", "03_merged")
  fig_d <- file.path(fig_d, "03_integrated")
}
for (d in c(fig_d, tab_d, log_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

PRIMARY_RES <- as.numeric(Sys.getenv("SC_RES", "0.6"))
RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0)
DIMS <- 1:20
PT <- 0.6   # reset after loading (depends on cell number)

lgf <- file(file.path(log_d, paste0("42_harmony", tag, ".log")), "wt")
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf); flush(lgf)
}
options(error = function() {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ERROR: ", geterrmessage(), "\n"),
      file = lgf, append = TRUE); flush(lgf)
  q(status = 1)
})

f_pre <- file.path(obj_d, paste0("03_preintegration", tag, ".rds"))
f_out <- file.path(obj_d, paste0("04_integrated", tag, ".rds"))

lg("=== 42 Harmony integration (seed=", SEED, ", dims=", max(DIMS), ") ===")
so <- readRDS(f_pre)
PT <- if (ncol(so) > 80000) 0.25 else 0.6
lg("loaded: ", ncol(so), " cells x ", nrow(so), " genes | reductions: ",
   paste(Reductions(so), collapse = ","))
lg("cells per sample range: ", min(table(so$sample)), "-", max(table(so$sample)))
lg("harmony ", as.character(packageVersion("harmony")), " | mclust ",
   as.character(packageVersion("mclust")))

# ------------------------------------------------------------ 1. Harmony ----
t0 <- Sys.time()
so <- RunHarmony(so, group.by.vars = "sample", reduction.use = "pca",
                 dims.use = DIMS, reduction.save = "harmony",
                 theta = 2, lambda = 1, sigma = 0.1, project.dim = FALSE,
                 max_iter = 20, early_stop = TRUE, verbose = FALSE)
lg(sprintf("RunHarmony done in %.1f min | harmony dims = %d",
           as.numeric(difftime(Sys.time(), t0, units = "mins")),
           ncol(Embeddings(so, "harmony"))))
invisible(gc(verbose = FALSE))

# ------------------------------------------- 2. graph / UMAP / clustering ---
so <- FindNeighbors(so, reduction = "harmony", dims = DIMS, k.param = 20,
                    graph.name = "harmony_snn", verbose = FALSE)
so <- RunUMAP(so, reduction = "harmony", dims = DIMS, n.neighbors = 30,
              min.dist = 0.3, seed.use = SEED, verbose = FALSE,
              reduction.name = "umap.harmony", reduction.key = "UMAPha_")
lg("integrated UMAP done")

t0 <- Sys.time()
so <- FindClusters(so, graph.name = "harmony_snn", resolution = RESOLUTIONS,
                   random.seed = SEED, verbose = FALSE)
lg(sprintf("FindClusters at %d resolutions done in %.1f min",
           length(RESOLUTIONS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res_cols <- paste0("harmony_snn_res.", RESOLUTIONS)
res_cols <- res_cols[res_cols %in% colnames(so@meta.data)]
if (!length(res_cols)) stop("resolution columns not found; check FindClusters naming")
lg("resolution columns: ", paste(res_cols, collapse = ", "))

# --------------------------------------------- 3. resolution stability -------
cl <- so@meta.data[, res_cols, drop = FALSE]
for (j in seq_along(res_cols)) cl[[j]] <- as.integer(as.character(cl[[j]]))

stab <- expand.grid(res1 = RESOLUTIONS, res2 = RESOLUTIONS)
stab$ARI <- mapply(function(a, b) adjustedRandIndex(cl[[paste0("harmony_snn_res.", a)]],
                                                    cl[[paste0("harmony_snn_res.", b)]]),
                   stab$res1, stab$res2)
write.csv(stab, file.path(tab_d, paste0("03_resolution_ARI", tag, ".csv")), row.names = FALSE)

sz <- lapply(seq_along(RESOLUTIONS), function(i) {
  k <- cl[[i]]
  data.frame(resolution = RESOLUTIONS[i], n_clusters = length(unique(k)),
             min_size = min(table(k)), median_size = median(table(k)),
             max_size = max(table(k)),
             n_clusters_lt20 = sum(table(k) < 20))
}) %>% bind_rows()
sz$ARI_vs_prev <- c(NA, sapply(2:nrow(sz), function(i)
  adjustedRandIndex(cl[[i]], cl[[i - 1]])))
sz$ARI_vs_primary <- sapply(seq_along(RESOLUTIONS), function(i)
  adjustedRandIndex(cl[[i]], cl[[which.min(abs(RESOLUTIONS - PRIMARY_RES))]]))
write.csv(sz, file.path(tab_d, paste0("03_resolution_stability", tag, ".csv")), row.names = FALSE)
lg("resolution summary:")
for (i in seq_len(nrow(sz))) lg(sprintf("  res=%.1f  clusters=%d  size %d/%d/%d  ARI_vs_prev=%s",
      sz$resolution[i], sz$n_clusters[i], sz$min_size[i], sz$median_size[i], sz$max_size[i],
      ifelse(is.na(sz$ARI_vs_prev[i]), "-", sprintf("%.3f", sz$ARI_vs_prev[i]))))

# cluster-level Jaccard between adjacent resolutions (max-overlap matching)
jac <- function(a, b) {
  tb <- table(a, b)
  inter <- as.matrix(tb)
  ta <- rowSums(inter); tbv <- colSums(inter)
  m <- inter / (ta + rep(tbv, each = nrow(inter)) - inter)
  data.frame(cluster_low = rownames(m),
             best_match = colnames(m)[apply(m, 1, which.max)],
             jaccard = round(apply(m, 1, max), 3),
             n_cells = as.integer(ta))
}
jac_list <- lapply(2:length(RESOLUTIONS), function(i) {
  d <- jac(cl[[i - 1]], cl[[i]])
  d$from <- RESOLUTIONS[i - 1]; d$to <- RESOLUTIONS[i]; d
})
jac_df <- bind_rows(jac_list)
write.csv(jac_df, file.path(tab_d, paste0("03_cluster_jaccard_adjacent", tag, ".csv")),
          row.names = FALSE)

# ------------------------------------------------ 4. composition metrics -----
comp_metrics <- function(cluster_vec, so, label) {
  if (is.data.frame(cluster_vec)) cluster_vec <- cluster_vec[[1]]
  cl_chr <- as.character(cluster_vec)
  cl_ord <- order(suppressWarnings(as.integer(cl_chr)), cl_chr, na.last = TRUE)
  df <- data.frame(cluster = cl_chr,
                   sample = as.character(so@meta.data[["sample"]]),
                   donor = as.character(so@meta.data[["donor"]]),
                   condition = as.character(so@meta.data[["condition"]]),
                   stringsAsFactors = FALSE)
  ns <- length(unique(df$sample))
  do.call(rbind, lapply(unique(df$cluster)[order(unique(suppressWarnings(
    as.integer(unique(df$cluster)))))], function(k) {
    d <- df[df$cluster == k, ]
    p <- table(d$sample) / nrow(d)
    data.frame(resolution = label, cluster = k, n_cells = nrow(d),
               n_samples = sum(table(d$sample) > 0),
               n_donors = length(unique(d$donor)),
               max_sample_pct = round(max(p) * 100, 2),
               sample_entropy_norm = round(-sum(p[p > 0] * log(p[p > 0])) / log(ns), 4),
               pct_disease = round(mean(d$condition == "Disease") * 100, 2))
  }))
}
cm_pre  <- comp_metrics(so@meta.data[["preint_cluster"]], so, "preint_0.5")
cm_prim <- comp_metrics(so@meta.data[[paste0("harmony_snn_res.", PRIMARY_RES)]], so,
                        paste0("harmony_", PRIMARY_RES))
write.csv(bind_rows(cm_pre, cm_prim),
          file.path(tab_d, paste0("03_cluster_composition_metrics", tag, ".csv")),
          row.names = FALSE)

# Cramer's V: how strongly clustering tracks sample vs condition (over-correction check)
cramers_v <- function(a, b) {
  tb <- table(a, b); n <- sum(tb)
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (n * (min(dim(tb)) - 1))))
}
cv <- data.frame(
  clustering = c("preint_0.5", paste0("harmony_", PRIMARY_RES)),
  V_vs_sample = c(cramers_v(so@meta.data[["preint_cluster"]], so@meta.data[["sample"]]),
                  cramers_v(so@meta.data[[paste0("harmony_snn_res.", PRIMARY_RES)]],
                            so@meta.data[["sample"]])),
  V_vs_donor = c(cramers_v(so@meta.data[["preint_cluster"]], so@meta.data[["donor"]]),
                 cramers_v(so@meta.data[[paste0("harmony_snn_res.", PRIMARY_RES)]],
                           so@meta.data[["donor"]])),
  V_vs_condition = c(cramers_v(so@meta.data[["preint_cluster"]], so@meta.data[["condition"]]),
                     cramers_v(so@meta.data[[paste0("harmony_snn_res.", PRIMARY_RES)]],
                               so@meta.data[["condition"]])))
write.csv(cv, file.path(tab_d, paste0("03_cramersV_before_after", tag, ".csv")), row.names = FALSE)
lg("Cramer's V (sample / donor / condition):")
for (i in seq_len(nrow(cv)))
  lg(sprintf("  %-14s %.3f / %.3f / %.3f", cv$clustering[i],
             cv$V_vs_sample[i], cv$V_vs_donor[i], cv$V_vs_condition[i]))

# per-sample cluster distribution (row = sample, col = cluster, %)
PCOL <- paste0("harmony_snn_res.", PRIMARY_RES)
pc <- as.data.frame.matrix(table(as.character(so@meta.data[["sample"]]),
                                 as.character(so@meta.data[[PCOL]])))
pc_pct <- round(sweep(as.matrix(pc), 1, rowSums(as.matrix(pc)), "/") * 100, 2)
write.csv(data.frame(sample = rownames(pc_pct), pc_pct, check.names = FALSE),
          file.path(tab_d, paste0("03_sample_by_cluster_pct", tag, ".csv")), row.names = FALSE)
write.csv(data.frame(cluster = colnames(pc), t(pc), check.names = FALSE),
          file.path(tab_d, paste0("03_cluster_by_sample_counts", tag, ".csv")), row.names = FALSE)

# per-cluster donor / condition composition (primary resolution)
pc2 <- as.data.frame.matrix(table(as.character(so@meta.data[[PCOL]]),
                                  as.character(so@meta.data[["donor"]])))
pc2 <- data.frame(cluster = rownames(pc2), pc2, check.names = FALSE)
write.csv(pc2, file.path(tab_d, paste0("03_cluster_by_donor_counts", tag, ".csv")), row.names = FALSE)
pc3 <- as.data.frame.matrix(table(as.character(so@meta.data[[PCOL]]),
                                  as.character(so@meta.data[["condition"]])))
write.csv(data.frame(cluster = rownames(pc3), pc3, check.names = FALSE),
          file.path(tab_d, paste0("03_cluster_by_condition_counts", tag, ".csv")), row.names = FALSE)

# ------------------------------------------------------------ 5. figures ----
Idents(so) <- PCOL
so$integrated_cluster <- factor(so@meta.data[[PCOL]])

pal25 <- tryCatch(DiscretePalette(n = nlevels(so$sample), palette = "glasbey"),
                  error = function(e) DiscretePalette(n = nlevels(so$sample), palette = "polychrome"))
names(pal25) <- levels(so$sample)
cols_cond <- c(Control = "#4C72B0", Disease = "#C44E52")

q1 <- DimPlot(so, reduction = "umap.harmony", group.by = "sample", cols = pal25,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Harmony UMAP - sample") + theme(legend.text = element_text(size = 6))
q2 <- DimPlot(so, reduction = "umap.harmony", group.by = "condition", cols = cols_cond,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Harmony UMAP - condition")
ggsave(file.path(fig_d, paste0("umap_harmony_sample_condition", tag, ".png")),
       q1 | q2, width = 15, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_harmony_sample_condition", tag, ".pdf")),
       q1 | q2, width = 15, height = 6.4)

q3 <- DimPlot(so, reduction = "umap.harmony", group.by = "integrated_cluster",
              label = TRUE, repel = TRUE, pt.size = PT, alpha = 0.8,
              raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle(paste0("Harmony UMAP - cluster (res ", PRIMARY_RES, ")")) + NoLegend()
ggsave(file.path(fig_d, paste0("umap_harmony_cluster", tag, ".png")), q3,
       width = 7.4, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_harmony_cluster", tag, ".pdf")), q3,
       width = 7.4, height = 6.4)

## before / after side by side
r1 <- DimPlot(so, reduction = "umap.unintegrated", group.by = "sample", cols = pal25,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Before integration (PCA) - sample") + theme(legend.text = element_text(size = 6))
r2 <- DimPlot(so, reduction = "umap.harmony", group.by = "sample", cols = pal25,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("After integration (Harmony) - sample") + theme(legend.text = element_text(size = 6))
ggsave(file.path(fig_d, paste0("umap_before_after_sample", tag, ".png")),
       r1 | r2, width = 15, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_before_after_sample", tag, ".pdf")),
       r1 | r2, width = 15, height = 6.4)

s1 <- DimPlot(so, reduction = "umap.unintegrated", group.by = "condition", cols = cols_cond,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Before integration (PCA) - condition")
s2 <- DimPlot(so, reduction = "umap.harmony", group.by = "condition", cols = cols_cond,
              pt.size = PT, alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("After integration (Harmony) - condition")
ggsave(file.path(fig_d, paste0("umap_before_after_condition", tag, ".png")),
       s1 | s2, width = 11, height = 5.6, dpi = 300)
ggsave(file.path(fig_d, paste0("umap_before_after_condition", tag, ".pdf")),
       s1 | s2, width = 11, height = 5.6)

## multi-resolution small multiples
pl <- lapply(res_cols, function(cc) {
  DimPlot(so, reduction = "umap.harmony", group.by = cc, label = TRUE, label.size = 2.2,
          pt.size = PT * 0.6, alpha = 0.8, raster = TRUE, raster.dpi = c(250, 250)) +
    ggtitle(sub("harmony_snn_res.", "res = ", cc)) + NoLegend() + theme(plot.title = element_text(size = 10))
})
ggsave(file.path(fig_d, paste0("umap_harmony_resolutions", tag, ".png")),
       wrap_plots(pl, nrow = 2), width = 11, height = 10, dpi = 300)

## cluster size + ARI overview
pa <- ggplot(sz, aes(resolution, n_clusters)) + geom_line() + geom_point() +
  geom_vline(xintercept = PRIMARY_RES, linetype = 2, colour = "red") +
  labs(x = "resolution", y = "number of clusters", title = "Clusters vs resolution") +
  theme_bw(base_size = 11)
pb <- ggplot(stab, aes(factor(res1), factor(res2), fill = ARI)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", ARI)), size = 3) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0.5,
                       limits = c(0, 1)) +
  labs(x = "resolution", y = "resolution", title = "ARI between resolutions") + theme_bw(base_size = 11)
ggsave(file.path(fig_d, paste0("resolution_stability", tag, ".png")), pa | pb,
       width = 10, height = 4.4, dpi = 300)

## sample composition stacked bar (primary resolution)
comp_df <- as.data.frame(table(cluster = so$integrated_cluster, sample = so$sample))
comp_df <- comp_df %>% group_by(cluster) %>% mutate(pct = Freq / sum(Freq) * 100) %>% ungroup()
pd <- ggplot(comp_df, aes(cluster, pct, fill = sample)) + geom_col(width = 0.85) +
  scale_fill_manual(values = pal25) +
  labs(x = "cluster", y = "% of cells in cluster", title = "Sample composition per cluster") +
  theme_bw(base_size = 10) + theme(legend.text = element_text(size = 5), legend.key.size = unit(0.3, "cm"))
ggsave(file.path(fig_d, paste0("cluster_sample_composition", tag, ".png")), pd,
       width = 11, height = 5.2, dpi = 300)

## sample x cluster heatmap
hm <- pc_pct
pheatmap::pheatmap(hm, cluster.rows = TRUE, cluster.cols = TRUE,
                   color = colorRampPalette(c("white", "#FEE0B6", "#D7191C"))(100),
                   fontsize_row = 6, fontsize_col = 7,
                   filename = file.path(fig_d, paste0("sample_cluster_heatmap", tag, ".png")),
                   width = 9, height = 6.5, main = "Row % of each sample per cluster")

saveRDS(so, f_out, compress = FALSE)
lg("saved: ", basename(f_out), sprintf(" (%.2f GB)", file.size(f_out) / 1024^3))
writeLines(capture.output(sessionInfo()),
           file.path(log_d, paste0("sessionInfo_42", tag, ".txt")))
lg("=== 42 END ===")
close(lgf)
