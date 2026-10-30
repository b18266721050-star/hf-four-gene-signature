# =============================================================================
# 43_sc_stage2_markers.R
# Stage 2 - Step 4: JoinLayers (official Seurat v5 flow) and marker discovery
# on the RNA assay, with per-cluster tables and DotPlot / FeaturePlot / heatmap.
#
# Env vars:
#   SC_SUBSET=n  -> validation mode
#   SC_RES=0.6   -> cluster column used (must match script 42)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork); library(dplyr)
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
  fig_d <- file.path(fig_d, "04_markers_validation")
  tab_d <- file.path(sdir, "tables", "validation")
} else {
  obj_d <- file.path(sdir, "objects", "03_merged")
  fig_d <- file.path(fig_d, "04_markers")
}
for (d in c(fig_d, tab_d, log_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

PRIMARY_RES <- as.numeric(Sys.getenv("SC_RES", "0.6"))
MAX_CELLS_PER_IDENT <- 2000

lgf <- file(file.path(log_d, paste0("43_markers", tag, ".log")), "wt")
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf); flush(lgf)
}
options(error = function() {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ERROR: ", geterrmessage(), "\n"),
      file = lgf, append = TRUE); flush(lgf)
  q(status = 1)
})

f_in  <- file.path(obj_d, paste0("04_integrated", tag, ".rds"))
f_out <- file.path(obj_d, paste0("05_markers", tag, ".rds"))

lg("=== 43 markers (seed=", SEED, ", res=", PRIMARY_RES, ") ===")
lg("presto available: ", requireNamespace("presto", quietly = TRUE),
   " (used automatically by Seurat wilcox for speed)")
so <- readRDS(f_in)
PT <- if (ncol(so) > 80000) 0.25 else 0.6
lg("loaded: ", ncol(so), " cells | layers: ", paste(Layers(so[["RNA"]]), collapse = ","))

# ------------------------------------------------ 1. JoinLayers (official) ---
## only per-sample layers (counts.GSM..., data.GSM...) need joining; a plain
## counts + data pair is already single-layer and joining it would copy the matrix
per_sample <- grep("^(counts|data)\\.GSM", Layers(so[["RNA"]]), value = TRUE)
if (length(per_sample) > 0) {
  so <- JoinLayers(so)
  lg("JoinLayers applied (", length(per_sample), " per-sample layers) -> ",
     paste(Layers(so[["RNA"]]), collapse = ","))
} else {
  lg("layers are already single (", paste(Layers(so[["RNA"]]), collapse = ","),
     "), JoinLayers is a no-op")
}

CL_COL <- paste0("harmony_snn_res.", PRIMARY_RES)
cl_chr <- as.character(so@meta.data[[CL_COL]])
lv <- as.character(sort(unique(as.integer(cl_chr))))
so$cluster <- factor(cl_chr, levels = lv)
Idents(so) <- "cluster"
lg("clusters: ", nlevels(so$cluster), " | sizes: ",
   paste(as.integer(table(so$cluster)), collapse = "/"))

## marker testing only needs the normalised data layer; dropping counts from the
## in-memory copy halves the footprint (the saved objects keep both layers)
if ("counts" %in% Layers(so[["RNA"]])) {
  ok <- tryCatch({
    so[["RNA"]] <- DietSeurat(so[["RNA"]], layers = "data"); TRUE
  }, error = function(e) FALSE)
  lg("counts layer dropped from memory: ", ok, " | layers now: ",
     paste(Layers(so[["RNA"]]), collapse = ","))
  invisible(gc(verbose = FALSE))
}

# ------------------------------------------------------- 2. FindAllMarkers ---
## loop explicitly instead of FindAllMarkers so each cluster is logged and the
## transient copies are released between clusters
t0 <- Sys.time()
lev <- levels(so$cluster)
res_list <- vector("list", length(lev))
for (i in seq_along(lev)) {
  k <- lev[i]
  m <- FindMarkers(so, ident.1 = k, assay = "RNA", test.use = "wilcox",
                   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                   max.cells.per.ident = MAX_CELLS_PER_IDENT,
                   random.seed = SEED, verbose = FALSE)
  if (nrow(m) > 0) { m$cluster <- k; m$gene <- rownames(m) }
  res_list[[i]] <- m
  lg(sprintf("  [%2d/%d] cluster %-2s -> %5d markers  (%.1f min elapsed)",
             i, length(lev), k, nrow(m),
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(gc(verbose = FALSE))
}
mk <- do.call(rbind, res_list)
rownames(mk) <- NULL
mk$cluster <- factor(mk$cluster, levels = lv)
lg(sprintf("markers done in %.1f min | %d rows, %d genes, %d clusters",
           as.numeric(difftime(Sys.time(), t0, units = "mins")),
           nrow(mk), length(unique(mk$gene)), length(unique(mk$cluster))))
lg("note: max.cells.per.ident=", MAX_CELLS_PER_IDENT,
   " (downsampling for speed; ranking is stable, p-values are approximate)")

ann <- read.csv(file.path(sdir, "tables", "feature_annotation.csv"), stringsAsFactors = FALSE)
mk$gene_symbol <- ann$gene_symbol[match(mk$gene, ann$ensembl_id)]
mk$pct_diff <- mk$pct.1 - mk$pct.2
mk <- mk[order(mk$cluster, -mk$avg_log2FC), ]
write.csv(mk, file.path(tab_d, paste0("04_markers_all", tag, ".csv")), row.names = FALSE)

top10 <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
top3  <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3) %>% ungroup()
write.csv(top10, file.path(tab_d, paste0("04_markers_top10_per_cluster", tag, ".csv")), row.names = FALSE)

summ <- mk %>% group_by(cluster) %>%
  summarise(n_markers = n(),
            n_pct50 = sum(pct.1 >= 0.5),
            top1 = first(gene_symbol), top1_lfc = round(first(avg_log2FC), 2),
            top1_pct = round(first(pct.1), 3), .groups = "drop")
write.csv(summ, file.path(tab_d, paste0("04_marker_summary_per_cluster", tag, ".csv")), row.names = FALSE)
lg("markers per cluster: ", paste(summ$n_markers, collapse = "/"))

# ------------------------------------------------------------- 3. plots -----
sym_map <- setNames(ifelse(is.na(mk$gene_symbol), mk$gene, mk$gene_symbol), mk$gene)
sym_map <- sym_map[!duplicated(names(sym_map))]

## DotPlot, split into pages of 5 clusters
lab_fun <- function(x) { y <- unname(sym_map[x]); ifelse(is.na(y), x, y) }

clusts <- levels(so$cluster)
pages <- split(clusts, ceiling(seq_along(clusts) / 5))
dp_last <- NULL
for (k in seq_along(pages)) {
  cc <- pages[[k]]
  feats <- top3$gene[top3$cluster %in% cc]
  dp <- DotPlot(so, features = feats, group.by = "cluster", assay = "RNA") +
    RotatedAxis() +
    scale_x_discrete(labels = lab_fun) +
    labs(title = paste0("Top 3 markers per cluster (res ", PRIMARY_RES, ") - clusters ",
                        paste(cc, collapse = ",")),
         x = NULL, y = "cluster") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(fig_d, sprintf("dotplot_top3_page%d%s.png", k, tag)), dp,
         width = 13, height = 4.6, dpi = 300)
  dp_last <- dp
}
ggsave(file.path(fig_d, paste0("dotplot_top3_page1", tag, ".pdf")),
       dp_last, width = 13, height = 4.6)

## FeaturePlot of the top 2 markers per cluster, 6 per page
top2 <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 2) %>% ungroup()
fpages <- split(top2$gene, ceiling(seq_along(top2$gene) / 6))
for (k in seq_along(fpages)) {
  fp <- FeaturePlot(so, features = fpages[[k]], ncol = 3, pt.size = PT,
                    alpha = 0.8, raster = TRUE, raster.dpi = c(300, 300), order = TRUE,
                    reduction = "umap.harmony")
  fp <- fp & scale_colour_gradientn(colours = c("grey90", "#FEE0B6", "#D7191C"))
  ggsave(file.path(fig_d, sprintf("featureplot_top2_page%d%s.png", k, tag)), fp,
         width = 11, height = 7, dpi = 250)
}

## heatmap: top 10 per cluster on downsampled cells (needs scale.data)
set.seed(SEED)
sub <- subset(so, downsample = 100)
hm_genes <- unique(top10$gene)
sub <- ScaleData(sub, features = hm_genes, verbose = FALSE)
ph <- DoHeatmap(sub, features = hm_genes, group.by = "cluster", size = 2.6,
                raster = FALSE, angle = 90) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C")
ph_lab <- tryCatch(ph + scale_y_discrete(labels = lab_fun),
                   error = function(e) { lg("heatmap gene relabel skipped: ", conditionMessage(e)); ph })
ggsave(file.path(fig_d, paste0("heatmap_top10", tag, ".png")), ph_lab,
       width = 15, height = 13, dpi = 250)
ggsave(file.path(fig_d, paste0("heatmap_top10", tag, ".pdf")), ph_lab, width = 15, height = 13)
rm(sub); invisible(gc(verbose = FALSE))

saveRDS(mk, f_out)
lg("saved marker table: ", basename(f_out))
writeLines(capture.output(sessionInfo()),
           file.path(log_d, paste0("sessionInfo_43", tag, ".txt")))
lg("=== 43 END ===")
close(lgf)
