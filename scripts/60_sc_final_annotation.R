# =============================================================================
# 60_sc_final_annotation.R
# Stage 4 - final major cell-type annotation.
#
# Decisions baked in (2026-10-29)
#   * cluster numbering is preserved: the original `cluster` /
#     `harmony_snn_res.0.6` columns are NOT touched, RenameIdents is NOT used
#   * a new `celltype` column carries the label; no existing column is overwritten
#   * C8 stays "Uncertain" (never auto-filled)
#   * C16 is removed from the object (151,007 -> 150,750 nuclei)
#   * `group` (NF / NICM / DoxCM) is added from the sample manifest;
#     `condition` (Control / Disease) is left as it was
#
# Outputs
#   results/Singlecel/objects/03_merged/05_final_annotated.rds
#   results/annotation/tables/final_celltype_mapping.csv
#   results/annotation/tables/celltype_counts_by_sample.csv
#   results/annotation/tables/celltype_proportions_by_sample.csv
#   results/annotation/tables/celltype_counts_by_group.csv
#   results/annotation/tables/celltype_proportions_by_group.csv
#   results/annotation/tables/celltype_overall.csv
#   results/annotation/tables/removed_cells_C16.csv
#   results/annotation/figures/final_umap_celltype.{png,pdf}
#   results/annotation/figures/final_umap_by_sample.{png,pdf}
#   results/annotation/figures/final_umap_sample_split.{png,pdf}
#   results/annotation/figures/final_umap_by_group.{png,pdf}
#   results/annotation/figures/final_composition_by_sample.{png,pdf}
#   results/annotation/figures/final_composition_by_group.{png,pdf}
#   results/annotation/figures/final_counts_by_sample.{png,pdf}
#   results/annotation/figures/final_dotplot_markers.{png,pdf}
#   results/annotation/figures/final_marker_heatmap.{png,pdf}
#   results/annotation/60_final_annotation.log
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork); library(scales)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
obj_d <- file.path(sdir, "objects", "03_merged")
adr   <- file.path(proj, "results", "annotation")
tab_d <- file.path(adr, "tables")
fig_d <- file.path(adr, "figures")
for (d in c(tab_d, fig_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file.path(adr, "60_final_annotation.log")
if (file.exists(lgf)) file.remove(lgf)
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = "")
  cat(m, "\n", sep = "", file = lgf, append = TRUE)
}
lg("=== 60 final annotation | seed ", SEED, " ===")

# ---------------------------------------------------------------- load ------
t0 <- Sys.time()
so <- readRDS(file.path(obj_d, "04_integrated.rds"))
lg(sprintf("loaded %d genes x %d cells in %.1f min", nrow(so), ncol(so),
           as.numeric(difftime(Sys.time(), t0, units = "mins"))))

CL_COL <- "harmony_snn_res.0.6"
stopifnot(CL_COL %in% colnames(so@meta.data))
# the object only carries the resolution column; add a plain `cluster` copy
# (new column, nothing overwritten) so downstream code has a stable name
if (!"cluster" %in% colnames(so@meta.data)) {
  so@meta.data[["cluster"]] <- as.character(so@meta.data[[CL_COL]])
  lg("added `cluster` column as a copy of ", CL_COL)
}
cl <- as.character(so@meta.data[[CL_COL]])
lg("original clusters: ", paste(sort(unique(as.integer(cl))), collapse = ","))

# ------------------------------------------------------- label mapping ------
# label, confidence and decision come from the reviewed evidence tables
map <- data.frame(
  cluster = as.character(0:16),
  celltype = c("Cardiomyocyte", "Fibroblast", "Endothelial cell", "Pericyte",
               "Fibroblast", "Macrophage", "Cardiomyocyte", "Endothelial cell",
               "Uncertain", "Monocyte/DC", "T/NK cell", "Smooth muscle cell",
               "Neural/glial cell", "Lymphatic endothelial cell", "Mast cell",
               "Proliferating cell", "Adipocyte"),
  confidence = c("High","High","High","High","High","High","High","High",
                 "Low","Moderate","High","High","High","High","High",
                 "Moderate","Excluded"),
  include = c(rep(TRUE, 16), FALSE),
  reason = c(
    "CM panel 0.89, sarcomere programme complete",
    "FB panel 0.59, PDGFRA/DCN/ABCA9 positive",
    "EC panel 0.57, BTNL9/CA4 capillary-arterial",
    "PC panel 0.67, RGS5/HIGD1B/EGFLAM positive",
    "FB panel 0.71, MFAP5/SERPINE2 activated",
    "MYE panel 0.52, MRC1/CD163/F13A1 resident",
    "CM panel 0.77, MYBPC3/RBM20/MYL2 ventricular",
    "EC panel 0.50, ACKR1 0.23 venous",
    "decision postponed by user; endocardial-like EC but not confirmed",
    "MYE panel 0.35, LYZ/FCN1/ITGAX non-resident",
    "TNK panel 0.38, THEMIS/BCL11B/ITK/SKAP1",
    "SMC panel 0.59, MYH11/MYOCD/LMOD1",
    "NEUR panel 0.65, SOX10/PLP1/PTPRZ1; CDH19 withdrawn (also 0.82-0.90 in fibroblasts)",
    "LEC panel 0.62, PROX1/CCL21/FLT4/MMRN1",
    "MAST panel 0.62, KIT/CPA3/TPSB2/MS4A2",
    "PROL panel 0.66, cycling myeloid-biased",
    "user decision: remove (257 nuclei, 18/25 samples)"),
  stringsAsFactors = FALSE)

stopifnot(all(cl %in% map$cluster))

# ------------------------------------------------- add columns (no overwrite) --
so$celltype <- factor(map$celltype[match(cl, map$cluster)],
                      levels = c("Cardiomyocyte","Fibroblast","Endothelial cell",
                                 "Pericyte","Smooth muscle cell","Macrophage",
                                 "Monocyte/DC","T/NK cell","Lymphatic endothelial cell",
                                 "Mast cell","Neural/glial cell","Proliferating cell",
                                 "Uncertain"))
so$cluster_id <- as.character(so@meta.data[[CL_COL]])   # explicit numeric cluster id
so$celltype_confidence <- map$confidence[match(cl, map$cluster)]

# group from the deposited manifest (NF / NICM / DoxCM); condition untouched
man <- read.csv(file.path(proj, "clean_data", "GSE292067_sample_manifest.csv"),
                stringsAsFactors = FALSE)
grp <- setNames(man$group, man$sample)
# unname(): a named vector would be read as cell names by AddMetaData
so$group <- factor(unname(grp[as.character(so$sample)]),
                   levels = c("NF", "NICM", "DoxCM"))
lg("group table:"); print(table(so$group, useNA = "ifany"))
stopifnot(!any(is.na(so$group)))

# ------------------------------------------------------- remove C16 ---------
rm_cells <- colnames(so)[cl == "16"]
write.csv(data.frame(cell = rm_cells,
                     cluster = "16",
                     reason = "adipocyte cluster removed by user decision 2026-10-29",
                     stringsAsFactors = FALSE),
          file.path(tab_d, "removed_cells_C16.csv"), row.names = FALSE)
lg("removing C16: ", length(rm_cells), " nuclei")
keep <- setdiff(colnames(so), rm_cells)
so <- subset(so, cells = keep)
lg("object now: ", ncol(so), " nuclei x ", nrow(so), " genes")

# drop heavy slots that are not needed downstream (keeps the saved object small)
so <- tryCatch(
  DietSeurat(so, assays = "RNA", layers = c("counts", "data"),
             dimreducs = c("pca", "umap.harmony")),
  error = function(e) { lg("DietSeurat failed: ", conditionMessage(e)); so })
lg("layers after diet: ", paste(Layers(so[["RNA"]]), collapse = ","),
   " | graphs: ", length(so@graphs),
   " | reductions: ", paste(Reductions(so), collapse = ","))
invisible(gc(verbose = FALSE))

# ------------------------------------------------------- composition --------
md <- so@meta.data
ct_levels <- levels(so$celltype)
md$celltype <- factor(as.character(md$celltype), levels = ct_levels)

cnt_by_sample <- as.data.frame.matrix(table(md$sample, md$celltype))
cnt_by_sample <- data.frame(sample = rownames(cnt_by_sample),
                            group = grp[rownames(cnt_by_sample)],
                            cnt_by_sample, check.names = FALSE, row.names = NULL)
prop_by_sample <- cnt_by_sample
prop_by_sample[, ct_levels] <- round(sweep(as.matrix(cnt_by_sample[, ct_levels]),
                                           1, rowSums(cnt_by_sample[, ct_levels]),
                                           "/") * 100, 3)

cnt_by_group <- as.data.frame.matrix(table(md$group, md$celltype))
cnt_by_group <- data.frame(group = rownames(cnt_by_group), cnt_by_group,
                           check.names = FALSE, row.names = NULL)
prop_by_group <- cnt_by_group
prop_by_group[, ct_levels] <- round(sweep(as.matrix(cnt_by_group[, ct_levels]),
                                          1, rowSums(cnt_by_group[, ct_levels]),
                                          "/") * 100, 3)

overall <- data.frame(celltype = ct_levels,
                      n_cells = as.integer(table(md$celltype)[ct_levels]),
                      stringsAsFactors = FALSE)
overall$pct_of_all <- round(overall$n_cells / ncol(so) * 100, 3)
overall$confidence <- map$confidence[match(overall$celltype,
                                           map$celltype[map$include])]
overall$n_samples <- sapply(ct_levels, function(k) length(unique(md$sample[md$celltype == k])))

write.csv(cnt_by_sample,  file.path(tab_d, "celltype_counts_by_sample.csv"), row.names = FALSE)
write.csv(prop_by_sample, file.path(tab_d, "celltype_proportions_by_sample.csv"), row.names = FALSE)
write.csv(cnt_by_group,   file.path(tab_d, "celltype_counts_by_group.csv"), row.names = FALSE)
write.csv(prop_by_group,  file.path(tab_d, "celltype_proportions_by_group.csv"), row.names = FALSE)
write.csv(overall,        file.path(tab_d, "celltype_overall.csv"), row.names = FALSE)
lg("composition tables written")
print(overall, row.names = FALSE)
print(prop_by_group, row.names = FALSE)

# mapping table (documentation of what was done to each cluster)
final_map <- data.frame(
  cluster = map$cluster,
  original_cluster_column = CL_COL,
  n_cells_before = as.integer(table(factor(cl, levels = map$cluster))),
  celltype = map$celltype,
  confidence = map$confidence,
  included = map$include,
  reason = map$reason,
  stringsAsFactors = FALSE)
write.csv(final_map, file.path(tab_d, "final_celltype_mapping.csv"), row.names = FALSE)

# ------------------------------------------------------------- save ---------
out_obj <- file.path(obj_d, "05_final_annotated.rds")
t1 <- Sys.time()
saveRDS(so, out_obj)
lg(sprintf("saved %s (%.1f MB) in %.1f min",
           basename(out_obj), file.size(out_obj) / 1024^2,
           as.numeric(difftime(Sys.time(), t1, units = "mins"))))

# ------------------------------------------------------------- plots --------
pal <- c("Cardiomyocyte"="#D1495B","Fibroblast"="#3E7CB1","Endothelial cell"="#2A9D8F",
         "Pericyte"="#8E6C8A","Smooth muscle cell"="#6A994E","Macrophage"="#E9A03B",
         "Monocyte/DC"="#C97B4A","T/NK cell"="#4C5FD7","Lymphatic endothelial cell"="#7BC8B0",
         "Mast cell"="#B5527A","Neural/glial cell"="#8C8C3A","Proliferating cell"="#9B5DE5",
         "Uncertain"="#BDBDBD")

um <- Embeddings(so, "umap.harmony")
pd <- data.frame(UMAP_1 = um[, 1], UMAP_2 = um[, 2],
                 celltype = md$celltype, sample = md$sample, group = md$group)
cent <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = pd, FUN = median)

p_ct <- ggplot(pd, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = celltype), size = 0.05, stroke = 0) +
  geom_label(data = cent, aes(label = celltype), size = 2.5, label.size = 0,
             fill = "white", alpha = 0.75, label.padding = unit(0.08, "lines")) +
  scale_colour_manual(values = pal, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "Final cell-type annotation (major types)",
       subtitle = sprintf("n = %s nuclei; C8 retained as Uncertain; C16 removed",
                          format(ncol(so), big.mark = ",")),
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 7.5),
        legend.key.height = unit(0.36, "cm"))
ggsave(file.path(fig_d, "final_umap_celltype.png"), p_ct, width = 9.6, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, "final_umap_celltype.pdf"), p_ct, width = 9.6, height = 6.4)

# by sample (single panel, 25 colours) - draw largest cell types first
ord <- pd[order(-as.integer(table(pd$celltype)[pd$celltype])), ]
p_s <- ggplot(ord, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = sample), size = 0.05, stroke = 0) +
  scale_colour_manual(values = scales::hue_pal()(length(unique(ord$sample))), name = "sample") +
  guides(colour = guide_legend(override.aes = list(size = 2.2), ncol = 2)) +
  labs(title = "Harmony UMAP coloured by sample",
       subtitle = "25 libraries, 151,007 -> 150,750 nuclei after C16 removal",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 6),
        legend.key.height = unit(0.28, "cm"))
ggsave(file.path(fig_d, "final_umap_by_sample.png"), p_s, width = 8.6, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, "final_umap_by_sample.pdf"), p_s, width = 8.6, height = 6.4)

# per-sample small multiples (that sample highlighted)
samp_lev <- unique(pd$sample[order(match(pd$group, c("NF","NICM","DoxCM")), pd$sample)])
pd$sample_f <- factor(pd$sample, levels = samp_lev)
bg <- pd[sample(nrow(pd), 20000), c("UMAP_1", "UMAP_2")]   # thin background for all facets
p_ss <- ggplot() +
  geom_point(data = bg, aes(UMAP_1, UMAP_2), colour = "grey88", size = 0.05, stroke = 0) +
  geom_point(data = pd, aes(UMAP_1, UMAP_2, colour = group), size = 0.05, stroke = 0) +
  facet_wrap(~ sample_f, ncol = 5) +
  scale_colour_manual(values = c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B"),
                      name = "group") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Per-sample contribution to the integrated UMAP",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.5))
ggsave(file.path(fig_d, "final_umap_sample_split.png"), p_ss, width = 11, height = 8.2,
       dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "final_umap_sample_split.pdf"), p_ss, width = 11, height = 8.2)

# by group
p_g <- ggplot(ord, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = group), size = 0.05, stroke = 0) +
  scale_colour_manual(values = c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B"),
                      name = "group") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Harmony UMAP coloured by group",
       subtitle = "NF = non-failing (8), NICM (6), DoxCM (11)",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_d, "final_umap_by_group.png"), p_g, width = 7.6, height = 6.0, dpi = 300)
ggsave(file.path(fig_d, "final_umap_by_group.pdf"), p_g, width = 7.6, height = 6.0)

# ------------------------------------------------- composition bars ---------
long_s <- do.call(rbind, lapply(ct_levels, function(k) {
  data.frame(sample = cnt_by_sample$sample, group = cnt_by_sample$group,
             celltype = k,
             n = cnt_by_sample[[k]],
             pct = prop_by_sample[[k]],
             stringsAsFactors = FALSE)
}))
long_s$celltype <- factor(long_s$celltype, levels = ct_levels)
long_s$sample <- factor(long_s$sample, levels = samp_lev)
write.csv(long_s, file.path(tab_d, "celltype_composition_long.csv"), row.names = FALSE)

p_cs <- ggplot(long_s, aes(sample, pct, fill = celltype)) +
  geom_col(width = 0.82) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Cell-type composition per sample",
       subtitle = "proportions are descriptive; nuclei are not biological replicates",
       x = NULL, y = "% of nuclei in sample") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 6.8),
        legend.key.height = unit(0.32, "cm"),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "final_composition_by_sample.png"), p_cs, width = 10.4, height = 6.6, dpi = 300)
ggsave(file.path(fig_d, "final_composition_by_sample.pdf"), p_cs, width = 10.4, height = 6.6)

long_g <- do.call(rbind, lapply(ct_levels, function(k) {
  data.frame(group = cnt_by_group$group, celltype = k,
             n = cnt_by_group[[k]], pct = prop_by_group[[k]], stringsAsFactors = FALSE)
}))
long_g$celltype <- factor(long_g$celltype, levels = ct_levels)
long_g$group <- factor(long_g$group, levels = c("NF","NICM","DoxCM"))

p_cg <- ggplot(long_g, aes(group, pct, fill = celltype)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.1f", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.4, colour = "white") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Cell-type composition per group",
       subtitle = "pooled nuclei per group; descriptive only, no test across groups",
       x = NULL, y = "% of nuclei in group") +
  theme_bw(base_size = 9.5) +
  theme(panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.34, "cm"))
ggsave(file.path(fig_d, "final_composition_by_group.png"), p_cg, width = 8.4, height = 5.8, dpi = 300)
ggsave(file.path(fig_d, "final_composition_by_group.pdf"), p_cg, width = 8.4, height = 5.8)

p_ns <- ggplot(long_s, aes(sample, n, fill = celltype)) +
  geom_col(width = 0.82) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = comma) +
  labs(title = "Nuclei per sample by cell type", x = NULL, y = "nuclei") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 6.8),
        legend.key.height = unit(0.32, "cm"),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "final_counts_by_sample.png"), p_ns, width = 10.4, height = 6.6, dpi = 300)
ggsave(file.path(fig_d, "final_counts_by_sample.pdf"), p_ns, width = 10.4, height = 6.6)

# ------------------------------------------------------ marker figures ------
dot_genes <- list(
  CM    = c("TTN","MYH7","TNNT2","MYBPC3","RBM20","CHRM2","ALPK2"),
  FB    = c("DCN","LUM","COL1A1","PDGFRA","MFAP5","SERPINE2"),
  EC    = c("PECAM1","CDH5","VWF","FLT1","EMCN","BTNL9","ACKR1"),
  LEC   = c("CCL21","PROX1","MMRN1","FLT4","PKHD1L1"),
  PC    = c("RGS5","KCNJ8","ABCC9","PDGFRB","HIGD1B","EGFLAM"),
  SMC   = c("MYH11","ACTA2","MYOCD","LMOD1","CARMN"),
  MYE   = c("PTPRC","LYZ","CD68","CSF1R","MRC1","CD163","F13A1","ITGAX","FCN1"),
  TNK   = c("CD3D","CD3E","CD2","IL7R","NKG7","SKAP1","THEMIS"),
  MAST  = c("CPA3","TPSB2","KIT","MS4A2","HDC"),
  NEUR  = c("SOX10","PLP1","MPZ","PTPRZ1","CDH19"),
  PROL  = c("MKI67","TOP2A","BUB1B"),
  ENDO  = c("NPR3","NRG1","NFATC1"))

gs_map <- so[["RNA"]][[]][["gene_symbol"]]
ens <- rownames(so)[match(unlist(dot_genes), gs_map)]
names(ens) <- unlist(dot_genes)
ens <- ens[!is.na(ens)]
lg("dotplot genes resolved: ", length(ens), " / ", length(unlist(dot_genes)))

dp <- DotPlot(so, features = names(ens), group.by = "celltype", assay = "RNA")
dd <- dp$data
dd$features.plot <- factor(as.character(dd$features.plot), levels = rev(names(ens)))
dd$id <- factor(as.character(dd$id), levels = ct_levels)
p_dot <- ggplot(dd, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027", midpoint = 0,
                         name = "scaled mean\nexpression") +
  scale_size_area(max_size = 4, name = "% expressing") +
  labs(title = "Canonical marker panel across final cell types",
       subtitle = "C8 shown as Uncertain; C16 removed; x-axis ordered by cell type",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.4),
        axis.text.y = element_text(size = 6.4, face = "italic"),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_d, "final_dotplot_markers.png"), p_dot, width = 7.6,
       height = 9.6, dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "final_dotplot_markers.pdf"), p_dot, width = 7.6, height = 9.6)

# marker heatmap: mean scaled expression per cell type (computed directly to
# avoid AverageExpression group-name sanitising)
dat <- GetAssayData(so, assay = "RNA", layer = "data")[names(ens), , drop = FALSE]
avg <- sapply(ct_levels, function(k) Matrix::rowMeans(dat[, md$celltype == k, drop = FALSE]))
colnames(avg) <- ct_levels
z <- t(scale(t(avg)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, 2), -2)
hz <- data.frame(gene = rep(rownames(z), ncol(z)),
                 celltype = rep(colnames(z), each = nrow(z)),
                 z = as.vector(z), stringsAsFactors = FALSE)
hz$gene <- factor(hz$gene, levels = rev(rownames(z)))
hz$celltype <- factor(hz$celltype, levels = ct_levels)
p_hm <- ggplot(hz, aes(celltype, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"),
                       limits = c(-2, 2), name = "scaled mean\nexpression") +
  labs(title = "Marker expression heatmap by final cell type", x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.6),
        axis.text.y = element_text(size = 6.2, face = "italic"),
        panel.grid = element_blank())
ggsave(file.path(fig_d, "final_marker_heatmap.png"), p_hm, width = 8.6, height = 7.8, dpi = 300)
ggsave(file.path(fig_d, "final_marker_heatmap.pdf"), p_hm, width = 8.6, height = 7.8)

# ------------------------------------------------------------- session ------
writeLines(capture.output(sessionInfo()), file.path(adr, "sessionInfo_60.txt"))
writeLines(c(sprintf("seed = %d", SEED),
             sprintf("date = %s", Sys.Date()),
             sprintf("nuclei_before = %d", ncol(so) + length(rm_cells)),
             sprintf("nuclei_after = %d", ncol(so)),
             sprintf("removed = %d (C16)", length(rm_cells))),
           file.path(adr, "60_run_meta.txt"))
lg("=== done ===")
