# =============================================================================
# 61_sc_final_figures.R
# Stage 4 - figures for the final annotation (reads the annotated object only).
#
# Input : results/Singlecel/objects/03_merged/05_final_annotated.rds
#         results/annotation/tables/celltype_composition_long.csv
# Output: results/annotation/figures/final_*.{png,pdf}
#
# Note: the RNA assay is keyed by Ensembl ID; gene symbols live in the
#       `gene_symbol` feature metadata, so every feature lookup goes through it.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork); library(scales)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
obj_f <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
adr   <- file.path(proj, "results", "annotation")
tab_d <- file.path(adr, "tables")
fig_d <- file.path(adr, "figures")

lgf <- file.path(adr, "61_final_figures.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

so <- readRDS(obj_f)
lg("loaded ", ncol(so), " nuclei x ", nrow(so), " genes")
md <- so@meta.data
ct_levels <- levels(so$celltype)
lg("cell types: ", paste(ct_levels, collapse = " | "))

pal <- c("Cardiomyocyte"="#D1495B","Fibroblast"="#3E7CB1","Endothelial cell"="#2A9D8F",
         "Pericyte"="#6A4C93","Smooth muscle cell"="#6A994E","Macrophage"="#E9A03B",
         "Monocyte/DC"="#C97B4A","T/NK cell"="#4C5FD7","Lymphatic endothelial cell"="#7BC8B0",
         "Mast cell"="#B5527A","Neural/glial cell"="#8C8C3A","Proliferating cell"="#9B5DE5",
         "Uncertain"="#9E9E9E")

um <- Embeddings(so, "umap.harmony")
pd <- data.frame(UMAP_1 = um[, 1], UMAP_2 = um[, 2],
                 celltype = md$celltype, sample = md$sample,
                 group = factor(md$group, levels = c("NF", "NICM", "DoxCM")))
cent <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = pd, FUN = median)
# only label cell types that occupy a visible area; rare types stay in the legend
cent$pct <- as.numeric(table(pd$celltype)[as.character(cent$celltype)]) / nrow(pd) * 100
cent_lab <- cent[cent$pct >= 2, ]

# ------------------------------------------------------ 1. celltype UMAP ----
p_ct <- ggplot(pd, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = celltype), size = 0.05, stroke = 0) +
  geom_label(data = cent_lab, aes(label = celltype), size = 2.3, linewidth = 0,
             fill = "white", alpha = 0.8, label.padding = unit(0.07, "lines")) +
  scale_colour_manual(values = pal, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "Final cell-type annotation",
       subtitle = sprintf("n = %s nuclei; 12 major types; C16 removed; labels shown for types >= 2%%",
                          format(ncol(so), big.mark = ",")),
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 7.5),
        legend.key.height = unit(0.36, "cm"))
ggsave(file.path(fig_d, "final_umap_celltype.png"), p_ct, width = 9.6, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, "final_umap_celltype.pdf"), p_ct, width = 9.6, height = 6.4)
lg("umap celltype done")

# -------------------------------------------------------- 2. by sample ------
ord <- pd[order(-as.integer(table(pd$celltype)[pd$celltype])), ]
p_s <- ggplot(ord, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = sample), size = 0.05, stroke = 0) +
  scale_colour_manual(values = scales::hue_pal()(length(unique(ord$sample))), name = "sample") +
  guides(colour = guide_legend(override.aes = list(size = 2.2), ncol = 2)) +
  labs(title = "Harmony UMAP coloured by sample",
       subtitle = "25 libraries; 151,007 -> 150,750 nuclei after C16 removal",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) +
  theme(legend.text = element_text(size = 6),
        legend.key.height = unit(0.28, "cm"))
ggsave(file.path(fig_d, "final_umap_by_sample.png"), p_s, width = 8.8, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, "final_umap_by_sample.pdf"), p_s, width = 8.8, height = 6.4)
lg("umap by sample done")

# per-sample small multiples
samp_lev <- unique(pd$sample[order(match(pd$group, c("NF","NICM","DoxCM")), pd$sample)])
pd$sample_f <- factor(pd$sample, levels = samp_lev)
bg <- pd[sample(nrow(pd), 20000), c("UMAP_1", "UMAP_2")]
p_ss <- ggplot() +
  geom_point(data = bg, aes(UMAP_1, UMAP_2), colour = "grey88", size = 0.05, stroke = 0) +
  geom_point(data = pd, aes(UMAP_1, UMAP_2, colour = group), size = 0.05, stroke = 0) +
  facet_wrap(~ sample_f, ncol = 5) +
  scale_colour_manual(values = c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B"),
                      name = "group") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Per-sample contribution to the integrated UMAP",
       subtitle = "grey = all nuclei; coloured = nuclei of the facet sample",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.4))
ggsave(file.path(fig_d, "final_umap_sample_split.png"), p_ss, width = 11, height = 8.4,
       dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "final_umap_sample_split.pdf"), p_ss, width = 11, height = 8.4)
lg("umap sample split done")

# ---------------------------------------------------------- 3. by group -----
p_g <- ggplot(ord, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = group), size = 0.05, stroke = 0) +
  scale_colour_manual(values = c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B"),
                      name = "group") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Harmony UMAP coloured by group",
       subtitle = "NF = non-failing (8 samples), NICM (6), DoxCM (11)",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_d, "final_umap_by_group.png"), p_g, width = 7.8, height = 6.0, dpi = 300)
ggsave(file.path(fig_d, "final_umap_by_group.pdf"), p_g, width = 7.8, height = 6.0)
lg("umap by group done")

# ------------------------------------------------------ 4. composition ------
long_s <- read.csv(file.path(tab_d, "celltype_composition_long.csv"), stringsAsFactors = FALSE)
long_s$celltype <- factor(long_s$celltype, levels = ct_levels)
long_s$sample <- factor(long_s$sample, levels = samp_lev)
long_s$group <- factor(long_s$group, levels = c("NF", "NICM", "DoxCM"))

p_cs <- ggplot(long_s, aes(sample, pct, fill = celltype)) +
  geom_col(width = 0.82) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Cell-type composition per sample",
       subtitle = "proportions are descriptive; nuclei are not treated as biological replicates",
       x = NULL, y = "% of nuclei in sample") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 6.8),
        legend.key.height = unit(0.32, "cm"),
        plot.margin = margin(4, 4, 16, 4),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "final_composition_by_sample.png"), p_cs, width = 10.4, height = 7.8, dpi = 300)
ggsave(file.path(fig_d, "final_composition_by_sample.pdf"), p_cs, width = 10.4, height = 7.8)
lg("composition by sample done")

cnt_g <- read.csv(file.path(tab_d, "celltype_counts_by_group.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
prop_g <- read.csv(file.path(tab_d, "celltype_proportions_by_group.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(ct_levels %in% colnames(cnt_g)))
long_g <- do.call(rbind, lapply(ct_levels, function(k) {
  data.frame(group = cnt_g$group, celltype = k, n = cnt_g[[k]], pct = prop_g[[k]],
             stringsAsFactors = FALSE)
}))
long_g$celltype <- factor(long_g$celltype, levels = ct_levels)
long_g$group <- factor(long_g$group, levels = c("NF", "NICM", "DoxCM"))

p_cg <- ggplot(long_g, aes(group, pct, fill = celltype)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.1f", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.4, colour = "white") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(title = "Cell-type composition per group",
       subtitle = "pooled nuclei per group; descriptive only, no significance testing",
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
        plot.margin = margin(4, 4, 16, 4),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "final_counts_by_sample.png"), p_ns, width = 10.4, height = 7.8, dpi = 300)
ggsave(file.path(fig_d, "final_counts_by_sample.pdf"), p_ns, width = 10.4, height = 7.8)
lg("counts by sample done")

# ---------------------------------------------------------- 5. markers ------
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
  NEUR  = c("SOX10","PLP1","MPZ","PTPRZ1"),
  PROL  = c("MKI67","TOP2A","BUB1B"))

gs_map <- so[["RNA"]][[]][["gene_symbol"]]
ens <- rownames(so)[match(unlist(dot_genes), gs_map)]
names(ens) <- unlist(dot_genes)
ens <- ens[!is.na(ens)]
sym_of <- setNames(names(ens), unname(ens))   # ensembl -> symbol
lg("marker genes resolved: ", length(ens))

dp <- DotPlot(so, features = unname(ens), group.by = "celltype", assay = "RNA")
dd <- dp$data
dd$symbol <- factor(sym_of[as.character(dd$features.plot)], levels = rev(names(ens)))
dd$id <- factor(as.character(dd$id), levels = ct_levels)
p_dot <- ggplot(dd, aes(x = id, y = symbol)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027", midpoint = 0,
                         name = "scaled mean\nexpression") +
  scale_size_area(max_size = 4, name = "% expressing") +
  labs(title = "Canonical marker panel across final cell types",
       subtitle = "C8 resolved as endocardial endothelial; C16 removed",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.4),
        axis.text.y = element_text(size = 6.4, face = "italic"),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_d, "final_dotplot_markers.png"), p_dot, width = 7.8,
       height = 9.8, dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "final_dotplot_markers.pdf"), p_dot, width = 7.8, height = 9.8)
lg("dotplot done")

# marker heatmap (computed from the data layer; means per cell type)
dat <- GetAssayData(so, assay = "RNA", layer = "data")[unname(ens), , drop = FALSE]
avg <- sapply(ct_levels, function(k) Matrix::rowMeans(dat[, md$celltype == k, drop = FALSE]))
colnames(avg) <- ct_levels
rownames(avg) <- names(ens)
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
ggsave(file.path(fig_d, "final_marker_heatmap.png"), p_hm, width = 8.8, height = 7.8, dpi = 300)
ggsave(file.path(fig_d, "final_marker_heatmap.pdf"), p_hm, width = 8.8, height = 7.8)
lg("marker heatmap done")

writeLines(capture.output(sessionInfo()), file.path(adr, "sessionInfo_61.txt"))
lg("=== done ===")
