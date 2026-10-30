# =============================================================================
# 52_sc_annotation_figures.R
# Stage 3 - Step 2: annotation review tables + evidence figures.
#
# Inputs : results/Singlecel/annotation/objects/annotation_small.rds
#          results/Singlecel/annotation/tables/panel_{avg,pct,z}score.csv
# Outputs: results/annotation/  (review CSV is written separately)
#          tables/module_mean_pct.csv, module_mean_avg.csv, cluster_labels.csv
#          figures/umap_candidate_annotation.{png,pdf}
#                  marker_evidence_dotplot.{png,pdf}
#                  module_detection_heatmap.{png,pdf}
#                  cluster_overview.{png,pdf}
# Notes  : module z-scores can be inflated for modules that are absent almost
#          everywhere (a cluster with 2% detection gets a high z), so module
#          ranking here uses the ABSOLUTE mean detection rate.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(Matrix); library(scales)
})

proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
in_d  <- file.path(sdir, "annotation")
out_d <- file.path(proj, "results", "annotation")
tab_d <- file.path(out_d, "tables")
fig_d <- file.path(out_d, "figures")
for (d in c(out_d, tab_d, fig_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

SM <- readRDS(file.path(in_d, "objects", "annotation_small.rds"))
avg <- SM$avg; pct <- SM$pct; z <- SM$z
meta <- SM$meta
clu  <- sub("^C", "", colnames(avg))            # c("0",...)
ccol <- colnames(avg)

# ------------------------------------------------------- candidate labels ----
lab <- c("0"="Cardiomyocyte", "1"="Fibroblast", "2"="Endothelial cell",
         "3"="Pericyte", "4"="Fibroblast", "5"="Macrophage", "6"="Cardiomyocyte",
         "7"="Endothelial cell", "8"="Endothelial cell", "9"="Monocyte/DC",
         "10"="T/NK cell", "11"="Smooth muscle cell", "12"="Neural/glial cell",
         "13"="Lymphatic endothelial cell", "14"="Mast cell",
         "15"="Proliferating cell", "16"="Adipocyte")
sub_hint <- c("0"="CM (atrial-like, MYH6/CHRM2 high)", "6"="CM (ventricular, MYBPC3/RBM20)",
              "1"="FB (quiescent, PDGFRA/ABCA9)", "4"="FB (activated, MFAP5/SERPINE2)",
              "2"="EC (capillary/arterial, BTNL9/CA4)", "7"="EC (venous, ACKR1/ABCB1)",
              "8"="EC (endocardial-like, NPR3/NRG1/POSTN)", "3"="Pericyte (RGS5/HIGD1B)",
              "5"="Macrophage (resident, MRC1/CD163)", "9"="Monocyte/DC (LYZ/ITGAX/FCN1)",
              "10"="T/NK (CD2/IL7R/SKAP1)", "11"="SMC (MYH11/MYOCD/LMOD1)",
              "12"="Schwann/neural glia (SOX10/PLP1/MPZ)", "13"="LEC (CCL21/PROX1/MMRN1)",
              "14"="Mast (CPA3/TPSB2/KIT)", "15"="Cycling (myeloid-biased)",
              "16"="Adipocyte (PLIN1/ADIPOQ/GPAM)")

write.csv(data.frame(cluster = clu, column = ccol, label = unname(lab[clu]),
                     sub_hint = unname(sub_hint[clu]), check.names = FALSE),
          file.path(tab_d, "cluster_labels.csv"), row.names = FALSE)

# ------------------------------------------------ module detection rates -----
panel <- SM$panel
mod_pct <- sapply(ccol, function(cc) sapply(names(panel), function(m) {
  g <- intersect(panel[[m]], rownames(pct))
  if (!length(g)) return(NA_real_)
  mean(pct[g, cc])
}))
mod_avg <- sapply(ccol, function(cc) sapply(names(panel), function(m) {
  g <- intersect(panel[[m]], rownames(avg))
  if (!length(g)) return(NA_real_)
  mean(avg[g, cc])
}))
mod_pct <- as.data.frame(mod_pct); mod_avg <- as.data.frame(mod_avg)
write.csv(data.frame(module = rownames(mod_pct), mod_pct, check.names = FALSE),
          file.path(tab_d, "module_mean_pct.csv"), row.names = FALSE)
write.csv(data.frame(module = rownames(mod_avg), mod_avg, check.names = FALSE),
          file.path(tab_d, "module_mean_avg.csv"), row.names = FALSE)

best <- sapply(ccol, function(cc) {
  v <- mod_pct[[cc]]; names(v) <- rownames(mod_pct)
  v <- v[!names(v) %in% c("AMBI")]
  v <- sort(v, decreasing = TRUE)
  sprintf("%s(%.2f) / 2nd %s(%.2f)", names(v)[1], v[1], names(v)[2], v[2])
})
write.csv(data.frame(cluster = clu, auto_module_rank = unname(best), check.names = FALSE),
          file.path(tab_d, "cluster_auto_module.csv"), row.names = FALSE)

# --------------------------------------------------------------- UMAP -------
meta$label <- factor(lab[meta$cluster],
                     levels = c("Cardiomyocyte","Fibroblast","Endothelial cell",
                                "Pericyte","Smooth muscle cell","Macrophage",
                                "Monocyte/DC","T/NK cell","Lymphatic endothelial cell",
                                "Mast cell","Neural/glial cell","Adipocyte",
                                "Proliferating cell"))
pal <- c("Cardiomyocyte"="#D1495B","Fibroblast"="#3E7CB1","Endothelial cell"="#2A9D8F",
         "Pericyte"="#8E6C8A","Smooth muscle cell"="#6A994E","Macrophage"="#E9A03B",
         "Monocyte/DC"="#C97B4A","T/NK cell"="#4C5FD7","Lymphatic endothelial cell"="#7BC8B0",
         "Mast cell"="#B5527A","Neural/glial cell"="#8C8C3A","Adipocyte"="#00A3C4",
         "Proliferating cell"="#9B5DE5")
meta$label <- droplevels(meta$label)

cent <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = meta, FUN = median)

p_id <- ggplot(meta, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = cluster), size = 0.05, stroke = 0) +
  geom_text(data = cent, aes(label = cluster), size = 3.6, fontface = "bold") +
  labs(title = "Harmony UMAP - cluster (res 0.6)", x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

# draw large clusters first so that small clusters stay visible on top
meta_ord <- meta[order(-as.integer(table(meta$cluster)[meta$cluster])), ]

p_lab <- ggplot(meta_ord, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = label), size = 0.05, stroke = 0) +
  scale_colour_manual(values = pal, name = "candidate label") +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "Candidate annotation (major types only, no renaming applied)",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 10) + theme(legend.text = element_text(size = 7.5),
                                   legend.key.height = unit(0.38, "cm"))

ggsave(file.path(fig_d, "umap_candidate_annotation.png"), p_id | p_lab,
       width = 13.2, height = 6.2, dpi = 300)
ggsave(file.path(fig_d, "umap_candidate_annotation.pdf"), p_id | p_lab,
       width = 13.2, height = 6.2)

# ------------------------------------------------------------ dot plot ------
show_genes <- list(
  CM    = c("TTN","MYH6","MYH7","TNNT2","NPPA","MYBPC3","RBM20","CHRM2"),
  FB    = c("DCN","LUM","COL1A1","PDGFRA","MFAP5","SERPINE2","POSTN","FBN1"),
  EC    = c("PECAM1","CDH5","VWF","FLT1","EMCN","BTNL9","ACKR1"),
  ENDO  = c("NPR3","NRG1"),
  LEC   = c("CCL21","PROX1","MMRN1","FLT4","PKHD1L1"),
  PC    = c("RGS5","KCNJ8","ABCC9","PDGFRB","HIGD1B"),
  SMC   = c("MYH11","ACTA2","MYOCD","LMOD1","CARMN"),
  MYE   = c("PTPRC","LYZ","CD68","CD14","CSF1R","MRC1","CD163","F13A1","ITGAX","FCN1"),
  TNK   = c("CD3D","CD3E","CD2","IL7R","NKG7","SKAP1","THEMIS"),
  BPL   = c("CD79A","MS4A1","MZB1"),
  MAST  = c("CPA3","TPSB2","KIT","MS4A2","HDC"),
  NEUR  = c("SOX10","PLP1","MPZ","NRXN1","PTPRZ1"),
  ADIPO = c("PLIN1","ADIPOQ","GPAM","CIDEA"),
  PROL  = c("MKI67","TOP2A","BUB1B"),
  MESO  = c("UPK3B","MSLN"),
  AMBI  = c("MUC3A","KRT8","EPCAM")
)
g_order <- unlist(show_genes, use.names = FALSE)
g_order <- g_order[g_order %in% rownames(z)]
g_mod   <- rep(names(show_genes), lengths(show_genes))[match(g_order, unlist(show_genes))]

dd <- expand.grid(gene = g_order, cluster = clu, stringsAsFactors = FALSE)
dd$expr  <- as.vector(z[cbind(match(dd$gene, rownames(z)), match(paste0("C", dd$cluster), ccol))])
dd$pct   <- as.vector(pct[cbind(match(dd$gene, rownames(pct)), match(paste0("C", dd$cluster), ccol))])
dd$expr  <- pmax(pmin(dd$expr, 2.5), -1.5)
dd$module <- factor(g_mod[match(dd$gene, g_order)], levels = names(show_genes))
dd$gene   <- factor(dd$gene, levels = rev(g_order))
xl <- paste0("C", clu, "\n", lab[clu])

p_dot <- ggplot(dd, aes(x = factor(paste0("C", cluster), levels = paste0("C", clu)),
                        y = gene)) +
  geom_point(aes(size = pct, colour = expr)) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027", midpoint = 0,
                         name = "scaled mean\nexpression") +
  scale_size_area(max_size = 3.6, name = "fraction\nexpressing", limits = c(0, 1)) +
  scale_x_discrete(labels = xl) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Marker evidence per cluster (canonical panel)",
       subtitle = "colour = mean log-normalised expression scaled per gene; size = detection rate",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 6.6, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6.6, face = "italic"),
        strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right", legend.key.height = unit(0.5, "cm"))

ggsave(file.path(fig_d, "marker_evidence_dotplot.png"), p_dot,
       width = 9.2, height = 15.5, dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "marker_evidence_dotplot.pdf"), p_dot,
       width = 9.2, height = 15.5, limitsize = FALSE)

# --------------------------------------------- module detection heatmap -----
hm <- expand.grid(module = rownames(mod_pct), cluster = clu, stringsAsFactors = FALSE)
hm$value <- as.vector(as.matrix(mod_pct)[cbind(match(hm$module, rownames(mod_pct)),
                                               match(paste0("C", hm$cluster), ccol))])
hm$cluster_lab <- factor(paste0("C", hm$cluster, " ", lab[hm$cluster]),
                         levels = paste0("C", clu, " ", lab[clu]))
hm$module <- factor(hm$module, levels = rev(rownames(mod_pct)))

p_hm <- ggplot(hm, aes(cluster_lab, module, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.3,
            colour = ifelse(hm$value > 0.55, "white", "grey15")) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       limits = c(0, 1), name = "mean detection\nrate of module") +
  labs(title = "Lineage module detection rate per cluster",
       subtitle = "mean fraction of module genes detected; absolute values avoid z-score inflation of near-absent modules",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid = element_blank())

ggsave(file.path(fig_d, "module_detection_heatmap.png"), p_hm,
       width = 9.4, height = 5.6, dpi = 300)
ggsave(file.path(fig_d, "module_detection_heatmap.pdf"), p_hm,
       width = 9.4, height = 5.6)

# ------------------------------------------------------- cluster overview ---
ov <- data.frame(cluster = clu,
                 n = SM$qc$n_cells,
                 pct_all = SM$qc$n_cells / nrow(meta) * 100,
                 pct_disease = SM$qc$pct_disease,
                 label = unname(lab[clu]),
                 stringsAsFactors = FALSE)
ov$cluster_lab <- factor(paste0("C", ov$cluster), levels = paste0("C", ov$cluster))
p_ov <- ggplot(ov, aes(cluster_lab, pct_all, fill = label)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%d\n%.1f%%", n, pct_all)), size = 2.4, vjust = -0.25,
            lineheight = 0.85) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Cluster size and candidate label", x = NULL,
       y = "% of all nuclei (n = 151,007)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

p_dis <- ggplot(ov, aes(cluster_lab, pct_disease, fill = label)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 65.8, linetype = 2, colour = "grey30") +
  geom_text(aes(label = sprintf("%.1f", pct_disease)), size = 2.3, vjust = -0.3) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Disease fraction per cluster (dashed = global 65.8%)", x = NULL,
       y = "% Disease nuclei") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_d, "cluster_overview.png"), p_ov / p_dis, width = 9.0, height = 6.4, dpi = 300)
ggsave(file.path(fig_d, "cluster_overview.pdf"), p_ov / p_dis, width = 9.0, height = 6.4)

cat("figures written to", fig_d, "\n")
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_52.txt"))
