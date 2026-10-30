# =============================================================================
# 68_sc_figure6_composite.R
# Assemble the four-gene single-cell figure (Figure 6) as one composite:
#   A  DotPlot      4 genes x 12 cell types
#   B  FeaturePlot  4 UMAPs, colour = expression relative to each gene's 99th pct
#   C  donor-level  boxplots: top row = all nuclei, bottom row = carrier cell type
#   D  donor-level  cell-type composition, key compartments (added 2026-10-29)
#
# Input : fourgene_slim.rds (64), 65_*.csv tables (65)
# Output: figures/Figure6_fourgene_composite.{pdf,png}
#         tables/68_sc_celltype_composition_{donor,stats}.csv (panel D numbers)
#         Figure6_fourgene_legend.txt
#
# Layout notes
#   * cell types are ordered by lineage (parenchyma -> stroma/vascular -> immune)
#     instead of alphabetically, which is what a reviewer expects to read
#   * panel B uses a relative 0-1 scale so the four UMAPs share one legend
#   * panel C keeps one point = one donor; nuclei are never replicates
#   * panel D is computed here from the slim object meta so that the composition
#     numbers in the figure are traceable without re-running the Composition
#     module; it is cross-checked against Composition/01_tables/T16 on the fly
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(ggplot2); library(patchwork); library(scales); library(ggrastr)
  library(ggrepel)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "fourgene")
tab_d <- file.path(out_d, "tables"); fig_d <- file.path(out_d, "figures")

GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")
CT_ORDER <- c("Cardiomyocyte", "Fibroblast", "Endothelial cell",
              "Lymphatic endothelial cell", "Pericyte", "Smooth muscle cell",
              "Macrophage", "Monocyte/DC", "T/NK cell", "Mast cell",
              "Neural/glial cell", "Proliferating cell")
GRP <- c("NF", "NICM", "DoxCM")
GCOL <- c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B")

fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                            sprintf("%.3f", p)))

slim <- readRDS(file.path(out_d, "objects", "fourgene_slim.rds"))
cnt <- slim$counts; dat <- slim$data; md <- slim$meta; um <- slim$umap
md$group <- factor(as.character(md$group), levels = GRP)
md$celltype <- factor(as.character(md$celltype), levels = CT_ORDER)

dot <- read.csv(file.path(tab_d, "65_dotplot_celltype.csv"), stringsAsFactors = FALSE)
pbc <- read.csv(file.path(tab_d, "65_pseudobulk_sample_x_celltype.csv"), stringsAsFactors = FALSE)
pbd <- read.csv(file.path(tab_d, "65_pseudobulk_donor_tissue.csv"), stringsAsFactors = FALSE)
st_ct <- read.csv(file.path(tab_d, "65_stats_celltype.csv"), stringsAsFactors = FALSE)
st_dn <- read.csv(file.path(tab_d, "65_stats_donor_whole_tissue.csv"), stringsAsFactors = FALSE)

# =============================================================== panel A =====
dot$celltype <- factor(dot$celltype, levels = CT_ORDER)
dot$gene <- factor(dot$gene, levels = rev(GENES))
dot$scaled <- ave(dot$avg, dot$gene, FUN = function(x) {
  s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })

pA <- ggplot(dot, aes(celltype, gene)) +
  geom_point(aes(size = pct, colour = scaled)) +
  scale_colour_gradient2(low = "#3B6EA5", mid = "#F2F2F2", high = "#C62828",
                         midpoint = 0, name = "scaled mean\nexpression",
                         guide = guide_colourbar(barheight = unit(1.5, "cm"),
                                                 barwidth = unit(0.3, "cm"),
                                                 title.theme = element_text(size = 6.4))) +
  scale_size_area(max_size = 5.2, breaks = c(1, 5, 10, 20), name = "% nuclei\nexpressing",
                  guide = guide_legend(keyheight = unit(0.34, "cm"),
                                       title.theme = element_text(size = 6.4),
                                       override.aes = list(colour = "grey30"))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.2),
        axis.text.y = element_text(face = "italic", size = 6.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        legend.key = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

# =============================================================== panel B =====
cen <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype,
                 data = data.frame(um, celltype = md$celltype), FUN = median)
cen$pct <- as.numeric(table(md$celltype)[as.character(cen$celltype)]) / nrow(md) * 100
cen <- cen[cen$pct >= 1.5, ]
# labels only once (first facet); repeating them four times just adds clutter
cen_lab <- transform(cen, gene = factor("AQP3", levels = GENES))

rel <- do.call(rbind, lapply(GENES, function(g) {
  v <- as.numeric(dat[g, ])
  cap <- as.numeric(quantile(v[v > 0], 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(v)
  r <- pmin(v / cap, 1)
  o <- order(r)
  data.frame(gene = g, UMAP_1 = um[o, 1], UMAP_2 = um[o, 2], rel = r[o],
             stringsAsFactors = FALSE)
}))
rel$gene <- factor(rel$gene, levels = GENES)

pB <- ggplot(rel, aes(UMAP_1, UMAP_2)) +
  rasterise(geom_point(data = rel[rel$rel == 0, ], colour = "grey88", size = 0.05, stroke = 0), dpi = 500) +
  rasterise(geom_point(data = rel[rel$rel > 0, ], aes(colour = rel), size = 0.09, stroke = 0), dpi = 500) +
  geom_text_repel(data = cen_lab, aes(label = celltype), size = 2.1, colour = "grey10",
                  min.segment.length = 0, segment.size = 0.15, segment.colour = "grey55",
                  box.padding = 0.10, point.padding = 0.06, max.overlaps = 30,
                  seed = SEED, show.legend = FALSE) +
  facet_wrap(~ gene, nrow = 1) +
  coord_fixed(ratio = 1, expand = TRUE) +
  scale_colour_gradientn(colours = c("#FFF3B0", "#FDBB5E", "#EF6A3C", "#A50F15"),
                         limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                         name = "relative expression
(per-gene 99th pct = 1)",
                         guide = guide_colourbar(barheight = unit(1.4, "cm"),
                                                 barwidth = unit(0.28, "cm"),
                                                 title.theme = element_text(size = 6.4))) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "italic", size = 7.2),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 6.4),
        plot.margin = margin(2, 2, 2, 2))

# =============================================================== panel C =====
main_ct <- sapply(GENES, function(g) {
  d <- dot[dot$gene == g, ]
  ok <- st_ct$celltype[st_ct$gene == g & !is.na(st_ct$p_wilcox_NF_vs_Disease)]
  d <- d[as.character(d$celltype) %in% ok, ]
  as.character(d$celltype[which.max(d$avg)])
})

C_long <- do.call(rbind, lapply(GENES, function(g) {
  a <- pbd[pbd$gene == g, c("sample", "group", "log1p_cpm")]
  a$gene <- g; a$panel <- paste0(g, " | all nuclei"); a$row <- 1
  b <- pbc[pbc$gene == g & pbc$celltype == main_ct[[g]] & pbc$usable,
           c("sample", "group", "log1p_cpm")]
  b$gene <- g; b$panel <- paste0(g, " | ", main_ct[[g]]); b$row <- 2
  rbind(a, b)
}))
lv <- c(paste0(GENES, " | all nuclei"), paste0(GENES, " | ", main_ct[GENES]))
C_long$panel <- factor(C_long$panel, levels = lv)
C_long$group <- factor(C_long$group, levels = GRP)

ann <- rbind(
  data.frame(gene = GENES, panel = paste0(GENES, " | all nuclei"),
             fdr = st_dn$fdr[match(GENES, st_dn$gene)],
             lfc = st_dn$log2FC_median[match(GENES, st_dn$gene)], stringsAsFactors = FALSE),
  data.frame(gene = GENES, panel = paste0(GENES, " | ", main_ct[GENES]),
             fdr = st_ct$fdr[match(paste(GENES, main_ct[GENES]),
                                   paste(st_ct$gene, st_ct$celltype))],
             lfc = st_ct$log2FC_median[match(paste(GENES, main_ct[GENES]),
                                             paste(st_ct$gene, st_ct$celltype))],
             stringsAsFactors = FALSE))
ann$panel <- factor(ann$panel, levels = lv)
ann$lab <- ifelse(is.na(ann$fdr), "", sprintf("FDR %s", fmt_p(ann$fdr)))
ann$lab2 <- ifelse(is.na(ann$lfc), "", sprintf("log2FC %+.2f", ann$lfc))

pC <- ggplot(C_long, aes(group, log1p_cpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.55, linewidth = 0.25) +
  geom_point(aes(colour = group), size = 0.85, stroke = 0.12,
             position = position_jitter(width = 0.13, seed = SEED)) +
  geom_text(data = ann, aes(x = 2, y = Inf, label = lab), vjust = 1.35,
            size = 1.95, inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = 2, y = Inf, label = lab2), vjust = 3.1,
            size = 1.8, colour = "grey35", inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.34))) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL, name = NULL) +
  labs(x = NULL, y = "log1p(CPM), donor level") +
  theme_bw(base_size = 7) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 5.9),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.1),
        axis.text.y = element_text(size = 5.9),
        legend.position = "none",
        panel.spacing = unit(0.28, "lines"),
        plot.margin = margin(2, 2, 2, 2))

# =============================================================== panel D =====
# Donor-level cell-type composition. A nucleus is never a replicate: the unit is
# "share of a donor's nuclei", so 25 values per cell type, NF n = 8 vs 17.
KEY_CT <- c("Cardiomyocyte", "Fibroblast", "Endothelial cell",
            "Pericyte", "Macrophage", "T/NK cell")

ct_tab <- as.data.frame(table(sample = md$sample, celltype = md$celltype),
                        stringsAsFactors = FALSE)
names(ct_tab) <- c("sample", "celltype", "n_nuclei")
ct_tab$n_total <- ave(ct_tab$n_nuclei, ct_tab$sample, FUN = sum)
ct_tab$pct     <- ct_tab$n_nuclei / ct_tab$n_total * 100
grp_map <- unique(md[, c("sample", "group")])
grp_map$sample <- as.character(grp_map$sample)
ct_tab <- merge(ct_tab, grp_map, by = "sample")
ct_tab$celltype <- factor(ct_tab$celltype, levels = CT_ORDER)
ct_tab$group    <- factor(as.character(ct_tab$group), levels = GRP)
ct_tab$grp2     <- factor(ifelse(ct_tab$group == "NF", "NF", "Disease"),
                          levels = c("NF", "Disease"))
ct_tab <- ct_tab[order(ct_tab$group, ct_tab$sample, ct_tab$celltype), ]

comp_stats <- do.call(rbind, lapply(CT_ORDER, function(ct) {
  d  <- ct_tab[ct_tab$celltype == ct, ]
  nf <- d$pct[d$grp2 == "NF"]; dis <- d$pct[d$grp2 == "Disease"]
  data.frame(celltype = ct,
             n_NF = length(nf), n_Disease = length(dis),
             median_NF_pct = median(nf), median_disease_pct = median(dis),
             delta_pp = median(dis) - median(nf),
             p = tryCatch(wilcox.test(dis, nf, exact = FALSE)$p.value,
                          error = function(e) NA_real_),
             stringsAsFactors = FALSE)
}))
comp_stats$fdr <- p.adjust(comp_stats$p, "BH")
comp_stats <- comp_stats[order(comp_stats$p), ]

write.csv(comp_stats,
          file.path(tab_d, "68_sc_celltype_composition_stats.csv"), row.names = FALSE)
write.csv(ct_tab[, c("sample", "group", "celltype", "n_nuclei", "n_total", "pct")],
          file.path(tab_d, "68_sc_celltype_composition_donor.csv"), row.names = FALSE)

# cross-check against the Composition module (same object, different script)
t16 <- file.path(proj, "results", "Composition", "01_tables", "T16_sc_group_comparison.csv")
if (file.exists(t16)) {
  ref <- read.csv(t16, stringsAsFactors = FALSE)
  m <- merge(comp_stats[, c("celltype", "p")], ref[, c("celltype", "p")], by = "celltype")
  cat(sprintf("panel D vs Composition T16: %d cell types, max |dp| = %.4f, max |dFDR| = %.4f\n",
              nrow(m), max(abs(m$p.x - m$p.y)),
              max(abs(p.adjust(m$p.x, "BH") - p.adjust(m$p.y, "BH")))))
}

dsub <- ct_tab[as.character(ct_tab$celltype) %in% KEY_CT, ]
dsub$celltype <- factor(as.character(dsub$celltype), levels = KEY_CT)

annD <- comp_stats[match(KEY_CT, comp_stats$celltype), ]
annD$x <- NA
annD$celltype <- factor(KEY_CT, levels = KEY_CT)
annD$lab <- sprintf("FDR %s", fmt_p(annD$fdr))
annD$lab2 <- sprintf("%+.1f pp", annD$delta_pp)

pD <- ggplot(dsub, aes(group, pct, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.55, linewidth = 0.25) +
  geom_point(aes(colour = group), size = 0.85, stroke = 0.12,
             position = position_jitter(width = 0.13, seed = SEED)) +
  geom_text(data = annD, aes(x = 2, y = Inf, label = lab), vjust = 1.35,
            size = 1.95, inherit.aes = FALSE) +
  geom_text(data = annD, aes(x = 2, y = Inf, label = lab2), vjust = 3.05,
            size = 1.8, colour = "grey35", inherit.aes = FALSE) +
  facet_wrap(~ celltype, ncol = 6, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.36))) +
  scale_fill_manual(values = GCOL, name = NULL) +
  scale_colour_manual(values = GCOL, name = NULL) +
  labs(x = NULL, y = "% of donor nuclei") +
  theme_bw(base_size = 7) +
  theme(strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 5.9),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.1),
        axis.text.y = element_text(size = 5.9),
        legend.position = "none",
        panel.spacing = unit(0.28, "lines"),
        plot.margin = margin(2, 2, 2, 2))

# ============================================================== assemble =====
fig <- (pA / pB / pC / pD) +
  plot_layout(heights = c(1.05, 1.05, 1.30, 0.78)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        plot.margin = margin(3, 6, 3, 3))

ggsave(file.path(fig_d, "Figure6_fourgene_composite.pdf"), fig,
       width = 7.1, height = 9.5, device = cairo_pdf)
ggsave(file.path(fig_d, "Figure6_fourgene_composite.png"), fig,
       width = 7.1, height = 9.5, dpi = 600)

# ---------------------------------------------------------------- legend -----
cap <- c(
  "Figure 6. Single-cell localisation of the four model genes, donor-level expression and cell-type composition.",
  "",
  paste0("(A) Dot plot of AQP3, GGT5, ISLR and LAD1 across the 12 annotated cell types of the ",
         "snRNA-seq atlas (GSE292067; 150,750 nuclei from 25 donors). Dot size is the percentage ",
         "of nuclei with at least one count in that cell type; colour is the mean log-normalised ",
         "expression scaled across cell types within each gene."),
  "",
  paste0("(B) UMAP feature plots of the same four genes. Colour is expression relative to the 99th ",
         "percentile of expressing nuclei for each gene (0-1), so that the four panels share one ",
         "colour scale; grey denotes nuclei without detected transcripts. Cell-type labels are drawn ",
         "once (AQP3 panel) at the median coordinate of each cluster and apply to all four panels; ",
         "the underlying atlas is the one described in Figure 1."),
  "",
  paste0("(C) Donor-level expression. Top row: pseudobulk over all nuclei of a donor. Bottom row: ",
         "pseudobulk restricted to the cell type that carries the gene most strongly ",
         "(AQP3 endothelium, GGT5 pericytes, ISLR fibroblasts, LAD1 cardiomyocytes). ",
         "Each point is one donor (NF n = 8, NICM n = 6, DoxCM n = 11); values are log1p(CPM) ",
         "computed on the summed counts of the same nuclei. p values are two-sided Wilcoxon ",
         "rank-sum tests (NF vs NICM + DoxCM) with Benjamini-Hochberg correction within the ",
         "gene x cell-type family; nuclei were never treated as replicates."),
  "",
  paste0("(D) Donor-level cell-type composition for the six compartments that carry the four genes or ",
         "the myeloid signal. Each point is one donor and each value is the percentage of that donor's ",
         "nuclei assigned to the cell type (NF n = 8, NICM n = 6, DoxCM n = 11). p values are two-sided ",
         "Wilcoxon rank-sum tests (NF vs NICM + DoxCM) with Benjamini-Hochberg correction across the ",
         "12 cell types; pp = percentage-point difference between group medians. Composition is shown ",
         "separately from per-cell expression because a shift in the share of a cell type can mimic or ",
         "mask a shift in expression, and the two are never conflated in the text.")
)
writeLines(cap, file.path(out_d, "Figure6_fourgene_legend.txt"))

# quick numeric sanity dump
cat("main cell types:", paste(GENES, main_ct[GENES], sep = "=", collapse = " | "), "\n")
cat("panels in C:", paste(levels(C_long$panel), collapse = " // "), "\n")
cat("figure saved:", file.path(fig_d, "Figure6_fourgene_composite.pdf"), "\n")
