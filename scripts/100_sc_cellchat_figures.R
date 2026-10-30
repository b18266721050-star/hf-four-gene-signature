# =============================================================================
# 100_sc_cellchat_figures.R
# Stage 6 (A4) - Figure 7: communication overview + focal pathway (NRG).
#
# Input : cellchat/objects/{cc_input.rds, cc_networks.rds, cc_donor.rds}
#         cellchat/tables/99_*.csv
# Output: cellchat/figures/Figure7_cellchat_composite.{pdf,png}
#         cellchat/figures/100_panel{A..F}_*.{pdf,png}
#         cellchat/tables/100_focal_pair_detail.csv, 100_focal_edge_matrix.csv
#         cellchat/Figure7_legend.txt
#
# Panels
#   A  total communication strength and number of detected LR pairs per group
#   B  relative information flow of the 20 most active pathways (row %)
#   C  differential pathway volcano (disease vs NF, donor-level Wilcoxon)
#   D  focal pathway, sender -> receiver log2FC heatmap (FDR-masked)
#   E  focal pathway donor-level score (one point = one donor)
#   F  focal pathway circle plots, NF vs disease
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(ggplot2); library(patchwork); library(scales); library(ggrepel)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "cellchat")
tab_d <- file.path(out_d, "tables"); fig_d <- file.path(out_d, "figures")
obj_d <- file.path(out_d, "objects")
if (!dir.exists(fig_d)) dir.create(fig_d, recursive = TRUE)

GRP   <- c("NF", "NICM", "DoxCM")
GCOL  <- c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B")
KH <- 0.5; EXPR_FLOOR <- 0.05; MIN_CELLS <- 20
FOCAL <- "NRG"

cc   <- readRDS(file.path(obj_d, "cc_input.rds"))
net  <- readRDS(file.path(obj_d, "cc_networks.rds"))
don  <- readRDS(file.path(obj_d, "cc_donor.rds"))
CT   <- don$params$ct
gsum <- read.csv(file.path(tab_d, "98_group_summary.csv"), stringsAsFactors = FALSE)
flow <- read.csv(file.path(tab_d, "98_pathway_flow.csv"), stringsAsFactors = FALSE)
tests <- read.csv(file.path(tab_d, "99_pathway_tests.csv"), stringsAsFactors = FALSE)
edge <- read.csv(file.path(tab_d, "99_focal_edge_tests.csv"), stringsAsFactors = FALSE)
stopifnot(FOCAL %in% tests$pathway, tests$fdr[tests$pathway == FOCAL] < 0.05)

save2 <- function(p, name, w, h) {
  ggsave(file.path(fig_d, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(fig_d, paste0(name, ".png")), p, width = w, height = h, dpi = 600)
}
fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))

# =========================================================== panel A =========
gA <- gsum[gsum$group %in% c("NF", "NICM", "DoxCM"), ]
gA$group <- factor(gA$group, levels = GRP)
rel <- data.frame(group = factor(c(GRP, GRP), levels = GRP),
                  metric = rep(c("Detected\nLR pairs", "Total\nstrength"), each = 3),
                  value = c(gA$n_interactions_detected, gA$total_strength))
pA <- ggplot(rel, aes(group, value, fill = group)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(metric == "Total\nstrength",
                               sprintf("%.0f", value), sprintf("%d", as.integer(value)))),
            vjust = -0.35, size = 2.2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = NULL, title = "A  Communication output per group") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 6.6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(3, 3, 3, 3))

# =========================================================== panel B =========
fl <- flow[flow$pathway %in% tests$pathway[!tests$neuro], ]
fl$All <- NULL
top <- head(fl[order(-fl$NF), ], 20)
tot <- colSums(fl[, GRP])
relB <- do.call(rbind, lapply(GRP, function(g)
  data.frame(pathway = top$pathway, group = g, pct = top[[g]] / tot[[g]] * 100,
             stringsAsFactors = FALSE)))
relB$group <- factor(relB$group, levels = GRP)
ord <- top$pathway[order(top$NF / tot["NF"])]
relB$pathway <- factor(relB$pathway, levels = ord)
relB$val_z <- ave(relB$pct, relB$pathway, FUN = function(x) {
  s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })
pB <- ggplot(relB, aes(group, pathway, fill = val_z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", pct)), size = 1.9, colour = "grey15") +
  scale_fill_gradient2(low = "#3B6EA5", mid = "#F7F7F7", high = "#C62828",
                       midpoint = 0, name = "row z") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "B  Relative information flow (% of non-neurotransmitter total)") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6.2),
        axis.text.x = element_text(size = 6.6),
        legend.key.width = unit(0.24, "cm"), legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(3, 3, 3, 3))

# =========================================================== panel C =========
tv <- tests[!tests$neuro & !is.na(tests$p) &
            tests$median_NF > 0.05 & tests$median_disease > 0.05, ]
tv$neglog10 <- -log10(pmax(tv$p, 1e-6))
tv$log2_ratio_raw <- tv$log2_ratio
tv$log2_ratio <- pmax(pmin(tv$log2_ratio, 4), -4)   # on/off pathways would stretch the axis
tv$capped <- abs(tv$log2_ratio_raw) > 4
tv$cls <- ifelse(tv$fdr < 0.05 & tv$log2_ratio > 0, "up in disease",
          ifelse(tv$fdr < 0.05 & tv$log2_ratio < 0, "down in disease", "n.s."))
tv$cls <- factor(tv$cls, levels = c("up in disease", "down in disease", "n.s."))
pick <- function(d, n) { d <- d[order(-d$neglog10), ]; d[seq_len(min(n, nrow(d))), ] }
lab <- rbind(pick(tv[tv$cls == "up in disease", ], 6),
             pick(tv[tv$cls == "down in disease", ], 6),
             tv[tv$pathway == FOCAL, ])
lab <- lab[!duplicated(lab$pathway), ]
pC <- ggplot(tv, aes(log2_ratio, neglog10)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.22, colour = "grey55") +
  geom_point(aes(colour = cls, size = median_NF + median_disease), alpha = 0.85, stroke = 0) +
  geom_point(data = tv[tv$pathway == FOCAL, ], shape = 21, size = 3.4, colour = "black", fill = NA, stroke = 0.5) +
  geom_text_repel(data = lab, aes(label = pathway), size = 2.1, max.overlaps = 40,
                  min.segment.length = 0, segment.size = 0.15, seed = SEED,
                  box.padding = 0.35, point.padding = 0.2) +
  scale_colour_manual(values = c("up in disease" = "#C62828", "down in disease" = "#2C7FB8",
                                 "n.s." = "grey72"), name = NULL) +
  scale_size_area(max_size = 3.2, name = "baseline\nflow") +
  scale_x_continuous(breaks = seq(-4, 4, 2)) +
  labs(x = expression(log[2]~"(disease / NF), capped at " %+-% 4),
       y = expression(-log[10]~italic(p)),
       title = "C  Differential pathway activity (donor level)") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 5.6),
        legend.title = element_text(size = 5.8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(3, 3, 3, 3))

# =========================================================== panel D =========
ed <- edge[edge$pathway == FOCAL, ]
m <- matrix(NA_real_, length(CT), length(CT), dimnames = list(CT, CT))
mf <- matrix(1, length(CT), length(CT), dimnames = list(CT, CT))
for (i in seq_len(nrow(ed))) {
  undetected <- ed$median_NF[i] <= 1e-9 & ed$median_disease[i] <= 1e-9
  m[ed$sender[i], ed$receiver[i]] <- if (undetected) NA_real_ else ed$log2_ratio[i]
  mf[ed$sender[i], ed$receiver[i]] <- ed$fdr[i]
}
m[!is.na(m) & abs(m) > 3] <- sign(m[!is.na(m) & abs(m) > 3]) * 3
dD <- expand.grid(sender = CT, receiver = CT, stringsAsFactors = FALSE)
dD$log2 <- as.vector(m); dD$fdr <- as.vector(mf)
dD$sender <- factor(dD$sender, levels = rev(CT)); dD$receiver <- factor(dD$receiver, levels = CT)
dD$star <- ifelse(is.na(dD$log2), "", ifelse(dD$fdr < 0.05, "*", ifelse(dD$fdr < 0.1, ".", "")))
write.csv(data.frame(sender = rep(CT, times = length(CT)), receiver = rep(CT, each = length(CT)),
                     log2_ratio = as.vector(m), fdr = as.vector(mf), pathway = FOCAL),
          file.path(tab_d, "100_focal_edge_matrix.csv"), row.names = FALSE)
pD <- ggplot(dD, aes(receiver, sender, fill = log2)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = star), size = 3.2, vjust = 0.78, colour = "grey10") +
  scale_fill_gradient2(low = "#3B6EA5", mid = "#F7F7F7", high = "#C62828",
                       midpoint = 0, limits = c(-3, 3), oob = squish, na.value = "grey88",
                       name = expression(log[2]~"FC")) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "receiver", y = "sender",
       title = paste0("D  ", FOCAL, " pathway: sender -> receiver")) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, size = 6),
        axis.text.y = element_text(size = 6),
        legend.key.width = unit(0.24, "cm"), legend.key.height = unit(0.55, "cm"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(3, 3, 3, 3))

# =========================================================== panel E =========
sc <- don$score
dg <- don$donor_grp
dv <- data.frame(sample = rownames(sc), score = sc[, FOCAL],
                 group = factor(as.character(dg$group), levels = GRP))
tt <- tests[tests$pathway == FOCAL, ]
pE <- ggplot(dv, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35, linewidth = 0.28) +
  geom_point(shape = 21, size = 1.5, stroke = 0.25, colour = "grey20",
             position = position_jitter(width = 0.07, seed = SEED)) +
  annotate("text", x = 1.5, y = max(dv$score) * 1.02,
           label = sprintf("NF vs disease: log2FC = %.2f, p = %s, FDR = %s",
                           tt$log2_ratio, fmt_p(tt$p), fmt_p(tt$fdr)), size = 2.1) +
  scale_fill_manual(values = GCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = paste0(FOCAL, " pathway score"), title = "E  Donor-level score") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(3, 3, 3, 3))

# =========================================================== panel F =========
# pair-level donor probabilities for the focal pathway (recomputed here: cheap)
cnt <- cc$counts; meta <- cc$meta; pairs <- cc$pairs
lib <- meta$nCount_RNA / 1e4
dat <- cnt %*% Matrix::Diagonal(x = 1 / lib); dat@x <- log1p(dat@x)
samples <- rownames(sc)
strata <- interaction(meta$sample, meta$celltype, drop = TRUE, lex.order = TRUE)
S2 <- Matrix::sparseMatrix(i = seq_len(nrow(meta)), j = as.integer(strata), x = 1,
                           dims = c(nrow(meta), nlevels(strata)))
n_str <- as.numeric(Matrix::colSums(S2))
mean_expr <- sweep(as.matrix(dat %*% S2), 2, pmax(n_str, 1), "/")
lv <- levels(strata); sp <- factor(sub("\\..*$", "", lv), levels = samples)
ctt <- factor(sub("^.*\\.", "", lv), levels = CT); usable <- n_str >= MIN_CELLS
E <- array(0, dim = c(nrow(mean_expr), length(CT), length(samples)),
           dimnames = list(rownames(mean_expr), CT, samples))
for (k in which(usable)) E[, as.integer(ctt[k]), as.integer(sp[k])] <- mean_expr[, k]

sel <- which(pairs$pathway == FOCAL)
net_edge <- array(0, dim = c(length(samples), length(CT), length(CT)),
                  dimnames = list(samples, CT, CT))
pair_rows <- list()
for (i in sel) {
  a <- pairs$ligand[[i]]; a <- a[a %in% dimnames(E)[[1]]]
  b <- pairs$receptor[[i]]; b <- b[b %in% dimnames(E)[[1]]]
  if (!length(a) || !length(b)) next
  Lm <- if (length(a) == 1) E[a, , ] else apply(E[a, , , drop = FALSE], c(2, 3), min)
  Rm <- if (length(b) == 1) E[b, , ] else apply(E[b, , , drop = FALSE], c(2, 3), min)
  for (d in seq_along(samples)) {
    x <- outer(Lm[, d], Rm[, d]); p <- x / (x + KH); p[x < EXPR_FLOOR] <- 0
    net_edge[d, , ] <- net_edge[d, , ] + p
  }
  prob_pair <- vapply(seq_along(samples), function(d) {
    x <- outer(Lm[, d], Rm[, d]); sum(x / (x + KH))
  }, numeric(1))
  pair_rows[[length(pair_rows) + 1]] <- data.frame(
    interaction = pairs$interaction_name[i], pathway = FOCAL,
    ligand = pairs$ligand_str[i], receptor = pairs$receptor_str[i],
    NF_median = median(prob_pair[dg$group == "NF"]),
    disease_median = median(prob_pair[dg$group != "NF"]),
    p = suppressWarnings(wilcox.test(prob_pair[dg$group != "NF"],
                                     prob_pair[dg$group == "NF"], exact = FALSE)$p.value),
    stringsAsFactors = FALSE)
}
pd <- do.call(rbind, pair_rows)
pd$log2_ratio <- log2((pd$disease_median + 1e-9) / (pd$NF_median + 1e-9))
pd$fdr <- p.adjust(pd$p, method = "BH")
pd <- pd[order(-pd$disease_median), ]
write.csv(pd, file.path(tab_d, "100_focal_pair_detail.csv"), row.names = FALSE)

circ_dat <- function(grp_lab, idx) {
  m <- apply(net_edge[idx, , , drop = FALSE], c(2, 3), sum)
  data.frame(sender = rep(rownames(m), times = ncol(m)),
             receiver = rep(colnames(m), each = nrow(m)),
             w = as.vector(m), group = grp_lab, stringsAsFactors = FALSE)
}
ci <- rbind(circ_dat("NF", which(dg$group == "NF")),
            circ_dat("Disease", which(dg$group != "NF")))
ci$group <- factor(ci$group, levels = c("NF", "Disease"))
wmax <- quantile(ci$w[ci$w > 0], 0.9)
ci <- ci[ci$w > wmax * 0.02, ]

node <- data.frame(celltype = CT, ang = seq(90, -270 + 360 / length(CT), length.out = length(CT)) * pi / 180)
node$x <- cos(node$ang); node$y <- sin(node$ang)
SHORT <- c(Cardiomyocyte = "CM", Fibroblast = "Fibroblast", `Endothelial cell` = "Endo-\nthelial",
           `Lymphatic endothelial cell` = "Lymphatic\nEC", Pericyte = "Pericyte",
           `Smooth muscle cell` = "Smooth\nmuscle", Macrophage = "Macro-\nphage",
           `Monocyte/DC` = "Mono/\nDC", `T/NK cell` = "T/NK", `Mast cell` = "Mast",
           `Neural/glial cell` = "Neural/\nglial", `Proliferating cell` = "Prolif.")
node$lab <- SHORT[as.character(node$celltype)]
ci$x <- node$x[match(ci$sender, node$celltype)]
ci$y <- node$y[match(ci$sender, node$celltype)]
ci$xend <- node$x[match(ci$receiver, node$celltype)]
ci$yend <- node$y[match(ci$receiver, node$celltype)]
self <- ci[ci$sender == ci$receiver, ]
cross <- ci[ci$sender != ci$receiver, ]

pF <- ggplot() +
  geom_curve(data = cross, aes(x = x, y = y, xend = xend, yend = yend, linewidth = w),
             curvature = 0.28, colour = "#B03A2E", alpha = 0.45,
             arrow = arrow(length = unit(0.05, "cm"), type = "closed")) +
  geom_point(data = self, aes(x = x * 1.12, y = y * 1.12, size = w),
             shape = 21, fill = "#B03A2E", colour = "grey20", stroke = 0.2, alpha = 0.6) +
  geom_point(data = node, aes(x * 1.0, y * 1.0), size = 2.4, shape = 21,
             fill = "grey96", colour = "grey25", stroke = 0.4) +
  geom_text(data = node, aes(x * 1.30, y * 1.30, label = lab), size = 1.75, lineheight = 0.85) +
  facet_wrap(~ group, nrow = 1) +
  coord_equal(clip = "off", xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55)) +
  scale_linewidth_continuous(range = c(0.1, 1.6), name = "strength") +
  scale_size_area(max_size = 4, name = "self-loop") +
  labs(title = paste0("F  ", FOCAL, " network, NF vs disease (edge width = summed probability)")) +
  theme_void(base_size = 7) +
  theme(strip.text = element_text(size = 7.4, face = "bold"),
        panel.spacing = unit(0.9, "lines"),
        legend.position = "right", legend.key.size = unit(0.32, "cm"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(4, 4, 4, 4))

# =========================================================== composite =======
comp <- (pA | pC) / (pB | pE) / (pD | pF) +
  plot_layout(heights = c(1, 1.15, 1.05))
save2(comp, "Figure7_cellchat_composite", 9.2, 9.6)
save2(pA, "100_panelA_group_totals", 4.2, 2.4)
save2(pB, "100_panelB_pathway_heatmap", 4.4, 3.6)
save2(pC, "100_panelC_pathway_volcano", 4.6, 3.4)
save2(pD, "100_panelD_focal_edge_heatmap", 4.4, 4.0)
save2(pE, "100_panelE_focal_donor_score", 3.6, 3.0)
save2(pF, "100_panelF_focal_circle", 7.6, 3.8)

legend <- c(
  "Figure 7. Cell-cell communication landscape of the failing human left ventricle (GSE292067, 25 donors).",
  "A  Total communication strength (sum of CellChat mass-action probabilities over all ligand-receptor pairs and all sender-receiver combinations) and the number of detected ligand-receptor pairs, per group. Group pools are descriptive only; no test is performed on pooled nuclei.",
  "B  Relative information flow of the 20 most active non-neurotransmitter pathways, expressed as percentage of each group's total non-neurotransmitter information flow and coloured by row z-score. Numbers are the percentages.",
  "C  Differential pathway activity between disease (NICM + DoxCM, n = 17) and non-failing donors (NF, n = 8). Each point is a pathway; x = log2 ratio of median donor-level pathway score, y = -log10 Wilcoxon p. Filled red/blue = FDR < 0.05; grey = not significant. The focal pathway NRG is circled.",
  "D  NRG pathway, sender -> receiver log2 fold change (disease / NF) of the summed probability. Asterisks mark donor-level FDR < 0.05, dots FDR < 0.10 (Wilcoxon, 8 vs 17 donors, BH within pathway edges). Values are capped at +/-6 for display.",
  "E  Donor-level NRG pathway score; one point = one donor. Box hinges = IQR, line = median. p and FDR from the donor-level Wilcoxon test.",
  "F  NRG network in NF and disease donors. Nodes are the 12 annotated cell types arranged on a circle; edge width is proportional to the summed communication probability, arrows point from sender to receiver, and filled circles next to a node are self-loops.",
  "Statistics: communication probabilities were computed per donor from donor x cell type pseudobulk profiles, so the unit of replication is the donor (NF 8 vs disease 17), not the nucleus. Pathways were tested with two-sided Wilcoxon rank-sum tests and BH correction across non-neurotransmitter pathways; effect size is Cliff's delta."
)
writeLines(legend, file.path(out_d, "Figure7_legend.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_100.txt"))
cat("done\n")
