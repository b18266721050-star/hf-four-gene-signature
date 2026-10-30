# =============================================================================
# 65_sc_fourgene_plots_donor.R
# Stage 5 (A3) - four ML genes in the snRNA-seq atlas + donor-level statistics.
#
# Input : results/Singlecel/fourgene/objects/fourgene_slim.rds  (from 64)
#         results/DEG/GSE57338_DEGall.csv                        (bulk log2FC)
# Output: results/Singlecel/fourgene/figures/65_*.{png,pdf}
#         results/Singlecel/fourgene/tables/65_*.csv / .xlsx
#         results/Singlecel/fourgene/65_plots_donor.log
#
# Statistical unit = DONOR (sample), never the nucleus.
#   * pseudobulk = sum of counts per sample (x cell type), CPM on the summed
#     nCount_RNA of the same nuclei -> removes the cell-number artefact
#   * primary contrast NF vs Disease (NICM + DoxCM), Wilcoxon rank-sum,
#     BH-FDR within the gene x cell-type family, Cliff's delta as effect size
#   * secondary: 3-group Kruskal-Wallis + pairwise vs NF, and a linear mixed
#     model with a donor random intercept as a repeated-measures sensitivity
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(ggplot2); library(patchwork); library(scales)
  library(lme4); library(car); library(openxlsx); library(ggrastr)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "fourgene")
fig_d <- file.path(out_d, "figures")
tab_d <- file.path(out_d, "tables")
obj_d <- file.path(out_d, "objects")
for (d in c(fig_d, tab_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file.path(out_d, "65_plots_donor.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

GRP <- c("NF", "NICM", "DoxCM")
GCOL <- c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B")
MIN_CELLS <- 50      # nuclei needed for a sample x cell-type pseudobulk
MIN_DONOR <- 4       # donors needed per group for a test

slim <- readRDS(file.path(obj_d, "fourgene_slim.rds"))
cnt  <- slim$counts; dat <- slim$data; md <- slim$meta; um <- slim$umap
GENES <- slim$genes
lg("slim loaded: ", paste(GENES, collapse = ","), " | ", ncol(cnt), " nuclei")

md$group    <- factor(as.character(md$group), levels = GRP)
md$celltype <- factor(as.character(md$celltype))
md$sample   <- factor(as.character(md$sample))
CT <- levels(md$celltype)
lg("cell types (", length(CT), "): ", paste(CT, collapse = " | "))
lg("donors per group: ", paste(sprintf("%s=%d", GRP,
    sapply(GRP, function(g) length(unique(md$sample[md$group == g])))), collapse = " "))

# bulk reference (direction check)
bulk <- read.csv(file.path(proj, "results", "DEG", "GSE57338_DEGall.csv"),
                 stringsAsFactors = FALSE)
bulk4 <- bulk[match(GENES, bulk$GeneSymbol), c("GeneSymbol", "log2FoldChange", "padj")]
names(bulk4) <- c("gene", "bulk_log2FC", "bulk_padj")
lg("bulk log2FC: ", paste(sprintf("%s %+.2f", bulk4$gene, bulk4$bulk_log2FC), collapse = " | "))

# =============================================================== helpers ====
agg_mean <- function(mat, f) {
  D <- sparse.model.matrix(~0 + f); colnames(D) <- levels(f)
  n <- setNames(as.numeric(colSums(D)), colnames(D))
  pos <- (mat > 0) * 1
  list(avg = sweep(as.matrix(mat %*% D), 2, n, "/"),
       pct = sweep(as.matrix(pos %*% D), 2, n, "/") * 100,
       n = n, D = D)
}
cliff <- function(x, y) {                     # x vs y
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  n1 <- length(x); n2 <- length(y)
  U <- sum(outer(x, y, ">")) + 0.5 * sum(outer(x, y, "=="))
  2 * U / (n1 * n2) - 1
}
pfmt <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
}
plab <- function(fdr) ifelse(is.na(fdr), "", ifelse(fdr < 0.001, "***",
                    ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "ns"))))

# ---------------------------------------------------- dot-plot matrices -----
ct_all <- agg_mean(dat, md$celltype)
dot <- expand.grid(gene = GENES, celltype = CT, stringsAsFactors = FALSE)
gi <- match(dot$gene, rownames(ct_all$avg))
ci <- match(dot$celltype, colnames(ct_all$avg))
dot$avg <- ct_all$avg[cbind(gi, ci)]
dot$pct <- ct_all$pct[cbind(gi, ci)]
dot$n_nuclei <- as.numeric(ct_all$n[ci])
dot$scaled <- ave(dot$avg, dot$gene, FUN = function(x) {
  s <- sd(x); if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })
dot$gene <- factor(dot$gene, levels = rev(GENES))
dot$celltype <- factor(dot$celltype, levels = CT)
write.csv(dot, file.path(tab_d, "65_dotplot_celltype.csv"), row.names = FALSE)

# same, split by group (raw mean expression on a shared colour scale)
gc_key <- interaction(md$group, md$celltype, drop = TRUE, sep = "|")
gc <- agg_mean(dat, gc_key)
kv <- data.frame(key = colnames(gc$avg), stringsAsFactors = FALSE)
kv$group <- sub("\\|.*$", "", kv$key)
kv$celltype <- sub("^.*\\|", "", kv$key)
dotg <- do.call(rbind, lapply(GENES, function(g) data.frame(
  gene = g, group = kv$group, celltype = kv$celltype,
  avg = as.vector(gc$avg[g, kv$key]), pct = as.vector(gc$pct[g, kv$key]),
  n = gc$n[kv$key], stringsAsFactors = FALSE)))
dotg$group <- factor(dotg$group, levels = GRP)
dotg$celltype <- factor(dotg$celltype, levels = CT)
dotg$gene <- factor(dotg$gene, levels = rev(GENES))
write.csv(dotg, file.path(tab_d, "65_dotplot_celltype_by_group.csv"), row.names = FALSE)

# ---------------------------------------------------------- pseudobulk ------
pb_one <- function(key) {
  f <- factor(key)
  D <- sparse.model.matrix(~0 + f); colnames(D) <- levels(f)
  cnts <- as.matrix(cnt %*% D)
  tot  <- setNames(as.numeric(crossprod(D, md$nCount_RNA)), colnames(D))
  n    <- setNames(as.numeric(colSums(D)), colnames(D))
  det  <- as.matrix(((cnt > 0) * 1) %*% D)
  list(counts = cnts, total = tot, n = n, det = det, keys = levels(f), D = D)
}

# (a) sample x cell type
k1 <- pb_one(paste(md$sample, md$celltype, sep = "|"))
p1 <- data.frame(key = k1$keys, stringsAsFactors = FALSE)
p1$sample <- sub("\\|.*$", "", p1$key)
p1$celltype <- sub("^.*\\|", "", p1$key)
p1$n_nuclei <- as.integer(k1$n[p1$key])
p1$total_counts <- as.numeric(k1$total[p1$key])
p1$group <- as.character(md$group[match(p1$sample, as.character(md$sample))])
pbc <- do.call(rbind, lapply(GENES, function(g) data.frame(
  gene = g, p1[, c("sample", "group", "celltype", "n_nuclei", "total_counts")],
  counts = as.numeric(k1$counts[g, p1$key]),
  pct_detect = as.numeric(k1$det[g, p1$key]) / p1$n_nuclei * 100,
  stringsAsFactors = FALSE)))
pbc$cpm <- pbc$counts / pbc$total_counts * 1e6
pbc$log1p_cpm <- log1p(pbc$cpm)
pbc$usable <- pbc$n_nuclei >= MIN_CELLS
pbc$group <- factor(pbc$group, levels = GRP)
write.csv(pbc, file.path(tab_d, "65_pseudobulk_sample_x_celltype.csv"), row.names = FALSE)

# (b) whole tissue, one value per donor
k2 <- pb_one(md$sample)
p2 <- data.frame(sample = k2$keys, stringsAsFactors = FALSE)
p2$n_nuclei <- as.integer(k2$n[p2$sample])
p2$total_counts <- as.numeric(k2$total[p2$sample])
p2$group <- as.character(md$group[match(p2$sample, as.character(md$sample))])
pbd <- do.call(rbind, lapply(GENES, function(g) data.frame(
  gene = g, p2[, c("sample", "group", "n_nuclei", "total_counts")],
  counts = as.numeric(k2$counts[g, p2$sample]),
  pct_detect = as.numeric(k2$det[g, p2$sample]) / p2$n_nuclei * 100,
  stringsAsFactors = FALSE)))
pbd$cpm <- pbd$counts / pbd$total_counts * 1e6
pbd$log1p_cpm <- log1p(pbd$cpm)
pbd$group <- factor(pbd$group, levels = GRP)
write.csv(pbd, file.path(tab_d, "65_pseudobulk_donor_tissue.csv"), row.names = FALSE)
lg("pseudobulk built: ", nrow(pbc), " gene x sample x celltype rows, ",
   sum(pbc$usable), " usable (>= ", MIN_CELLS, " nuclei)")

# ============================================================= statistics ===
test_block <- function(df, label) {
  # df: gene, sample, group, value, extra columns already filtered
  do.call(rbind, lapply(split(df, df$gene), function(d) {
    nf  <- d$value[d$group == "NF"]
    dis <- d$value[d$group != "NF"]
    ni  <- d$value[d$group == "NICM"]
    dx  <- d$value[d$group == "DoxCM"]
    ok  <- length(nf) >= MIN_DONOR && length(dis) >= MIN_DONOR
    pw  <- if (ok) tryCatch(wilcox.test(dis, nf, exact = FALSE)$p.value,
                            error = function(e) NA_real_) else NA_real_
    kw  <- tryCatch(kruskal.test(value ~ group, data = d)$p.value,
                    error = function(e) NA_real_)
    pni <- if (length(ni) >= MIN_DONOR && length(nf) >= MIN_DONOR)
      tryCatch(wilcox.test(ni, nf, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
    pdx <- if (length(dx) >= MIN_DONOR && length(nf) >= MIN_DONOR)
      tryCatch(wilcox.test(dx, nf, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
    data.frame(
      gene = d$gene[1], unit = label,
      n_NF = length(nf), n_NICM = length(ni), n_DoxCM = length(dx),
      median_NF = median(nf), median_Disease = median(dis),
      log2FC_median = log2((median(dis) + 0.05) / (median(nf) + 0.05)),
      p_wilcox_NF_vs_Disease = pw,
      cliffs_delta = cliff(dis, nf),
      p_kruskal_3group = kw,
      p_NICM_vs_NF = pni, p_DoxCM_vs_NF = pdx,
      stringsAsFactors = FALSE)
  }))
}

tab_donor <- test_block(
  data.frame(gene = pbd$gene, sample = pbd$sample, group = pbd$group,
             value = pbd$log1p_cpm, stringsAsFactors = FALSE), "donor_whole_tissue")
tab_donor$fdr <- p.adjust(tab_donor$p_wilcox_NF_vs_Disease, "BH")
tab_donor$direction <- ifelse(tab_donor$cliffs_delta > 0, "up_in_disease",
                       ifelse(tab_donor$cliffs_delta < 0, "down_in_disease", "none"))
tab_donor <- merge(tab_donor, bulk4, by = "gene", all.x = TRUE)
tab_donor$bulk_consistent <- with(tab_donor, sign(bulk_log2FC) == sign(cliffs_delta))
write.csv(tab_donor, file.path(tab_d, "65_stats_donor_whole_tissue.csv"), row.names = FALSE)
lg("--- donor-level (whole tissue) NF vs Disease ---")
for (i in seq_len(nrow(tab_donor)))
  lg(sprintf("  %s: log2FC %+.2f (bulk %+.2f) p=%.3g FDR=%.3g delta=%+.2f",
             tab_donor$gene[i], tab_donor$log2FC_median[i], tab_donor$bulk_log2FC[i],
             tab_donor$p_wilcox_NF_vs_Disease[i], tab_donor$fdr[i], tab_donor$cliffs_delta[i]))

# per cell type
pbc_u <- pbc[pbc$usable, ]
tab_ct <- do.call(rbind, lapply(split(pbc_u, list(pbc_u$gene, pbc_u$celltype), drop = TRUE),
  function(d) {
    nf <- d$log1p_cpm[d$group == "NF"]; dis <- d$log1p_cpm[d$group != "NF"]
    ni <- d$log1p_cpm[d$group == "NICM"]; dx <- d$log1p_cpm[d$group == "DoxCM"]
    ok <- length(nf) >= MIN_DONOR && length(dis) >= MIN_DONOR
    data.frame(
      gene = d$gene[1], celltype = d$celltype[1], unit = "sample_x_celltype",
      n_nuclei_total = sum(d$n_nuclei),
      n_NF = length(nf), n_NICM = length(ni), n_DoxCM = length(dx),
      pct_detect_NF = mean(d$pct_detect[d$group == "NF"]),
      pct_detect_Disease = mean(d$pct_detect[d$group != "NF"]),
      median_NF = median(nf), median_Disease = median(dis),
      log2FC_median = log2((median(dis) + 0.05) / (median(nf) + 0.05)),
      p_wilcox_NF_vs_Disease = if (ok)
        tryCatch(wilcox.test(dis, nf, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_,
      cliffs_delta = cliff(dis, nf),
      p_kruskal_3group = tryCatch(kruskal.test(log1p_cpm ~ group, data = d)$p.value,
                                  error = function(e) NA_real_),
      p_NICM_vs_NF = if (length(ni) >= MIN_DONOR && length(nf) >= MIN_DONOR)
        tryCatch(wilcox.test(ni, nf, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_,
      p_DoxCM_vs_NF = if (length(dx) >= MIN_DONOR && length(nf) >= MIN_DONOR)
        tryCatch(wilcox.test(dx, nf, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_,
      stringsAsFactors = FALSE)
  }))
tab_ct$fdr <- p.adjust(tab_ct$p_wilcox_NF_vs_Disease, "BH")
tab_ct$direction <- ifelse(tab_ct$cliffs_delta > 0, "up_in_disease",
                    ifelse(tab_ct$cliffs_delta < 0, "down_in_disease", "none"))
lg("tab_ct cols: ", paste(colnames(tab_ct), collapse = ","), " rows=", nrow(tab_ct))
tab_ct <- merge(tab_ct, bulk4, by = "gene", all.x = TRUE)
tab_ct$bulk_consistent <- sign(tab_ct$bulk_log2FC) == sign(tab_ct$cliffs_delta)
tab_ct <- tab_ct[order(tab_ct$gene, -abs(tab_ct$cliffs_delta)), ]
write.csv(tab_ct, file.path(tab_d, "65_stats_celltype.csv"), row.names = FALSE)
lg("celltype tests: ", sum(!is.na(tab_ct$p_wilcox_NF_vs_Disease)), " testable rows; ",
   sum(tab_ct$fdr < 0.05, na.rm = TRUE), " at FDR<0.05")

# mixed model sensitivity: value ~ group + celltype + (1|donor)
mm <- do.call(rbind, lapply(GENES, function(g) {
  d <- pbc_u[pbc_u$gene == g, ]
  if (length(unique(d$group)) < 2 || nrow(d) < 12) return(NULL)
  fit <- tryCatch(lmer(log1p_cpm ~ group + celltype + (1 | sample), data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  a <- tryCatch(car::Anova(fit, type = 2, test.statistic = "Chisq"),
                error = function(e) NULL)
  if (is.null(a)) return(NULL)
  data.frame(gene = g, term = rownames(a), chisq = a$Chisq, df = a$Df,
             p = a$`Pr(>Chisq)`, n_obs = nrow(d),
             n_donor = length(unique(d$sample)), stringsAsFactors = FALSE)
}))
if (!is.null(mm)) {
  write.csv(mm, file.path(tab_d, "65_lmm_group_effect.csv"), row.names = FALSE)
  lg("LMM group effect: ", paste(sprintf("%s p=%.3g", mm$gene[mm$term == "group"],
                                        mm$p[mm$term == "group"]), collapse = " | "))
}

# one workbook for convenience
wb <- createWorkbook()
for (nm in c("65_stats_donor_whole_tissue", "65_stats_celltype", "65_lmm_group_effect",
             "65_pseudobulk_donor_tissue", "65_dotplot_celltype")) {
  f <- file.path(tab_d, paste0(nm, ".csv"))
  if (file.exists(f)) addWorksheet(wb, substr(nm, 1, 31)) else next
  writeData(wb, substr(nm, 1, 31), read.csv(f, stringsAsFactors = FALSE))
}
saveWorkbook(wb, file.path(tab_d, "65_fourgene_donor_stats.xlsx"), overwrite = TRUE)

# =================================================================== figures =
# ---- F1: dot plot, 4 genes x cell types ------------------------------------
ymax <- max(dot$pct)
p1f <- ggplot(dot, aes(celltype, gene)) +
  geom_point(aes(size = pct, colour = scaled)) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027",
                         midpoint = 0, name = "scaled mean\nexpression",
                         guide = guide_colourbar(barheight = unit(2.0, "cm"),
                                                 barwidth = unit(0.36, "cm"))) +
  scale_size_area(max_size = 7, breaks = c(1, 5, 10, 20), name = "% nuclei\nexpressing",
                  guide = guide_legend(keyheight = unit(0.42, "cm"),
                                       override.aes = list(colour = "grey35"))) +
  labs(title = "Four model genes across the 12 annotated cell types",
       subtitle = sprintf("snRNA-seq, %s nuclei, 25 donors; sizes are raw detection rates (max %.1f%%)",
                          format(ncol(cnt), big.mark = ","), ymax),
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 9),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.62, "cm"),
        legend.key.width = unit(0.42, "cm"),
        plot.margin = margin(10, 12, 8, 8))
ggsave(file.path(fig_d, "65_dotplot_celltype.png"), p1f, width = 7.6, height = 3.4, dpi = 300)
ggsave(file.path(fig_d, "65_dotplot_celltype.pdf"), p1f, width = 7.6, height = 3.4)

# ---- F2: dot plot split by group -------------------------------------------
p2f <- ggplot(dotg, aes(celltype, gene)) +
  geom_point(aes(size = pct, colour = avg)) +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_gradientn(colours = c("#F7F7F7", "#FDAE61", "#D73027", "#7F0000"),
                         trans = "sqrt",
                         name = "mean log-\nnormalised\nexpression\n(sqrt scale)") +
  scale_size_area(max_size = 5.5, name = "% nuclei\nexpressing") +
  labs(title = "Same panel split by group (raw scale, comparable across panels)",
       subtitle = "cell-level descriptive view; composition differences are not controlled here - use the donor-level panels for inference",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(face = "italic", size = 8),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "65_dotplot_celltype_by_group.png"), p2f, width = 10.4, height = 3.6, dpi = 300)
ggsave(file.path(fig_d, "65_dotplot_celltype_by_group.pdf"), p2f, width = 10.4, height = 3.6)

# ---- F3: FeaturePlot, 4 panels ---------------------------------------------
cen <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = data.frame(um, celltype = md$celltype),
                 FUN = median)
cen$pct <- as.numeric(table(md$celltype)[as.character(cen$celltype)]) / nrow(md) * 100
cen <- cen[cen$pct >= 1.5, ]

fp_panel <- function(g) {
  v <- as.numeric(dat[g, ])
  o <- order(v)
  d <- data.frame(UMAP_1 = um[o, 1], UMAP_2 = um[o, 2], v = v[o])
  cap <- as.numeric(quantile(v[v > 0], 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(v)
  d$v[d$v > cap] <- cap
  ggplot(d, aes(UMAP_1, UMAP_2)) +
    rasterise(geom_point(data = d[d$v == 0, ], colour = "grey88", size = 0.05, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = d[d$v > 0, ], aes(colour = v), size = 0.09, stroke = 0), dpi = 300) +
    geom_text(data = cen, aes(label = celltype), size = 1.9, colour = "grey25") +
    scale_colour_gradientn(colours = c("#FFEDA0", "#FEB24C", "#F03B20", "#7F0000"),
                           name = "log-norm\nexpression") +
    labs(title = g, subtitle = sprintf("%.2f%% nuclei > 0", mean(v > 0) * 100),
         x = "UMAP_1", y = "UMAP_2") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "italic"),
          legend.key.height = unit(0.55, "cm"),
          legend.key.width = unit(0.34, "cm"))
}
p3f <- wrap_plots(lapply(GENES, fp_panel), nrow = 2) +
  plot_annotation(title = "Feature plots of the four model genes",
                  subtitle = "all 150,750 nuclei; grey = not detected; colour capped at the 99th percentile of expressing nuclei")
ggsave(file.path(fig_d, "65_featureplot_fourgene.png"), p3f, width = 9.2, height = 7.4, dpi = 300)
ggsave(file.path(fig_d, "65_featureplot_fourgene.pdf"), p3f, width = 9.2, height = 7.4)

# ---- F4: FeaturePlot split by group (positive cells only) ------------------
grp_of <- as.character(md$group)
mk <- do.call(rbind, lapply(GENES, function(g) {
  v <- as.numeric(dat[g, ])
  do.call(rbind, lapply(GRP, function(gr) {
    idx <- which(grp_of == gr & v > 0)
    if (!length(idx)) return(NULL)
    data.frame(gene = g, group = gr, UMAP_1 = um[idx, 1], UMAP_2 = um[idx, 2],
               v = v[idx], stringsAsFactors = FALSE)
  }))
}))
mk$gene <- factor(mk$gene, levels = GENES); mk$group <- factor(mk$group, levels = GRP)
mkb <- data.frame(UMAP_1 = um[, 1], UMAP_2 = um[, 2])
mkb <- mkb[sample(nrow(mkb), 30000), ]
p4f <- ggplot(mk, aes(UMAP_1, UMAP_2)) +
  rasterise(geom_point(data = mkb, colour = "grey93", size = 0.05, stroke = 0), dpi = 300) +
  rasterise(geom_point(aes(colour = v), size = 0.12, stroke = 0), dpi = 300) +
  facet_grid(gene ~ group) +
  scale_colour_gradientn(colours = c("#FFEDA0", "#FEB24C", "#F03B20", "#7F0000"),
                         name = "log-norm\nexpression") +
  labs(title = "Detecting nuclei per group (positive nuclei only)",
       subtitle = "grey = all nuclei; facets not equalised for nuclei number - compare shapes, not density",
       x = "UMAP_1", y = "UMAP_2") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(face = "italic", size = 7),
        strip.background = element_rect(fill = "grey94", colour = NA))
ggsave(file.path(fig_d, "65_featureplot_by_group.png"), p4f, width = 7.4, height = 9.0, dpi = 300)
ggsave(file.path(fig_d, "65_featureplot_by_group.pdf"), p4f, width = 7.4, height = 9.0)

# ---- F5: donor-level boxplots, whole tissue --------------------------------
pbd$gene <- factor(pbd$gene, levels = GENES)
lab_d <- data.frame(gene = factor(tab_donor$gene, levels = GENES),
                    y = Inf,
                    lab = sprintf("FDR %s", pfmt(tab_donor$fdr)))
p5f <- ggplot(pbd, aes(group, log1p_cpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.55, linewidth = 0.3) +
  geom_point(aes(colour = group), size = 1.5, stroke = 0.2,
             position = position_jitter(width = 0.12, seed = SEED)) +
  geom_text(data = lab_d, aes(x = 2, y = y, label = lab), vjust = 1.4,
            size = 2.9, inherit.aes = FALSE) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  labs(title = "Donor-level pseudobulk expression (whole tissue)",
       subtitle = "one point = one donor; CPM on summed nuclei counts, log1p; FDR = NF vs Disease (NICM + DoxCM)",
       x = NULL, y = "log1p(CPM)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_d, "65_donor_boxplot_whole_tissue.png"), p5f, width = 9.4, height = 3.4, dpi = 300)
ggsave(file.path(fig_d, "65_donor_boxplot_whole_tissue.pdf"), p5f, width = 9.4, height = 3.4)

# ---- F6: donor-level boxplots per cell type (top 3 per gene + all nuclei) --
top3 <- do.call(rbind, lapply(GENES, function(g) {
  # only cell types with at least MIN_DONOR donors in both arms can carry a test
  okct <- tab_ct$celltype[tab_ct$gene == g & !is.na(tab_ct$p_wilcox_NF_vs_Disease)]
  d <- dot[dot$gene == g & as.character(dot$celltype) %in% okct, ]
  d <- d[order(-d$avg), ]
  n <- min(3, nrow(d))
  data.frame(gene = g, celltype = as.character(d$celltype[1:n]), rank = seq_len(n),
             stringsAsFactors = FALSE)
}))
rows <- do.call(rbind, lapply(GENES, function(g) {
  t3 <- top3$celltype[top3$gene == g]
  a <- pbd[pbd$gene == g, c("sample", "group", "log1p_cpm")]
  a$gene <- g; a$celltype <- "All nuclei"; a$rank <- 0
  b <- pbc[pbc$gene == g & pbc$celltype %in% t3 & pbc$usable,
           c("sample", "group", "log1p_cpm", "celltype")]
  b$gene <- g; b$rank <- match(b$celltype, t3)
  rbind(a, b)
}))
lev <- unname(unlist(lapply(GENES, function(g)
  paste0(g, " - ", c("All nuclei", top3$celltype[top3$gene == g])))))
rows$panel <- factor(paste0(rows$gene, " - ", rows$celltype), levels = lev)
ann <- do.call(rbind, lapply(GENES, function(g) {
  t3 <- top3$celltype[top3$gene == g]
  data.frame(gene = g, celltype = c("All nuclei", t3),
             fdr = c(tab_donor$fdr[tab_donor$gene == g],
                     tab_ct$fdr[match(paste(g, t3), paste(tab_ct$gene, tab_ct$celltype))]),
             stringsAsFactors = FALSE)
}))
ann$panel <- factor(paste0(ann$gene, " - ", ann$celltype), levels = lev)
ann <- ann[ann$panel %in% rows$panel, ]
ann$lab <- ifelse(is.na(ann$fdr), "", paste0("FDR ", pfmt(ann$fdr)))
p6f <- ggplot(rows, aes(group, log1p_cpm, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.5, linewidth = 0.28) +
  geom_point(aes(colour = group), size = 1.15, stroke = 0.15,
             position = position_jitter(width = 0.13, seed = SEED)) +
  geom_text(data = ann, aes(x = 2, y = Inf, label = lab), vjust = 1.4,
            size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~ panel, scales = "free_y", ncol = 4, drop = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  labs(title = "Donor-level expression in the cell types that actually carry each gene",
       subtitle = "one point = one donor (>= 50 nuclei per sample x cell type); FDR = NF vs Disease; row = gene, leftmost panel = all nuclei",
       x = NULL, y = "log1p(CPM)") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6.8),
        strip.background = element_rect(fill = "grey94", colour = NA),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        panel.spacing = unit(0.26, "lines"))
ggsave(file.path(fig_d, "65_donor_boxplot_by_celltype.png"), p6f, width = 11.2,
       height = 7.0, dpi = 300, limitsize = FALSE)
ggsave(file.path(fig_d, "65_donor_boxplot_by_celltype.pdf"), p6f, width = 11.2, height = 7.0)

# ---- F7: log2FC heatmap, gene x cell type ----------------------------------
hm <- tab_ct
hm$celltype <- factor(hm$celltype, levels = rev(CT))
hm$gene <- factor(hm$gene, levels = GENES)
hm$star <- ifelse(is.na(hm$fdr) & pmax(hm$pct_detect_NF, hm$pct_detect_Disease) >= 0.05,
                  "n.t.", plab(hm$fdr))
hm$lfc <- pmax(pmin(hm$log2FC_median, 1.5), -1.5)
# cells where the gene is essentially never detected carry no fold change
hm$lfc[pmax(hm$pct_detect_NF, hm$pct_detect_Disease) < 0.05] <- NA
hm$star[is.na(hm$lfc)] <- ""
p7f <- ggplot(hm, aes(gene, celltype, fill = lfc)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = star), size = 2.8, colour = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1.5, 1.5), na.value = "grey92",
                       name = "log2 FC\nDisease / NF\n(median pseudobulk)") +
  labs(title = "Donor-level fold change per cell type",
       subtitle = "* FDR < 0.05, ** < 0.01, *** < 0.001 (NF vs Disease)\ngrey = never detected; n.t. = not testable (< 4 donors per group)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(face = "italic", size = 9),
        panel.grid = element_blank())
ggsave(file.path(fig_d, "65_heatmap_log2fc_celltype.png"), p7f, width = 5.6, height = 4.6, dpi = 300)
ggsave(file.path(fig_d, "65_heatmap_log2fc_celltype.pdf"), p7f, width = 5.6, height = 4.6)

# ---- F8: detection rate per donor ------------------------------------------
p8f <- ggplot(pbd, aes(group, pct_detect, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, linewidth = 0.3) +
  geom_point(aes(colour = group), size = 1.4,
             position = position_jitter(width = 0.12, seed = SEED)) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  labs(title = "Per-donor detection rate (fraction of nuclei with > 0 counts)",
       subtitle = "normalisation-free readout; a shift here cannot be explained by CPM scaling",
       x = NULL, y = "% nuclei expressing") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_d, "65_detection_rate_donor.png"), p8f, width = 9.4, height = 3.2, dpi = 300)
ggsave(file.path(fig_d, "65_detection_rate_donor.pdf"), p8f, width = 9.4, height = 3.2)

lg("figures written to ", fig_d)
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_65.txt"))
lg("=== done ===")
