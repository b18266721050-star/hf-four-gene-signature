###############################################################################
## 96d_axis_replication.R
## Does the 70-gene myeloid axis of GSE57338 replicate outside GSE57338?
##
## Three independent layers, because "the genes correlate with macrophages"
## can be produced by three very different things:
##   L1  module cohesion      do the 70 genes co-vary with EACH OTHER?
##   L2  myeloid coupling     do they track myeloid abundance?
##   L3  disease direction    is the axis lower in failing myocardium?
##
## Every layer is computed twice:
##   *_raw     across all samples   (same estimand as the published rho = 0.59)
##   *_within  after removing each group's mean
##             (group fixed effect removed -- this is what separates true
##              co-regulation from "everything is lower in disease", which by
##              itself manufactures correlation between any two down-regulated
##              genes)
##
## Myeloid abundance is measured two ways:
##   xCell Macrophages   (same tool as the discovery cohort)
##   marker score        (mean z of pan-macrophage markers, 70 genes EXCLUDED
##                        from the marker set to avoid circularity)
##
## Run: Rscript scripts/96d_axis_replication.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

proj <- "E:/another data file/GEObulk analysis project"
outr <- file.path(proj, "results", "AxisExternal")
d_tab <- file.path(outr, "01_tables"); d_fig <- file.path(outr, "02_figures")
d_log <- file.path(outr, "03_logs");   d_dat <- file.path(outr, "04_data")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

con <- file(file.path(d_log, "96d_replication_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); writeLines(s, con); cat(s, "\n"); flush.console() }
hr  <- function() say(strrep("-", 74))

###############################################################################
## 0. definitions
###############################################################################
cand <- fread(file.path(proj, "results", "HubGene", "GSE57338_candidate_genes.csv"))
g70  <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
say("candidate70 n = ", length(g70))

## pan-macrophage markers; intersection with g70 is removed below so that the
## marker score can never be correlated with the axis merely by sharing genes
mac_markers <- c("CD68", "CD14", "CSF1R", "ITGAM", "TYROBP", "AIF1", "FCGR1A",
                 "FCGR2A", "MSR1", "MRC1", "C1QA", "C1QB", "C1QC", "LYZ",
                 "ITGB2", "CYBB", "TLR2", "FCER1G", "CTSS", "CTSB", "MS4A7",
                 "LILRB1", "MARCO", "SIGLEC1")
overlap <- intersect(mac_markers, g70)
mac_markers <- setdiff(mac_markers, g70)
say("macrophage markers: ", length(mac_markers), " used; dropped as axis members: ",
    ifelse(length(overlap) == 0, "none", paste(overlap, collapse = ",")))

key_genes <- c("AQP3", "GGT5", "LAD1", "CD163", "VSIG4", "SERPINA3", "S100A9")

###############################################################################
## helpers
###############################################################################
auc_disease_vs_nf <- function(dis, nf) {          # AUC with disease as positive
  dis <- dis[!is.na(dis)]; nf <- nf[!is.na(nf)]
  if (length(dis) < 3 || length(nf) < 3) return(NA_real_)
  r <- rank(c(dis, nf)); n1 <- length(dis); n2 <- length(nf)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

## remove each group's mean from every gene -> group fixed effect gone
center_within_group <- function(e, grp) {
  out <- e
  for (lv in unique(grp)) {
    idx <- which(grp == lv)
    if (length(idx) > 1) out[, idx] <- out[, idx] - rowMeans(out[, idx, drop = FALSE])
  }
  out
}

## median pairwise Spearman between the axis genes themselves
module_cohesion <- function(e, g) {
  if (length(g) < 3) return(NA_real_)
  cm <- cor(t(e[g, , drop = FALSE]), method = "spearman")
  median(cm[upper.tri(cm)])
}

pc1_var <- function(e, g) {
  if (length(g) < 3) return(NA_real_)
  x <- t(e[g, , drop = FALSE])
  x <- scale(x)
  s <- svd(x, nu = 0, nv = 0)$d
  s[1]^2 / sum(s^2) * 100
}

rho_to_score <- function(e, g, score) {
  sapply(g, function(x) {
    if (sd(e[x, ]) == 0) return(NA_real_)
    suppressWarnings(cor(e[x, ], score, method = "spearman"))
  })
}

zscore_axis <- function(e, g) {
  z <- t(scale(t(e[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

###############################################################################
## 1. discovery cohort, recomputed here so all three cohorts use one code path
###############################################################################
hr(); say("[1] discovery cohort GSE57338")
tr <- fread(file.path(proj, "clean_data", "GSE57338_clean_expression_matrix.csv"),
            data.table = FALSE)
et <- as.matrix(tr[, -1, drop = FALSE]); rownames(et) <- tr[[1]]
rm(tr); invisible(gc())
comp <- fread(file.path(proj, "results", "Composition", "01_tables",
                        "T04_composition_merged.csv"))
et <- et[, comp$sample, drop = FALSE]
grp_t <- comp$group
say("  ", nrow(et), " genes x ", ncol(et), " samples ; groups: ",
    paste(names(table(grp_t)), table(grp_t), sep = "=", collapse = ", "))

g70_t <- intersect(g70, rownames(et))
say("  70-gene coverage: ", length(g70_t), "/", length(g70))

mac_t <- comp$xcell_Macrophages
mk_t  <- colMeans(t(scale(t(et[intersect(mac_markers, rownames(et)), , drop = FALSE]))),
                  na.rm = TRUE)
say("  markers used: ", length(intersect(mac_markers, rownames(et))),
    " ; Macrophages zero% = ", round(100 * mean(mac_t == 0), 1))

et_w <- center_within_group(et, grp_t)

dis <- data.table(
  cohort = "GSE57338",
  n_sample = ncol(et), n_gene = length(g70_t),
  mac_zero_pct = round(100 * mean(mac_t == 0), 1),
  rho_mac_raw = median(rho_to_score(et, g70_t, mac_t), na.rm = TRUE),
  rho_mac_within = median(rho_to_score(et_w, g70_t, mac_t), na.rm = TRUE),
  rho_marker_raw = median(rho_to_score(et, g70_t, mk_t), na.rm = TRUE),
  rho_marker_within = median(rho_to_score(et_w, g70_t, mk_t), na.rm = TRUE),
  cohesion_raw = module_cohesion(et, g70_t),
  cohesion_within = module_cohesion(et_w, g70_t),
  pc1_var_pct = pc1_var(et, g70_t)
)

axis_t <- zscore_axis(et, g70_t)
dis[, axis_auc := auc_disease_vs_nf(axis_t[grp_t != "NF"], axis_t[grp_t == "NF"])]
dis[, axis_delta := median(axis_t[grp_t != "NF"]) - median(axis_t[grp_t == "NF"])]
dis[, wilcox_p := signif(suppressWarnings(wilcox.test(
  axis_t[grp_t != "NF"], axis_t[grp_t == "NF"]))$p.value, 3)]
say("  median rho_mac raw = ", round(dis$rho_mac_raw, 3),
    " ; within = ", round(dis$rho_mac_within, 3))
say("  cohesion raw = ", round(dis$cohesion_raw, 3),
    " ; within = ", round(dis$cohesion_within, 3))

rho_store <- list(data.table(cohort = "GSE57338", gene = g70_t,
                             rho_mac_raw = rho_to_score(et, g70_t, mac_t),
                             rho_mac_within = rho_to_score(et_w, g70_t, mac_t),
                             rho_marker_raw = rho_to_score(et, g70_t, mk_t),
                             rho_marker_within = rho_to_score(et_w, g70_t, mk_t)))
grp_store <- list(data.table(cohort = "GSE57338", sample = comp$sample,
                             group = grp_t, axis = axis_t))

###############################################################################
## 2. external cohorts
###############################################################################
run_cohort <- function(cohort, xcell_tag, id_col = "sample") {
  hr(); say("[2] ", cohort)
  e <- readRDS(file.path(d_dat, sprintf("expr_%s_symbol.rds", cohort)))
  info <- fread(file.path(d_dat, sprintf("sample_group_%s.csv", cohort)))
  idx <- match(colnames(e), info[[id_col]])
  info <- info[idx]
  stopifnot(!anyNA(info$sample))
  grp <- info$group

  xc <- fread(file.path(d_dat, sprintf("xcell_%s.csv", xcell_tag)))
  xm <- as.matrix(xc[, -1]); rownames(xm) <- xc$cell_type
  xm <- xm[, colnames(e), drop = FALSE]
  mac <- xm["Macrophages", ]

  mk_genes <- setdiff(intersect(mac_markers, rownames(e)), g70)
  mk <- colMeans(t(scale(t(e[mk_genes, , drop = FALSE]))), na.rm = TRUE)
  say("  ", ncol(e), " samples ; markers used: ", length(mk_genes),
      " ; Macrophages zero% = ", round(100 * mean(mac == 0), 1))

  g <- intersect(g70, rownames(e))
  say("  70-gene coverage: ", length(g), "/", length(g70),
      " ; missing: ", paste(setdiff(g70, g), collapse = ","))

  ew <- center_within_group(e, grp)

  row <- data.table(
    cohort = cohort, n_sample = ncol(e), n_gene = length(g),
    mac_zero_pct = round(100 * mean(mac == 0), 1),
    rho_mac_raw = median(rho_to_score(e, g, mac), na.rm = TRUE),
    rho_mac_within = median(rho_to_score(ew, g, mac), na.rm = TRUE),
    rho_marker_raw = median(rho_to_score(e, g, mk), na.rm = TRUE),
    rho_marker_within = median(rho_to_score(ew, g, mk), na.rm = TRUE),
    cohesion_raw = module_cohesion(e, g),
    cohesion_within = module_cohesion(ew, g),
    pc1_var_pct = pc1_var(e, g)
  )
  ax <- zscore_axis(e, g)
  row[, axis_auc := auc_disease_vs_nf(ax[grp != "NF"], ax[grp == "NF"])]
  row[, axis_delta := median(ax[grp != "NF"]) - median(ax[grp == "NF"])]

  say("  rho_mac   raw = ", round(row$rho_mac_raw, 3),
      " ; within = ", round(row$rho_mac_within, 3))
  say("  rho_marker raw = ", round(row$rho_marker_raw, 3),
      " ; within = ", round(row$rho_marker_within, 3))
  say("  cohesion  raw = ", round(row$cohesion_raw, 3),
      " ; within = ", round(row$cohesion_within, 3))
  say("  axis AUC (disease vs NF) = ", round(row$axis_auc, 3),
      " ; delta = ", round(row$axis_delta, 3))
  for (gg in c("NF", "DCM", "ICM")) {
    if (any(grp == gg)) say("    ", gg, " (n=", sum(grp == gg), ") axis median = ",
                            round(median(ax[grp == gg]), 3))
  }
  ## Wilcoxon disease vs NF
  wt <- suppressWarnings(wilcox.test(ax[grp != "NF"], ax[grp == "NF"]))
  say("    Wilcoxon disease vs NF: P = ", signif(wt$p.value, 3))
  row[, wilcox_p := signif(wt$p.value, 3)]

  list(row = row,
       gene = data.table(cohort = cohort, gene = g,
                         rho_mac_raw = rho_to_score(e, g, mac),
                         rho_mac_within = rho_to_score(ew, g, mac),
                         rho_marker_raw = rho_to_score(e, g, mk),
                         rho_marker_within = rho_to_score(ew, g, mk)),
       grp = data.table(cohort = cohort, sample = colnames(e), group = grp, axis = ax),
       pval = wt$p.value)
}

r5406 <- run_cohort("GSE5406", "GSE5406_array", id_col = "sample")
r116  <- run_cohort("GSE116250", "GSE116250_rna", id_col = "title")

summary_tab <- rbind(dis, r5406$row, r116$row, fill = TRUE)
rho_store <- rbindlist(list(rho_store[[1]], r5406$gene, r116$gene), fill = TRUE)
grp_store <- rbindlist(list(grp_store[[1]], r5406$grp, r116$grp), fill = TRUE)

###############################################################################
## 3. concordance with the discovery cohort
###############################################################################
hr(); say("[3] concordance of gene-level rho with GSE57338")
base <- rho_store[cohort == "GSE57338"]
conc <- rbindlist(lapply(c("GSE5406", "GSE116250"), function(co) {
  x <- rho_store[cohort == co]
  mm <- merge(base[, .(gene, d_mac = rho_mac_raw, d_mac_w = rho_mac_within,
                       d_mk = rho_marker_raw)],
              x[, .(gene, e_mac = rho_mac_raw, e_mac_w = rho_mac_within,
                    e_mk = rho_marker_raw)], by = "gene")
  data.table(cohort = co, n_gene = nrow(mm),
             spearman_rho_mac = cor(mm$d_mac, mm$e_mac, method = "spearman"),
             spearman_rho_marker = cor(mm$d_mk, mm$e_mk, method = "spearman"),
             sign_concordance_mac = mean(sign(mm$d_mac) == sign(mm$e_mac)),
             sign_concordance_marker = mean(sign(mm$d_mk) == sign(mm$e_mk)))
}))
print(conc)
for (i in seq_len(nrow(conc))) {
  say(sprintf("  %-9s n=%2d  rho-mac vs discovery: Spearman %.3f, sign %.0f%%; marker: %.3f, sign %.0f%%",
              conc$cohort[i], conc$n_gene[i], conc$spearman_rho_mac[i],
              100 * conc$sign_concordance_mac[i], conc$spearman_rho_marker[i],
              100 * conc$sign_concordance_marker[i]))
}

## key genes
hr(); say("[4] key gene rho")
keytab <- rho_store[gene %in% key_genes]
keytab <- dcast(keytab, gene ~ cohort, value.var = "rho_mac_raw")
print(keytab)

###############################################################################
## 4. write tables
###############################################################################
fwrite(summary_tab, file.path(d_tab, "T01_axis_replication_summary.csv"))
fwrite(conc,        file.path(d_tab, "T02_concordance_with_discovery.csv"))
fwrite(rho_store,   file.path(d_tab, "T03_gene_level_rho.csv"))
fwrite(grp_store,   file.path(d_tab, "T04_axis_score_by_sample.csv"))
fwrite(keytab,      file.path(d_tab, "T05_key_gene_rho.csv"))
say(""); say("tables written to ", d_tab)

###############################################################################
## 5. figures
###############################################################################
pal <- c("GSE57338" = "#2c3e50", "GSE5406" = "#c0392b", "GSE116250" = "#2471a3")

## F01 gene-level rho, discovery vs external
pd <- merge(rho_store[cohort == "GSE57338", .(gene, discovery = rho_mac_raw)],
            rho_store[cohort != "GSE57338", .(gene, cohort, external = rho_mac_raw)],
            by = "gene")
pd[, key := gene %in% key_genes]
pd[, cohort := factor(cohort, levels = c("GSE5406", "GSE116250"))]
p1 <- ggplot(pd, aes(discovery, external)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(aes(colour = key), size = 1.9, alpha = .85) +
  ggrepel::geom_text_repel(data = pd[key == TRUE], aes(label = gene), size = 2.5,
                           min.segment.length = 0, box.padding = .3,
                           segment.colour = "grey60", segment.size = .25,
                           max.overlaps = 20) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = "#7f8c8d", `TRUE` = "#c0392b"),
                      guide = "none") +
  labs(x = "Spearman rho with xCell Macrophages, GSE57338",
       y = "Spearman rho with xCell Macrophages, external cohort",
       title = "Gene-level coupling to the myeloid compartment") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(d_fig, "F01_rho_concordance.png"), p1, width = 7.2, height = 3.7, dpi = 300)
ggsave(file.path(d_fig, "F01_rho_concordance.pdf"), p1, width = 7.2, height = 3.7)

## F02 axis score by group
gd <- copy(grp_store)
gd[, group := factor(group, levels = c("NF", "DCM", "ICM"))]
gd[, cohort := factor(cohort, levels = c("GSE57338", "GSE5406", "GSE116250"))]
p2 <- ggplot(gd, aes(group, axis, fill = group)) +
  geom_boxplot(outlier.size = .6, width = .6, alpha = .75, linewidth = .3) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50", linewidth = .3) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = c(NF = "#5d6d7e", DCM = "#c0392b", ICM = "#e59866"),
                    guide = "none") +
  labs(x = NULL, y = "70-gene axis score (mean z)",
       title = "Axis score by group in discovery and external cohorts") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(d_fig, "F02_axis_score_by_group.png"), p2, width = 7.2, height = 3.4, dpi = 300)
ggsave(file.path(d_fig, "F02_axis_score_by_group.pdf"), p2, width = 7.2, height = 3.4)

## F03 summary bars: raw vs within-group
sdl <- melt(summary_tab[, .(cohort, rho_mac_raw, rho_mac_within,
                            rho_marker_raw, rho_marker_within,
                            cohesion_raw, cohesion_within)],
            id.vars = "cohort")
metric_lab <- c(rho_mac_raw = "Macrophages", rho_mac_within = "Macrophages",
                rho_marker_raw = "marker score", rho_marker_within = "marker score",
                cohesion_raw = "cohesion", cohesion_within = "cohesion")
sdl[, metric := metric_lab[as.character(variable)]]
sdl[, adjust := fifelse(grepl("within", as.character(variable)), "within group", "raw")]
sdl[, cohort := factor(cohort, levels = c("GSE57338", "GSE5406", "GSE116250"))]
p3 <- ggplot(sdl, aes(metric, value, fill = cohort)) +
  geom_col(position = position_dodge(.8), width = .72, colour = "grey25", linewidth = .2) +
  geom_hline(yintercept = 0, linewidth = .3) +
  facet_wrap(~adjust) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "median Spearman rho", fill = NULL,
       title = "Effect of removing the group mean") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(d_fig, "F03_raw_vs_within.png"), p3, width = 7.0, height = 3.3, dpi = 300)
ggsave(file.path(d_fig, "F03_raw_vs_within.pdf"), p3, width = 7.0, height = 3.3)

say(""); say("figures written to ", d_fig)
hr(); say("DONE"); close(con)
cat("REPLICATION DONE\n")
