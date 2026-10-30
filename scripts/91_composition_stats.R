###############################################################################
## 91_composition_stats.R
## GSE57338 -- composition statistics: group differences, single-feature
## discrimination, correlation with the four-gene score, and a composition-only
## classifier evaluated on the SAME train / internal-validation split as M2.
##
## Question answered here
##   1. Which compartments differ between non-failing and failing myocardium,
##      and by how much (effect size, not only p)?
##   2. Can composition alone reproduce the AUC of the four-gene model? If yes,
##      the near-perfect discrimination is a tissue-composition contrast, not a
##      property of the four genes.
##   3. Are the four genes (and the PPI hubs) correlated with the compartment
##      they are supposed to mark? (ISLR vs fibroblasts is the critical one.)
##
## Inputs
##   results/Composition/01_tables/T04_composition_merged.csv
##   results/Nomogram/04_data/modelpredictions_{train,validation}.csv
##   clean_data/GSE57338_clean_expression_matrix.csv
##
## Outputs (results/Composition/)
##   01_tables/T05_xcell_group_comparison.csv
##   01_tables/T06_epic_group_comparison.csv
##   01_tables/T07_single_feature_auc.csv
##   01_tables/T08_composition_vs_genes_correlation.csv
##   01_tables/T09_composition_classifier_auc.csv
##   01_tables/T10_composition_pca.csv
##   02_figures/F01_key_compartment_boxplots.{pdf,png}
##   02_figures/F02_composition_gene_correlation.{pdf,png}
##   02_figures/F03_composition_pca.{pdf,png}
##   02_figures/F04_single_feature_auc.{pdf,png}
##   03_logs/91_run_log.txt , sessionInfo_91.txt
##
## Run:  Rscript scripts/91_composition_stats.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(glmnet)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(20261029)

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
out_root <- file.path(proj_dir, "results", "Composition")
d_tab <- file.path(out_root, "01_tables")
d_fig <- file.path(out_root, "02_figures")
d_log <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = ""); flush.console()
}
hr <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}
save_fig <- function(p, stem, w, h) {
  ggsave(file.path(d_fig, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(d_fig, paste0(stem, ".png")), p, width = w, height = h, dpi = 600)
}

COL_GROUP <- c(NF = "#4C78A8", DCM = "#E45756", ICM = "#F2A93B")
COL_GROUP2 <- c(Control = "#4C78A8", Disease = "#E45756")
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey95", colour = NA)))

## ===========================================================================
## 1. Data
## ===========================================================================
hr("=")
log_msg("Composition statistics (GSE57338, 313 samples)")
hr("=")
log_msg("[1] Loading ...")

comp <- data.table::fread(file.path(d_tab, "T04_composition_merged.csv"))
pred <- rbind(
  data.table::fread(file.path(proj_dir, "results", "Nomogram", "04_data", "modelpredictions_train.csv")),
  data.table::fread(file.path(proj_dir, "results", "Nomogram", "04_data", "modelpredictions_validation.csv"))
)
pred <- unique(pred, by = "sample")
comp <- merge(comp, pred[, .(sample, set, lp_M2, prob_M2_gene)], by = "sample", all.x = TRUE)
if (anyNA(comp$lp_M2)) stop("model predictions missing for some samples")
comp[, group2 := factor(ifelse(group == "NF", "Control", "Disease"),
                        levels = c("Control", "Disease"))]
log_msg("  ", nrow(comp), " samples; train = ", sum(comp$set == "train"),
        ", validation = ", sum(comp$set == "validation"))

xc_cols <- grep("^xcell_", names(comp), value = TRUE)
ep_cols <- grep("^epic_", names(comp), value = TRUE)
log_msg("  xCell features: ", length(xc_cols), " | EPIC features: ", length(ep_cols))

## expression of the four genes + PPI hubs, for the correlation analysis
expr_dt <- data.table::fread(file.path(proj_dir, "clean_data",
                                       "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")
HUBS  <- c("CD163", "VSIG4", "S100A9", "FPR1", "STAT3", "ANPEP", "IL1R2")
stopifnot(all(c(GENES, HUBS) %in% rownames(expr)))
gmat <- t(expr[c(GENES, HUBS), comp$sample, drop = FALSE])
colnames(gmat) <- c(GENES, HUBS)
comp <- cbind(comp, as.data.table(gmat))

## ===========================================================================
## 2. Group comparisons
## ===========================================================================
hr()
log_msg("[2] Group comparisons ...")

cliff_delta <- function(x, y) {           # x = Disease, y = Control
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  u <- as.numeric(wilcox.test(x, y, exact = FALSE)$statistic)
  2 * u / (length(x) * length(y)) - 1
}

group_compare <- function(cols, tag) {
  rows <- lapply(cols, function(cc) {
    d <- comp[[cc]]
    ax <- d[comp$group2 == "Disease"]; ay <- d[comp$group2 == "Control"]
    p_all <- tryCatch(wilcox.test(ax, ay, exact = FALSE)$p.value, error = function(e) NA_real_)
    p_dcm <- tryCatch(wilcox.test(d[comp$group == "DCM"], d[comp$group == "NF"],
                                  exact = FALSE)$p.value, error = function(e) NA_real_)
    p_icm <- tryCatch(wilcox.test(d[comp$group == "ICM"], d[comp$group == "NF"],
                                  exact = FALSE)$p.value, error = function(e) NA_real_)
    auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(comp$group2, d, quiet = TRUE,
                                                   direction = "<"))),
                    error = function(e) NA_real_)
    med <- tapply(d, comp$group, median, na.rm = TRUE)
    data.table(feature = sub(paste0("^", tag, "_"), "", cc),
               median_NF = round(med[["NF"]], 4),
               median_DCM = round(med[["DCM"]], 4),
               median_ICM = round(med[["ICM"]], 4),
               delta_disease_vs_NF = round(cliff_delta(ax, ay), 3),
               p_disease_vs_NF = signif(p_all, 4),
               p_DCM_vs_NF = signif(p_dcm, 4),
               p_ICM_vs_NF = signif(p_icm, 4),
               auc_disease_vs_NF = round(auc, 3))
  })
  out <- rbindlist(rows)
  out[, fdr_disease_vs_NF := signif(p.adjust(p_disease_vs_NF, "BH"), 4)]
  out[, fdr_DCM_vs_NF := signif(p.adjust(p_DCM_vs_NF, "BH"), 4)]
  out[, fdr_ICM_vs_NF := signif(p.adjust(p_ICM_vs_NF, "BH"), 4)]
  out[, direction := ifelse(delta_disease_vs_NF > 0, "higher in disease",
                            "lower in disease")]
  setorder(out, p_disease_vs_NF)
  out
}

cmp_x <- group_compare(xc_cols, "xcell")
cmp_e <- group_compare(ep_cols, "epic")
write_csv(cmp_x, file.path(d_tab, "T05_xcell_group_comparison.csv"))
write_csv(cmp_e, file.path(d_tab, "T06_epic_group_comparison.csv"))

log_msg("  xCell: ", sum(cmp_x$fdr_disease_vs_NF < 0.05), " / ", nrow(cmp_x),
        " features differ at FDR<0.05")
log_msg("  top 12 by p (delta = Cliff's delta, + means higher in disease):")
for (i in seq_len(min(12, nrow(cmp_x)))) {
  log_msg(sprintf("    %-22s NF=%7.3f  DCM=%7.3f  ICM=%7.3f  delta=%+.2f  p=%.2e  FDR=%.2e  AUC=%.2f",
                  cmp_x$feature[i], cmp_x$median_NF[i], cmp_x$median_DCM[i],
                  cmp_x$median_ICM[i], cmp_x$delta_disease_vs_NF[i],
                  cmp_x$p_disease_vs_NF[i], cmp_x$fdr_disease_vs_NF[i],
                  cmp_x$auc_disease_vs_NF[i]))
}
log_msg("  EPIC top 8:")
for (i in seq_len(min(8, nrow(cmp_e)))) {
  log_msg(sprintf("    %-18s NF=%6.3f  DCM=%6.3f  ICM=%6.3f  delta=%+.2f  FDR=%.2e  AUC=%.2f",
                  cmp_e$feature[i], cmp_e$median_NF[i], cmp_e$median_DCM[i],
                  cmp_e$median_ICM[i], cmp_e$delta_disease_vs_NF[i],
                  cmp_e$fdr_disease_vs_NF[i], cmp_e$auc_disease_vs_NF[i]))
}

## ===========================================================================
## 3. Single-feature discrimination
## ===========================================================================
hr()
log_msg("[3] Single-feature AUC (Disease vs Control) ...")
auc_tab <- rbind(
  cmp_x[, .(method = "xCell", feature, auc = auc_disease_vs_NF,
            delta = delta_disease_vs_NF, fdr = fdr_disease_vs_NF)],
  cmp_e[, .(method = "EPIC", feature, auc = auc_disease_vs_NF,
            delta = delta_disease_vs_NF, fdr = fdr_disease_vs_NF)]
)
auc_tab[, auc_dist_from_half := round(abs(auc - 0.5) + 0.5, 3)]
setorder(auc_tab, -auc_dist_from_half)
write_csv(auc_tab, file.path(d_tab, "T07_single_feature_auc.csv"))
log_msg("  best single features:")
for (i in seq_len(min(8, nrow(auc_tab)))) {
  log_msg(sprintf("    %-6s %-22s AUC=%.3f  (%s)", auc_tab$method[i], auc_tab$feature[i],
                  auc_tab$auc_dist_from_half[i],
                  ifelse(auc_tab$auc[i] > 0.5, "higher in disease", "lower in disease")))
}

## ===========================================================================
## 4. Correlation with the four genes / hubs / model score
## ===========================================================================
hr()
log_msg("[4] Composition vs gene expression ...")
targets <- c(GENES, HUBS, "lp_M2")
feat_use <- c(xc_cols, ep_cols)
cors <- rbindlist(lapply(feat_use, function(f) {
  rbindlist(lapply(targets, function(g) {
    ct <- suppressWarnings(cor.test(comp[[f]], comp[[g]], method = "spearman",
                                    exact = FALSE))
    data.table(feature = sub("^(xcell|epic)_", "", f), method = sub("_.*$", "", f),
               target = g, rho = round(unname(ct$estimate), 3),
               p = signif(ct$p.value, 3))
  }))
}))
cors[, fdr := signif(p.adjust(p, "BH"), 3)]
## xCell and EPIC share a few compartment names (Macrophages, CAFs, ...): keep an
## explicit method-labelled id so the heatmap does not overlay two rows
cors[, feature_lab := paste0(ifelse(method == "xcell", "xCell", "EPIC"), ": ",
                             gsub("_", " ", feature))]
write_csv(cors, file.path(d_tab, "T08_composition_vs_genes_correlation.csv"))

key_feats <- c("xCell: Myocytes", "xCell: Fibroblasts", "xCell: Endothelial cells",
               "xCell: Macrophages", "xCell: Macrophages M2", "xCell: Pericytes",
               "xCell: CD8 T cells", "xCell: ImmuneScore", "xCell: StromaScore",
               "EPIC: CAFs", "EPIC: Macrophages")
for (g in c(GENES, "lp_M2")) {
  sub <- cors[target == g & feature_lab %in% key_feats][order(-abs(rho))]
  log_msg("  ", g, ": ", paste(sprintf("%s rho=%+.2f", sub$feature_lab, sub$rho),
                               collapse = "; "))
}

## ===========================================================================
## 5. Composition-only classifier (same split as M2)
## ===========================================================================
hr()
log_msg("[5] Composition-only classifier ...")

X <- scale(as.matrix(comp[, ..xc_cols]))
X[is.na(X)] <- 0
y <- as.integer(comp$group2 == "Disease")
tr <- comp$set == "train"; va <- comp$set == "validation"

ridge <- function(Xtr, ytr, Xva, yva, tag) {
  cvfit <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0, nfolds = 5,
                     type.measure = "auc")
  pr_va <- as.numeric(predict(cvfit, newx = Xva, s = "lambda.min", type = "response"))
  pr_tr <- as.numeric(predict(cvfit, newx = Xtr, s = "lambda.min", type = "response"))
  a_va <- as.numeric(pROC::auc(pROC::roc(yva, pr_va, quiet = TRUE, direction = "<")))
  a_tr <- as.numeric(pROC::auc(pROC::roc(ytr, pr_tr, quiet = TRUE, direction = "<")))
  ci <- as.numeric(pROC::ci.auc(pROC::roc(yva, pr_va, quiet = TRUE, direction = "<")))
  log_msg(sprintf("    %-28s train AUC=%.3f  validation AUC=%.3f (%.3f-%.3f)  lambda=%.4f",
                  tag, a_tr, a_va, ci[1], ci[3], cvfit$lambda.min))
  data.table(model = tag, train_auc = round(a_tr, 4), validation_auc = round(a_va, 4),
             ci_lower = round(ci[1], 4), ci_upper = round(ci[3], 4),
             lambda_min = signif(cvfit$lambda.min, 3),
             n_features = ncol(Xtr))
}

class_rows <- list()
class_rows[["xCell 67 features (ridge)"]] <- ridge(X[tr, ], y[tr], X[va, ], y[va],
                                                   "xCell 67 features (ridge)")
## single best feature, and the three-compartment model (myocyte / fibroblast / endothelial)
single <- function(col, tag) {
  dtr <- data.frame(score = X[tr, col], y = y[tr])
  dva <- data.frame(score = X[va, col])
  f <- glm(y ~ score, data = dtr, family = binomial)
  pr <- predict(f, newdata = dva, type = "response")
  a <- as.numeric(pROC::auc(pROC::roc(y[va], pr, quiet = TRUE, direction = "<")))
  ci <- as.numeric(pROC::ci.auc(pROC::roc(y[va], pr, quiet = TRUE, direction = "<")))
  log_msg(sprintf("    %-28s validation AUC=%.3f (%.3f-%.3f)", tag, a, ci[1], ci[3]))
  data.table(model = tag, train_auc = NA_real_, validation_auc = round(a, 4),
             ci_lower = round(ci[1], 4), ci_upper = round(ci[3], 4),
             lambda_min = NA_real_, n_features = 1)
}
i_myo <- which(xc_cols == "xcell_Myocytes")
i_fib <- which(xc_cols == "xcell_Fibroblasts")
i_end <- which(xc_cols == "xcell_Endothelial_cells")
class_rows[["xCell Myocytes only"]] <- single(i_myo, "xCell Myocytes only")
three <- c(i_myo, i_fib, i_end)
Xtr3 <- data.frame(score_myo = X[tr, i_myo], score_fib = X[tr, i_fib],
                   score_endo = X[tr, i_end], y = y[tr])
Xva3 <- data.frame(score_myo = X[va, i_myo], score_fib = X[va, i_fib],
                   score_endo = X[va, i_end])
f3 <- glm(y ~ ., data = Xtr3, family = binomial)
pr3 <- predict(f3, newdata = Xva3, type = "response")
a3 <- as.numeric(pROC::auc(pROC::roc(y[va], pr3, quiet = TRUE, direction = "<")))
ci3 <- as.numeric(pROC::ci.auc(pROC::roc(y[va], pr3, quiet = TRUE, direction = "<")))
log_msg(sprintf("    %-28s validation AUC=%.3f (%.3f-%.3f)",
                "Myocytes+Fibroblasts+Endothelial", a3, ci3[1], ci3[3]))
class_rows[["Myocytes+Fibro+Endo (3 features)"]] <-
  data.table(model = "Myocytes+Fibro+Endo (3 features)", train_auc = NA_real_,
             validation_auc = round(a3, 4), ci_lower = round(ci3[1], 4),
             ci_upper = round(ci3[3], 4), lambda_min = NA_real_, n_features = 3)
## reference rows (already published elsewhere, repeated here for the comparison table)
class_rows[["M2 four-gene model (reference)"]] <-
  data.table(model = "M2 four-gene model (reference)", train_auc = 0.9856,
             validation_auc = 0.9926, ci_lower = 0.9815, ci_upper = 1,
             lambda_min = NA_real_, n_features = 4)
class_rows[["GSM-order block only (reference)"]] <-
  data.table(model = "GSM-order block only (reference)", train_auc = NA_real_,
             validation_auc = 0.9618, ci_lower = 0.921, ci_upper = 1,
             lambda_min = NA_real_, n_features = 16)
class_tab <- rbindlist(class_rows)
write_csv(class_tab, file.path(d_tab, "T09_composition_classifier_auc.csv"))

## ===========================================================================
## 6. Figures
## ===========================================================================
hr()
log_msg("[6] Figures ...")

plot_feats <- c("xcell_Myocytes", "xcell_Fibroblasts", "xcell_Endothelial_cells",
                "xcell_Macrophages", "xcell_CD8_T_cells", "xcell_ImmuneScore")
lab_map <- c(xcell_Myocytes = "Myocytes", xcell_Fibroblasts = "Fibroblasts",
             xcell_Endothelial_cells = "Endothelial cells",
             xcell_Macrophages = "Macrophages", xcell_CD8_T_cells = "CD8+ T cells",
             xcell_ImmuneScore = "Immune score")
pl <- lapply(plot_feats, function(f) {
  dt <- comp[, .(group, value = get(f))]
  ggplot(dt, aes(group, value, fill = group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.62, alpha = 0.85) +
    scale_fill_manual(values = COL_GROUP, guide = "none") +
    labs(x = NULL, y = lab_map[[f]],
         subtitle = sprintf("p = %.1e", cmp_x[feature == sub("^xcell_", "", f)]$p_disease_vs_NF[1])) +
    theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
})
p1 <- wrap_plots(pl, ncol = 3)
save_fig(p1, "F01_key_compartment_boxplots", 9.6, 5.6)

## correlation heatmap: key features x (4 genes + hubs + score)
## xCell rows only (EPIC duplicates the names of a few compartments); EPIC is
## added once, explicitly, for the stromal / myeloid compartments it measures
## independently.
hm_xc <- c("xCell: Myocytes", "xCell: Fibroblasts", "xCell: Endothelial cells",
           "xCell: Pericytes", "xCell: Macrophages", "xCell: Macrophages M1",
           "xCell: Macrophages M2", "xCell: Monocytes", "xCell: CD8 T cells",
           "xCell: CD4 T cells", "xCell: Plasma cells", "xCell: ImmuneScore",
           "xCell: StromaScore")
hm_ep <- c("EPIC: CAFs", "EPIC: Macrophages")
hm_feats <- c(hm_xc, hm_ep)
hm_targets <- c(GENES, HUBS, "lp_M2")
hm <- cors[feature_lab %in% hm_feats & target %in% hm_targets]
stopifnot(nrow(hm) == length(hm_feats) * length(hm_targets))
hm[, feature_lab := factor(feature_lab, levels = rev(hm_feats))]
hm[, target := factor(target, levels = hm_targets)]
p2 <- ggplot(hm, aes(target, feature_lab, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(rho) >= 0.3, sprintf("%.2f", rho), "")),
            size = 2.6, colour = "grey15") +
  scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030",
                       midpoint = 0, limits = c(-0.8, 0.8), name = "Spearman rho") +
  labs(x = NULL, y = NULL,
       title = "Cell composition vs four-gene panel, PPI hubs and model score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p2, "F02_composition_gene_correlation", 8.4, 6.0)

## composition PCA
pca <- prcomp(scale(as.matrix(comp[, ..xc_cols])), center = FALSE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_dt <- data.table(sample = comp$sample, group = comp$group,
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2])
write_csv(pca_dt, file.path(d_tab, "T10_composition_pca.csv"))
eta <- function(v, g) {
  summary(lm(v ~ g))$r.squared
}
log_msg(sprintf("  composition PC1 %.1f%% (eta2 group = %.2f), PC2 %.1f%% (eta2 group = %.2f)",
                ve[1], eta(pca_dt$PC1, comp$group), ve[2], eta(pca_dt$PC2, comp$group)))
p3 <- ggplot(pca_dt, aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.5, alpha = 0.85) +
  stat_ellipse(level = 0.68, linewidth = 0.5) +
  scale_colour_manual(values = COL_GROUP) +
  labs(x = sprintf("PC1 (%.1f%% of composition variance)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "PCA of xCell composition scores",
       subtitle = sprintf("eta2 of group on PC1 = %.2f", eta(pca_dt$PC1, comp$group)))
save_fig(p3, "F03_composition_pca", 6.4, 5.2)

## single-feature AUC
top_auc <- rbind(auc_tab[method == "xCell"][order(-auc_dist_from_half)][1:12],
                 auc_tab[method == "EPIC"][order(-auc_dist_from_half)][1:6])
top_auc[, label := paste0(method, ": ", gsub("_", " ", feature))]
top_auc[, label := factor(label, levels = rev(label[order(auc_dist_from_half)]))]
p4 <- ggplot(top_auc, aes(auc_dist_from_half, label, fill = method)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c(xCell = "#4C78A8", EPIC = "#8C6BB1"), name = NULL) +
  coord_cartesian(xlim = c(0.5, 1)) +
  labs(x = "AUC, disease vs non-failing (direction: larger = higher in disease)",
       y = NULL, title = "Single composition feature discrimination")
save_fig(p4, "F04_single_feature_auc", 7.6, 5.4)

hr("=")
log_msg("Done.")
hr("=")
writeLines(run_log, file.path(d_log, "91_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_91.txt"))
