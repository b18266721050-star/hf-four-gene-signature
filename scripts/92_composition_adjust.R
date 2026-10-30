###############################################################################
## 92_composition_adjust.R
## GSE57338 -- does the "everything myeloid goes down" pattern survive
## adjustment for cell composition?
##
## Logic
##   91 showed that failing myocardium has markedly LOWER xCell macrophage and
##   endothelial scores and HIGHER fibroblast scores, and that a composition-only
##   classifier reaches validation AUC 0.979 (four-gene model: 0.993). If the
##   uniform down-regulation of the 70 candidate genes is carried by that shift,
##   adding the compartment scores as covariates should attenuate their log2FC.
##   If log2FC survives, the change is per-cell (transcriptional), not
##   compositional -- and the Discussion has to say so.
##
## Caveat stated in the output
##   The compartment scores are strongly associated with group, so adjusting for
##   them is over-adjustment and the attenuation is an UPPER bound on the
##   compositional contribution. Reported as a sensitivity analysis, never as the
##   primary model.
##
## Inputs
##   clean_data/GSE57338_clean_expression_matrix.csv
##   results/Composition/01_tables/T04_composition_merged.csv
##   results/HubGene/GSE57338_candidate_genes.csv   (the 70-gene list)
##   results/Nomogram/04_data/modelpredictions_{train,validation}.csv
##
## Outputs (results/Composition/)
##   01_tables/T11_adjusted_deg_summary.csv
##   01_tables/T12_adjusted_70_genes.csv
##   01_tables/T13_adjusted_key_genes.csv
##   01_tables/T14_incremental_value_auc.csv
##   02_figures/F05_log2fc_before_after_adjustment.{pdf,png}
##   02_figures/F06_incremental_auc.{pdf,png}
##   03_logs/92_run_log.txt , sessionInfo_92.txt
##
## Run:  Rscript scripts/92_composition_adjust.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(pROC)
  library(glmnet)
  library(ggplot2)
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
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank()))

LFC_CUT  <- 0.585
PADJ_CUT <- 0.05
GENES    <- c("AQP3", "GGT5", "ISLR", "LAD1")
HUBS     <- c("CD163", "VSIG4", "S100A9", "FPR1", "STAT3", "ANPEP", "IL1R2")

## ===========================================================================
## 1. Data
## ===========================================================================
hr("=")
log_msg("Composition-adjusted re-analysis (GSE57338)")
hr("=")
log_msg("[1] Loading ...")

comp <- data.table::fread(file.path(d_tab, "T04_composition_merged.csv"))
expr_dt <- data.table::fread(file.path(proj_dir, "clean_data",
                                       "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

## align
stopifnot(setequal(comp$sample, colnames(expr)))
comp <- comp[match(colnames(expr), sample)]
comp[, group2 := factor(ifelse(group == "NF", "Control", "Disease"),
                        levels = c("Control", "Disease"))]
comp[, sex := factor(sex, levels = c("Female", "Male"))]
if (anyNA(comp$age_num)) {
  comp[is.na(age_num), age_num := median(comp$age_num, na.rm = TRUE)]
}

## compartment covariates (xCell scores, z-scaled)
COMPART <- c(Macrophages = "xcell_Macrophages",
             Fibroblasts = "xcell_Fibroblasts",
             Endothelial = "xcell_Endothelial_cells",
             Myocytes = "xcell_Myocytes")
for (nm in names(COMPART)) {
  v <- comp[[COMPART[[nm]]]]
  comp[, (paste0("z_", nm)) := as.numeric(scale(v))]
}
zcols <- paste0("z_", names(COMPART))
log_msg("  covariates: ", paste(zcols, collapse = ", "))

cand <- data.table::fread(file.path(proj_dir, "results", "HubGene",
                                    "GSE57338_candidate_genes.csv"))
cand70 <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
log_msg("  candidate genes (WGCNA x DEG): ", length(cand70))

## ===========================================================================
## 2. limma: unadjusted vs composition-adjusted
## ===========================================================================
hr()
log_msg("[2] limma, unadjusted vs adjusted ...")

fit <- function(rhs, tag) {
  form <- as.formula(paste("~", rhs))
  design <- model.matrix(form, data = comp)
  f <- lmFit(expr, design)
  f <- eBayes(f, trend = TRUE, robust = TRUE)
  tt <- topTable(f, coef = "group2Disease", number = Inf, sort.by = "none",
                 adjust.method = "BH")
  d <- data.table(GeneSymbol = rownames(tt), log2FC = round(tt$logFC, 4),
                  padj = tt$adj.P.Val, AveExpr = tt$AveExpr)
  d[, sig := !is.na(padj) & padj < PADJ_CUT & abs(log2FC) > LFC_CUT]
  log_msg(sprintf("  %-14s sig=%3d (up %3d / down %3d)", tag, sum(d$sig),
                  sum(d$sig & d$log2FC > 0), sum(d$sig & d$log2FC < 0)))
  d
}

d_raw  <- fit("group2 + sex + age_num", "unadjusted")
d_adj  <- fit(paste("group2 + sex + age_num +", paste(zcols, collapse = " + ")),
              "adjusted")
setnames(d_raw, c("log2FC", "padj", "sig"), c("log2FC_raw", "padj_raw", "sig_raw"))
setnames(d_adj, c("log2FC", "padj", "sig"), c("log2FC_adj", "padj_adj", "sig_adj"))
m <- merge(d_raw[, .(GeneSymbol, log2FC_raw, padj_raw, sig_raw)],
           d_adj[, .(GeneSymbol, log2FC_adj, padj_adj, sig_adj)], by = "GeneSymbol")
m[, delta_log2FC := round(log2FC_adj - log2FC_raw, 4)]
m[, sign_flip := sign(log2FC_raw) != sign(log2FC_adj)]

summary_tab <- data.table(
  metric = c("significant genes (unadjusted)",
             "significant genes (composition-adjusted)",
             "retained among unadjusted DEGs (%)",
             "direction flipped among unadjusted DEGs (n)",
             "70 candidate genes significant (unadjusted)",
             "70 candidate genes significant (adjusted)",
             "70 candidate genes direction flipped (n)",
             "median |log2FC| of DEGs (unadjusted)",
             "median |log2FC| of DEGs (adjusted)",
             "median log2FC of 70 candidates (unadjusted)",
             "median log2FC of 70 candidates (adjusted)"),
  value = c(sum(m$sig_raw),
            sum(m$sig_adj),
            round(100 * sum(m$sig_raw & m$sig_adj) / sum(m$sig_raw), 1),
            sum(m$sig_raw & m$sign_flip),
            sum(m[GeneSymbol %in% cand70]$sig_raw),
            sum(m[GeneSymbol %in% cand70]$sig_adj),
            sum(m[GeneSymbol %in% cand70]$sign_flip),
            round(median(abs(m[sig_raw == TRUE]$log2FC_raw)), 4),
            round(median(abs(m[sig_raw == TRUE]$log2FC_adj)), 4),
            round(median(m[GeneSymbol %in% cand70]$log2FC_raw), 4),
            round(median(m[GeneSymbol %in% cand70]$log2FC_adj), 4))
)
write_csv(summary_tab, file.path(d_tab, "T11_adjusted_deg_summary.csv"))
log_msg("  summary:")
for (i in seq_len(nrow(summary_tab))) {
  log_msg(sprintf("    %-52s %s", summary_tab$metric[i], summary_tab$value[i]))
}

## 70 candidate genes in detail
det70 <- merge(data.table(GeneSymbol = cand70), m, by = "GeneSymbol", all.x = TRUE)
det70 <- merge(det70, cand[, .(GeneSymbol, Regulation, module, MM_kME)],
               by = "GeneSymbol", all.x = TRUE)
det70[, attenuation_pct := round(100 * (1 - abs(log2FC_adj) / abs(log2FC_raw)), 1)]
det70[, abs_raw := abs(log2FC_raw)]
setorder(det70, -abs_raw)
write_csv(det70, file.path(d_tab, "T12_adjusted_70_genes.csv"))
log_msg("  70 genes: median attenuation ",
        round(median(det70$attenuation_pct, na.rm = TRUE), 1), "%",
        "; genes with >50% attenuation: ",
        sum(det70$attenuation_pct > 50, na.rm = TRUE), "/", nrow(det70),
        "; becoming non-significant: ",
        sum(!det70$sig_adj, na.rm = TRUE))

key <- m[GeneSymbol %in% c(GENES, HUBS)]
key[, role := ifelse(GeneSymbol %in% GENES, "model gene", "PPI hub")]
key[, attenuation_pct := round(100 * (1 - abs(log2FC_adj) / abs(log2FC_raw)), 1)]
key[, abs_raw := abs(log2FC_raw)]
setorder(key, role, -abs_raw)
write_csv(key, file.path(d_tab, "T13_adjusted_key_genes.csv"))
for (i in seq_len(nrow(key))) {
  log_msg(sprintf("    %-8s %-10s raw=%+.3f adj=%+.3f (%+5.1f%%)  padj_adj=%.2e",
                  key$GeneSymbol[i], key$role[i], key$log2FC_raw[i],
                  key$log2FC_adj[i], key$attenuation_pct[i], key$padj_adj[i]))
}

## ===========================================================================
## 3. Incremental value: does the four-gene panel add anything to composition?
## ===========================================================================
hr()
log_msg("[3] Incremental value of the four genes over composition ...")

pred <- rbind(
  data.table::fread(file.path(proj_dir, "results", "Nomogram", "04_data", "modelpredictions_train.csv")),
  data.table::fread(file.path(proj_dir, "results", "Nomogram", "04_data", "modelpredictions_validation.csv"))
)
pred <- unique(pred, by = "sample")
comp <- merge(comp, pred[, .(sample, set, lp_M2, AQP3_z, GGT5_z, ISLR_z, LAD1_z)],
              by = "sample", all.x = TRUE)
stopifnot(!anyNA(comp$lp_M2))

y <- as.integer(comp$group2 == "Disease")
tr <- comp$set == "train"; va <- comp$set == "validation"
Xcomp <- scale(as.matrix(comp[, c("xcell_Macrophages", "xcell_Fibroblasts",
                                  "xcell_Endothelial_cells", "xcell_Myocytes",
                                  "xcell_CD8_T_cells", "xcell_ImmuneScore"),
                            with = FALSE]))
Xgene <- as.matrix(comp[, .(AQP3_z, GGT5_z, ISLR_z, LAD1_z)])

auc_va <- function(pr) {
  r <- pROC::roc(y[va], pr, quiet = TRUE, direction = "<")
  c(auc = as.numeric(pROC::auc(r)), low = as.numeric(pROC::ci.auc(r))[1],
    high = as.numeric(pROC::ci.auc(r))[3])
}
ridge_fit <- function(Xtr, Xva) {
  cv <- cv.glmnet(Xtr, y[tr], family = "binomial", alpha = 0, nfolds = 5,
                  type.measure = "auc")
  as.numeric(predict(cv, newx = Xva, s = "lambda.min", type = "response"))
}
glm_fit <- function(Xtr, Xva) {
  dtr <- data.frame(Xtr, y = y[tr]); dva <- data.frame(Xva)
  f <- glm(y ~ ., data = dtr, family = binomial)
  predict(f, newdata = dva, type = "response")
}

rows <- list()
pr_comp <- ridge_fit(Xcomp[tr, ], Xcomp[va, ]); a <- auc_va(pr_comp)
rows[["composition only (ridge)"]] <- data.table(model = "Composition only (ridge)",
                                                 n_terms = ncol(Xcomp), auc = a["auc"],
                                                 ci_lower = a["low"], ci_upper = a["high"])
pr <- glm_fit(Xgene[tr, ], Xgene[va, ]); a <- auc_va(pr)
rows[["four genes only (glm)"]] <- data.table(model = "Four genes only (glm)",
                                              n_terms = 4, auc = a["auc"],
                                              ci_lower = a["low"], ci_upper = a["high"])
pr_both <- ridge_fit(cbind(Xcomp, Xgene)[tr, ], cbind(Xcomp, Xgene)[va, ]); a <- auc_va(pr_both)
rows[["composition + four genes (ridge)"]] <-
  data.table(model = "Composition + four genes (ridge)", n_terms = ncol(Xcomp) + 4,
             auc = a["auc"], ci_lower = a["low"], ci_upper = a["high"])
## DeLong test: composition vs composition + genes (paired ROC on the same samples)
r1 <- pROC::roc(y[va], pr_comp, quiet = TRUE, direction = "<")
r2 <- pROC::roc(y[va], pr_both, quiet = TRUE, direction = "<")
dl <- pROC::roc.test(r1, r2, method = "delong")
log_msg(sprintf("  DeLong, composition vs composition + 4 genes: p = %.3f", dl$p.value))
inc <- rbindlist(rows)
inc[, auc := round(as.numeric(auc), 4)]
inc[, ci_lower := round(as.numeric(ci_lower), 4)]
inc[, ci_upper := round(as.numeric(ci_upper), 4)]
write_csv(inc, file.path(d_tab, "T14_incremental_value_auc.csv"))
for (i in seq_len(nrow(inc))) {
  log_msg(sprintf("    %-34s terms=%2d  validation AUC=%.3f (%.3f-%.3f)",
                  inc$model[i], inc$n_terms[i], inc$auc[i], inc$ci_lower[i], inc$ci_upper[i]))
}

## ===========================================================================
## 4. Figures
## ===========================================================================
hr()
log_msg("[4] Figures ...")

sc <- m[sig_raw == TRUE]
sc[, group70 := ifelse(GeneSymbol %in% cand70, "70 candidate genes", "other DEGs")]
p5 <- ggplot(sc, aes(log2FC_raw, log2FC_adj)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(data = sc[group70 == "other DEGs"], colour = "grey70", size = 1, alpha = 0.6) +
  geom_point(data = sc[group70 == "70 candidate genes"], colour = "#C53030", size = 1.4,
             alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, colour = "#2B6CB0", linewidth = 0.6,
              formula = y ~ x) +
  labs(x = "log2FC, unadjusted", y = "log2FC, adjusted for composition",
       title = "Effect of composition adjustment on differential expression",
       subtitle = sprintf("%d DEGs; red = 70 candidate genes; dashed = identity", nrow(sc)))
save_fig(p5, "F05_log2fc_before_after_adjustment", 6.2, 5.4)

inc[, label := factor(model, levels = rev(model))]
p6 <- ggplot(inc, aes(auc, label, fill = grepl("four genes", model))) +
  geom_col(width = 0.62) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.18, orientation = "y") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = "#4C78A8"), guide = "none") +
  coord_cartesian(xlim = c(0.5, 1.02)) +
  labs(x = "Validation AUC (95% CI)", y = NULL,
       title = "Does the four-gene panel add anything to composition?")
save_fig(p6, "F06_incremental_auc", 7.2, 3.6)

hr("=")
log_msg("Done.")
hr("=")
writeLines(run_log, file.path(d_log, "92_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_92.txt"))
