###############################################################################
## 137_trainonly_panel.R
##
## R10 of the pre-submission audit: what the four-gene panel would have looked
## like if feature selection had been done inside the training split only.
##
## Background
##   09b_ml_trainonly.R re-ran the three feature-selection algorithms on the
##   219 training samples. The three-method intersection then contained EIGHT
##   genes (AQP3, ASPN, CREB5, FCN3, GGT5, LAD1, RARRES1, SMOC2) and ISLR was
##   recovered by random forest only. The published panel (AQP3, GGT5, ISLR,
##   LAD1) therefore depends, for one of its four members, on selection over the
##   full 313 samples.
##
## What this script quantifies
##   1. Held-out discrimination of the eight-gene model fitted on the 219
##      training samples, evaluated on the same 94 held-out samples used for the
##      published internal validation.
##   2. Paired DeLong test against the published four-gene model on those 94.
##   3. External discrimination of the eight-gene model in GSE116250, the only
##      external cohort in which all eight genes are measurable (GSE5406 is a
##      HG-U133A array and lacks CREB5 and SMOC2).
##
## Outputs (results/ExternalValidation/01_tables)
##   T15_trainonly_panel_internal.csv
##   T16_trainonly_panel_external.csv
##   03_logs/137_run_log.txt , sessionInfo_137.txt
##
## Run: Rscript scripts/137_trainonly_panel.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

set.seed(20261030)
proj <- "E:/another data file/GEObulk analysis project"
d_tab <- file.path(proj, "results", "ExternalValidation", "01_tables")
d_log <- file.path(proj, "results", "ExternalValidation", "03_logs")
dir.create(d_log, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = ""); flush.console()
}
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

PANEL8 <- c("AQP3", "ASPN", "CREB5", "FCN3", "GGT5", "LAD1", "RARRES1", "SMOC2")
PANEL4 <- c("AQP3", "GGT5", "ISLR", "LAD1")

## ---------------------------------------------------------------------------
## 1. training data
## ---------------------------------------------------------------------------
log_msg("R10  train-only eight-gene panel")
expr_dt <- data.table::fread(file.path(proj, "clean_data", "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

tr <- data.table::fread(file.path(proj, "results", "Nomogram", "04_data", "modelpredictions_train.csv"))
va <- data.table::fread(file.path(proj, "results", "Nomogram", "04_data", "modelpredictions_validation.csv"))
tr_ids <- tr$sample; va_ids <- va$sample
log_msg("  training n = ", length(tr_ids), " ; held-out n = ", length(va_ids))
stopifnot(all(c(tr_ids, va_ids) %in% colnames(expr)))
stopifnot(all(PANEL8 %in% rownames(expr)))

Xtr <- t(expr[PANEL8, tr_ids, drop = FALSE])
Xva <- t(expr[PANEL8, va_ids, drop = FALSE])
mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd)
Ztr <- scale(Xtr, center = mu, scale = sdv)
Zva <- scale(Xva, center = mu, scale = sdv)

ytr <- tr$HF; yva <- va$HF
dtr <- as.data.frame(Ztr); names(dtr) <- PANEL8; dtr$y <- ytr
dva <- as.data.frame(Zva); names(dva) <- PANEL8
fit8 <- glm(y ~ ., data = dtr, family = binomial)
lp8 <- as.numeric(predict(fit8, newdata = dva))

## ---------------------------------------------------------------------------
## 2. internal held-out discrimination
## ---------------------------------------------------------------------------
r8 <- roc(yva, lp8, quiet = TRUE, direction = "<")
ci8 <- as.numeric(ci.auc(r8))
## published four-gene model on the same samples
r4 <- roc(yva, va$prob_M2_gene, quiet = TRUE, direction = "<")
ci4 <- as.numeric(ci.auc(r4))
dl <- roc.test(r8, r4, method = "delong")

log_msg(sprintf("  eight-gene (train-only selection), held-out AUC = %.4f (%.4f-%.4f)",
                as.numeric(auc(r8)), ci8[1], ci8[3]))
log_msg(sprintf("  four-gene (published, full-cohort selection), AUC = %.4f (%.4f-%.4f)",
                as.numeric(auc(r4)), ci4[1], ci4[3]))
log_msg(sprintf("  paired DeLong: Z = %.3f, p = %.3f", dl$statistic, dl$p.value))

## per-gene coefficients for the record
co <- summary(fit8)$coefficients
coef_tab <- data.table(gene = rownames(co)[-1],
                       beta = round(co[-1, 1], 4),
                       se = round(co[-1, 2], 4),
                       p = signif(co[-1, 4], 4))

int_tab <- data.table(
  model = c("Eight-gene, selected in the 219 training samples",
            "Four-gene, selected in all 313 samples (published)"),
  n_genes = c(8, 4),
  auc_heldout = round(c(as.numeric(auc(r8)), as.numeric(auc(r4))), 4),
  ci_low = round(c(ci8[1], ci4[1]), 4),
  ci_high = round(c(ci8[3], ci4[3]), 4),
  delong_p_vs_four_gene = c(signif(dl$p.value, 4), NA)
)
write_csv(int_tab, file.path(d_tab, "T15_trainonly_panel_internal.csv"))
log_msg("")
for (i in seq_len(nrow(coef_tab))) {
  log_msg(sprintf("    %-8s beta = %+7.4f  p = %s", coef_tab$gene[i], coef_tab$beta[i],
                  format(coef_tab$p[i], scientific = TRUE, digits = 3)))
}

## ---------------------------------------------------------------------------
## 3. external: GSE116250 (all eight genes measurable)
## ---------------------------------------------------------------------------
log_msg("")
log_msg("[3] GSE116250 (all eight genes present)")
e2 <- readRDS(file.path(proj, "results", "AxisExternal", "04_data", "expr_GSE116250_symbol.rds"))
g2 <- data.table::fread(file.path(proj, "results", "AxisExternal", "04_data", "sample_group_GSE116250.csv"))
## columns are titles; map to sample then group
grp <- g2$group[match(colnames(e2), g2$title)]
keep <- !is.na(grp)
Xe <- t(e2[PANEL8, keep, drop = FALSE])
## LAD1 carried the placeholder 999999999 in two control samples (NF11, NF12);
## exclude them so that both models are scored on the same 62 samples
Xe[Xe > 25] <- NA
Xe[colnames(e2)[keep] %in% c("NF11", "NF12"), ] <- NA
ge <- grp[keep]
ye <- as.integer(ge %in% c("DCM", "ICM"))
log_msg("  samples with group: ", length(ye), " ; cases ", sum(ye == 1), " controls ", sum(ye == 0))

## quantile-map every gene onto the training distribution, then apply the
## frozen coefficients (same rule as the published external validation)
qp <- function(x, ref) {
  out <- rep(NA_real_, length(x)); ok <- !is.na(x)
  q <- rank(x[ok], ties.method = "average") / (sum(ok) + 1)
  out[ok] <- as.numeric(quantile(ref, probs = q, type = 7, names = FALSE))
  out
}
Ze <- sapply(seq_along(PANEL8), function(j) qp(Xe[, j], Xtr[, j]))
colnames(Ze) <- PANEL8
de <- as.data.frame(scale(Ze, center = mu, scale = sdv)); names(de) <- PANEL8
lp_e <- as.numeric(predict(fit8, newdata = de))

ok <- !is.na(lp_e)
r_e <- roc(ye[ok], lp_e[ok], quiet = TRUE, direction = "<")
ci_e <- as.numeric(ci.auc(r_e))
log_msg(sprintf("  eight-gene external AUC = %.4f (%.4f-%.4f) on %d samples",
                as.numeric(auc(r_e)), ci_e[1], ci_e[3], sum(ok)))

## published four-gene predictions on the same cohort
pe <- data.table::fread(file.path(proj, "results", "ExternalValidation", "04_data",
                                  "external_predictions_all_versions.csv"))
pe <- pe[cohort == "GSE116250" & version == "V2"]
## match by title (the external matrix columns are titles)
idx <- match(colnames(e2)[keep][ok], pe$title)
y4 <- pe$y[idx]; p4 <- pe$prob[idx]
good <- !is.na(p4) & !is.na(y4)
r4e <- roc(y4[good], p4[good], quiet = TRUE, direction = "<")
dle <- roc.test(roc(ye[ok][good], lp_e[ok][good], quiet = TRUE, direction = "<"),
                r4e, method = "delong")
log_msg(sprintf("  four-gene external AUC (same samples) = %.4f ; paired DeLong p = %.3f",
                as.numeric(auc(r4e)), dle$p.value))

ext_tab <- data.table(
  cohort = "GSE116250",
  model = c("Eight-gene, selected in training only", "Four-gene, published"),
  n = c(sum(ok), sum(good)),
  n_case = c(sum(ye[ok] == 1), sum(y4[good] == 1)),
  n_control = c(sum(ye[ok] == 0), sum(y4[good] == 0)),
  auc = round(c(as.numeric(auc(r_e)), as.numeric(auc(r4e))), 4),
  ci_low = round(c(ci_e[1], as.numeric(ci.auc(r4e))[1]), 4),
  ci_high = round(c(ci_e[3], as.numeric(ci.auc(r4e))[3]), 4),
  delong_p_vs_four_gene = c(signif(dle$p.value, 4), NA)
)
write_csv(ext_tab, file.path(d_tab, "T16_trainonly_panel_external.csv"))
log_msg("  GSE5406 not evaluable: CREB5 and SMOC2 are absent from GPL96 after probe collapsing")

log_msg("")
log_msg("Done.")
writeLines(run_log, file.path(d_log, "137_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_137.txt"))
