###############################################################################
## 71_combat_reanalysis.R
## GSE57338 -- batch correction (ComBat) and re-analysis, before/after contrast
##
## Conditions compared
##   C0  original            : limma ~ group + sex + age           (as published)
##   C1  ComBat(order block) : ComBat batch = GSM-order block B01-B16, mod = ~group+sex+age
##   C2  ComBat(ID block)    : ComBat batch = patient-ID quartile within prefix
##   C3  batch as covariate  : limma ~ group + sex + age + order block
##
## For each condition the script recomputes
##   - DE gene counts (up / down) at |log2FC| > 0.585 & padj < 0.05
##   - overlap with the original DEG list (Jaccard), log2FC Spearman correlation
##   - direction of the 70 WGCNA-DEG candidate genes and of the 7 PPI hubs
##   - M2 (4-gene) logistic model: coefficients, training / internal-validation AUC
##     with DeLong 95% CI, using the ORIGINAL train/validation split
##   - two cross-validation schemes: random stratified 5-fold (naive) and
##     batch-aware 5-fold (folds defined by GSM-order block)
##   - PCA before/after, group and batch eta-squared on PC1
##
## Inputs
##   clean_data/GSE57338_clean_expression_matrix.csv
##   results/Batch/01_tables/T01_sample_batch_metadata.csv
##   results/Nomogram/04_data/modelpredictions_{train,validation}.csv  (split definition)
##   results/HubGene/GSE57338_candidate_genes.csv
##   results/PPI/05_node_degree_table.csv
##
## Outputs (results/Batch/)
##   01_tables/T08_deg_summary_before_after.csv
##   01_tables/T09_key_gene_log2fc_before_after.csv
##   01_tables/T10_candidate70_direction_before_after.csv
##   01_tables/T11_model_auc_before_after.csv
##   01_tables/T12_cv_scheme_comparison.csv
##   01_tables/T13_pca_variance_before_after.csv
##   01_tables/T14_comparison_table_main.csv     <- the before/after table
##   02_figures/F04_pca_before_after.{pdf,png}
##   02_figures/F05_roc_before_after.{pdf,png}
##   02_figures/F06_key_gene_log2fc_shift.{pdf,png}
##   02_figures/F07_four_genes_boxplot_before_after.{pdf,png}
##   03_logs/71_run_log.txt , sessionInfo_71.txt
##
## Run:  Rscript scripts/71_combat_reanalysis.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
  library(pROC)
  library(rms)
})

options(stringsAsFactors = FALSE)
set.seed(20261029)

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
clean_dir <- file.path(proj_dir, "clean_data")
res_dir   <- file.path(proj_dir, "results")
out_root  <- file.path(res_dir, "Batch")
d_tab     <- file.path(out_root, "01_tables")
d_fig     <- file.path(out_root, "02_figures")
d_log     <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) { txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = "") }
hr <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

LFC_CUT  <- 0.585
PADJ_CUT <- 0.05
GENES    <- c("AQP3", "GGT5", "ISLR", "LAD1")

## ===========================================================================
## 1. Load and align
## ===========================================================================
log_msg("[1] Loading expression matrix and batch metadata ...")

expr_dt <- data.table::fread(file.path(clean_dir, "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

bmeta <- data.table::fread(file.path(d_tab, "T01_sample_batch_metadata.csv"))
bmeta <- bmeta[match(colnames(expr), sample)]
if (anyNA(bmeta$sample)) stop("batch metadata does not cover all matrix columns")
bmeta[, group2 := factor(group2, levels = c("Control", "Disease"))]
bmeta[, sex := factor(sex, levels = c("Female", "Male"))]
bmeta[, HF := as.integer(group2 == "Disease")]
bmeta[, batch_order := factor(batch_order)]
bmeta[, batch_idnum := factor(batch_idnum)]
bmeta[, batch_prefix := factor(batch_prefix)]
## age: 0 is a recording error in adult hearts; impute with the median (same rule as 04_deg.R)
bmeta[!is.na(age_num) & age_num == 0, age_num := NA_real_]
if (anyNA(bmeta$age_num)) {
  med_age <- median(bmeta$age_num, na.rm = TRUE)
  log_msg("  age missing for ", sum(is.na(bmeta$age_num)),
          " sample(s); imputed with the median ", med_age)
  bmeta[is.na(age_num), age_num := med_age]
}
log_msg("  ", nrow(expr), " genes x ", ncol(expr), " samples; blocks: ",
        nlevels(bmeta$batch_order))

## drop zero-variance genes (ComBat cannot handle them)
gv <- apply(expr, 1, var)
keep <- is.finite(gv) & gv > 0
if (sum(!keep) > 0) log_msg("  dropped ", sum(!keep), " zero-variance genes")
expr <- expr[keep, , drop = FALSE]

## ===========================================================================
## 2. Batch correction
## ===========================================================================
log_msg("")
log_msg("[2] ComBat ...")

mod <- model.matrix(~ group2 + sex + age_num, data = bmeta)
combat <- function(batchvar, tag) {
  x <- tryCatch(
    sva::ComBat(dat = expr, batch = bmeta[[batchvar]], mod = mod,
                par.prior = TRUE, prior.plots = FALSE),
    error = function(e) { log_msg("  [", tag, "] ComBat failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(x)) log_msg("  [", tag, "] done; range ", sprintf("%.2f-%.2f", min(x), max(x)))
  x
}

expr_c1 <- combat("batch_order", "C1 order-block")
expr_c2 <- combat("batch_idnum", "C2 ID-block")
if (is.null(expr_c1) || is.null(expr_c2)) stop("ComBat failed")

## ===========================================================================
## 3. Differential expression under each condition
## ===========================================================================
log_msg("")
log_msg("[3] limma under each condition ...")

fit_limma <- function(mat, form, coef_name) {
  design <- model.matrix(form, data = bmeta)
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
  data.table(GeneSymbol = rownames(tt), log2FoldChange = tt$logFC,
             AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value, padj = tt$adj.P.Val)
}

f_main <- ~ group2 + sex + age_num
f_bcov <- ~ group2 + sex + age_num + batch_order

res <- list()
res$C0_original <- fit_limma(expr, f_main, "group2Disease")
res$C1_combat_order <- fit_limma(expr_c1, f_main, "group2Disease")
res$C2_combat_idblock <- fit_limma(expr_c2, f_main, "group2Disease")
res$C3_batch_covariate <- fit_limma(expr, f_bcov, "group2Disease")

for (nm in names(res)) {
  d <- res[[nm]]
  ## round log2FC to 4 decimals BEFORE flagging, exactly as 04_deg.R does
  d[, log2FoldChange := round(log2FoldChange, 4)]
  d[, sig := !is.na(padj) & padj < PADJ_CUT & abs(log2FoldChange) > LFC_CUT]
  d[, direction := ifelse(!sig, "NS", ifelse(log2FoldChange > 0, "Up", "Down"))]
  res[[nm]] <- d
  log_msg(sprintf("  %-22s sig=%4d (Up %3d / Down %3d)  padj<0.05: %5d  |log2FC|>1 & padj<0.05: %3d",
                  nm, sum(d$sig), sum(d$direction == "Up"), sum(d$direction == "Down"),
                  sum(d$padj < PADJ_CUT, na.rm = TRUE),
                  sum(abs(d$log2FoldChange) > 1 & d$padj < PADJ_CUT, na.rm = TRUE)))
}

ref <- res$C0_original
ref_sig <- ref[sig == TRUE, GeneSymbol]

deg_summary <- rbindlist(lapply(names(res), function(nm) {
  d <- res[[nm]]
  s <- d[sig == TRUE]
  common <- intersect(d$GeneSymbol, ref$GeneSymbol)
  a <- d[match(common, GeneSymbol), log2FoldChange]
  b <- ref[match(common, GeneSymbol), log2FoldChange]
  inter <- length(intersect(s$GeneSymbol, ref_sig))
  uni   <- length(union(s$GeneSymbol, ref_sig))
  new_g <- setdiff(s$GeneSymbol, ref_sig)
  data.table(
    condition = nm,
    n_sig = nrow(s),
    n_up = sum(s$direction == "Up"),
    n_down = sum(s$direction == "Down"),
    n_padj005 = sum(d$padj < PADJ_CUT, na.rm = TRUE),
    n_absfc1 = sum(abs(d$log2FoldChange) > 1 & d$padj < PADJ_CUT, na.rm = TRUE),
    overlap_with_original = inter,
    jaccard_vs_original = round(inter / max(1, uni), 3),
    pct_original_retained = round(100 * inter / length(ref_sig), 1),
    n_new_vs_original = length(new_g),
    log2fc_spearman_vs_original = round(cor(a, b, method = "spearman"), 4),
    log2fc_sign_concordance = round(mean(sign(a) == sign(b)), 4),
    median_abs_log2fc_sig = round(median(abs(s$log2FoldChange)), 3),
    attenuation_ratio = round(median(abs(s$log2FoldChange)) /
                              median(abs(ref[sig == TRUE, log2FoldChange])), 3)
  )
}))
write_csv(deg_summary, file.path(d_tab, "T08_deg_summary_before_after.csv"))
log_msg("")
for (i in seq_len(nrow(deg_summary))) {
  log_msg(sprintf("  %-22s sig=%4d  overlap=%3d (%.1f%% of original)  new=%3d  rho(logFC)=%.3f  attenuation=%.2f",
                  deg_summary$condition[i], deg_summary$n_sig[i],
                  deg_summary$overlap_with_original[i], deg_summary$pct_original_retained[i],
                  deg_summary$n_new_vs_original[i],
                  deg_summary$log2fc_spearman_vs_original[i],
                  deg_summary$attenuation_ratio[i]))
}

## ===========================================================================
## 4. Key genes: 4 model genes + 7 PPI hubs
## ===========================================================================
log_msg("")
log_msg("[4] key genes ...")

ppi <- data.table::fread(file.path(res_dir, "PPI", "05_node_degree_table.csv"))
hubs <- head(ppi[order(-Degree), GeneSymbol], 7)
key_genes <- unique(c(GENES, hubs))
log_msg("  key genes: ", paste(key_genes, collapse = ", "))

vals <- rbindlist(lapply(names(res), function(nm) {
  d <- res[[nm]][match(key_genes, GeneSymbol)]
  data.table(condition = nm, gene = key_genes,
             log2FoldChange = round(d$log2FoldChange, 3),
             padj = signif(d$padj, 4),
             significant = d$sig,
             direction = d$direction)
}))
## long-format key gene table
key_view <- vals[, .(gene, condition, log2FoldChange, padj, significant)]
write_csv(key_view, file.path(d_tab, "T09_key_gene_log2fc_long.csv"))

key_wide <- merge(
  vals[condition == "C0_original", .(gene, log2FC_before = log2FoldChange,
                                     padj_before = padj, sig_before = significant)],
  vals[condition == "C1_combat_order", .(gene, log2FC_after_order = log2FoldChange,
                                         padj_after_order = padj, sig_after_order = significant)],
  by = "gene")
key_wide <- merge(key_wide,
  vals[condition == "C2_combat_idblock", .(gene, log2FC_after_id = log2FoldChange,
                                           padj_after_id = padj)],
  by = "gene", all.x = TRUE)
key_wide <- merge(key_wide,
  vals[condition == "C3_batch_covariate", .(gene, log2FC_covariate = log2FoldChange,
                                            padj_covariate = padj)],
  by = "gene", all.x = TRUE)
key_wide[, direction_change := sign(log2FC_before) != sign(log2FC_after_order)]
key_wide[, role := ifelse(gene %in% GENES, "M2 model gene", "PPI hub")]
write_csv(key_wide, file.path(d_tab, "T09_key_gene_log2fc_before_after.csv"))
print(key_wide[, .(gene, role, log2FC_before, log2FC_after_order, padj_before, padj_after_order)])

## ===========================================================================
## 5. The 70 WGCNA-DEG candidate genes: direction stability
## ===========================================================================
log_msg("")
log_msg("[5] 70 candidate genes ...")

cand <- data.table::fread(file.path(res_dir, "HubGene", "GSE57338_candidate_genes.csv"))
cand_genes <- cand$GeneSymbol

cand_dir <- rbindlist(lapply(names(res), function(nm) {
  d <- res[[nm]][match(cand_genes, GeneSymbol)]
  data.table(condition = nm,
             n_tested = sum(!is.na(d$GeneSymbol)),
             n_still_sig = sum(d$sig, na.rm = TRUE),
             n_same_direction = sum(sign(d$log2FoldChange) == sign(cand$log2FoldChange),
                                    na.rm = TRUE),
             n_flipped = sum(sign(d$log2FoldChange) != sign(cand$log2FoldChange),
                             na.rm = TRUE),
             median_abs_log2fc = round(median(abs(d$log2FoldChange), na.rm = TRUE), 3))
}))
write_csv(cand_dir, file.path(d_tab, "T10_candidate70_direction_before_after.csv"))
for (i in seq_len(nrow(cand_dir))) {
  log_msg(sprintf("  %-22s still significant: %2d/70   same direction: %2d   flipped: %d   median|logFC|=%.2f",
                  cand_dir$condition[i], cand_dir$n_still_sig[i],
                  cand_dir$n_same_direction[i], cand_dir$n_flipped[i],
                  cand_dir$median_abs_log2fc[i]))
}

## ===========================================================================
## 6. The M2 four-gene model, re-fitted under each condition
## ===========================================================================
log_msg("")
log_msg("[6] models, original train/validation split ...")

tr_file <- file.path(res_dir, "Nomogram", "04_data", "modelpredictions_train.csv")
va_file <- file.path(res_dir, "Nomogram", "04_data", "modelpredictions_validation.csv")
if (!file.exists(tr_file)) stop("original split file not found: ", tr_file)
train_samples <- data.table::fread(tr_file)$sample
valid_samples <- data.table::fread(va_file)$sample
log_msg("  original split reused: train ", length(train_samples),
        " / validation ", length(valid_samples))
stopifnot(length(intersect(train_samples, valid_samples)) == 0)

## common data frame builder: genes + covariates, split preserved
make_df <- function(mat, with_batch = FALSE) {
  g <- t(mat[GENES, , drop = FALSE]); colnames(g) <- GENES
  df <- data.table(sample = bmeta$sample, HF = bmeta$HF,
                   as.data.table(g),
                   block = bmeta$batch_order,
                   set = ifelse(bmeta$sample %in% train_samples, "train", "validation"))
  stopifnot(!anyNA(df$HF))
  df
}

## fit on train, evaluate on validation; z parameters from train only
run_model <- function(df, rhs, tag, ztransform = TRUE) {
  tr <- df[set == "train"]; va <- df[set == "validation"]
  if (ztransform) {
    for (x in GENES) {
      m <- mean(tr[[x]]); s <- sd(tr[[x]])
      tr[[paste0(x, "_z")]] <- (tr[[x]] - m) / s
      va[[paste0(x, "_z")]] <- (va[[x]] - m) / s
    }
  }
  fit <- suppressWarnings(glm(as.formula(paste("HF ~", rhs)), data = tr, family = binomial))
  pr_tr <- predict(fit, type = "response")
  pr_va <- predict(fit, newdata = va, type = "response")
  a_tr <- pROC::roc(tr$HF, pr_tr, quiet = TRUE)
  a_va <- pROC::roc(va$HF, pr_va, quiet = TRUE)
  ci_tr <- as.numeric(pROC::ci.auc(a_tr)); ci_va <- as.numeric(pROC::ci.auc(a_va))
  list(tag = tag, coef = coef(fit),
       auc_train = as.numeric(a_tr$auc), train_l = ci_tr[1], train_u = ci_tr[3],
       auc_valid = as.numeric(a_va$auc), valid_l = ci_va[1], valid_u = ci_va[3],
       pr_valid = pr_va, y_valid = va$HF, pr_train = pr_tr, y_train = tr$HF,
       brier_valid = mean((pr_va - va$HF)^2))
}

df0 <- make_df(expr)
df1 <- make_df(expr_c1)
df2 <- make_df(expr_c2)

m2 <- list(
  C0_original       = run_model(df0, "AQP3_z + GGT5_z + ISLR_z + LAD1_z", "C0_original"),
  C1_combat_order   = run_model(df1, "AQP3_z + GGT5_z + ISLR_z + LAD1_z", "C1_combat_order"),
  C2_combat_idblock = run_model(df2, "AQP3_z + GGT5_z + ISLR_z + LAD1_z", "C2_combat_idblock"),
  C3_genes_plus_batch = run_model(df0, "AQP3_z + GGT5_z + ISLR_z + LAD1_z + block",
                                  "C3_genes_plus_batch"),
  C4_batch_only     = run_model(df0, "block", "C4_batch_only", ztransform = FALSE),
  C5_combat_plus_batch = run_model(df1, "AQP3_z + GGT5_z + ISLR_z + LAD1_z + block",
                                   "C5_combat_plus_batch")
)
log_msg("  models: 4-gene (C0-C2), 4-gene + batch block (C3, C5), batch block only (C4)")

## frozen original coefficients applied to ComBat-corrected validation data
frozen <- local({
  zp <- data.table(gene = GENES,
                   m = sapply(GENES, function(x) mean(df0[set == "train"][[x]])),
                   s = sapply(GENES, function(x) sd(df0[set == "train"][[x]])))
  b <- m2$C0_original$coef
  va <- df1[set == "validation"]
  z <- sapply(GENES, function(x) (va[[x]] - zp[gene == x, m]) / zp[gene == x, s])
  lp <- b["(Intercept)"] + as.numeric(z %*% b[paste0(GENES, "_z")])
  pr <- 1 / (1 + exp(-lp))
  a <- pROC::roc(va$HF, pr, quiet = TRUE); ci <- as.numeric(pROC::ci.auc(a))
  list(tag = "C1_combat_order_frozen_model", auc = as.numeric(a$auc),
       l = ci[1], u = ci[3], pr = pr, y = va$HF)
})
log_msg(sprintf("  frozen-model AUC on ComBat-corrected validation set = %.4f (%.4f-%.4f)",
                frozen$auc, frozen$l, frozen$u))

auc_tab <- rbindlist(lapply(m2, function(x) {
  data.table(condition = x$tag,
             auc_train = round(x$auc_train, 4), train_l = round(x$train_l, 4),
             train_u = round(x$train_u, 4),
             auc_valid = round(x$auc_valid, 4), valid_l = round(x$valid_l, 4),
             valid_u = round(x$valid_u, 4),
             brier_valid = round(x$brier_valid, 4),
             coef_AQP3 = round(x$coef["AQP3_z"], 4),
             coef_GGT5 = round(x$coef["GGT5_z"], 4),
             coef_ISLR = round(x$coef["ISLR_z"], 4),
             coef_LAD1 = round(x$coef["LAD1_z"], 4))
}))
auc_tab <- rbind(auc_tab, data.table(
  condition = frozen$tag, auc_train = NA_real_, train_l = NA_real_, train_u = NA_real_,
  auc_valid = round(frozen$auc, 4), valid_l = round(frozen$l, 4), valid_u = round(frozen$u, 4),
  brier_valid = round(mean((frozen$pr - frozen$y)^2), 4),
  coef_AQP3 = NA_real_, coef_GGT5 = NA_real_, coef_ISLR = NA_real_, coef_LAD1 = NA_real_))
write_csv(auc_tab, file.path(d_tab, "T11_model_auc_before_after.csv"))
for (i in seq_len(nrow(auc_tab))) {
  log_msg(sprintf("  %-32s train AUC=%.4f  valid AUC=%.4f (%.3f-%.3f)  Brier=%.4f",
                  auc_tab$condition[i], auc_tab$auc_train[i], auc_tab$auc_valid[i],
                  auc_tab$valid_l[i], auc_tab$valid_u[i], auc_tab$brier_valid[i]))
}
## 7. Cross-validation schemes: naive vs batch-aware
## ===========================================================================
log_msg("")
log_msg("[7] naive vs batch-aware cross-validation ...")

blocks <- as.character(bmeta$batch_order)

assign_folds <- function(blocks, K = 5) {
  tab <- table(blocks)
  ord <- names(tab)[order(-as.integer(tab))]
  fold_of_block <- integer(length(ord)); names(fold_of_block) <- ord
  load <- numeric(K)
  for (b in ord) {
    k <- which.min(load)
    fold_of_block[b] <- k
    load[k] <- load[k] + as.integer(tab[b])
  }
  unname(fold_of_block[blocks])
}

cv_auc_for <- function(mat, folds, stratified = FALSE, reps = 1) {
  g <- t(mat[GENES, , drop = FALSE]); colnames(g) <- GENES
  df <- data.table(sample = bmeta$sample, HF = bmeta$HF, as.data.table(g))
  out <- numeric(0)
  for (r in seq_len(reps)) {
    fold <- folds
    if (stratified) {
      fold <- integer(nrow(df))
      for (lv in 0:1) {
        idx <- which(df$HF == lv)
        set.seed(2000 + r)
        fold[idx] <- sample(rep(1:5, length.out = length(idx)))
      }
    }
    pr <- rep(NA_real_, nrow(df))
    for (k in sort(unique(fold))) {
      tr <- fold != k; te <- fold == k
      zp <- sapply(GENES, function(x) c(mean(df[[x]][tr]), sd(df[[x]][tr])))
      ztr <- sapply(GENES, function(x) (df[[x]][tr] - zp[1, x]) / zp[2, x])
      zte <- sapply(GENES, function(x) (df[[x]][te] - zp[1, x]) / zp[2, x])
      dtr <- data.table(HF = df$HF[tr], as.data.table(ztr))
      fit <- glm(HF ~ ., data = dtr, family = binomial)
      pr[te] <- predict(fit, newdata = data.table(as.data.table(zte)), type = "response")
    }
    out <- c(out, as.numeric(pROC::auc(pROC::roc(df$HF, pr, quiet = TRUE))))
  }
  out
}

bm_folds <- assign_folds(blocks, 5)
log_msg("  batch-aware folds: ", paste(table(bm_folds), collapse = " / "),
        " samples per fold; blocks per fold: ",
        paste(tapply(blocks, bm_folds, function(x) length(unique(x))), collapse = " / "))

cv_tab <- rbindlist(list(
  data.table(condition = "C0_original", scheme = "random stratified 5-fold CV (20 reps)",
             auc_mean = round(mean(cv_auc_for(expr, NULL, TRUE, 20)), 4),
             auc_sd = round(sd(cv_auc_for(expr, NULL, TRUE, 20)), 4)),
  data.table(condition = "C0_original", scheme = "batch-aware 5-fold CV (GSM-order block)",
             auc_mean = round(mean(cv_auc_for(expr, bm_folds, FALSE, 1)), 4),
             auc_sd = NA_real_),
  data.table(condition = "C1_combat_order", scheme = "random stratified 5-fold CV (20 reps)",
             auc_mean = round(mean(cv_auc_for(expr_c1, NULL, TRUE, 20)), 4),
             auc_sd = round(sd(cv_auc_for(expr_c1, NULL, TRUE, 20)), 4)),
  data.table(condition = "C1_combat_order", scheme = "batch-aware 5-fold CV (GSM-order block)",
             auc_mean = round(mean(cv_auc_for(expr_c1, bm_folds, FALSE, 1)), 4),
             auc_sd = NA_real_)
))
cv_tab[, optimism := round(max(auc_mean[grepl("random", scheme)]) -
                             auc_mean[grepl("batch-aware", scheme)], 4), by = condition]
write_csv(cv_tab, file.path(d_tab, "T12_cv_scheme_comparison.csv"))
for (i in seq_len(nrow(cv_tab))) log_msg(sprintf("  %-18s %-45s AUC=%.4f", cv_tab$condition[i],
                                                 cv_tab$scheme[i], cv_tab$auc_mean[i]))

## ===========================================================================
## 8. PCA / variance decomposition before and after
## ===========================================================================
log_msg("")
log_msg("[8] PCA before/after ...")

eta2 <- function(x, g) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- factor(g[ok])
  ss <- summary(aov(x ~ g))[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
}
gv <- apply(expr, 1, var)
topg <- names(sort(gv, decreasing = TRUE))[1:2000]

pca_cmp <- rbindlist(lapply(c("C0_original", "C1_combat_order"), function(nm) {
  m <- if (nm == "C0_original") expr else expr_c1
  p <- prcomp(t(m[topg, , drop = FALSE]), center = TRUE, scale. = FALSE)
  ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  data.table(condition = nm,
             pc1_var_pct = round(ve[1], 2), pc2_var_pct = round(ve[2], 2),
             pc1_eta2_group = round(eta2(p$x[, 1], bmeta$group2), 4),
             pc1_eta2_block = round(eta2(p$x[, 1], bmeta$batch_order), 4),
             pc2_eta2_group = round(eta2(p$x[, 2], bmeta$group2), 4),
             pc2_eta2_block = round(eta2(p$x[, 2], bmeta$batch_order), 4),
             mean_pergene_eta2_group = round(mean(apply(m[topg, ], 1, function(y) eta2(y, bmeta$group2))), 4),
             mean_pergene_eta2_block = round(mean(apply(m[topg, ], 1, function(y) eta2(y, bmeta$batch_order))), 4))
}))
write_csv(pca_cmp, file.path(d_tab, "T13_pca_variance_before_after.csv"))
for (i in seq_len(nrow(pca_cmp))) {
  log_msg(sprintf("  %-18s PC1 %.1f%%  eta2(group)=%.3f  eta2(block)=%.3f  | mean per-gene eta2 group=%.3f block=%.3f",
                  pca_cmp$condition[i], pca_cmp$pc1_var_pct[i], pca_cmp$pc1_eta2_group[i],
                  pca_cmp$pc1_eta2_block[i], pca_cmp$mean_pergene_eta2_group[i],
                  pca_cmp$mean_pergene_eta2_block[i]))
}

## ===========================================================================
## 9. The before/after comparison table
## ===========================================================================
log_msg("")
log_msg("[9] assembling the comparison table ...")

gv_row <- function(df, cond, col) df[condition == cond, get(col)]
mk <- function(metric, before, c1, c2, c3, note = "") {
  data.table(metric = metric, before = before, after_combat_order = c1,
             after_combat_idblock = c2, batch_as_covariate = c3, note = note)
}
cmp <- rbindlist(list(
  mk("Significant DE genes (n)", gv_row(deg_summary, "C0_original", "n_sig"),
     gv_row(deg_summary, "C1_combat_order", "n_sig"),
     gv_row(deg_summary, "C2_combat_idblock", "n_sig"),
     gv_row(deg_summary, "C3_batch_covariate", "n_sig"),
     "|log2FC|>0.585 & padj<0.05"),
  mk("  up-regulated (n)", gv_row(deg_summary, "C0_original", "n_up"),
     gv_row(deg_summary, "C1_combat_order", "n_up"),
     gv_row(deg_summary, "C2_combat_idblock", "n_up"),
     gv_row(deg_summary, "C3_batch_covariate", "n_up"), ""),
  mk("  down-regulated (n)", gv_row(deg_summary, "C0_original", "n_down"),
     gv_row(deg_summary, "C1_combat_order", "n_down"),
     gv_row(deg_summary, "C2_combat_idblock", "n_down"),
     gv_row(deg_summary, "C3_batch_covariate", "n_down"), ""),
  mk("Genes with padj<0.05 (n)", gv_row(deg_summary, "C0_original", "n_padj005"),
     gv_row(deg_summary, "C1_combat_order", "n_padj005"),
     gv_row(deg_summary, "C2_combat_idblock", "n_padj005"),
     gv_row(deg_summary, "C3_batch_covariate", "n_padj005"), "no log2FC filter"),
  mk("Original DEGs retained (%)", 100,
     gv_row(deg_summary, "C1_combat_order", "pct_original_retained"),
     gv_row(deg_summary, "C2_combat_idblock", "pct_original_retained"),
     gv_row(deg_summary, "C3_batch_covariate", "pct_original_retained"), "of 269"),
  mk("New DEGs not in original list (n)", 0,
     gv_row(deg_summary, "C1_combat_order", "n_new_vs_original"),
     gv_row(deg_summary, "C2_combat_idblock", "n_new_vs_original"),
     gv_row(deg_summary, "C3_batch_covariate", "n_new_vs_original"), "over-correction signal"),
  mk("log2FC Spearman vs original", 1,
     gv_row(deg_summary, "C1_combat_order", "log2fc_spearman_vs_original"),
     gv_row(deg_summary, "C2_combat_idblock", "log2fc_spearman_vs_original"),
     gv_row(deg_summary, "C3_batch_covariate", "log2fc_spearman_vs_original"), "all genes"),
  mk("log2FC sign concordance", 1,
     gv_row(deg_summary, "C1_combat_order", "log2fc_sign_concordance"),
     gv_row(deg_summary, "C2_combat_idblock", "log2fc_sign_concordance"),
     gv_row(deg_summary, "C3_batch_covariate", "log2fc_sign_concordance"), "all genes"),
  mk("Median |log2FC| of DEGs", gv_row(deg_summary, "C0_original", "median_abs_log2fc_sig"),
     gv_row(deg_summary, "C1_combat_order", "median_abs_log2fc_sig"),
     gv_row(deg_summary, "C2_combat_idblock", "median_abs_log2fc_sig"),
     gv_row(deg_summary, "C3_batch_covariate", "median_abs_log2fc_sig"),
     "attenuation ratio in T08"),
  mk("70 candidate genes still significant", gv_row(cand_dir, "C0_original", "n_still_sig"),
     gv_row(cand_dir, "C1_combat_order", "n_still_sig"),
     gv_row(cand_dir, "C2_combat_idblock", "n_still_sig"),
     gv_row(cand_dir, "C3_batch_covariate", "n_still_sig"), "of 70"),
  mk("70 candidate genes same direction", gv_row(cand_dir, "C0_original", "n_same_direction"),
     gv_row(cand_dir, "C1_combat_order", "n_same_direction"),
     gv_row(cand_dir, "C2_combat_idblock", "n_same_direction"),
     gv_row(cand_dir, "C3_batch_covariate", "n_same_direction"), "of 70"),
  mk("M2 training AUC", round(m2$C0_original$auc_train, 4),
     round(m2$C1_combat_order$auc_train, 4),
     round(m2$C2_combat_idblock$auc_train, 4),
     round(m2$C3_genes_plus_batch$auc_train, 4), "4-gene logistic, n=219"),
  mk("M2 internal-validation AUC", round(m2$C0_original$auc_valid, 4),
     round(m2$C1_combat_order$auc_valid, 4),
     round(m2$C2_combat_idblock$auc_valid, 4),
     round(m2$C3_genes_plus_batch$auc_valid, 4), "n=94"),
  mk("M2 validation AUC 95% CI lower", round(m2$C0_original$valid_l, 4),
     round(m2$C1_combat_order$valid_l, 4),
     round(m2$C2_combat_idblock$valid_l, 4),
     round(m2$C3_genes_plus_batch$valid_l, 4), "DeLong"),
  mk("M2 validation AUC 95% CI upper", round(m2$C0_original$valid_u, 4),
     round(m2$C1_combat_order$valid_u, 4),
     round(m2$C2_combat_idblock$valid_u, 4),
     round(m2$C3_genes_plus_batch$valid_u, 4), "DeLong"),
  mk("M2 validation Brier score", round(m2$C0_original$brier_valid, 4),
     round(m2$C1_combat_order$brier_valid, 4),
     round(m2$C2_combat_idblock$brier_valid, 4),
     round(m2$C3_genes_plus_batch$brier_valid, 4), "lower is better"),
  mk("Validation AUC: batch block alone", round(m2$C4_batch_only$auc_valid, 4), NA, NA, NA,
     "16-level block as the only predictor"),
  mk("Validation AUC: 4 genes + batch block", round(m2$C3_genes_plus_batch$auc_valid, 4),
     round(m2$C5_combat_plus_batch$auc_valid, 4), NA, NA,
     "how much the 4 genes add on top of batch"),
  mk("Validation AUC: frozen model on ComBat data", NA, round(frozen$auc, 4), NA, NA,
     "original coefficients, no refitting"),
  mk("Random stratified CV AUC", cv_tab[condition == "C0_original" & grepl("random", scheme), auc_mean],
     cv_tab[condition == "C1_combat_order" & grepl("random", scheme), auc_mean],
     NA, NA, "20 repeats, 5-fold, z and coefficients refitted inside folds"),
  mk("Batch-aware CV AUC", cv_tab[condition == "C0_original" & grepl("batch-aware", scheme), auc_mean],
     cv_tab[condition == "C1_combat_order" & grepl("batch-aware", scheme), auc_mean],
     NA, NA, "folds defined by GSM-order block"),
  mk("PC1 eta2: group", pca_cmp[condition == "C0_original", pc1_eta2_group],
     pca_cmp[condition == "C1_combat_order", pc1_eta2_group], NA, NA, "share of PC1 variance"),
  mk("PC1 eta2: batch block", pca_cmp[condition == "C0_original", pc1_eta2_block],
     pca_cmp[condition == "C1_combat_order", pc1_eta2_block], NA, NA, "share of PC1 variance")
))
write_csv(cmp, file.path(d_tab, "T14_comparison_table_main.csv"))
log_msg("  comparison table written: T14_comparison_table_main.csv")

## ===========================================================================
## 10. Figures
## ===========================================================================
log_msg("")
log_msg("[10] figures ...")

save_both <- function(stem, expr_code, width = 8, height = 6) {
  e <- substitute(expr_code); env <- parent.frame()
  grDevices::pdf(file.path(d_fig, paste0(stem, ".pdf")), width = width, height = height)
  eval(e, env); grDevices::dev.off()
  grDevices::png(file.path(d_fig, paste0(stem, ".png")), width = width, height = height,
                 units = "in", res = 300)
  eval(e, env); grDevices::dev.off()
}
grp_col <- c(Control = "#2C7FB8", Disease = "#D7301F")

pca_of <- function(m) prcomp(t(m[topg, , drop = FALSE]), center = TRUE, scale. = FALSE)
p0 <- pca_of(expr); p1 <- pca_of(expr_c1)
ve0 <- 100 * p0$sdev^2 / sum(p0$sdev^2); ve1 <- 100 * p1$sdev^2 / sum(p1$sdev^2)

save_both("F04_pca_before_after", {
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1), cex = 0.85)
  plot(p0$x[, 1], p0$x[, 2], col = grp_col[as.character(bmeta$group2)], pch = 19, cex = 0.75,
       xlab = sprintf("PC1 (%.1f%%)", ve0[1]), ylab = sprintf("PC2 (%.1f%%)", ve0[2]),
       main = "Before: raw RMA expression")
  legend("topright", legend = names(grp_col), col = grp_col, pch = 19, bty = "n", cex = 0.8)
  plot(p1$x[, 1], p1$x[, 2], col = grp_col[as.character(bmeta$group2)], pch = 19, cex = 0.75,
       xlab = sprintf("PC1 (%.1f%%)", ve1[1]), ylab = sprintf("PC2 (%.1f%%)", ve1[2]),
       main = "After: ComBat (GSM-order block)")
  legend("topright", legend = names(grp_col), col = grp_col, pch = 19, bty = "n", cex = 0.8)
  par(op)
}, width = 11, height = 5)

save_both("F05_roc_before_after", {
  op <- par(mar = c(4.2, 4.2, 3, 1))
  cols <- c("C0 original" = "#333333", "C1 ComBat(order)" = "#D7301F",
            "C2 ComBat(ID)" = "#EF8A62", "C1 frozen model" = "#2C7FB8")
  plot(1, type = "n", xlim = c(1, 0), ylim = c(0, 1), xlab = "Specificity",
       ylab = "Sensitivity", main = "M2 gene model: internal validation set (n = 94)")
  abline(0, 1, lty = 2, col = "grey70")
  add <- function(y, pr, lab, col) {
    r <- pROC::roc(y, pr, quiet = TRUE)
    lines(1 - r$specificities, r$sensitivities, col = col, lwd = 2)
  }
  add(m2$C0_original$y_valid, m2$C0_original$pr_valid,
      "C0 original", cols[["C0 original"]])
  add(m2$C1_combat_order$y_valid, m2$C1_combat_order$pr_valid,
      "C1 ComBat", cols[["C1 ComBat(order)"]])
  add(m2$C2_combat_idblock$y_valid, m2$C2_combat_idblock$pr_valid,
      "C2 ComBat ID", cols[["C2 ComBat(ID)"]])
  add(frozen$y, frozen$pr, "frozen", cols[["C1 frozen model"]])
  legend("bottomright", bty = "n", lwd = 2, cex = 0.85, col = unname(cols),
         legend = sprintf("%s (AUC %.3f)", names(cols),
                          c(m2$C0_original$auc_valid, m2$C1_combat_order$auc_valid,
                            m2$C2_combat_idblock$auc_valid, frozen$auc)))
  par(op)
}, width = 6.2, height = 6)

save_both("F06_key_gene_log2fc_shift", {
  kd <- key_wide[order(log2FC_before)]
  op <- par(mar = c(4.2, 8, 3, 1))
  y <- seq_len(nrow(kd))
  plot(kd$log2FC_before, y, pch = 19, col = "#333333", cex = 0.9,
       xlim = range(c(kd$log2FC_before, kd$log2FC_after_order, 0)),
       yaxt = "n", ylab = "", xlab = "log2 fold change (disease vs control)",
       main = "Key genes: before vs after ComBat")
  axis(2, at = y, labels = kd$gene, las = 1, cex.axis = 0.8)
  segments(kd$log2FC_before, y, kd$log2FC_after_order, y, col = "grey60")
  points(kd$log2FC_after_order, y, pch = 17, col = "#D7301F", cex = 0.9)
  abline(v = 0, lty = 3)
  legend("bottomright", pch = c(19, 17), col = c("#333333", "#D7301F"),
         legend = c("before (raw)", "after ComBat"), bty = "n", cex = 0.8)
  par(op)
}, width = 7, height = 6)

save_both("F07_four_genes_boxplot_before_after", {
  op <- par(mfrow = c(2, 4), mar = c(2.5, 3.6, 2.2, 0.6), cex = 0.8)
  for (gg in GENES) {
    boxplot(expr[gg, ] ~ bmeta$group2, col = c("#9ECAE1", "#FC9272"), outline = FALSE,
            main = paste0(gg, " (raw)"), ylab = "log2 expression", xlab = "")
    boxplot(expr_c1[gg, ] ~ bmeta$group2, col = c("#9ECAE1", "#FC9272"), outline = FALSE,
            main = paste0(gg, " (ComBat)"), ylab = "log2 expression", xlab = "")
  }
  par(op)
}, width = 10, height = 6)

## F08: validation AUC by model (dot + DeLong 95% CI)
save_both("F08_validation_auc_comparison", {
  d <- auc_tab[!is.na(auc_valid)][order(auc_valid)]
  lab <- c(C0_original = "4 genes, raw (C0)",
           C1_combat_order = "4 genes, ComBat order (C1)",
           C2_combat_idblock = "4 genes, ComBat ID (C2)",
           C3_genes_plus_batch = "4 genes + batch block (C3)",
           C4_batch_only = "batch block only (C4)",
           C5_combat_plus_batch = "4 genes + block, ComBat (C5)",
           C1_combat_order_frozen_model = "frozen model on ComBat (C1)")
  d[, label := lab[condition]]
  d[, col := ifelse(grepl("only", label), "#7F7F7F",
             ifelse(grepl("ComBat", label), "#D7301F", "#2C7FB8"))]
  op <- par(mar = c(4.2, 18, 3, 4))
  y <- seq_len(nrow(d))
  xl <- range(c(d$valid_l, d$valid_u), na.rm = TRUE)
  plot(d$auc_valid, y, xlim = c(xl[1] - 0.01, xl[2] + 0.012), ylim = c(0.4, nrow(d) + 0.6),
       pch = 19, cex = 1.1, col = d$col, yaxt = "n", ylab = "",
       xlab = "AUC on the internal validation set (n = 94), DeLong 95% CI",
       main = "Discrimination before and after batch handling")
  axis(2, at = y, labels = d$label, las = 1, cex.axis = 0.78, tick = FALSE)
  segments(d$valid_l, y, d$valid_u, y, col = d$col, lwd = 2)
  segments(d$valid_l, y - 0.12, d$valid_l, y + 0.12, col = d$col, lwd = 2)
  segments(d$valid_u, y - 0.12, d$valid_u, y + 0.12, col = d$col, lwd = 2)
  text(d$auc_valid, y, sprintf("%.3f", d$auc_valid), pos = 3, cex = 0.72, col = "grey20")
  abline(v = 0.5, lty = 2, col = "grey75")
  par(op)
}, width = 9, height = 5.2)

## ===========================================================================
## 11. Logs
## ===========================================================================
writeLines(run_log, file.path(d_log, "71_run_log.txt"), useBytes = TRUE)
con <- file(file.path(d_log, "sessionInfo_71.txt"), open = "w"); writeLines(capture.output(sessionInfo()), con); close(con)
saveRDS(list(m2 = m2, frozen = frozen, deg = deg_summary), file.path(d_log, "71_objects.rds"))
log_msg("")
hr("=")
log_msg("Done. Comparison table: ", file.path(d_tab, "T14_comparison_table_main.csv"))
hr("=")
