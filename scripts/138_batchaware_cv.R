###############################################################################
## 138_batchaware_cv.R
##
## R9 of the pre-submission audit: make the batch-aware cross-validation the
## primary estimate of internal performance, with uncertainty.
##
## The published table (Batch/T12) reports one fixed block-to-fold assignment for
## the batch-aware scheme and 20 repeats for the random scheme, so the two are
## not directly comparable and the batch-aware number has no interval. Here both
## schemes are repeated 20 times with the same fold-flexibility, and the paired
## difference is reported.
##
## Outputs (results/Batch/01_tables)
##   T24_batchaware_cv_uncertainty.csv
##   03_logs/138_run_log.txt , sessionInfo_138.txt
##
## Run: Rscript scripts/138_batchaware_cv.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

proj <- "E:/another data file/GEObulk analysis project"
d_tab <- file.path(proj, "results", "Batch", "01_tables")
d_log <- file.path(proj, "results", "Batch", "03_logs")
dir.create(d_log, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = ""); flush.console()
}
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")
REPS <- 20

expr_dt <- data.table::fread(file.path(proj, "clean_data", "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())

bmeta <- data.table::fread(file.path(d_tab, "T01_sample_batch_metadata.csv"))
bmeta <- bmeta[match(colnames(expr), sample)]
stopifnot(!anyNA(bmeta$sample))
df <- data.table(sample = bmeta$sample, HF = as.integer(bmeta$group2 == "Disease"),
                 block = bmeta$batch_order)
g <- t(expr[GENES, , drop = FALSE]); colnames(g) <- GENES
df <- cbind(df, as.data.table(g))
log_msg("samples: ", nrow(df), " ; cases ", sum(df$HF), " ; blocks ", length(unique(df$block)))

## ---- fold assignment -------------------------------------------------------
## stratified random folds
folds_random <- function(y, K = 5, seed) {
  set.seed(seed)
  f <- integer(length(y))
  for (lv in 0:1) {
    idx <- which(y == lv)
    f[idx] <- sample(rep(1:K, length.out = length(idx)))
  }
  f
}
## block-to-fold assignment: blocks kept intact, random greedy balancing
folds_block <- function(block, K = 5, seed) {
  set.seed(seed)
  tab <- table(block)
  ord <- sample(names(tab))                 # random order, then greedy balancing
  fold_of_block <- integer(length(ord)); names(fold_of_block) <- ord
  load <- numeric(K)
  for (b in ord) {
    k <- which.min(load)
    fold_of_block[b] <- k
    load[k] <- load[k] + as.integer(tab[b])
  }
  unname(fold_of_block[block])
}

cv_once <- function(fold) {
  pr <- rep(NA_real_, nrow(df))
  for (k in sort(unique(fold))) {
    tr <- fold != k; te <- fold == k
    if (sum(df$HF[tr]) == 0 || sum(df$HF[tr] == 0) == 0) return(NA_real_)
    zp <- sapply(GENES, function(x) c(mean(df[[x]][tr]), sd(df[[x]][tr])))
    ztr <- sapply(GENES, function(x) (df[[x]][tr] - zp[1, x]) / zp[2, x])
    zte <- sapply(GENES, function(x) (df[[x]][te] - zp[1, x]) / zp[2, x])
    dtr <- data.table(HF = df$HF[tr], as.data.table(ztr))
    fit <- glm(HF ~ ., data = dtr, family = binomial)
    pr[te] <- predict(fit, newdata = data.table(as.data.table(zte)), type = "response")
  }
  as.numeric(pROC::auc(pROC::roc(df$HF, pr, quiet = TRUE)))
}

log_msg("[1] ", REPS, " repeats of each scheme ...")
rnd <- sapply(1:REPS, function(r) cv_once(folds_random(df$HF, 5, 3000 + r)))
bat <- sapply(1:REPS, function(r) cv_once(folds_block(df$block, 5, 5000 + r)))
log_msg("  random stratified : mean ", round(mean(rnd), 4), " sd ", round(sd(rnd), 4))
log_msg("  batch-aware       : mean ", round(mean(bat), 4), " sd ", round(sd(bat), 4))

d <- bat - rnd
ci <- function(x) as.numeric(quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE))
ci_t <- function(tt) as.numeric(tt$conf.int)
tt <- t.test(bat, rnd, paired = TRUE)
out <- rbind(
  data.table(scheme = "random stratified 5-fold CV", repeats = REPS,
             auc_mean = round(mean(rnd), 4), auc_sd = round(sd(rnd), 4),
             ci_low = round(ci(rnd)[1], 4), ci_high = round(ci(rnd)[2], 4),
             paired_diff = NA_real_, diff_ci_low = NA_real_, diff_ci_high = NA_real_),
  data.table(scheme = "batch-aware 5-fold CV (GSM-order blocks intact)", repeats = REPS,
             auc_mean = round(mean(bat), 4), auc_sd = round(sd(bat), 4),
             ci_low = round(ci(bat)[1], 4), ci_high = round(ci(bat)[2], 4),
             paired_diff = round(mean(d), 4),
             diff_ci_low = round(ci_t(tt)[1], 4), diff_ci_high = round(ci_t(tt)[2], 4))
)
write_csv(out, file.path(d_tab, "T24_batchaware_cv_uncertainty.csv"))
log_msg("  paired difference (batch-aware - random): ", round(mean(d), 4),
        " (", round(ci(d)[1], 4), " to ", round(ci(d)[2], 4), ")")

## paired t test on the 20 repeats
log_msg("  paired t test p = ", signif(tt$p.value, 3))
log_msg("Done.")
writeLines(run_log, file.path(d_log, "138_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_138.txt"))
