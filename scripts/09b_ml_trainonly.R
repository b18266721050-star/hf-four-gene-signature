###############################################################################
## 09b_ml_trainonly.R
## Sensitivity analysis: repeat LASSO / SVM-RFE / random forest with feature
## selection restricted to the 219 training samples of the 7:3 split.
##
## Inputs (all already produced by the main pipeline):
##   results/ML/01_matrix/GSE57338_ML_matrix_sample_by_gene.csv  (313 x 212)
##   results/Nomogram/04_data/modelpredictions_train.csv         (train IDs)
## Output: results/ML_trainonly/
##
## Two runs:
##   RUN A (all 313)  -> sanity check of this re-implementation
##   RUN B (219 only) -> the actual sensitivity analysis
## Pure ASCII. Run with R_LIBS_USER set.
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(glmnet); library(randomForest); library(e1071)
})

proj <- "E:/another data file/GEObulk analysis project"
outd <- file.path(proj, "results", "ML_trainonly")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(outd, "run_log.txt")
con <- file(logf, open = "wt")
msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  writeLines(txt, con)
  flush(con)
}

SEED <- 20261027L
set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")

msg("=== 09b: feature selection restricted to the training split ===")
msg("R ", getRversion(), " | seed ", SEED)

## ---- data ------------------------------------------------------------------
mat <- read.csv(file.path(proj, "results", "ML", "01_matrix",
                          "GSE57338_ML_matrix_sample_by_gene.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
gene_cols <- setdiff(colnames(mat), c("sample", "group", "group_binary", "y"))
msg("matrix: ", nrow(mat), " samples x ", length(gene_cols), " genes")

tr <- read.csv(file.path(proj, "results", "Nomogram", "04_data",
                         "modelpredictions_train.csv"),
               stringsAsFactors = FALSE, encoding = "UTF-8", fileEncoding = "UTF-8-BOM")
names(tr)[1] <- "sample"
train_ids <- tr$sample[tr$set == "train"]
msg("training IDs from modelpredictions_train.csv: ", length(train_ids))
cat("overlap with matrix: ", sum(train_ids %in% mat$sample), "\n", sep = "")
msg("overlap with matrix: ", sum(train_ids %in% mat$sample))

## ---- helpers ---------------------------------------------------------------
make_folds <- function(y, k = 5L, seed = SEED) {
  set.seed(seed)
  idx <- list()
  for (cls in sort(unique(y))) {
    ii <- which(y == cls)
    ii <- sample(ii)
    f <- rep(seq_len(k), length.out = length(ii))
    for (j in seq_len(k)) idx[[length(idx) + 1L]] <- ii[f == j]
  }
  fold <- integer(length(y))
  for (j in seq_along(idx)) fold[idx[[j]]] <- ((j - 1L) %% k) + 1L
  fold
}

cv_error <- function(x, y, fold, cost = 1) {
  k <- max(fold); err <- numeric(k)
  for (j in seq_len(k)) {
    te <- which(fold == j); trn <- setdiff(seq_along(y), te)
    m <- svm(x[trn, , drop = FALSE], y = factor(y[trn]), kernel = "linear",
             cost = cost, scale = TRUE)
    p <- predict(m, x[te, , drop = FALSE])
    err[j] <- mean(p != factor(y[te]))
  }
  mean(err)
}

run_lasso <- function(x, y) {
  set.seed(SEED)
  fit <- cv.glmnet(as.matrix(x), y, family = "binomial", alpha = 1,
                   nfolds = 5L, type.measure = "deviance", standardize = TRUE)
  cf <- as.matrix(coef(fit, s = "lambda.min"))
  g <- rownames(cf)[cf[, 1] != 0]
  g <- setdiff(g, "(Intercept)")
  list(genes = g, lambda = fit$lambda.min)
}

run_svmrfe <- function(x, y, drop_prop = 0.10, min_feat = 2L, cost = 1) {
  set.seed(SEED)
  fold <- make_folds(y, 5L, SEED)
  keep <- colnames(x)
  res <- data.frame(n = integer(), err = numeric())
  path <- list()
  while (length(keep) >= min_feat) {
    e <- cv_error(x[, keep, drop = FALSE], y, fold, cost)
    res <- rbind(res, data.frame(n = length(keep), err = e))
    path[[as.character(length(keep))]] <- keep
    if (length(keep) == min_feat) break
    m <- svm(x[, keep, drop = FALSE], y = factor(y), kernel = "linear",
             cost = cost, scale = TRUE)
    w <- as.numeric(t(m$coefs) %*% m$SV)
    names(w) <- keep
    nd <- max(1L, floor(length(keep) * drop_prop))
    worst <- names(sort(abs(w)))[seq_len(nd)]
    keep <- setdiff(keep, worst)
  }
  best_err <- min(res$err)
  ## most features among the tied minima: standard "first minimum from the top"
  best_n <- max(res$n[res$err <= best_err + 1e-9])
  list(genes = path[[as.character(best_n)]], best_n = best_n,
       best_err = best_err, curve = res)
}

run_rf <- function(x, y, ntree = 500L) {
  set.seed(SEED)
  fit <- randomForest(x = x, y = factor(y), ntree = ntree,
                      importance = TRUE, scale = FALSE)
  imp <- importance(fit, type = 1)[, 1]
  thr <- mean(imp) + sd(imp)
  g <- names(imp)[imp > thr]
  list(genes = g, thr = thr, imp = imp)
}

run_all <- function(dat, tag) {
  x <- as.matrix(dat[, gene_cols, drop = FALSE])
  y <- dat$y
  msg("")
  msg("---- ", tag, ": n = ", nrow(x), " (cases ", sum(y == 1),
      " / controls ", sum(y == 0), ") ----")
  l <- run_lasso(x, y)
  msg(tag, " LASSO  : ", length(l$genes), " genes (lambda.min = ",
      signif(l$lambda, 4), ")")
  s <- run_svmrfe(x, y)
  msg(tag, " SVM-RFE: ", length(s$genes), " genes at best n = ", s$best_n,
      " (CV error ", signif(s$best_err, 4), ")")
  r <- run_rf(x, y)
  msg(tag, " RF     : ", length(r$genes), " genes (thr = ", signif(r$thr, 4), ")")
  inter <- Reduce(intersect, list(l$genes, s$genes, r$genes))
  msg(tag, " INTERSECTION (", length(inter), "): ", paste(inter, collapse = ", "))
  panel <- c("AQP3", "GGT5", "ISLR", "LAD1")
  msg(tag, " panel genes retained -> LASSO: ", paste(panel[panel %in% l$genes], collapse = ","),
      " | SVM: ", paste(panel[panel %in% s$genes], collapse = ","),
      " | RF: ", paste(panel[panel %in% r$genes], collapse = ","))
  write.csv(data.frame(gene = l$genes), file.path(outd, paste0(tag, "_lasso.csv")), row.names = FALSE)
  write.csv(data.frame(gene = s$genes), file.path(outd, paste0(tag, "_svmrfe.csv")), row.names = FALSE)
  write.csv(data.frame(gene = names(r$imp), mda = as.numeric(r$imp)),
            file.path(outd, paste0(tag, "_rf_importance.csv")), row.names = FALSE)
  write.csv(s$curve, file.path(outd, paste0(tag, "_svmrfe_curve.csv")), row.names = FALSE)
  list(lasso = l$genes, svm = s$genes, rf = r$genes, inter = inter)
}

## ---- RUN A: all 313 (sanity check) ----------------------------------------
a <- run_all(mat, "all313")

## ---- RUN B: training split only -------------------------------------------
sub <- mat[mat$sample %in% train_ids, , drop = FALSE]
b <- run_all(sub, "train219")

msg("")
msg("=== comparison ===")
msg("LASSO  overlap: ", length(intersect(a$lasso, b$lasso)), " genes common")
msg("SVM    overlap: ", length(intersect(a$svm, b$svm)))
msg("RF     overlap: ", length(intersect(a$rf, b$rf)))
msg("panel in RUN A intersection: ", paste(intersect(a$inter, c("AQP3","GGT5","ISLR","LAD1")), collapse = ", "))
msg("panel in RUN B intersection: ", paste(intersect(b$inter, c("AQP3","GGT5","ISLR","LAD1")), collapse = ", "))
msg("done")
close(con)
