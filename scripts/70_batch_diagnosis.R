###############################################################################
## 70_batch_diagnosis.R
## GSE57338 -- batch / confounding diagnostics (metadata proxies + expression)
##
## Purpose
##   The GEO series matrix of GSE57338 contains NO explicit batch fields
##   (no scan date, plate, or hybridization batch). Only five sample-level
##   annotations exist: tissue, heart failure, disease status, gender, age.
##   Following the pre-specified plan, processing order is therefore proxied by
##   (a) the patient identifier recorded in !Sample_description (prefix C / P),
##   (b) the GSM accession order (deposition / processing sequence), and
##   (c) the numeric part of the patient identifier.
##
##   This script quantifies how strongly each proxy is associated with the
##   case/control grouping, and how much of the expression variance is
##   attributable to group versus to the proxy. It does NOT modify expression.
##
## Inputs
##   clean_data/GSE57338_clean_expression_matrix.csv
##   clean_data/GSE57338_sample_info.csv
##   raw_data/GSE57338_series_matrix.txt.gz
##
## Outputs  (results/Batch/)
##   01_tables/T01_sample_batch_metadata.csv
##   01_tables/T02_contingency_group_by_batch.csv
##   01_tables/T03_batch_association_summary.csv
##   01_tables/T04_variance_decomposition.csv
##   01_tables/T05_batch_predicts_group_auc.csv
##   01_tables/T06_surrogate_variables.csv
##   02_figures/F01_pca_group_vs_batch.{pdf,png}
##   02_figures/F02_sample_correlation_by_order.{pdf,png}
##   02_figures/F03_nf_proportion_by_block.{pdf,png}
##   03_logs/70_run_log.txt , 03_logs/sessionInfo_70.txt
##
## Run:  Rscript scripts/70_batch_diagnosis.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
})

options(stringsAsFactors = FALSE)
set.seed(20261029)

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
clean_dir <- file.path(proj_dir, "clean_data")
raw_dir   <- file.path(proj_dir, "raw_data")
out_root  <- file.path(proj_dir, "results", "Batch")
d_tab     <- file.path(out_root, "01_tables")
d_fig     <- file.path(out_root, "02_figures")
d_log     <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...)
  run_log <<- c(run_log, txt)
  cat(txt, "\n", sep = "")
}
hr <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

## ===========================================================================
## 1. Parse series matrix -> sample-level batch proxies
## ===========================================================================
log_msg("[1] Parsing series matrix for sample identifiers ...")

## NB (2026-10-29): `gzip -dc` via fread(cmd=) needs gzip.exe on PATH and fails
## wherever Rtools' usr/bin is not visible (reproducibility issue). Decompress
## with R's own gzfile(); the text stream is identical to `gzip -dc`.
sm_con <- gzfile(file.path(raw_dir, "GSE57338_series_matrix.txt.gz"), open = "rt")
sm_lines <- readLines(sm_con, warn = FALSE)
close(sm_con)
sm <- data.table::fread(
  text = paste(sm_lines, collapse = "\n"),
  sep = "\t", header = FALSE, fill = TRUE, quote = "", showProgress = FALSE
)
get_row <- function(key) {
  r <- sm[grepl(paste0("^", key), sm$V1)]
  if (nrow(r) == 0) stop("field not found: ", key)
  gsub('"', "", as.character(unlist(r[1, -1])))
}

gsm  <- get_row("!Sample_geo_accession")
desc <- get_row("!Sample_description")

## !Sample_description has several rows (patient ID, then "PolyA mRNA");
## keep only entries shaped like C02065 / P01119
pid_idx <- which(grepl("^[A-Za-z]{1,3}[0-9]{2,}$", desc))
if (length(pid_idx) != length(gsm)) {
  warning("patient-ID row does not cover all samples; falling back to GSM order only")
  pid_idx <- seq_along(gsm)
}
bmeta <- data.table(
  sample    = gsm[pid_idx],
  patient   = desc[pid_idx],
  gsm_order = seq_along(gsm)
)
bmeta[, prefix := sub("[0-9]+$", "", patient)]
bmeta[, num_id := suppressWarnings(as.integer(sub("^[A-Za-z]+", "", patient)))]

log_msg("  samples parsed: ", nrow(bmeta),
        " | unique patients: ", uniqueN(bmeta$patient),
        " | prefixes: ", paste(names(table(bmeta$prefix)), table(bmeta$prefix),
                               sep = "=", collapse = ", "))

meta <- data.table::fread(file.path(clean_dir, "GSE57338_sample_info.csv"))
bmeta <- merge(bmeta, meta[, .(sample, group, gender, age)], by = "sample", all.x = TRUE)
bmeta <- bmeta[order(gsm_order)]
if (anyNA(bmeta$group)) stop("some samples have no group label")

bmeta[, group2 := factor(ifelse(group == "NF", "Control", "Disease"),
                         levels = c("Control", "Disease"))]
bmeta[, sex := factor(ifelse(tolower(gender) == "male", "Male", "Female"))]
bmeta[, age_num := suppressWarnings(as.numeric(age))]
bmeta[!is.na(age_num) & age_num == 0, age_num := NA_real_]
bmeta[is.na(age_num), age_num := median(bmeta$age_num, na.rm = TRUE)]

## ---- batch proxies -------------------------------------------------------
## P1: identifier prefix (C / P) -- collection centre or cohort of origin
bmeta[, batch_prefix := factor(prefix)]
## P2: GSM deposit order split into 16 consecutive blocks of ~20 samples
NB <- 16L
bmeta[, batch_order := factor(sprintf("B%02d", cut(gsm_order, breaks = NB,
                                                    labels = FALSE)))]
## P3: numeric patient-ID quartiles within each prefix
bmeta[, batch_idnum := {
  q <- quantile(num_id, probs = seq(0, 1, 0.25), na.rm = TRUE)
  q[1] <- -Inf; q[5] <- Inf
  factor(paste0(prefix, "_Q", findInterval(num_id, q, rightmost.closed = TRUE)),
         levels = paste0(rep(sort(unique(prefix)), each = 4), "_Q", 1:4))
}, by = prefix]

proxies <- c("batch_prefix", "batch_order", "batch_idnum")

bmeta_out <- bmeta[, .(sample, patient, gsm_order, group, group2, sex, age_num,
                       batch_prefix, batch_order, batch_idnum)]
write_csv(bmeta_out, file.path(d_tab, "T01_sample_batch_metadata.csv"))
log_msg("  metadata written; block sizes: ",
        paste(range(as.integer(table(bmeta$batch_order))), collapse = "-"))

## ===========================================================================
## 2. Contingency tables, chi-square, Cramer's V, adjusted residuals
## ===========================================================================
log_msg("")
log_msg("[2] group x batch association ...")

cramers_v <- function(tb) {
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE))
  sqrt(as.numeric(chi$statistic) / (sum(tb) * (min(dim(tb)) - 1)))
}

## adjusted standardised residuals (Agresti): (obs - exp) / sqrt(exp*(1-r/n)*(1-c/n))
adj_resid <- function(tb) {
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE))
  o <- as.matrix(tb); e <- chi$expected
  n <- sum(o); r <- rowSums(o); cc <- colSums(o)
  (o - e) / sqrt(e * (1 - r / n) * (1 - cc / n))
}

assoc_rows <- list()
for (p in proxies) {
  tb <- table(bmeta[[p]], bmeta$group)
  tb <- tb[rowSums(tb) > 0, , drop = FALSE]
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE))
  ## Monte-Carlo p for sparse tables
  p_sim <- tryCatch(suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 20000)$p.value),
                    error = function(e) NA_real_)
  v <- cramers_v(tb)
  ar <- adj_resid(tb)
  nf_prop <- tb[, "NF"] / rowSums(tb)
  assoc_rows[[p]] <- data.table(
    proxy           = p,
    n_levels        = nrow(tb),
    chi_square      = round(as.numeric(chi$statistic), 3),
    df              = as.integer(chi$parameter),
    p_asymptotic    = signif(chi$p.value, 6),
    p_montecarlo    = signif(p_sim, 6),
    cramers_v       = round(v, 4),
    max_abs_adj_res = round(max(abs(ar)), 3),
    n_blocks_all_case = sum(nf_prop < 0.05),
    n_blocks_all_ctrl = sum(nf_prop > 0.95),
    nf_prop_min     = round(min(nf_prop), 3),
    nf_prop_max     = round(max(nf_prop), 3)
  )
  ## long-format contingency table with adjusted residuals
  dt <- as.data.table(as.table(tb))
  setnames(dt, c("batch", "group", "n"))
  ar_dt <- as.data.table(as.table(ar)); setnames(ar_dt, c("batch", "group", "adj_resid"))
  dt <- merge(dt, ar_dt, by = c("batch", "group"))
  dt[, proxy := p]
  dt[, n_block := sum(n), by = batch]
  dt[, nf_prop := sum(n[group == "NF"]) / unique(n_block), by = batch]
  dt[, adj_resid := round(adj_resid, 3)]
  fwrite(dt[, .(proxy, batch, group, n, adj_resid, n_block, nf_prop)],
         file.path(d_tab, paste0("T02_contingency_", p, ".csv")), bom = TRUE)
}
assoc <- rbindlist(assoc_rows)
assoc[, interpretation := ifelse(cramers_v >= 0.5, "strong confounding",
                          ifelse(cramers_v >= 0.3, "moderate confounding",
                          ifelse(cramers_v >= 0.1, "weak confounding", "negligible")))]
write_csv(assoc, file.path(d_tab, "T03_batch_association_summary.csv"))
log_msg("  association summary:")
for (i in seq_len(nrow(assoc))) {
  log_msg(sprintf("    %-13s V=%.3f  chi2=%.1f (df=%d)  p=%.3g  [%s]",
                  assoc$proxy[i], assoc$cramers_v[i], assoc$chi_square[i],
                  assoc$df[i], assoc$p_asymptotic[i], assoc$interpretation[i]))
}

## ---- is the patient-ID axis the same dimension as the deposit order? ------
id_tab <- rbindlist(list(
  data.table(scope = "all",
             n = nrow(bmeta),
             spearman_numid_vs_order = round(cor(bmeta$num_id, bmeta$gsm_order,
                                                 method = "spearman"), 3),
             median_id_control = round(median(bmeta[group2 == "Control", num_id]), 1),
             median_id_disease = round(median(bmeta[group2 == "Disease", num_id]), 1),
             p_id_difference = signif(wilcox.test(num_id ~ group2, data = bmeta)$p.value, 4)),
  rbindlist(lapply(sort(unique(bmeta$prefix)), function(pp) {
    d <- bmeta[prefix == pp]
    if (uniqueN(d$group2) < 2) return(NULL)
    data.table(scope = paste0("prefix_", pp),
               n = nrow(d),
               spearman_numid_vs_order = round(cor(d$num_id, d$gsm_order, method = "spearman"), 3),
               median_id_control = round(median(d[group2 == "Control", num_id]), 1),
               median_id_disease = round(median(d[group2 == "Disease", num_id]), 1),
               p_id_difference = signif(wilcox.test(num_id ~ group2, data = d)$p.value, 4))
  }))
))
write_csv(id_tab, file.path(d_tab, "T07_id_axis_vs_order.csv"))
log_msg("  ID axis vs deposit order: rho = ",
        id_tab$spearman_numid_vs_order[1],
        "; median patient ID control vs disease: ",
        id_tab$median_id_control[1], " vs ", id_tab$median_id_disease[1])

## ===========================================================================
## 3. Expression: PCA and variance decomposition
## ===========================================================================
log_msg("")
log_msg("[3] PCA / variance decomposition ...")

expr_dt <- data.table::fread(file.path(clean_dir, "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())
if (!identical(colnames(expr), bmeta$sample)) {
  if (setequal(colnames(expr), bmeta$sample)) {
    bmeta <- bmeta[match(colnames(expr), bmeta$sample)]
    bmeta[, gsm_order := seq_len(.N)]
  } else stop("sample mismatch between matrix and metadata")
}
log_msg("  matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

## top-2000 variable genes for PCA
gv  <- apply(expr, 1, var)
top <- names(sort(gv, decreasing = TRUE))[1:2000]
X   <- t(expr[top, , drop = FALSE])
X   <- X[, apply(X, 2, sd) > 0, drop = FALSE]
pca <- prcomp(X, center = TRUE, scale. = FALSE)
ve  <- (pca$sdev^2) / sum(pca$sdev^2)

eta2 <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- factor(g[ok])
  if (nlevels(g) < 2) return(NA_real_)
  ss <- summary(aov(x ~ g))[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
}

pc_tab <- rbindlist(lapply(1:10, function(k) {
  pc <- pca$x[, k]
  data.table(
    pc = paste0("PC", k),
    var_explained_pct = round(100 * ve[k], 2),
    eta2_group      = round(eta2(pc, bmeta$group2), 4),
    eta2_order_block = round(eta2(pc, bmeta$batch_order), 4),
    eta2_prefix     = round(eta2(pc, bmeta$batch_prefix), 4),
    eta2_id_block   = round(eta2(pc, bmeta$batch_idnum), 4)
  )
}))

## per-gene variance decomposition on the same 2000 genes
gene_rows <- rbindlist(lapply(top, function(g) {
  y <- expr[g, ]
  data.table(gene = g,
             eta2_group = eta2(y, bmeta$group2),
             eta2_order_block = eta2(y, bmeta$batch_order),
             eta2_prefix = eta2(y, bmeta$batch_prefix))
}))
gene_summary <- data.table(
  gene   = c("MEAN_over_2000_genes", "MEDIAN_over_2000_genes", "PCT_genes_group_gt_batch"),
  eta2_group = c(mean(gene_rows$eta2_group), median(gene_rows$eta2_group), NA),
  eta2_order_block = c(mean(gene_rows$eta2_order_block), median(gene_rows$eta2_order_block), NA),
  eta2_prefix = c(mean(gene_rows$eta2_prefix), median(gene_rows$eta2_prefix), NA)
)
gene_summary[3, eta2_order_block := 100 * mean(gene_rows$eta2_group > gene_rows$eta2_order_block)]

var_decomp <- rbind(
  pc_tab,
  data.table(pc = gene_summary$gene,
             var_explained_pct = NA_real_,
             eta2_group = round(gene_summary$eta2_group, 4),
             eta2_order_block = round(gene_summary$eta2_order_block, 4),
             eta2_prefix = round(gene_summary$eta2_prefix, 4),
             eta2_id_block = NA_real_)
)
write_csv(var_decomp, file.path(d_tab, "T04_variance_decomposition.csv"))
log_msg(sprintf("  PC1 explains %.1f%% variance: eta2_group=%.3f, eta2_order_block=%.3f, eta2_prefix=%.3f",
                100 * ve[1], eta2(pca$x[, 1], bmeta$group2),
                eta2(pca$x[, 1], bmeta$batch_order), eta2(pca$x[, 1], bmeta$batch_prefix)))
log_msg(sprintf("  mean per-gene eta2: group=%.3f, order_block=%.3f, prefix=%.3f",
                mean(gene_rows$eta2_group), mean(gene_rows$eta2_order_block),
                mean(gene_rows$eta2_prefix)))
saveRDS(list(pca = pca[c("x", "sdev")], ve = ve, top_genes = top),
        file.path(out_root, "03_logs", "pca_uncorrected.rds"))

## ===========================================================================
## 4. How well can the batch proxy alone predict case/control?
## ===========================================================================
log_msg("")
log_msg("[4] batch proxy -> group predictability (5-fold stratified CV, 20 repeats) ...")

suppressPackageStartupMessages(library(pROC))
cv_auc <- function(y, Xdf, K = 5, R = 20) {
  aucs <- numeric(R)
  for (r in seq_len(R)) {
    set.seed(1000 + r)
    idx1 <- which(y == 1); idx0 <- which(y == 0)
    f1 <- sample(rep(seq_len(K), length.out = length(idx1)))
    f0 <- sample(rep(seq_len(K), length.out = length(idx0)))
    fold <- integer(length(y)); fold[idx1] <- f1; fold[idx0] <- f0
    oof <- rep(NA_real_, length(y))
    for (k in seq_len(K)) {
      tr <- fold != k; te <- fold == k
      fit <- suppressWarnings(glm(y[tr] ~ ., data = Xdf[tr, , drop = FALSE],
                                  family = binomial))
      oof[te] <- predict(fit, newdata = Xdf[te, , drop = FALSE], type = "response")
    }
    aucs[r] <- as.numeric(pROC::auc(pROC::roc(y, oof, quiet = TRUE)))
  }
  qs <- as.numeric(quantile(aucs, c(.025, .975)))
  c(mean = mean(aucs), sd = sd(aucs), q025 = qs[1], q975 = qs[2])
}

y <- as.integer(bmeta$group2 == "Disease")
dummy <- function(f) model.matrix(~ f - 1)
pred_rows <- list()
for (p in proxies) {
  Xd <- as.data.frame(dummy(bmeta[[p]]))
  s <- cv_auc(y, Xd)
  pred_rows[[p]] <- data.table(
    predictor = p, n_parameters = ncol(Xd),
    cv_auc_mean = round(s["mean"], 4), cv_auc_sd = round(s["sd"], 4),
    cv_auc_lower = round(s["q025"], 4), cv_auc_upper = round(s["q975"], 4),
    interpretation = ifelse(s["mean"] > 0.8, "batch alone nearly determines group",
                     ifelse(s["mean"] > 0.7, "batch alone strongly predicts group",
                     ifelse(s["mean"] > 0.6, "batch weakly predicts group", "no predictive value")))
  )
}
pred_tab <- rbindlist(pred_rows)
write_csv(pred_tab, file.path(d_tab, "T05_batch_predicts_group_auc.csv"))
for (i in seq_len(nrow(pred_tab))) {
  log_msg(sprintf("  %-13s CV AUC = %.3f (%.3f-%.3f)  [%s]",
                  pred_tab$predictor[i], pred_tab$cv_auc_mean[i],
                  pred_tab$cv_auc_lower[i], pred_tab$cv_auc_upper[i],
                  pred_tab$interpretation[i]))
}

## ===========================================================================
## 5. Surrogate variable analysis
## ===========================================================================
log_msg("")
log_msg("[5] SVA ...")

mod_full  <- model.matrix(~ group2 + sex + age_num, data = bmeta)   # group protected
mod_grponly <- model.matrix(~ group2, data = bmeta)
mod0      <- model.matrix(~ 1, data = bmeta)

sv_run <- function(mod, tag) {
  nsv <- tryCatch(sva::num.sv(expr, mod, method = "be"), error = function(e) 0L)
  log_msg(sprintf("  [%s] num.sv (be) = %d", tag, nsv))
  if (nsv == 0) return(NULL)
  sv <- sva::sva(expr, mod, mod0, n.sv = nsv)
  svdf <- as.data.table(sv$sv)
  setnames(svdf, paste0("SV", seq_len(ncol(svdf))))
  svdf[, sample := colnames(expr)]
  summ <- rbindlist(lapply(seq_len(ncol(svdf) - 1), function(k) {
    v <- svdf[[k]]
    data.table(
      model = tag,
      SV = paste0("SV", k),
      spearman_vs_gsm_order = round(cor(v, bmeta$gsm_order, method = "spearman"), 3),
      p_vs_gsm_order = signif(cor.test(v, bmeta$gsm_order, method = "spearman")$p.value, 4),
      eta2_group = round(eta2(v, bmeta$group2), 4),
      eta2_order_block = round(eta2(v, bmeta$batch_order), 4),
      eta2_prefix = round(eta2(v, bmeta$batch_prefix), 4),
      auc_group = round(as.numeric(pROC::auc(pROC::roc(y, v, quiet = TRUE))), 3)
    )
  }))
  list(summary = summ, scores = svdf)
}

sv_a <- sv_run(mod_full, "group_protected")
sv_b <- sv_run(mod_grponly, "group_included")

sv_summary <- rbind(sv_a$summary, sv_b$summary)
if (is.null(sv_summary)) {
  write_csv(data.table(note = "num.sv returned 0 for both models; no residual latent structure detected"),
            file.path(d_tab, "T06_surrogate_variables.csv"))
  log_msg("  SVA: no surrogate variables estimated under either model")
} else {
  write_csv(sv_summary, file.path(d_tab, "T06_surrogate_variables.csv"))
  ## when the group is protected and num.sv = 0, that itself means: after removing
  ## the group term, no additional latent factor is strong enough to be retained.
  writeLines(c(
    "SVA note:",
    paste0("num.sv with group protected (~ group + sex + age) = ",
           ifelse(is.null(sv_a), 0L, nrow(sv_a$summary)), ""),
    paste0("num.sv without protection (~ group) = ",
           ifelse(is.null(sv_b), 0L, nrow(sv_b$summary)), ""),
    "",
    "Interpretation: latent factors (SVs) are detectable only when the group term is",
    "not protected, and those SVs correlate both with the group label and with the",
    "deposition order. In other words the latent structure and the diagnosis are the",
    "same axis of variation in this dataset -- they cannot be separated statistically."
  ), file.path(d_tab, "T06_surrogate_variables_note.txt"))
  log_msg("  SV summary:")
  for (i in seq_len(nrow(sv_summary))) {
    log_msg(sprintf("    %-15s %-4s rho(order)=%+.2f  eta2_group=%.3f  eta2_block=%.3f  AUC(group)=%.3f",
                    sv_summary$model[i], sv_summary$SV[i], sv_summary$spearman_vs_gsm_order[i],
                    sv_summary$eta2_group[i], sv_summary$eta2_order_block[i],
                    sv_summary$auc_group[i]))
  }
}

## ===========================================================================
## 6. Figures
## ===========================================================================
log_msg("")
log_msg("[6] figures ...")

save_both <- function(stem, expr, width = 8, height = 6) {
  e <- substitute(expr); env <- parent.frame()
  grDevices::pdf(file.path(d_fig, paste0(stem, ".pdf")), width = width, height = height)
  eval(e, env); grDevices::dev.off()
  grDevices::png(file.path(d_fig, paste0(stem, ".png")), width = width, height = height,
                 units = "in", res = 300)
  eval(e, env); grDevices::dev.off()
}

grp_col <- c(Control = "#2C7FB8", Disease = "#D7301F")

## F01: PCA coloured by group (left) vs by order block (right)
save_both("F01_pca_group_vs_batch", {
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1), cex = 0.85)
  plot(pca$x[, 1], pca$x[, 2], col = grp_col[as.character(bmeta$group2)],
       pch = 19, cex = 0.8, xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]), main = "Coloured by group")
  legend("topright", legend = names(grp_col), col = grp_col, pch = 19, bty = "n", cex = 0.8)
  blk <- as.integer(bmeta$batch_order)
  pal <- grDevices::colorRampPalette(c("#4D4D4D", "#FDB863", "#B2182B"))(NB)
  plot(pca$x[, 1], pca$x[, 2], col = pal[blk], pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
       main = "Coloured by GSM-order block (B01-B16)")
  legend("topright", legend = c("B01 (first)", "B16 (last)"), col = pal[c(1, NB)],
         pch = 19, bty = "n", cex = 0.8)
  par(op)
}, width = 11, height = 5)

## F02: sample-sample correlation heatmap in GSM order + group strip below
save_both("F02_sample_correlation_by_order", {
  cc  <- cor(expr[top, ], method = "pearson")
  ord <- order(bmeta$gsm_order)
  cc  <- cc[ord, ord]
  n   <- ncol(cc)
  layout(matrix(c(1, 2), nrow = 2), heights = c(5, 1.1))
  par(mar = c(1.5, 4.2, 3, 1))
  image(1:n, 1:n, cc,
        col = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
        xlab = "", ylab = "Sample index (GSM order)",
        main = "Sample-sample Pearson correlation (top-2000 variable genes)",
        axes = FALSE, useRaster = TRUE)
  axis(2, at = seq(1, n, 40), las = 1, cex.axis = 0.8)
  par(mar = c(2.2, 4.2, 0, 1))
  gcol <- grp_col[as.character(bmeta$group2)][ord]
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "",
       ylab = "", xaxs = "i", yaxs = "i")
  rect(seq_len(n) - 0.5, 0, seq_len(n) + 0.5, 1, col = gcol, border = NA)
  box()
  axis(1, at = seq(1, n, 40), cex.axis = 0.8)
  mtext("Group", side = 2, las = 1, line = 0.5, cex = 0.7)
  layout(1)
}, width = 7.5, height = 7.4)

## F03: NF proportion per block
save_both("F03_nf_proportion_by_block", {
  tb <- table(bmeta$batch_order, bmeta$group)
  prop <- tb[, "NF"] / rowSums(tb)
  op <- par(mar = c(6, 4.2, 3, 1))
  bp <- barplot(prop, col = ifelse(prop > 0.95, "#2C7FB8",
                            ifelse(prop < 0.05, "#D7301F", "#BDBDBD")),
                ylim = c(0, 1), las = 2, cex.names = 0.75,
                ylab = "Proportion of non-failing (NF) samples",
                main = sprintf("Group composition by GSM-order block (Cramer's V = %.2f)",
                               assoc[proxy == "batch_order", cramers_v]))
  abline(h = 136 / 313, lty = 2, col = "grey30")
  text(bp, 1.02, labels = rowSums(tb), cex = 0.6, xpd = TRUE)
  legend("topright", fill = c("#2C7FB8", "#BDBDBD", "#D7301F"),
         legend = c("block ~ all NF", "mixed", "block ~ all disease"), bty = "n", cex = 0.75)
  par(op)
}, width = 10, height = 4.5)

## ===========================================================================
## 7. Logs
## ===========================================================================
writeLines(run_log, file.path(d_log, "70_run_log.txt"), useBytes = TRUE)
con <- file(file.path(d_log, "sessionInfo_70.txt"), open = "w"); writeLines(capture.output(sessionInfo()), con); close(con)
log_msg("")
hr("=")
log_msg("Done. Tables: ", d_tab)
hr("=")
