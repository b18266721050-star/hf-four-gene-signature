###############################################################################
# 03_immune_stats.R
# Group comparison (NF / DCM / ICM) + key-gene correlation
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(car)
})

ROOT   <- "E:/another data file/GEObulk analysis project"
OUT    <- file.path(ROOT, "results", "Immune")
DIR_S  <- file.path(OUT, "02_ssgsea")
DIR_T  <- file.path(OUT, "04_stats")
DIR_LOG<- file.path(OUT, "00_logs")
dir.create(DIR_T, recursive = TRUE, showWarnings = FALSE)

KEY_GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")

sink(file.path(DIR_LOG, "03_stats_console.log"), split = TRUE)
cat("=== 03_immune_stats.R started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n\n")

obj  <- readRDS(file.path(OUT, "01_qc", "qc_input.rds"))
expr <- obj$expr
es   <- readRDS(file.path(DIR_S, "GSE57338_immune_ssGSEA.rds"))$es
si   <- obj$si

stopifnot(identical(colnames(es), si$sample))

grp     <- factor(si$group, levels = c("NF", "DCM", "ICM"))
cells   <- rownames(es)
cat("[input]", length(cells), "cell types x", ncol(es), "samples\n")
cat("[groups]", paste(levels(grp), as.integer(table(grp)), sep = "=", collapse = "; "), "\n\n")

## ---------------------------------------------------------------------------
## 1. Distribution diagnostics -> test selection
## ---------------------------------------------------------------------------
norm_tab <- data.table()
for (cl in cells) {
  for (g in levels(grp)) {
    x <- as.numeric(es[cl, grp == g])
    sw <- tryCatch(shapiro.test(x), error = function(e) NULL)
    norm_tab <- rbind(norm_tab, data.table(
      cell_type = cl, group = g, n = length(x),
      mean = mean(x), sd = sd(x), median = median(x), IQR = IQR(x),
      skewness = as.numeric(moments_skew <- (sum((x - mean(x))^3) / length(x)) /
                              (sd(x)^3)),
      shapiro_W = if (is.null(sw)) NA_real_ else as.numeric(sw$statistic),
      shapiro_P = if (is.null(sw)) NA_real_ else sw$p.value
    ))
  }
}
norm_tab[, normal_0.05 := shapiro_P > 0.05]

## variance homogeneity (Levene, median-centred)
lev <- sapply(cells, function(cl) {
  x <- as.numeric(es[cl, ])
  tryCatch({ r <- leveneTest(x ~ grp, center = median)
             c(F = r$`F value`[1], P = r$`Pr(>F)`[1]) },
           error = function(e) c(F = NA, P = NA))
})
var_tab <- data.table(cell_type = cells, levene_F = lev["F", ], levene_P = lev["P", ])

n_nonnormal <- norm_tab[, .(n_group_nonnormal = sum(!normal_0.05), n_group = .N), by = cell_type]
n_nonnormal[, all_normal := n_group_nonnormal == 0]

cat("[normality] cell types with all 3 groups normal:",
    sum(n_nonnormal$all_normal), "/", nrow(n_nonnormal), "\n")
cat("[normality] Shapiro P<0.05 in at least one group:",
    sum(!n_nonnormal$all_normal), "cell types\n")
cat("[normality] median Shapiro P =", signif(median(norm_tab$shapiro_P, na.rm = TRUE), 4), "\n")
cat("[variance] Levene P<0.05 in:", sum(var_tab$levene_P < 0.05, na.rm = TRUE),
    "/", nrow(var_tab), "cell types\n")

## decision rule: any non-normal group OR heteroscedastic -> non-parametric
need_nonparam <- any(!n_nonnormal$all_normal) || any(var_tab$levene_P < 0.05, na.rm = TRUE)
cat("[decision] non-parametric pipeline used: ", need_nonparam, "\n\n", sep = "")

write.csv(norm_tab, file.path(DIR_T, "immune_normality_by_group.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(var_tab, file.path(DIR_T, "immune_variance_homogeneity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------------------------------------------------------
## 2. Global test: Kruskal-Wallis across 3 groups
## ---------------------------------------------------------------------------
kw <- rbindlist(lapply(cells, function(cl) {
  x  <- as.numeric(es[cl, ])
  kt <- kruskal.test(x ~ grp)
  data.table(cell_type = cl,
             kw_chi2 = as.numeric(kt$statistic),
             df = as.numeric(kt$parameter),
             p_value = kt$p.value,
             epsilon2 = as.numeric(kt$statistic) / (ncol(es) - 1))
}))
kw[, FDR := p.adjust(p_value, method = "BH")]
setorder(kw, p_value)

## ---------------------------------------------------------------------------
## 3. Pairwise tests: Wilcoxon rank-sum + Cliff's delta
## ---------------------------------------------------------------------------
cliff_delta <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  r <- rank(c(x, y))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2) * 2 - 1
}

contrasts <- list(c("NF", "DCM"), c("NF", "ICM"), c("DCM", "ICM"))
pw <- rbindlist(lapply(contrasts, function(cn) {
  g1 <- cn[1]; g2 <- cn[2]
  rbindlist(lapply(cells, function(cl) {
    x <- as.numeric(es[cl, grp == g1]); y <- as.numeric(es[cl, grp == g2])
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = TRUE))
    d  <- cliff_delta(x, y)
    data.table(
      contrast = paste0(g1, "_vs_", g2),
      cell_type = cl,
      n_1 = length(x), n_2 = length(y),
      median_1 = median(x), median_2 = median(y),
      mean_1 = mean(x), mean_2 = mean(y),
      median_diff = median(y) - median(x),
      W = as.numeric(wt$statistic),
      p_value = wt$p.value,
      cliffs_delta = d,
      effect_magnitude = ifelse(abs(d) < 0.147, "negligible",
                         ifelse(abs(d) < 0.33, "small",
                         ifelse(abs(d) < 0.474, "medium", "large"))),
      direction = ifelse(d > 0, paste0("higher_in_", g1),
                  ifelse(d < 0, paste0("higher_in_", g2), "no_direction"))
    )
  }))
}))
pw[, FDR := p.adjust(p_value, method = "BH"), by = contrast]
pw[, FDR_all := p.adjust(p_value, method = "BH")]
setorder(pw, contrast, p_value)

write.csv(kw, file.path(DIR_T, "immune_kruskal_global.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pw, file.path(DIR_T, "immune_pairwise_wilcoxon.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## significant summary
sig <- pw[FDR < 0.05]
cat("[global KW] FDR<0.05:", sum(kw$FDR < 0.05), "/", nrow(kw), "cell types\n")
cat("[pairwise] FDR<0.05:", nrow(sig), "of", nrow(pw), "tests\n")
cat("[pairwise] by contrast:\n"); print(sig[, .N, by = contrast])

## wide effect table for plotting
eff_wide <- dcast(pw, cell_type ~ contrast, value.var = c("p_value", "FDR", "cliffs_delta", "median_diff"))
write.csv(eff_wide, file.path(DIR_T, "immune_effect_summary_wide.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------------------------------------------------------
## 4. Key gene vs immune cell correlation (Spearman + BH)
## ---------------------------------------------------------------------------
genes_present <- intersect(KEY_GENES, rownames(expr))
cat("\n[key genes] used:", paste(genes_present, collapse = ", "), "\n")

cor_tab <- rbindlist(lapply(genes_present, function(g) {
  gv <- as.numeric(expr[g, ])
  rbindlist(lapply(cells, function(cl) {
    cv <- as.numeric(es[cl, ])
    ct <- suppressWarnings(cor.test(gv, cv, method = "spearman", exact = FALSE))
    data.table(gene = g, cell_type = cl,
               rho = as.numeric(ct$estimate),
               p_value = ct$p.value,
               n = length(gv))
  }))
}))
cor_tab[, FDR := p.adjust(p_value, method = "BH"), by = gene]
cor_tab[, FDR_all := p.adjust(p_value, method = "BH")]
cor_tab[, strength := ifelse(abs(rho) < 0.2, "very_weak",
                      ifelse(abs(rho) < 0.4, "weak",
                      ifelse(abs(rho) < 0.6, "moderate", "strong")))]
setorder(cor_tab, gene, p_value)

write.csv(cor_tab, file.path(DIR_T, "immune_keygene_spearman.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("[correlation] FDR<0.05:", sum(cor_tab$FDR < 0.05), "of", nrow(cor_tab), "gene-cell pairs\n")
print(cor_tab[FDR < 0.05, .(gene, cell_type, rho = round(rho, 3),
                             P = signif(p_value, 3), FDR = signif(FDR, 3))][order(gene, -abs(rho))])

## correlation among immune cells (for heatmap)
cm <- cor(t(es), method = "spearman")
write.csv(data.frame(cell_type = rownames(cm), cm, check.names = FALSE),
          file.path(DIR_T, "immune_cell_correlation_spearman.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## p-values for the correlation heatmap
cmp <- matrix(NA_real_, nrow(es), nrow(es),
              dimnames = list(cells, cells))
for (i in seq_along(cells)) for (j in seq_along(cells)) {
  if (i < j) {
    p <- suppressWarnings(cor.test(es[i, ], es[j, ], method = "spearman",
                                   exact = FALSE)$p.value)
    cmp[i, j] <- p; cmp[j, i] <- p
  }
}
diag(cmp) <- 0
write.csv(data.frame(cell_type = rownames(cmp), cmp, check.names = FALSE),
          file.path(DIR_T, "immune_cell_correlation_pvalue.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("[correlation] immune cell pairs FDR<0.05:",
    sum(p.adjust(cmp[upper.tri(cmp)], "BH") < 0.05), "of", sum(upper.tri(cmp)), "\n")

## ---------------------------------------------------------------------------
## 5. Save
## ---------------------------------------------------------------------------
saveRDS(list(kw = kw, pw = pw, cor_tab = cor_tab, cm = cm, cmp = cmp,
             norm = norm_tab, var = var_tab),
        file.path(DIR_T, "immune_stats.rds"))

writeLines(c(
  "GSE57338 immune infiltration - statistical testing summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Distribution assessment:",
  paste0("  Shapiro-Wilk (per cell type x group): ",
         sum(!n_nonnormal$all_normal), "/", nrow(n_nonnormal),
         " cell types had >=1 non-normal group (P<0.05)"),
  paste0("  Levene (median-centred, 3 groups): ",
         sum(var_tab$levene_P < 0.05, na.rm = TRUE), "/", nrow(var_tab),
         " cell types heteroscedastic (P<0.05)"),
  "  => Non-parametric pipeline selected.",
  "",
  "Tests:",
  "  Global: Kruskal-Wallis rank-sum test across NF/DCM/ICM; epsilon^2 effect size",
  "  Post-hoc: pairwise Wilcoxon rank-sum (normal approximation, continuity correction)",
  "  Effect size: Cliff's delta (>0 = first group of the contrast is higher;",
  "               <0 = second group is higher)",
  "  Multiplicity: Benjamini-Hochberg FDR within each contrast (28 tests);",
  "                FDR_all additionally across all 84 pairwise tests",
  "  Correlation: Spearman with BH FDR within each gene",
  "",
  paste0("Significant (FDR<0.05): global ", sum(kw$FDR < 0.05), "/", nrow(kw),
         "; pairwise ", nrow(sig), "/", nrow(pw)),
  paste0("Key gene x immune cell: ", sum(cor_tab$FDR < 0.05), "/", nrow(cor_tab),
         " pairs FDR<0.05")
), file.path(DIR_T, "immune_stats_summary.txt"))

cat("\n=== 03_immune_stats.R finished ===\n")
sink()
