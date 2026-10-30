## recompute Cramer's V from the saved contingency tables (no big object needed)
tab_dir <- "E:/another data file/GEObulk analysis project/results/Singlecel/tables"
out <- file.path(tab_dir, "03_cramersV_before_after.csv")

cramers_v <- function(tb) {
  tb <- as.matrix(tb)
  tb <- tb[rowSums(tb) > 0, colSums(tb) > 0, drop = FALSE]
  if (min(dim(tb)) < 2) return(NA_real_)
  n <- sum(tb)
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  if (!is.finite(chi)) return(NA_real_)
  as.numeric(sqrt(chi / (n * (min(dim(tb)) - 1))))
}

read_ct <- function(f, drop_cols = NULL) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(drop_cols)) d <- d[, setdiff(colnames(d), drop_cols), drop = FALSE]
  m <- as.matrix(d[, -1, drop = FALSE])
  rownames(m) <- d[[1]]
  m
}

## pre-integration (res 0.5)
p_s <- read_ct(file.path(tab_dir, "02_preint_cluster_by_sample.csv"), drop_cols = "n_cells")
p_c <- read_ct(file.path(tab_dir, "02_preint_cluster_by_condition.csv"))

## integrated (res 0.6)
h_s <- read_ct(file.path(tab_dir, "03_cluster_by_sample_counts.csv"))
h_c <- read_ct(file.path(tab_dir, "03_cluster_by_condition_counts.csv"))
h_d <- read_ct(file.path(tab_dir, "03_cluster_by_donor_counts.csv"))

res <- data.frame(
  clustering    = c("preint_0.5", "harmony_0.6"),
  V_vs_sample   = c(cramers_v(p_s), cramers_v(h_s)),
  V_vs_donor    = c(NA_real_, cramers_v(h_d)),
  V_vs_condition= c(cramers_v(p_c), cramers_v(h_c)),
  n_clusters    = c(nrow(p_s), nrow(h_s)),
  stringsAsFactors = FALSE)

cat("pre-int table dims:", paste(dim(p_s), collapse = "x"),
    "| harmony table dims:", paste(dim(h_s), collapse = "x"), "\n")
print(res)
write.csv(res, out, row.names = FALSE)
cat("written:", out, "\n")
