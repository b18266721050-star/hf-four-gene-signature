## Rebuild cluster composition metrics from the saved contingency tables.
## (the in-pipeline version silently expanded rows for the pre-integration block)
proj    <- "E:/another data file/GEObulk analysis project"
tab_dir <- file.path(proj, "results", "Singlecel", "tables")

man <- read.csv(file.path(proj, "clean_data", "GSE292067_sample_manifest.csv"),
                stringsAsFactors = FALSE)
s2d <- setNames(as.character(man$subject_id), man$sample)
N_SAMPLES <- length(unique(man$sample))

read_ct <- function(f) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE])
  rownames(m) <- as.character(d[[1]])
  m
}

metrics <- function(ct_s, ct_c, label) {
  stopifnot(nrow(ct_s) == nrow(ct_c))
  tot <- rowSums(ct_s)
  out <- do.call(rbind, lapply(seq_len(nrow(ct_s)), function(i) {
    p <- ct_s[i, ] / tot[i]
    p <- p[p > 0]
    smp <- colnames(ct_s)[ct_s[i, ] > 0]
    data.frame(resolution = label,
               cluster = rownames(ct_s)[i],
               n_cells = as.integer(tot[i]),
               n_samples = length(smp),
               n_donors = length(unique(s2d[smp])),
               max_sample_pct = round(max(p) * 100, 2),
               sample_entropy_norm = round(-sum(p * log(p)) / log(N_SAMPLES), 4),
               pct_disease = round(ct_c[i, "Disease"] / tot[i] * 100, 2),
               stringsAsFactors = FALSE)
  }))
  out[order(as.integer(out$cluster)), ]
}

p_s <- read_ct(file.path(tab_dir, "02_preint_cluster_by_sample.csv"))
p_s <- p_s[, setdiff(colnames(p_s), "n_cells"), drop = FALSE]
p_c <- read_ct(file.path(tab_dir, "02_preint_cluster_by_condition.csv"))

h_s <- read_ct(file.path(tab_dir, "03_cluster_by_sample_counts.csv"))
h_c <- read_ct(file.path(tab_dir, "03_cluster_by_condition_counts.csv"))

pre  <- metrics(p_s, p_c, "preint_0.5")
post <- metrics(h_s, h_c, "harmony_0.6")

stopifnot(sum(pre$n_cells) == sum(post$n_cells))
out <- rbind(pre, post)
write.csv(out, file.path(tab_dir, "03_cluster_composition_metrics.csv"), row.names = FALSE)

cat("total cells:", sum(pre$n_cells), "\n\n")
cat("== pre-integration summary ==\n")
cat("  clusters:", nrow(pre), " median entropy:", median(pre$sample_entropy_norm),
    " median max_sample_pct:", median(pre$max_sample_pct),
    " disease% range:", paste(range(pre$pct_disease), collapse = "-"), "\n")
cat("== integrated summary ==\n")
cat("  clusters:", nrow(post), " median entropy:", median(post$sample_entropy_norm),
    " median max_sample_pct:", median(post$max_sample_pct),
    " disease% range:", paste(range(post$pct_disease), collapse = "-"), "\n")
cat("\nclusters with <50% disease after integration:",
    sum(post$pct_disease < 50), "of", nrow(post), "\n")
