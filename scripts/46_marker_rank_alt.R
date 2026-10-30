## alternative marker ranking that punishes lowly-expressed genes:
## score = avg_log2FC * (pct.1 - pct.2), plus a "detected in >=50% of the cluster" view
proj <- "E:/another data file/GEObulk analysis project"
tab  <- file.path(proj, "results", "Singlecel", "tables")
d <- read.csv(file.path(tab, "04_markers_all.csv"), stringsAsFactors = FALSE)

d$score <- d$avg_log2FC * d$pct_diff
d <- d[order(d$cluster, -d$score), ]
top <- do.call(rbind, lapply(split(d, d$cluster), function(x) head(x, 10)))
write.csv(top, file.path(tab, "04_markers_top10_by_score.csv"), row.names = FALSE)

## one-line summary per cluster: top3 by score and top3 among genes detected in >=50% cells
s <- do.call(rbind, lapply(split(d, d$cluster), function(x) {
  hi <- x[x$pct.1 >= 0.5, ]
  data.frame(cluster = x$cluster[1],
             n_markers = nrow(x),
             top3_score = paste(head(x$gene_symbol, 3), collapse = "/"),
             top_lfc = round(x$avg_log2FC[1], 2),
             top_pct1 = round(x$pct.1[1], 2),
             top3_common = if (nrow(hi)) paste(head(hi$gene_symbol, 3), collapse = "/") else "",
             common_max_lfc = if (nrow(hi)) round(max(hi$avg_log2FC), 2) else NA,
             stringsAsFactors = FALSE)
}))
s <- s[order(as.integer(s$cluster)), ]
write.csv(s, file.path(tab, "04_marker_cluster_profile.csv"), row.names = FALSE)

print(s[, c("cluster", "n_markers", "top3_score", "top_lfc", "top_pct1", "top3_common")])
cat("\nrows written:", nrow(top), "\n")
