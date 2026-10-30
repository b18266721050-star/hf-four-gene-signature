###############################################################################
## 96e_supp_stats.R
## Why does the raw rho differ so much between cohorts? Quantify how much of
## each score's variance is between-group (disease vs NF) rather than within.
## A cohort dominated by one group (GSE5406: 7.6% NF) simply cannot produce the
## strong "everything moves together" correlation that a balanced cohort does.
##
## Also: per-cohort top/bottom genes, and the 52 genes shared by all cohorts.
##
## Run: Rscript scripts/96e_supp_stats.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

proj <- "E:/another data file/GEObulk analysis project"
outr <- file.path(proj, "results", "AxisExternal")
d_tab <- file.path(outr, "01_tables"); d_dat <- file.path(outr, "04_data")

say <- function(...) cat(paste0(...), "\n")

grp <- fread(file.path(d_tab, "T04_axis_score_by_sample.csv"))
rho <- fread(file.path(d_tab, "T03_gene_level_rho.csv"))
summ <- fread(file.path(d_tab, "T01_axis_replication_summary.csv"))

eta2 <- function(y, g) {
  y <- suppressWarnings(as.numeric(unlist(y)))
  g <- as.character(unlist(g))
  ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]
  grand <- mean(y)
  gm <- tapply(y, g, mean)
  gn <- tapply(y, g, length)
  ssb <- sum(gn * (gm - grand)^2)
  sst <- sum((y - grand)^2)
  as.numeric(ssb / sst)
}

say("== between-group share of variance (eta squared) ==")
res <- rbindlist(lapply(unique(grp$cohort), function(co) {
  d <- grp[cohort == co]
  e <- eta2(d$axis, d$group)
  nf <- mean(d$group == "NF")
  data.table(cohort = co, n = nrow(d), nf_pct = round(100 * nf, 1),
             eta2_axis = round(e, 3),
             n_groups = length(unique(d$group)))
}))
print(res)

## macrophage score eta2 per cohort
mac_eta <- rbindlist(lapply(c("GSE5406_array", "GSE116250_rna"), function(tag) {
  co <- sub("_(array|rna)$", "", tag)
  xc <- fread(file.path(d_dat, sprintf("xcell_%s.csv", tag)))
  m <- as.matrix(xc[, -1]); rownames(m) <- xc$cell_type
  g <- grp[cohort == co]
  g <- g[match(colnames(m), g$sample)]
  data.table(cohort = co, eta2_mac = round(eta2(m["Macrophages", ], g$group), 3))
}))
## discovery macrophage eta2
comp <- fread(file.path(proj, "results", "Composition", "01_tables",
                        "T04_composition_merged.csv"))
mac_eta <- rbind(data.table(cohort = "GSE57338",
                            eta2_mac = round(eta2(comp$xcell_Macrophages, comp$group), 3)),
                 mac_eta, fill = TRUE)
say(""); say("== eta2 of the Macrophages score ==")
print(mac_eta)

merged <- merge(res, mac_eta, by = "cohort", all.x = TRUE)
merged <- merge(merged, summ[, .(cohort, rho_mac_raw, rho_mac_within)], by = "cohort")
say(""); say("== combined ==")
print(merged)
fwrite(merged, file.path(d_tab, "T06_group_variance_share.csv"))

say(""); say("== top / bottom genes by rho with Macrophages ==")
for (co in c("GSE57338", "GSE5406", "GSE116250")) {
  d <- rho[cohort == co][order(-rho_mac_raw)]
  say("-- ", co)
  say("   top5   : ", paste(sprintf("%s(%.2f)", head(d$gene, 5), head(d$rho_mac_raw, 5)), collapse = ", "))
  say("   bottom5: ", paste(sprintf("%s(%.2f)", tail(d$gene, 5), tail(d$rho_mac_raw, 5)), collapse = ", "))
}

say(""); say("== genes measured in all three cohorts ==")
common <- Reduce(intersect, lapply(c("GSE57338", "GSE5406", "GSE116250"),
                                   function(co) rho[cohort == co, gene]))
say("n = ", length(common))
wide <- dcast(rho[gene %in% common], gene ~ cohort, value.var = "rho_mac_raw")
wide[, min_rho := pmin(GSE57338, GSE5406, GSE116250)]
wide[, all_positive := GSE57338 > 0 & GSE5406 > 0 & GSE116250 > 0]
say("genes positive in all three: ", sum(wide$all_positive), "/", nrow(wide))
say("genes above 0.3 in all three: ", sum(wide$min_rho > 0.3))
say("  -> ", paste(wide[min_rho > 0.3, gene], collapse = ", "))
fwrite(wide[order(-min_rho)], file.path(d_tab, "T07_common_genes_rho.csv"))
say("")
say("DONE")
