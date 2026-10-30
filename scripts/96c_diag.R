## 96c_diag.R -- quick diagnostic before the main correlation run
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

proj <- "E:/another data file/GEObulk analysis project"
ax   <- file.path(proj, "results", "AxisExternal", "04_data")
comp <- file.path(proj, "results", "Composition", "01_tables")

say <- function(...) cat(paste0(...), "\n")

q <- c(0, .25, .5, .75, 1)
summ <- function(x, nm) {
  say(sprintf("%-28s n=%3d  zero=%3d (%.0f%%)  q: %s",
              nm, length(x), sum(x == 0, na.rm = TRUE),
              100 * mean(x == 0, na.rm = TRUE),
              paste(sprintf("%.3f", quantile(x, q, na.rm = TRUE)), collapse = " ")))
}

say("== discovery cohort (GSE57338) ==")
t04 <- fread(file.path(comp, "T04_composition_merged.csv"))
summ(t04$xcell_Macrophages, "Macrophages")
summ(t04$xcell_Macrophages_M1, "Macrophages M1")
summ(t04$xcell_Macrophages_M2, "Macrophages M2")
summ(t04$xcell_Monocytes, "Monocytes")
summ(t04$xcell_ImmuneScore, "ImmuneScore")
say("")

for (co in c("GSE5406_array", "GSE116250_rna")) {
  say("== ", co, " ==")
  xc <- fread(file.path(ax, sprintf("xcell_%s.csv", co)))
  m <- as.matrix(xc[, -1]); rownames(m) <- xc$cell_type
  summ(m["Macrophages", ], "Macrophages")
  summ(m["Macrophages M1", ], "Macrophages M1")
  summ(m["Macrophages M2", ], "Macrophages M2")
  summ(m["Monocytes", ], "Monocytes")
  summ(m["ImmuneScore", ], "ImmuneScore")
  say("")
}

say("== expression spread of a few axis genes ==")
cand <- fread(file.path(proj, "results", "HubGene", "GSE57338_candidate_genes.csv"))
g70 <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
for (co in c("GSE5406", "GSE116250")) {
  e <- readRDS(file.path(ax, sprintf("expr_%s_symbol.rds", co)))
  ok <- intersect(g70, rownames(e))
  sds <- apply(e[ok, , drop = FALSE], 1, sd)
  say(co, ": covered=", length(ok), "  median gene SD=", round(median(sds), 3),
      "  min SD=", round(min(sds), 3), "  n SD<0.05: ", sum(sds < 0.05))
}
