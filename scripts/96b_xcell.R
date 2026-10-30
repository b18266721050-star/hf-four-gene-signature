###############################################################################
## 96b_xcell.R  -- run xCell on one external cohort
##
## Usage: Rscript scripts/96b_xcell.R <cohort> <rnaseq>
##   cohort  : GSE5406 | GSE116250
##   rnaseq  : TRUE | FALSE   (selects xCell's RNA-seq vs microarray
##             spillover calibration: xCell.data$spill vs $spill.array)
##
## Output: results/AxisExternal/04_data/xcell_<cohort>_<rna|array>.csv
##
## Note: xCell scores are RELATIVE enrichment scores on a cohort-internal
## scale. They are never compared across cohorts in absolute terms -- only the
## rank correlation structure and the direction of group differences are.
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(xCell) })

a <- commandArgs(trailingOnly = TRUE)
cohort <- a[1]; rnaseq <- as.logical(a[2])
stopifnot(cohort %in% c("GSE5406", "GSE116250"), !is.na(rnaseq))

proj <- "E:/another data file/GEObulk analysis project"
outr <- file.path(proj, "results", "AxisExternal")
d_data <- file.path(outr, "04_data"); d_log <- file.path(outr, "03_logs")
dir.create(d_log, recursive = TRUE, showWarnings = FALSE)

tag <- if (rnaseq) "rna" else "array"
logf <- file.path(d_log, sprintf("96b_xcell_%s_%s.txt", cohort, tag))
con <- file(logf, open = "wt")
say <- function(...) { s <- paste0(...); writeLines(s, con); cat(s, "\n"); flush.console() }

say("cohort = ", cohort, " ; rnaseq = ", rnaseq,
    " ; spill = ", ifelse(rnaseq, "xCell.data$spill", "xCell.data$spill.array"))

expr <- readRDS(file.path(d_data, sprintf("expr_%s_symbol.rds", cohort)))
say("input: ", nrow(expr), " genes x ", ncol(expr), " samples; range ",
    sprintf("%.2f-%.2f", min(expr), max(expr)))
say("duplicated symbols: ", sum(duplicated(rownames(expr))))
say("NA cells: ", sum(is.na(expr)))
if (anyNA(expr)) stop("expression matrix still contains NA")

t0 <- Sys.time()
xc <- xCellAnalysis(expr, rnaseq = rnaseq, parallel.sz = 4, parallel.type = "SOCK")
el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
say("xCell done in ", el, " s : ", nrow(xc), " cell types/scores x ", ncol(xc), " samples")
say("Macrophages present: ", "Macrophages" %in% rownames(xc))
say("Macrophages range: ", sprintf("%.3f-%.3f", min(xc["Macrophages", ]), max(xc["Macrophages", ])))

out <- file.path(d_data, sprintf("xcell_%s_%s.csv", cohort, tag))
dt <- as.data.table(as.data.frame(xc), keep.rownames = "cell_type")
fwrite(dt, out)
say("written: ", basename(out))

say("DONE"); close(con)
cat("XCELL DONE ", cohort, " ", tag, "\n")
