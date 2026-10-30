###############################################################################
## 90_composition_deconv.R
## GSE57338 -- bulk cell-composition deconvolution (three angles of the same
## question: is the "everything myeloid goes down" pattern a composition effect?)
##
## Why
##   The 70 WGCNA x DEG candidate genes are 100% down-regulated and the PPI hubs
##   (CD163, VSIG4, S100A9, FPR1, STAT3) are all down as well, which runs against
##   the mainstream view that inflammation is increased in heart failure. Two
##   explanations remain: (a) the failing samples simply contain fewer immune
##   cells (composition), or (b) immune cells are present but transcribe less
##   (per-cell change, or a donor-versus-recipient / RNA-quality artefact).
##   Deconvolution separates (a) from (b).
##
## Methods
##   xCell 1.1.0   ssGSEA enrichment scores, 64 cell types + 3 composite scores,
##                 spillover-compensated. Primary read-out (microarray-calibrated).
##   EPIC          constrained least squares, absolute fractions (7 compartments
##                 + uncharacterised). Sensitivity read-out for the main
##                 compartments (myocyte fraction is NOT in EPIC's reference, so
##                 EPIC cannot answer the cardiomyocyte question -- stated openly).
##   MCPcounter    NOT used: the package downloads its signature from GitHub at
##                 run time and this machine has no outbound access.
##
## Inputs
##   clean_data/GSE57338_clean_expression_matrix.csv   (18,726 x 313, RMA log2)
##   results/Batch/01_tables/T01_sample_batch_metadata.csv
##
## Outputs (results/Composition/)
##   01_tables/T01_xcell_raw.csv
##   01_tables/T02_xcell_adjusted.csv
##   01_tables/T03_epic_fractions.csv
##   01_tables/T04_composition_merged.csv
##   03_logs/90_run_log.txt , sessionInfo_90.txt
##
## Run:  Rscript scripts/90_composition_deconv.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(xCell)
  library(EPIC)
})

options(stringsAsFactors = FALSE)
set.seed(20261029)

proj_dir <- "E:/another data file/GEObulk analysis project"
if (!dir.exists(proj_dir)) proj_dir <- getwd()
clean_dir <- file.path(proj_dir, "clean_data")
out_root  <- file.path(proj_dir, "results", "Composition")
d_tab     <- file.path(out_root, "01_tables")
d_fig     <- file.path(out_root, "02_figures")
d_log     <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...)
  run_log <<- c(run_log, txt)
  cat(txt, "\n", sep = "")
  flush.console()
}
hr <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

## ===========================================================================
## 1. Load
## ===========================================================================
hr("=")
log_msg("GSE57338 cell-composition deconvolution")
hr("=")
log_msg("[1] Loading expression matrix ...")

expr_dt <- data.table::fread(file.path(clean_dir, "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE])
rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())
log_msg("  ", nrow(expr), " genes x ", ncol(expr), " samples; range ",
        sprintf("%.2f-%.2f", min(expr), max(expr)))
log_msg("  duplicated symbols: ", sum(duplicated(rownames(expr))))

bmeta <- data.table::fread(file.path(proj_dir, "results", "Batch", "01_tables",
                                     "T01_sample_batch_metadata.csv"))
bmeta <- bmeta[match(colnames(expr), sample)]
if (anyNA(bmeta$sample)) stop("metadata does not cover all samples")
bmeta[, group := factor(group, levels = c("NF", "DCM", "ICM"))]
log_msg("  groups: ", paste(names(table(bmeta$group)), table(bmeta$group),
                            sep = "=", collapse = ", "))

## ===========================================================================
## 2. xCell (primary)
## ===========================================================================
hr()
log_msg("[2] xCell (ssGSEA, spillover-compensated) ...")
t0 <- Sys.time()
xcell_raw <- xCellAnalysis(expr, rnaseq = FALSE, parallel.sz = 4,
                           parallel.type = "SOCK")
log_msg("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " s; ", nrow(xcell_raw), " cell types / scores x ", ncol(xcell_raw), " samples")

## xCell returns cell types in rows and samples in columns; the values are
## spillover-compensated (alpha = 0.5, package default) ssGSEA scores.
## These are enrichment scores on a relative scale, NOT absolute fractions --
## all downstream wording has to respect that (see the report, Limitations).
dt_x <- as.data.table(as.data.frame(xcell_raw), keep.rownames = "cell_type")
write_csv(dt_x, file.path(d_tab, "T02_xcell_scores.csv"))

comp_scores <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore")
cell_types  <- setdiff(rownames(xcell_raw), comp_scores)
log_msg("  cell types (", length(cell_types), "): ",
        paste(head(cell_types, 12), collapse = ", "), ", ...")
log_msg("  composite scores: ", paste(comp_scores, collapse = ", "))
log_msg("  key compartments present: ",
        paste(intersect(c("Myocytes", "Fibroblasts", "Endothelial cells", "Macrophages",
                          "Macrophages M1", "Macrophages M2", "Monocytes", "Pericytes",
                          "Smooth muscle", "CD8+ T-cells", "CD4+ T-cells", "B-cells",
                          "Plasma cells", "Mast cells", "Neutrophils", "pDC", "cDC"),
                        rownames(xcell_raw)), collapse = ", "))

## ===========================================================================
## 3. EPIC (sensitivity, absolute fractions)
## ===========================================================================
hr()
log_msg("[3] EPIC (absolute fractions) ...")
## EPIC expects TPM / linear scale -> back-transform the RMA log2 values
t0 <- Sys.time()
epic <- tryCatch(
  EPIC(2^expr, withOtherCells = TRUE),
  error = function(e) { log_msg("  EPIC failed: ", conditionMessage(e)); NULL }
)
if (!is.null(epic)) {
  ## EPIC returns samples x compartments
  frac <- epic$cellFractions
  log_msg("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
          " s; ", ncol(frac), " compartments: ", paste(colnames(frac), collapse = ", "))
  log_msg("  fraction sum per sample: ",
          sprintf("%.3f-%.3f", min(rowSums(frac)), max(rowSums(frac))))
  dt_e <- as.data.table(as.data.frame(frac), keep.rownames = "sample")
  write_csv(dt_e, file.path(d_tab, "T03_epic_fractions.csv"))
}

## ===========================================================================
## 4. Merge into one analysis-ready table
## ===========================================================================
hr()
log_msg("[4] Merging ...")
comp <- data.table(sample = colnames(xcell_raw))
comp <- cbind(comp, as.data.table(t(xcell_raw)))
setnames(comp, old = names(comp)[-1], new = paste0("xcell_", gsub("[^A-Za-z0-9]+", "_",
                                                                  names(comp)[-1])))
if (!is.null(epic)) {
  epi_dt <- as.data.table(as.data.frame(frac), keep.rownames = "sample")
  setnames(epi_dt, old = names(epi_dt)[-1], new = paste0("epic_",
                                                         gsub("[^A-Za-z0-9]+", "_",
                                                              names(epi_dt)[-1])))
  comp <- merge(comp, epi_dt, by = "sample", all.x = TRUE)
}
comp <- merge(bmeta[, .(sample, patient, gsm_order, group, group2, sex, age_num,
                        batch_prefix, batch_order, batch_idnum)],
              comp, by = "sample")
comp <- comp[order(gsm_order)]
log_msg("  merged table: ", nrow(comp), " samples x ", ncol(comp), " columns")
write_csv(comp, file.path(d_tab, "T04_composition_merged.csv"))

hr("=")
log_msg("Done. Tables: ", d_tab)
hr("=")

writeLines(run_log, file.path(d_log, "90_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_90.txt"))
