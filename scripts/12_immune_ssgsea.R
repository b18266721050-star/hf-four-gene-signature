###############################################################################
# 02_immune_ssgsea.R
# Immune infiltration estimation by ssGSEA (GSVA)
# Gene set: Charoentong et al. 2017, Cell Reports 18:248-262 (Table S6)
#           28 human immune cell types, 782 marker genes
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(GSEABase)
  library(parallel)
})

ROOT    <- "E:/another data file/GEObulk analysis project"
OUT     <- file.path(ROOT, "results", "Immune")
DIR_QC  <- file.path(OUT, "01_qc")
DIR_S   <- file.path(OUT, "02_ssgsea")
DIR_REF <- file.path(OUT, "05_reference")
DIR_LOG <- file.path(OUT, "00_logs")
for (d in c(DIR_S, DIR_REF, DIR_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

REF_FILE <- "E:/another data file/_tmp/gsref/Charoentong2017_tableS6_Immune.txt"

sink(file.path(DIR_LOG, "02_ssgsea_console.log"), split = TRUE)
cat("=== 02_immune_ssgsea.R started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n\n")

## ---------------------------------------------------------------------------
## 1. Load QC objects
## ---------------------------------------------------------------------------
obj  <- readRDS(file.path(DIR_QC, "qc_input.rds"))
expr <- obj$expr
si   <- obj$si
cat("[input] expression:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("[input] groups:", paste(names(table(si$group)), as.integer(table(si$group)),
                              sep = "=", collapse = "; "), "\n\n")

## ---------------------------------------------------------------------------
## 2. Build gene sets
## ---------------------------------------------------------------------------
ct <- fread(REF_FILE)
setnames(ct, c("gene", "cell_type", "immunity_class"))

## clean cosmetic issues present in the published table
ct$cell_type <- gsub("\\s+", " ", ct$cell_type)
ct$cell_type <- gsub("memeory", "memory", ct$cell_type)          # typo in source table
ct$cell_type <- trimws(ct$cell_type)

ct <- unique(ct[, .(gene, cell_type)])
ct <- ct[gene != "" & !is.na(gene)]

gs_list <- split(ct$gene, ct$cell_type)
cat("[geneset] cell types:", length(gs_list), " unique genes:", length(unique(ct$gene)), "\n")

## ---------------------------------------------------------------------------
## 3. Gene coverage QC (critical for ssGSEA validity)
## ---------------------------------------------------------------------------
cov <- data.frame(
  cell_type      = names(gs_list),
  n_marker       = lengths(gs_list),
  n_matched      = sapply(gs_list, function(g) sum(g %in% rownames(expr))),
  stringsAsFactors = FALSE
)
cov$match_rate    <- round(cov$n_matched / cov$n_marker, 4)
cov$missing_genes <- sapply(gs_list, function(g) paste(setdiff(g, rownames(expr)), collapse = ";"))
cov$status        <- ifelse(cov$n_matched >= 10, "PASS",
                     ifelse(cov$n_matched >= 5, "LOW_COVERAGE", "FAIL"))

write.csv(cov, file.path(DIR_S, "GSE57338_immune_gene_coverage.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n[coverage] per cell type:\n")
print(cov[, c("cell_type", "n_marker", "n_matched", "match_rate", "status")], row.names = FALSE)
cat("\n[coverage] overall matched:", sum(cov$n_matched), "/", sum(cov$n_marker),
    sprintf("(%.2f%%)\n", 100 * sum(cov$n_matched) / sum(cov$n_marker)))
cat("[coverage] status tally:", paste(names(table(cov$status)), table(cov$status),
                                      sep = "=", collapse = "; "), "\n")

## hard stop if any set is unusable
if (any(cov$status == "FAIL")) {
  sink(); stop("One or more immune gene sets have <5 matched genes; ssGSEA not reliable.")
}

## ---------------------------------------------------------------------------
## 4. Restrict gene sets to genes present in the matrix
## ---------------------------------------------------------------------------
gs_use <- lapply(gs_list, function(g) intersect(g, rownames(expr)))
gs_use <- gs_use[lengths(gs_use) > 0]

## export GMT for reproducibility
gmt_lines <- sapply(names(gs_use), function(n) {
  paste(c(n, "Charoentong2017_ssGSEA", gs_use[[n]]), collapse = "\t")
})
writeLines(gmt_lines, file.path(DIR_REF, "Charoentong2017_28immune_celltypes_matched.gmt"),
           useBytes = TRUE)
write.csv(cov, file.path(DIR_REF, "Charoentong2017_28immune_celltypes_manifest.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
file.copy(REF_FILE, file.path(DIR_REF, "Charoentong2017_tableS6_Immune_source.txt"),
          overwrite = TRUE)

## ---------------------------------------------------------------------------
## 5. ssGSEA
## ---------------------------------------------------------------------------
PARAM <- list(
  method      = "ssGSEA",
  kcdf        = "Gaussian",
  tau         = 1,
  maxDiff     = TRUE,
  absRanking  = FALSE,
  minSize     = 1,
  maxSize     = Inf,
  sparse      = FALSE,
  checkNA     = "yes",
  filterRows  = TRUE,
  n_cores     = max(1, parallel::detectCores() - 1)
)
cat("\n[ssGSEA] parameters:\n"); print(unlist(PARAM))

set.seed(20261027)
t0 <- Sys.time()
par <- gsvaParam(
  exprData  = expr,
  geneSets  = gs_use,
  minSize   = PARAM$minSize,
  maxSize   = PARAM$maxSize,
  kcdf      = PARAM$kcdf,
  tau       = PARAM$tau,
  maxDiff   = PARAM$maxDiff,
  absRanking= PARAM$absRanking,
  sparse    = PARAM$sparse,
  checkNA  = PARAM$checkNA,
  verbose   = TRUE
)
es <- gsva(par, verbose = TRUE)
t1 <- Sys.time()
cat("[ssGSEA] done in", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "min\n")
cat("[ssGSEA] output dim:", nrow(es), "x", ncol(es), "\n")
cat("[ssGSEA] score range:", round(min(es), 4), "~", round(max(es), 4), "\n")

## ---------------------------------------------------------------------------
## 6. Export matrices
## ---------------------------------------------------------------------------
es <- es[, colnames(expr), drop = FALSE]   # enforce sample order

res <- data.frame(cell_type = rownames(es), es, check.names = FALSE)
write.csv(res, file.path(DIR_S, "GSE57338_immune_ssGSEA_matrix_cell_by_sample.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## sample x cell (analysis-ready)
es_t <- as.data.frame(t(es), check.names = FALSE)
es_out <- data.frame(sample = rownames(es_t), es_t, check.names = FALSE)
es_out <- cbind(es_out, group = si$group[match(es_out$sample, si$sample)])
write.csv(es_out, file.path(DIR_S, "GSE57338_immune_ssGSEA_matrix_sample_by_cell.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## long format
long <- data.table(
  sample    = rep(colnames(es), each = nrow(es)),
  group     = rep(si$group[match(colnames(es), si$sample)], each = nrow(es)),
  cell_type = rep(rownames(es), times = ncol(es)),
  ssgsea    = as.vector(es)
)
long[, rank_within_sample := frank(-ssgsea, ties.method = "average"), by = sample]
setorder(long, sample, -ssgsea)
write.csv(long, file.path(DIR_S, "GSE57338_immune_ssgSEA_long.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## descriptive per cell type
desc <- long[, .(n = .N,
                 mean = mean(ssgsea), sd = sd(ssgsea),
                 median = median(ssgsea),
                 min = min(ssgsea), max = max(ssgsea)), by = cell_type]
setorder(desc, -mean)
write.csv(desc, file.path(DIR_S, "GSE57338_immune_ssgSEA_descriptive.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(es = es, gs_use = gs_use, cov = cov, param = PARAM),
        file.path(DIR_S, "GSE57338_immune_ssgSEA.rds"))

cat("\n[export] files written to:", DIR_S, "\n")
cat("=== 02_immune_ssgsea.R finished ===\n")
sink()
