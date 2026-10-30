###############################################################################
## 96a_build_external_matrices.R
## Build analysis-ready, symbol-level expression matrices for the two external
## cohorts so that the 70-gene myeloid axis can be tested outside GSE57338.
##
## Inputs  (raw_data/external)
##   GSE5406_series_matrix.txt.gz   GPL96 / HG-U133A / RMA log2
##   GPL96.annot.gz                 platform annotation
##   GSE116250_rpkm.txt.gz          RNA-seq RPKM (transcript level)
##
## Outputs (results/AxisExternal/04_data)
##   expr_GSE5406_symbol.rds        genes x samples, log2 RMA
##   expr_GSE116250_symbol.rds      genes x samples, log2(RPKM+1)
##   sample_group_GSE5406.csv
##   sample_group_GSE116250.csv
##
## Run: Rscript scripts/96a_build_external_matrices.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

proj <- "E:/another data file/GEObulk analysis project"
raw  <- file.path(proj, "raw_data", "external")
outr <- file.path(proj, "results", "AxisExternal")
d_data <- file.path(outr, "04_data"); d_log <- file.path(outr, "03_logs")
for (d in c(d_data, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

con <- file(file.path(d_log, "96a_build_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); writeLines(s, con); cat(s, "\n"); flush.console() }
hr  <- function() say(strrep("-", 70))

###############################################################################
## helpers
###############################################################################
read_series_matrix <- function(gz) {
  con2 <- gzfile(gz, "r"); ln <- readLines(con2, warn = FALSE); close(con2)
  bidx <- grep("^!series_matrix_table_begin", ln); eidx <- grep("^!series_matrix_table_end", ln)
  if (!length(bidx)) stop("no table block in ", gz)
  if (!length(eidx)) eidx <- length(ln) + 1
  body <- ln[(bidx + 1):(eidx - 1)]; body <- body[nzchar(body)]
  read.delim(text = body, header = TRUE, check.names = FALSE, quote = "\"",
             stringsAsFactors = FALSE)
}

## collapse multiple probes/transcripts per symbol: keep the one with the
## highest mean expression (same rule as 80_ext_parse.R, so numbers stay
## comparable with the published external validation)
collapse_to_symbol <- function(mat, sym) {
  keep <- !is.na(sym) & sym != "" & sym != "---"
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  mu  <- rowMeans(mat, na.rm = TRUE)
  ord <- order(sym, -mu)                 # within a symbol, best probe first
  sym_s <- sym[ord]; mat_s <- mat[ord, , drop = FALSE]
  first <- !duplicated(sym_s)
  out <- mat_s[first, , drop = FALSE]
  rownames(out) <- sym_s[first]
  list(expr = out, n_in = length(unique(sym)), n_probe_dropped = sum(!first))
}

###############################################################################
## 1. GSE5406
###############################################################################
hr(); say("[1] GSE5406 (GPL96, RMA log2)")
m <- read_series_matrix(file.path(raw, "GSE5406_series_matrix.txt.gz"))
say("  raw: ", nrow(m), " probes x ", ncol(m) - 1, " samples")
expr <- as.matrix(m[, -1, drop = FALSE]); rownames(expr) <- m[[1]]
storage.mode(expr) <- "double"
say("  value range: ", sprintf("%.2f - %.2f", min(expr, na.rm = TRUE), max(expr, na.rm = TRUE)))

ann <- fread(file.path(raw, "GPL96.annot.gz"), sep = "\t", header = TRUE,
             comment.char = "", quote = "", skip = "ID\t", data.table = FALSE)
sym_col <- grep("^Gene symbol$", colnames(ann), value = TRUE)[1]
ps <- setNames(ann[[sym_col]], ann$ID)[rownames(expr)]
ps <- sub(" ?///.*$", "", ps)                      # multi-symbol -> first

cg <- collapse_to_symbol(expr, ps)
e5406 <- cg$expr
say("  probes with symbol: ", sum(!is.na(ps) & ps != ""), " -> ",
    nrow(e5406), " unique symbols (", cg$n_probe_dropped, " redundant probes dropped)")

## drop genes that are flat / undetected: interquartile range == 0 adds nothing
## to a rank-based score and can distort ssGSEA ties
iqr <- apply(e5406, 1, IQR, na.rm = TRUE)
say("  IQR==0 genes: ", sum(iqr == 0), " ; genes with any NA: ", sum(rowSums(is.na(e5406)) > 0))
e5406 <- e5406[iqr > 0 & rowSums(is.na(e5406)) == 0, , drop = FALSE]
say("  retained: ", nrow(e5406), " genes x ", ncol(e5406), " samples")

info5406 <- fread(file.path(proj, "results", "ExternalValidation", "01_data",
                            "GSE5406_sample_info.csv"))
info5406 <- info5406[match(colnames(e5406), sample)]
say("  group counts: ", paste(names(table(info5406$group)), table(info5406$group),
                              sep = "=", collapse = ", "))
stopifnot(!anyNA(info5406$sample))

saveRDS(e5406, file.path(d_data, "expr_GSE5406_symbol.rds"))
fwrite(info5406[, .(sample, group)], file.path(d_data, "sample_group_GSE5406.csv"))

###############################################################################
## 2. GSE116250
###############################################################################
hr(); say("[2] GSE116250 (RNA-seq RPKM, transcript level)")
rpkm <- fread(file.path(raw, "GSE116250_rpkm.txt.gz"), sep = "\t", header = TRUE,
              data.table = FALSE, check.names = FALSE)
smp <- setdiff(colnames(rpkm), c("Gene", "Common_name"))
say("  raw: ", nrow(rpkm), " transcripts x ", length(smp), " samples")
rp <- as.matrix(rpkm[, smp, drop = FALSE]); storage.mode(rp) <- "double"
say("  raw range: ", sprintf("%.2f - %.2f", min(rp, na.rm = TRUE), max(rp, na.rm = TRUE)))

## sentinel / overflow values (1e9) are not real expression; treat as missing
n_sentinel <- sum(rp >= 1e6, na.rm = TRUE)
say("  sentinel values (>=1e6): ", n_sentinel, " in ",
    length(unique(which(rp >= 1e6, arr.ind = TRUE)[, 1])), " transcripts")
rp[rp >= 1e6] <- NA

sym11 <- rpkm$Common_name
cg11 <- collapse_to_symbol(rp, sym11)
e11 <- cg11$expr
say("  transcripts with symbol: ", sum(!is.na(sym11) & sym11 != ""), " -> ",
    nrow(e11), " unique symbols (", cg11$n_probe_dropped, " redundant dropped)")

## genes that are mostly sentinel/missing are unusable
frac_na <- rowMeans(is.na(e11))
say("  genes with >50% missing: ", sum(frac_na > 0.5))
e11 <- e11[frac_na <= 0.5, , drop = FALSE]
if (anyNA(e11)) {                                   # residual gaps -> row median
  for (i in which(rowSums(is.na(e11)) > 0)) {
    e11[i, is.na(e11[i, ])] <- median(e11[i, ], na.rm = TRUE)
  }
}
say("  after gap filling, remaining NA: ", sum(is.na(e11)))

## expression filter BEFORE log: a rank-based score is dominated by the
## thousands of all-zero transcripts, so keep genes expressed (>1 RPKM) in at
## least 20% of samples -- the standard RNA-seq 'expressed gene' cut
expressed <- rowMeans(e11 > 1) >= 0.2
say("  expressed in >=20% samples: ", sum(expressed), " / ", nrow(e11))
e11 <- e11[expressed, , drop = FALSE]

e11 <- log2(e11 + 1)
iqr11 <- apply(e11, 1, IQR, na.rm = TRUE)
say("  IQR==0 genes: ", sum(iqr11 == 0))
e11 <- e11[iqr11 > 0, , drop = FALSE]
say("  retained: ", nrow(e11), " genes x ", ncol(e11), " samples")
say("  log2 range: ", sprintf("%.2f - %.2f", min(e11), max(e11)))

info11 <- fread(file.path(proj, "results", "ExternalValidation", "01_data",
                          "GSE116250_sample_info.csv"))
## RPKM column order vs series-matrix order can differ -> match by title
info11 <- info11[match(colnames(e11), title)]
say("  group counts: ", paste(names(table(info11$group)), table(info11$group),
                              sep = "=", collapse = ", "))
stopifnot(!anyNA(info11$sample))

saveRDS(e11, file.path(d_data, "expr_GSE116250_symbol.rds"))
fwrite(info11[, .(sample, title, group, age, sex)],
       file.path(d_data, "sample_group_GSE116250.csv"))

###############################################################################
## 3. coverage of the 70-gene axis
###############################################################################
hr(); say("[3] 70-gene coverage")
cand <- fread(file.path(proj, "results", "HubGene", "GSE57338_candidate_genes.csv"))
g70  <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
say("  candidate70 n = ", length(g70))
cov5406 <- g70[g70 %in% rownames(e5406)]
cov11   <- g70[g70 %in% rownames(e11)]
say("  GSE5406  : ", length(cov5406), "/", length(g70),
    "  missing: ", paste(setdiff(g70, cov5406), collapse = ","))
say("  GSE116250: ", length(cov11), "/", length(g70),
    "  missing: ", paste(setdiff(g70, cov11), collapse = ","))
say("  covered by both: ", length(intersect(cov5406, cov11)))

hr(); say("DONE"); close(con)
cat("BUILD DONE\n")
