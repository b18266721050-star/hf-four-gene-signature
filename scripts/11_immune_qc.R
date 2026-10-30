###############################################################################
# 01_qc.R  ---  Immune infiltration: input QC (requirements 2 & 3)
# Project: GSE57338 (human left ventricle, GPL11532)
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

ROOT   <- "E:/another data file/GEObulk analysis project"
CLEAN  <- file.path(ROOT, "clean_data")
OUT    <- file.path(ROOT, "results", "Immune")
DIR_QC <- file.path(OUT, "01_qc")
DIR_LOG<- file.path(OUT, "00_logs")
dir.create(DIR_QC,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)

KEY_GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")

sink(file.path(DIR_LOG, "01_qc_console.log"), split = TRUE)
cat("=== 01_qc.R started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "===\n")
cat("R:", R.version.string, "\n")
cat("native encoding:", l10n_info()$`UTF-8`, " locale:", Sys.getlocale("LC_CTYPE"), "\n\n")

## ---------------------------------------------------------------------------
## 1. Read expression matrix
## ---------------------------------------------------------------------------
f_expr <- file.path(CLEAN, "GSE57338_clean_expression_matrix.csv")
expr_dt <- fread(f_expr)
cat("[read] expression matrix:", nrow(expr_dt), "rows x", ncol(expr_dt), "cols\n")

gene_col <- names(expr_dt)[1]
cat("[read] first column name:", gene_col, "\n")
cat("[read] first 3 gene ids:", paste(head(expr_dt[[1]], 3), collapse = ", "), "\n")
cat("[read] last  3 gene ids:", paste(tail(expr_dt[[1]], 3), collapse = ", "), "\n")

expr <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr) <- expr_dt[[1]]
samples <- colnames(expr)

## ---------------------------------------------------------------------------
## 2. Structure / orientation checks
## ---------------------------------------------------------------------------
qc <- list()
qc$n_gene_row    <- nrow(expr)
qc$n_sample_col  <- ncol(expr)
qc$gene_id_type  <- "HGNC gene symbol (header: GeneSymbol)"
qc$row_is_gene   <- all(grepl("^[A-Za-z][A-Za-z0-9.@_-]*$", head(rownames(expr), 500)))
qc$col_is_sample <- all(grepl("^GSM[0-9]+$", samples))
qc$matrix_numeric<- is.numeric(expr)
qc$storage_mode  <- storage.mode(expr)

## orientation sanity: number of genes >> number of samples
qc$orientation   <- if (nrow(expr) > ncol(expr)) "genes x samples (correct)" else "TRANSPOSED?"

## ---------------------------------------------------------------------------
## 3. Missing / invalid values
## ---------------------------------------------------------------------------
qc$n_na      <- sum(is.na(expr))
qc$n_nan     <- sum(is.nan(expr))
qc$n_inf     <- sum(is.infinite(expr))
qc$min_value <- min(expr, na.rm = TRUE)
qc$max_value <- max(expr, na.rm = TRUE)
qc$n_negative<- sum(expr < 0, na.rm = TRUE)
qc$n_zero    <- sum(expr == 0, na.rm = TRUE)

## integer / count check
is_int <- all(abs(expr[!is.na(expr)] - round(expr[!is.na(expr)])) < 1e-8)
qc$is_integer_like <- is_int

## log2 scale diagnosis
q <- quantile(expr, c(0, .01, .25, .5, .75, .99, 1), na.rm = TRUE)
qc$q000 <- q[1]; qc$q025 <- q[3]; qc$q500 <- q[4]; qc$q975 <- q[6]; qc$q100 <- q[7]

log2_like <- (qc$max_value < 30) && (qc$q500 > 1) && (qc$q500 < 15) && !is_int
qc$log2_like <- log2_like
qc$scale_call <- if (is_int && qc$max_value > 100) {
  "raw count / linear intensity"
} else if (log2_like) {
  "log2-transformed, continuous (microarray RMA)"
} else {
  "unclear linear scale"
}

## ---------------------------------------------------------------------------
## 4. Duplicate genes / zero variance
## ---------------------------------------------------------------------------
dup <- rownames(expr)[duplicated(rownames(expr))]
qc$n_dup_gene <- length(dup)
dup_tab <- if (length(dup)) {
  data.frame(gene = dup,
             n = as.integer(table(rownames(expr))[dup]),
             stringsAsFactors = FALSE)
} else {
  data.frame(gene = character(0), n = integer(0))
}
write.csv(dup_tab, file.path(DIR_QC, "duplicated_genes.csv"), row.names = FALSE)

v <- apply(expr, 1, var, na.rm = TRUE)
qc$n_zero_var <- sum(v == 0, na.rm = TRUE)
qc$n_low_var  <- sum(v < 1e-8, na.rm = TRUE)

## ---------------------------------------------------------------------------
## 5. Per-sample summary (distribution comparability)
## ---------------------------------------------------------------------------
samp_stat <- data.frame(
  sample    = samples,
  n_na      = colSums(is.na(expr)),
  mean      = colMeans(expr, na.rm = TRUE),
  median    = apply(expr, 2, median, na.rm = TRUE),
  sd        = apply(expr, 2, sd, na.rm = TRUE),
  min       = apply(expr, 2, min, na.rm = TRUE),
  max       = apply(expr, 2, max, na.rm = TRUE),
  q75_minus_q25 = apply(expr, 2, IQR, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(samp_stat, file.path(DIR_QC, "per_sample_distribution.csv"), row.names = FALSE)

qc$sample_median_range <- paste0(round(min(samp_stat$median), 3), " ~ ", round(max(samp_stat$median), 3))
qc$sample_iqr_range    <- paste0(round(min(samp_stat$q75_minus_q25), 3), " ~ ", round(max(samp_stat$q75_minus_q25), 3))

## ---------------------------------------------------------------------------
## 6. Sample info & matching
## ---------------------------------------------------------------------------
si <- fread(file.path(CLEAN, "GSE57338_sample_info.csv"))
cat("\n[read] sample_info:", nrow(si), "rows x", ncol(si), "cols\n")
cat("[read] columns:", paste(names(si), collapse = ", "), "\n")

qc$n_sample_info      <- nrow(si)
qc$sample_info_dup    <- sum(duplicated(si$sample))
qc$sample_info_na     <- sum(is.na(si$group)) + sum(si$group == "")
qc$group_table        <- paste(names(table(si$group)), table(si$group), sep = "=", collapse = "; ")

## exact match in both directions
m_in_matrix <- si$sample %in% samples
m_in_info   <- samples %in% si$sample
qc$n_unmatched_in_matrix <- sum(!m_in_matrix)
qc$n_unmatched_in_info   <- sum(!m_in_info)

unmatched <- data.frame(sample = character(0), present_in = character(0),
                        issue = character(0), stringsAsFactors = FALSE)
if (sum(!m_in_matrix) > 0) {
  unmatched <- rbind(unmatched, data.frame(
    sample = si$sample[!m_in_matrix], present_in = "sample_info",
    issue = "not found in expression matrix", stringsAsFactors = FALSE))
}
if (sum(!m_in_info) > 0) {
  unmatched <- rbind(unmatched, data.frame(
    sample = samples[!m_in_info], present_in = "expression matrix",
    issue = "not found in sample_info", stringsAsFactors = FALSE))
}
write.csv(unmatched, file.path(DIR_QC, "unmatched_samples.csv"), row.names = FALSE)

## order check
common_si   <- si$sample[si$sample %in% samples]
qc$order_identical <- identical(common_si, samples)
qc$set_identical   <- setequal(common_si, samples)

## case / whitespace look-alikes
norm_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
qc$n_case_only_mismatch <- sum(!(samples %in% si$sample) & (norm_id(samples) %in% norm_id(si$sample)))

## align (keep matrix order as reference)
si_m <- si[match(samples, si$sample), ]
stopifnot(all(si_m$sample == samples))
qc$align_ok <- TRUE

## group consistency with known totals
grp <- table(si_m$group)
qc$group_counts <- paste(names(grp), as.integer(grp), sep = "=", collapse = "; ")

## ---------------------------------------------------------------------------
## 7. Key gene presence
## ---------------------------------------------------------------------------
key_tab <- data.frame(
  gene      = KEY_GENES,
  in_matrix = KEY_GENES %in% rownames(expr),
  stringsAsFactors = FALSE
)
if (any(key_tab$in_matrix)) {
  key_tab$mean_expr <- NA_real_
  key_tab$median_expr <- NA_real_
  idx <- which(key_tab$in_matrix)
  key_tab$mean_expr[idx]   <- rowMeans(expr[key_tab$gene[idx], , drop = FALSE])
  key_tab$median_expr[idx] <- apply(expr[key_tab$gene[idx], , drop = FALSE], 1, median)
}
write.csv(key_tab, file.path(DIR_QC, "key_gene_presence.csv"), row.names = FALSE)
cat("\n[key genes]\n"); print(key_tab)

## alternative symbol hints for missing key genes
missing_key <- key_tab$gene[!key_tab$in_matrix]
if (length(missing_key)) {
  hint <- sapply(missing_key, function(g) {
    hit <- grep(g, rownames(expr), ignore.case = TRUE, value = TRUE)
    paste(head(hit, 8), collapse = "|")
  })
  write.csv(data.frame(gene = missing_key, fuzzy_hits = hint, stringsAsFactors = FALSE),
            file.path(DIR_QC, "key_gene_fuzzy_lookup.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------------------
## 8. Write QC tables
## ---------------------------------------------------------------------------
qc_df <- data.frame(
  item  = names(qc),
  value = sapply(qc, function(x) paste(as.character(x), collapse = " | ")),
  stringsAsFactors = FALSE
)
write.csv(qc_df, file.path(DIR_QC, "qc_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(expr = expr, si = si_m, qc = qc),
        file.path(DIR_QC, "qc_input.rds"))

cat("\n=== QC summary ===\n")
print(qc_df, row.names = FALSE)
cat("\n[order] sample_info reordered to matrix:", qc$order_identical,
    " set equal:", qc$set_identical, "\n")
cat("=== 01_qc.R finished ===\n")
sink()
