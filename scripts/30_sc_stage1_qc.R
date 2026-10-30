# =============================================================================
# 30_sc_stage1_qc.R
# Reproducible stage-1 Seurat QC for GSE292067 snRNA-seq.
#
# Modes:
#   --mode=validate
#   --mode=sample --sample=GSM8848306
#   --mode=finalize
#
# Stage 1 creates unfiltered per-sample Seurat objects, QC plots, risk summaries,
# and data-driven threshold recommendations. It never subsets cells or genes.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}
mode <- arg_value("--mode", "validate")
sample_arg <- arg_value("--sample", NA_character_)

proj_dir <- "E:/another data file/GEObulk analysis project"
raw_dir <- file.path(proj_dir, "raw_data", "GSE292067")
manifest_path <- file.path(proj_dir, "clean_data", "GSE292067_sample_manifest.csv")
out_dir <- file.path(proj_dir, "results", "Singlecel")
obj_dir <- file.path(out_dir, "objects")
table_dir <- file.path(out_dir, "tables")
per_sample_dir <- file.path(table_dir, "per_sample")
cell_qc_dir <- file.path(table_dir, "cell_qc")
violin_dir <- file.path(out_dir, "qc_plots", "violin")
scatter_dir <- file.path(out_dir, "qc_plots", "scatter")
log_dir <- file.path(out_dir, "logs")
report_path <- file.path(proj_dir, "report", "ac review.md")

required_dirs <- c(out_dir, obj_dir, table_dir, per_sample_dir, cell_qc_dir,
                   violin_dir, scatter_dir, log_dir, dirname(report_path))
for (d in required_dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path)
if (!dir.exists(raw_dir)) stop("Raw directory not found: ", raw_dir)

manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
required_manifest <- c("sample", "sample_title", "subject_id", "condition",
                       "tissue", "assay", "file_barcodes", "file_features",
                       "file_matrix", "file_doublet")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) stop("Manifest is missing columns: ", paste(missing_manifest, collapse = ", "))

sample_paths <- function(row) {
  list(
    barcodes = file.path(raw_dir, row$file_barcodes),
    features = file.path(raw_dir, row$file_features),
    matrix = file.path(raw_dir, row$file_matrix),
    doublet = file.path(raw_dir, row$file_doublet)
  )
}

read_features <- function(path) {
  x <- read.delim(gzfile(path), header = FALSE, sep = "\t", quote = "",
                  comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
  if (ncol(x) != 3L) stop("features file must contain exactly 3 columns: ", path)
  names(x) <- c("ensembl_id", "gene_symbol", "feature_type")
  x
}

read_barcodes <- function(path) {
  x <- read.delim(gzfile(path), header = FALSE, sep = "\t", quote = "",
                  comment.char = "", fill = TRUE, stringsAsFactors = FALSE)
  if (ncol(x) != 1L) stop("barcodes file must contain exactly 1 column: ", path)
  as.character(x[[1]])
}

read_doublet <- function(path) {
  x <- read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("barcodes", "score", "prediction")
  if (!all(required %in% names(x))) stop("doublet file missing required columns: ", path)
  x[, required, drop = FALSE]
}

read_mtx_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  if (!length(first)) stop("Empty matrix file: ", path)
  dims <- NULL
  for (i in seq_len(50L)) {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    if (!grepl("^%", ln) && nzchar(trimws(ln))) {
      vals <- suppressWarnings(as.numeric(strsplit(trimws(ln), "[[:space:]]+")[[1]]))
      if (length(vals) >= 3L && all(is.finite(vals[1:3]))) {
        dims <- vals[1:3]
        break
      }
    }
  }
  if (is.null(dims)) stop("Could not parse MatrixMarket dimensions: ", path)
  list(header = first, n_features = as.integer(dims[1]),
       n_cells = as.integer(dims[2]), nnz = as.numeric(dims[3]))
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, na = "")
}

prediction_to_logical <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "doublet")] <- TRUE
  out[y %in% c("false", "f", "0", "singlet")] <- FALSE
  out
}

robust_log_bounds <- function(x, nmads = 3) {
  z <- log1p(as.numeric(x))
  med <- median(z, na.rm = TRUE)
  s <- mad(z, center = med, constant = 1.4826, na.rm = TRUE)
  method <- paste0(nmads, " MAD on log1p scale")
  if (!is.finite(s) || s <= .Machine$double.eps) {
    q <- as.numeric(quantile(x, c(0.005, 0.995), na.rm = TRUE, names = FALSE, type = 8))
    return(c(lower = max(0, q[1]), upper = q[2], center = expm1(med), scale = s,
             method_code = 2))
  }
  c(lower = max(0, expm1(med - nmads * s)),
    upper = expm1(med + nmads * s), center = expm1(med), scale = s,
    method_code = 1)
}

robust_upper <- function(x, nmads = 3) {
  z <- log1p(as.numeric(x))
  med <- median(z, na.rm = TRUE)
  s <- mad(z, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= .Machine$double.eps) {
    q <- as.numeric(quantile(x, 0.995, na.rm = TRUE, names = FALSE, type = 8))
    return(c(upper = q, center = expm1(med), scale = s, method_code = 2))
  }
  c(upper = expm1(med + nmads * s), center = expm1(med), scale = s,
    method_code = 1)
}

safe_pct <- function(num, den) if (is.finite(den) && den > 0) 100 * num / den else NA_real_

if (identical(mode, "validate")) {
  audit <- vector("list", nrow(manifest))
  reference_features <- NULL
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    p <- sample_paths(row)
    exists_vec <- vapply(p, file.exists, logical(1))
    if (!all(exists_vec)) stop("Missing input for ", row$sample, ": ", paste(names(exists_vec)[!exists_vec], collapse = ", "))

    feat <- read_features(p$features)
    bc <- read_barcodes(p$barcodes)
    dbl <- read_doublet(p$doublet)
    mh <- read_mtx_header(p$matrix)
    pred <- prediction_to_logical(dbl$prediction)

    if (is.null(reference_features)) reference_features <- feat
    features_identical <- identical(feat, reference_features)
    file_prefix_ok <- all(startsWith(basename(unlist(p)), paste0(row$sample_title, "_")))
    ensembl_nonempty <- nzchar(feat$ensembl_id)
    symbol_nonempty <- nzchar(feat$gene_symbol)
    feature_type_ok <- all(feat$feature_type == "Gene Expression")
    barcode_pattern_ok <- grepl("^[ACGTN]+-[0-9]+$", bc)
    barcode_set_equal <- setequal(bc, as.character(dbl$barcodes))
    score_ok <- is.numeric(dbl$score) && all(is.finite(dbl$score)) && all(dbl$score >= 0 & dbl$score <= 1)
    prediction_ok <- !anyNA(pred)
    matrix_header_ok <- grepl("^%%MatrixMarket[[:space:]]+matrix[[:space:]]+coordinate", mh$header)

    audit[[i]] <- data.frame(
      sample = row$sample,
      donor = row$subject_id,
      condition = row$condition,
      all_files_present = all(exists_vec),
      file_prefix_matches_manifest = file_prefix_ok,
      matrixmarket_header_ok = matrix_header_ok,
      matrix_n_features = mh$n_features,
      features_rows = nrow(feat),
      matrix_n_cells = mh$n_cells,
      barcode_rows = length(bc),
      matrix_nnz = mh$nnz,
      matrix_feature_dimension_ok = mh$n_features == nrow(feat),
      matrix_cell_dimension_ok = mh$n_cells == length(bc),
      feature_table_three_columns = ncol(feat) == 3L,
      all_gene_ids_nonempty = all(ensembl_nonempty),
      all_gene_ids_ensembl_like = all(grepl("^ENSG[0-9]+([.][0-9]+)?$", feat$ensembl_id)),
      duplicate_ensembl_ids = sum(duplicated(feat$ensembl_id)),
      duplicate_gene_symbol_rows = sum(duplicated(feat$gene_symbol) | duplicated(feat$gene_symbol, fromLast = TRUE)),
      duplicate_gene_symbol_names = length(unique(feat$gene_symbol[duplicated(feat$gene_symbol)])),
      empty_gene_symbols = sum(!symbol_nonempty),
      feature_type_gene_expression_only = feature_type_ok,
      duplicate_barcodes = sum(duplicated(bc)),
      valid_barcode_format_rows = sum(barcode_pattern_ok),
      invalid_barcode_format_rows = sum(!barcode_pattern_ok),
      doublet_rows = nrow(dbl),
      duplicate_doublet_barcodes = sum(duplicated(dbl$barcodes)),
      barcode_doublet_row_count_ok = length(bc) == nrow(dbl),
      barcode_doublet_set_equal = barcode_set_equal,
      doublet_score_range_ok = score_ok,
      doublet_prediction_values_ok = prediction_ok,
      features_identical_to_first_sample = features_identical,
      stringsAsFactors = FALSE
    )
  }
  audit <- do.call(rbind, audit)
  write_csv(audit, file.path(table_dir, "input_file_audit.csv"))

  feature_key <- data.frame(
    ensembl_id = reference_features$ensembl_id,
    gene_symbol = reference_features$gene_symbol,
    feature_type = reference_features$feature_type,
    gene_symbol_duplicated = duplicated(reference_features$gene_symbol) |
      duplicated(reference_features$gene_symbol, fromLast = TRUE),
    stringsAsFactors = FALSE
  )
  write_csv(feature_key, file.path(table_dir, "feature_annotation.csv"))
  dup_symbols <- feature_key[feature_key$gene_symbol_duplicated, , drop = FALSE]
  write_csv(dup_symbols, file.path(table_dir, "duplicate_gene_symbols.csv"))

  checks <- c("all_files_present", "file_prefix_matches_manifest", "matrixmarket_header_ok",
              "matrix_feature_dimension_ok", "matrix_cell_dimension_ok",
              "feature_table_three_columns", "all_gene_ids_nonempty",
              "all_gene_ids_ensembl_like", "feature_type_gene_expression_only",
              "barcode_doublet_row_count_ok", "barcode_doublet_set_equal",
              "doublet_score_range_ok", "doublet_prediction_values_ok",
              "features_identical_to_first_sample")
  if (!all(vapply(audit[checks], all, logical(1)))) {
    stop("One or more input validation checks failed. See input_file_audit.csv")
  }
  if (any(audit$duplicate_ensembl_ids > 0) || any(audit$duplicate_barcodes > 0) ||
      any(audit$duplicate_doublet_barcodes > 0)) {
    stop("Duplicate stable IDs or barcodes found. See input_file_audit.csv")
  }

  writeLines(c(
    paste0("validated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("samples=", nrow(audit)),
    paste0("total_cells=", sum(audit$matrix_n_cells)),
    paste0("genes_per_sample=", paste(unique(audit$matrix_n_features), collapse = ",")),
    paste0("duplicate_ensembl_ids=", max(audit$duplicate_ensembl_ids)),
    paste0("duplicate_gene_symbol_rows=", max(audit$duplicate_gene_symbol_rows)),
    paste0("duplicate_gene_symbol_names=", max(audit$duplicate_gene_symbol_names)),
    paste0("features_identical_all_samples=", all(audit$features_identical_to_first_sample))
  ), file.path(log_dir, "validation_summary.txt"), useBytes = TRUE)
  message("Validation completed for ", nrow(audit), " samples.")
}

if (identical(mode, "sample")) {
  if (is.na(sample_arg) || !nzchar(sample_arg)) stop("--sample is required in sample mode")
  hit <- which(manifest$sample == sample_arg)
  if (length(hit) != 1L) stop("Sample not uniquely found in manifest: ", sample_arg)
  row <- manifest[hit, , drop = FALSE]
  p <- sample_paths(row)
  if (!all(vapply(p, file.exists, logical(1)))) stop("One or more sample files are missing")

  feature_map <- read_features(p$features)
  barcodes <- read_barcodes(p$barcodes)
  dbl <- read_doublet(p$doublet)
  pred <- prediction_to_logical(dbl$prediction)
  if (anyNA(pred)) stop("Unrecognized doublet prediction values for ", sample_arg)
  if (!setequal(barcodes, dbl$barcodes)) stop("Barcode/doublet mismatch for ", sample_arg)

  counts <- ReadMtx(mtx = p$matrix, cells = p$barcodes, features = p$features,
                    cell.column = 1, feature.column = 1, unique.features = TRUE,
                    strip.suffix = FALSE)
  if (nrow(counts) != nrow(feature_map) || ncol(counts) != length(barcodes))
    stop("Loaded matrix dimensions do not match sidecar files for ", sample_arg)
  if (!identical(rownames(counts), feature_map$ensembl_id))
    stop("Loaded Ensembl feature order changed for ", sample_arg)
  if (!identical(colnames(counts), barcodes))
    stop("Loaded barcode order changed for ", sample_arg)

  obj <- CreateSeuratObject(counts = counts, project = sample_arg,
                            assay = "RNA", min.cells = 0, min.features = 0)

  batch_value <- "not_reported"
  if ("batch" %in% names(row) && !is.na(row$batch) && nzchar(row$batch)) batch_value <- row$batch
  meta <- data.frame(
    sample = sample_arg,
    donor = as.character(row$subject_id),
    condition = as.character(row$condition),
    batch = batch_value,
    tissue = as.character(row$tissue),
    technology = as.character(row$assay),
    group = if ("group" %in% names(row)) as.character(row$group) else NA_character_,
    barcode_raw = barcodes,
    stringsAsFactors = FALSE,
    row.names = barcodes
  )
  obj <- AddMetaData(obj, metadata = meta)

  symbol_dup <- duplicated(feature_map$gene_symbol) |
    duplicated(feature_map$gene_symbol, fromLast = TRUE)
  feature_meta <- data.frame(
    gene_symbol = feature_map$gene_symbol,
    feature_type = feature_map$feature_type,
    gene_symbol_duplicated = symbol_dup,
    stringsAsFactors = FALSE,
    row.names = feature_map$ensembl_id
  )
  obj[["RNA"]] <- AddMetaData(obj[["RNA"]], metadata = feature_meta)

  mt_ids <- feature_map$ensembl_id[grepl("^MT-", feature_map$gene_symbol)]
  mt_ids <- intersect(mt_ids, rownames(obj))
  if (!length(mt_ids)) stop("No mitochondrial features matched gene symbols for ", sample_arg)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_ids, assay = "RNA")

  dbl_idx <- match(barcodes, dbl$barcodes)
  obj$doublet_score_source <- as.numeric(dbl$score[dbl_idx])
  obj$doublet_prediction_source <- factor(ifelse(pred[dbl_idx], "doublet", "singlet"),
                                          levels = c("singlet", "doublet"))

  fb <- robust_log_bounds(obj$nFeature_RNA)
  cb <- robust_log_bounds(obj$nCount_RNA)
  mtb <- robust_upper(obj$percent.mt)

  obj$qc_low_nFeature <- obj$nFeature_RNA < fb["lower"]
  obj$qc_high_nFeature <- obj$nFeature_RNA > fb["upper"]
  obj$qc_low_nCount <- obj$nCount_RNA < cb["lower"]
  obj$qc_high_nCount <- obj$nCount_RNA > cb["upper"]
  obj$qc_high_percent_mt <- obj$percent.mt > mtb["upper"]
  obj$qc_candidate_remove <- obj$qc_low_nFeature | obj$qc_high_nFeature |
    obj$qc_low_nCount | obj$qc_high_nCount | obj$qc_high_percent_mt
  obj$doublet_candidate_remove <- obj$doublet_prediction_source == "doublet"
  obj$union_candidate_remove <- obj$qc_candidate_remove | obj$doublet_candidate_remove

  n <- ncol(obj)
  n_qc_remove <- sum(obj$qc_candidate_remove)
  n_doublet <- sum(obj$doublet_candidate_remove)
  n_union <- sum(obj$union_candidate_remove)
  n_keep <- n - n_union
  expected_doublet_rate <- 0.008 * n / 1000
  expected_doublet_count <- round(expected_doublet_rate * n)

  threshold_row <- data.frame(
    sample = sample_arg,
    n_cells = n,
    nFeature_lower = as.numeric(fb["lower"]),
    nFeature_upper = as.numeric(fb["upper"]),
    nCount_lower = as.numeric(cb["lower"]),
    nCount_upper = as.numeric(cb["upper"]),
    percent_mt_upper = as.numeric(mtb["upper"]),
    nFeature_method = ifelse(fb["method_code"] == 1, "3_MAD_log1p", "q0.5_q99.5_fallback"),
    nCount_method = ifelse(cb["method_code"] == 1, "3_MAD_log1p", "q0.5_q99.5_fallback"),
    percent_mt_method = ifelse(mtb["method_code"] == 1, "3_MAD_log1p_upper", "q99.5_fallback"),
    qc_candidate_remove = n_qc_remove,
    qc_candidate_remove_pct = safe_pct(n_qc_remove, n),
    source_doublet_predicted = n_doublet,
    source_doublet_predicted_pct = safe_pct(n_doublet, n),
    expected_doublet_rate_pct_10x_rule = 100 * expected_doublet_rate,
    expected_doublet_count_10x_rule = expected_doublet_count,
    union_candidate_remove = n_union,
    union_candidate_remove_pct = safe_pct(n_union, n),
    proposed_retain = n_keep,
    proposed_retain_pct = safe_pct(n_keep, n),
    stringsAsFactors = FALSE
  )
  write_csv(threshold_row, file.path(per_sample_dir, paste0(sample_arg, "_thresholds.csv")))

  summary_row <- data.frame(
    sample = sample_arg,
    donor = as.character(row$subject_id),
    condition = as.character(row$condition),
    batch = batch_value,
    tissue = as.character(row$tissue),
    technology = as.character(row$assay),
    n_cells = n,
    n_genes = nrow(obj),
    n_mito_features = length(mt_ids),
    median_nFeature = median(obj$nFeature_RNA),
    q01_nFeature = as.numeric(quantile(obj$nFeature_RNA, 0.01, type = 8)),
    q99_nFeature = as.numeric(quantile(obj$nFeature_RNA, 0.99, type = 8)),
    median_nCount = median(obj$nCount_RNA),
    q01_nCount = as.numeric(quantile(obj$nCount_RNA, 0.01, type = 8)),
    q99_nCount = as.numeric(quantile(obj$nCount_RNA, 0.99, type = 8)),
    median_percent_mt = median(obj$percent.mt),
    q95_percent_mt = as.numeric(quantile(obj$percent.mt, 0.95, type = 8)),
    q99_percent_mt = as.numeric(quantile(obj$percent.mt, 0.99, type = 8)),
    stringsAsFactors = FALSE
  )
  write_csv(summary_row, file.path(per_sample_dir, paste0(sample_arg, "_summary.csv")))

  gene_totals <- Matrix::rowSums(counts)
  ord <- order(gene_totals, decreasing = TRUE)
  top_n <- min(20L, length(ord))
  top_idx <- ord[seq_len(top_n)]
  total_umis <- sum(gene_totals)
  top_genes <- data.frame(
    sample = sample_arg,
    rank = seq_len(top_n),
    ensembl_id = feature_map$ensembl_id[top_idx],
    gene_symbol = feature_map$gene_symbol[top_idx],
    total_counts = as.numeric(gene_totals[top_idx]),
    pct_sample_counts = safe_pct(as.numeric(gene_totals[top_idx]), total_umis),
    stringsAsFactors = FALSE
  )
  write_csv(top_genes, file.path(per_sample_dir, paste0(sample_arg, "_top_genes.csv")))

  cardiac_ambient_symbols <- c("TTN", "MYH6", "MYH7", "TNNT2", "TNNI3", "ACTC1",
                               "MYL2", "MYL3", "NPPA", "NPPB", "MB")
  panel_ids <- feature_map$ensembl_id[feature_map$gene_symbol %in% cardiac_ambient_symbols]
  panel_ids <- intersect(panel_ids, rownames(counts))
  low_cut <- as.numeric(quantile(obj$nCount_RNA, 0.10, type = 8))
  low_idx <- obj$nCount_RNA <= low_cut
  nonlow_idx <- !low_idx
  panel_counts_all <- if (length(panel_ids)) sum(counts[panel_ids, , drop = FALSE]) else 0
  panel_counts_low <- if (length(panel_ids) && any(low_idx)) sum(counts[panel_ids, low_idx, drop = FALSE]) else 0
  panel_counts_nonlow <- if (length(panel_ids) && any(nonlow_idx)) sum(counts[panel_ids, nonlow_idx, drop = FALSE]) else 0
  all_counts_low <- if (any(low_idx)) sum(counts[, low_idx, drop = FALSE]) else 0
  all_counts_nonlow <- if (any(nonlow_idx)) sum(counts[, nonlow_idx, drop = FALSE]) else 0
  panel_detect_any <- if (length(panel_ids)) Matrix::colSums(counts[panel_ids, , drop = FALSE] > 0) > 0 else rep(FALSE, n)
  top20_counts_low <- if (any(low_idx)) sum(counts[top_idx, low_idx, drop = FALSE]) else 0
  low_panel_pct <- safe_pct(panel_counts_low, all_counts_low)
  nonlow_panel_pct <- safe_pct(panel_counts_nonlow, all_counts_nonlow)
  ambient_row <- data.frame(
    sample = sample_arg,
    low_count_proxy_quantile = 0.10,
    low_count_proxy_cutoff_nCount = low_cut,
    low_count_proxy_cells = sum(low_idx),
    cardiac_panel_features_present = length(panel_ids),
    cardiac_panel_pct_all_counts = safe_pct(panel_counts_all, total_umis),
    cardiac_panel_pct_low_count_proxy = low_panel_pct,
    cardiac_panel_pct_other_cells = nonlow_panel_pct,
    cardiac_panel_low_vs_other_ratio = (low_panel_pct + 1e-8) / (nonlow_panel_pct + 1e-8),
    cells_with_any_cardiac_panel_pct = safe_pct(sum(panel_detect_any), n),
    top20_genes_pct_low_count_proxy = safe_pct(top20_counts_low, all_counts_low),
    raw_unfiltered_droplets_available = FALSE,
    stringsAsFactors = FALSE
  )
  write_csv(ambient_row, file.path(per_sample_dir, paste0(sample_arg, "_ambient_proxy.csv")))

  obj <- RenameCells(obj, add.cell.id = sample_arg)
  cell_qc <- data.frame(
    cell = colnames(obj),
    barcode_raw = obj$barcode_raw,
    sample = obj$sample,
    nFeature_RNA = obj$nFeature_RNA,
    nCount_RNA = obj$nCount_RNA,
    percent.mt = obj$percent.mt,
    doublet_score_source = obj$doublet_score_source,
    doublet_prediction_source = as.character(obj$doublet_prediction_source),
    qc_candidate_remove = obj$qc_candidate_remove,
    union_candidate_remove = obj$union_candidate_remove,
    stringsAsFactors = FALSE
  )
  write_csv_gz(cell_qc, file.path(cell_qc_dir, paste0(sample_arg, "_cell_qc.csv.gz")))

  p_vln <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0) +
    plot_annotation(title = paste0(sample_arg, " | unfiltered QC distributions")) &
    theme_minimal(base_size = 11)
  ggsave(file.path(violin_dir, paste0(sample_arg, "_QC_violin.png")), p_vln,
         width = 12, height = 4.5, dpi = 180, bg = "white")

  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by = "doublet_prediction_source", pt.size = 0.45,
                       cols = c("#4C78A8", "#D62728"), raster = TRUE) +
    scale_x_log10() + scale_y_log10() +
    geom_vline(xintercept = c(cb["lower"], cb["upper"]), linetype = 2, color = "#E45756") +
    geom_hline(yintercept = c(fb["lower"], fb["upper"]), linetype = 2, color = "#E45756") +
    ggtitle("Counts vs features (dashed: proposed robust bounds)") +
    theme_minimal(base_size = 10)
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt",
                       group.by = "doublet_prediction_source", pt.size = 0.45,
                       cols = c("#4C78A8", "#D62728"), raster = TRUE) +
    scale_x_log10() +
    geom_vline(xintercept = c(cb["lower"], cb["upper"]), linetype = 2, color = "#E45756") +
    geom_hline(yintercept = mtb["upper"], linetype = 2, color = "#E45756") +
    ggtitle("Counts vs mitochondrial percentage") +
    theme_minimal(base_size = 10)
  p_scatter <- (p1 + p2) + plot_annotation(title = paste0(sample_arg, " | unfiltered QC scatter"))
  ggsave(file.path(scatter_dir, paste0(sample_arg, "_QC_scatter.png")), p_scatter,
         width = 12, height = 5.5, dpi = 180, bg = "white")

  saveRDS(obj, file.path(obj_dir, paste0(sample_arg, "_unfiltered_qc.rds")), compress = "gzip")
  writeLines(c(
    paste0("sample=", sample_arg),
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("cells=", n),
    paste0("genes=", nrow(obj)),
    paste0("qc_candidate_remove=", n_qc_remove),
    paste0("source_doublet_predicted=", n_doublet),
    paste0("union_candidate_remove=", n_union),
    "filtering_executed=FALSE"
  ), file.path(log_dir, paste0(sample_arg, ".log")), useBytes = TRUE)
  rm(obj, counts)
  invisible(gc())
  message("Completed unfiltered QC for ", sample_arg)
}

if (identical(mode, "finalize")) {
  collect_files <- function(pattern) {
    f <- list.files(per_sample_dir, pattern = pattern, full.names = TRUE)
    if (length(f) != nrow(manifest)) stop("Expected ", nrow(manifest), " files for pattern ", pattern, ", found ", length(f))
    do.call(rbind, lapply(f, function(x) read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)))
  }
  thresholds <- collect_files("_thresholds[.]csv$")
  summaries <- collect_files("_summary[.]csv$")
  ambient <- collect_files("_ambient_proxy[.]csv$")
  top_genes <- collect_files("_top_genes[.]csv$")

  thresholds <- thresholds[match(manifest$sample, thresholds$sample), , drop = FALSE]
  summaries <- summaries[match(manifest$sample, summaries$sample), , drop = FALSE]
  ambient <- ambient[match(manifest$sample, ambient$sample), , drop = FALSE]
  ambient$relative_risk_rank <- rank(-ambient$cardiac_panel_pct_low_count_proxy, ties.method = "min")
  q75 <- as.numeric(quantile(ambient$cardiac_panel_pct_low_count_proxy, 0.75, na.rm = TRUE, type = 8))
  q25 <- as.numeric(quantile(ambient$cardiac_panel_pct_low_count_proxy, 0.25, na.rm = TRUE, type = 8))
  ambient$relative_risk_class <- ifelse(ambient$cardiac_panel_pct_low_count_proxy >= q75, "higher_relative",
                                        ifelse(ambient$cardiac_panel_pct_low_count_proxy <= q25,
                                               "lower_relative", "middle_relative"))

  write_csv(thresholds, file.path(table_dir, "qc_threshold_recommendations.csv"))
  write_csv(summaries, file.path(table_dir, "sample_qc_summary.csv"))
  write_csv(ambient, file.path(table_dir, "ambient_rna_risk_proxy.csv"))
  write_csv(top_genes, file.path(table_dir, "top_expressed_genes_by_sample.csv"))
  write_csv(thresholds[, c("sample", "n_cells", "source_doublet_predicted",
                           "source_doublet_predicted_pct", "expected_doublet_rate_pct_10x_rule",
                           "expected_doublet_count_10x_rule")],
            file.path(table_dir, "doublet_risk_summary.csv"))

  # ---- cross-sample overview figures (pre-filter, descriptive) ----
  summ <- summaries
  summ <- summ[match(thresholds$sample, summ$sample), , drop = FALSE]
  keep_samples <- thresholds$sample
  cell_list <- lapply(keep_samples, function(s) {
    x <- read.csv(gzfile(file.path(cell_qc_dir, paste0(s, "_cell_qc.csv.gz"))),
                  stringsAsFactors = FALSE, check.names = FALSE)
    x[, c("sample", "nFeature_RNA", "nCount_RNA", "percent.mt",
          "doublet_prediction_source")]
  })
  all_cells <- do.call(rbind, cell_list)
  cond_map <- summ$condition[match(keep_samples, summ$sample)]
  group_levels <- keep_samples[order(match(cond_map, c("Control", "Disease")), keep_samples)]
  all_cells$sample <- factor(all_cells$sample, levels = group_levels)

  metrics_long <- rbind(
    data.frame(sample = all_cells$sample, metric = "nFeature_RNA", value = all_cells$nFeature_RNA),
    data.frame(sample = all_cells$sample, metric = "nCount_RNA", value = all_cells$nCount_RNA),
    data.frame(sample = all_cells$sample, metric = "percent.mt", value = all_cells$percent.mt)
  )
  p_vln_all <- ggplot(metrics_long, aes(x = sample, y = value)) +
    geom_violin(fill = "#7FB3D5", color = "grey30", scale = "width") +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("All samples | unfiltered QC distributions") +
    xlab(NULL)
  ggsave(file.path(out_dir, "qc_plots", "all_samples_QC_violin.png"), p_vln_all,
         width = 12, height = 10, dpi = 160, bg = "white")

  sc1 <- ggplot(all_cells, aes(x = nCount_RNA, y = nFeature_RNA,
                               color = doublet_prediction_source)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_manual(values = c(singlet = "#4C78A8", doublet = "#D62728")) +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~sample, ncol = 5, scales = "free") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom") +
    ggtitle("All samples | nCount vs nFeature (log-log)")
  ggsave(file.path(out_dir, "qc_plots", "all_samples_counts_features.png"), sc1,
         width = 14, height = 10, dpi = 160, bg = "white")

  sc2 <- ggplot(all_cells, aes(x = nCount_RNA, y = percent.mt,
                               color = doublet_prediction_source)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_manual(values = c(singlet = "#4C78A8", doublet = "#D62728")) +
    scale_x_log10() +
    facet_wrap(~sample, ncol = 5, scales = "free") +
    theme_bw(base_size = 8) +
    theme(legend.position = "bottom") +
    ggtitle("All samples | nCount vs percent.mt")
  ggsave(file.path(out_dir, "qc_plots", "all_samples_counts_mt.png"), sc2,
         width = 14, height = 10, dpi = 160, bg = "white")

  session_lines <- capture.output(sessionInfo())
  writeLines(session_lines, file.path(log_dir, "sessionInfo.txt"), useBytes = TRUE)

  total_cells <- sum(thresholds$n_cells)
  total_qc <- sum(thresholds$qc_candidate_remove)
  total_doublet <- sum(thresholds$source_doublet_predicted)
  total_union <- sum(thresholds$union_candidate_remove)
  total_keep <- sum(thresholds$proposed_retain)
  audit <- read.csv(file.path(table_dir, "input_file_audit.csv"), stringsAsFactors = FALSE,
                    check.names = FALSE)
  feature_key <- read.csv(file.path(table_dir, "feature_annotation.csv"), stringsAsFactors = FALSE)
  n_dup_symbol_rows <- sum(feature_key$gene_symbol_duplicated)
  n_dup_symbol_names <- length(unique(feature_key$gene_symbol[feature_key$gene_symbol_duplicated]))
  top_ambient <- ambient[order(ambient$relative_risk_rank), , drop = FALSE]
  top_ambient <- head(top_ambient, 5)

  md <- c(
    "# GSE292067 Single-cell Stage-1 QC Review",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    "> STOP POINT: no cell or gene filtering was executed. All Seurat objects are unfiltered and await approval.",
    "",
    "## 1. Input and format audit",
    "",
    paste0("- Manifest used: `", manifest_path, "` (the exact filename `sample_manifest.csv` was not present)."),
    paste0("- Raw directory: `", raw_dir, "`."),
    paste0("- Samples: ", nrow(manifest), "; nuclei: ", format(total_cells, big.mark = ","), "."),
    paste0("- Every sample has one MatrixMarket matrix, one 3-column features file, one barcode file, and one source doublet prediction file."),
    paste0("- Matrix dimensions match feature/barcode rows for all samples: ", all(audit$matrix_feature_dimension_ok & audit$matrix_cell_dimension_ok), "."),
    paste0("- Gene namespace: Ensembl IDs in column 1; gene symbols in column 2; feature type in column 3."),
    paste0("- Duplicate Ensembl IDs: ", max(audit$duplicate_ensembl_ids), "; duplicate barcodes within sample: ", max(audit$duplicate_barcodes), "."),
    paste0("- Duplicate gene-symbol rows: ", n_dup_symbol_rows, " across ", n_dup_symbol_names,
           " symbols. Objects therefore use stable Ensembl IDs as feature names and preserve symbols in RNA feature metadata."),
    paste0("- Feature order is identical across all samples: ", all(audit$features_identical_to_first_sample), "."),
    "",
    "## 2. Seurat object construction",
    "",
    "- One independent Seurat object was created per sample with `min.cells = 0` and `min.features = 0`.",
    "- Required cell metadata: sample, donor, condition, batch, tissue, technology.",
    "- `donor` comes from subject_id; `technology` comes from assay.",
    "- No technical batch field exists in the manifest, so `batch = not_reported` was written explicitly rather than inventing a batch label.",
    "- Cell names were prefixed with the sample ID; the original 10x barcode is retained as `barcode_raw`.",
    "",
    "## 3. QC metrics and plots",
    "",
    "- `nFeature_RNA` and `nCount_RNA` are generated by CreateSeuratObject.",
    "- `percent.mt` is computed with PercentageFeatureSet using Ensembl IDs whose gene symbols match `^MT-`.",
    "- Each sample has one violin PNG and one scatter PNG under `results/Singlecel/qc_plots/`.",
    "",
    "## 4. Data-driven threshold proposal (not applied)",
    "",
    "Per sample, lower/upper nFeature_RNA and nCount_RNA bounds use median +/- 3 MAD on the log1p scale. percent.mt uses a one-sided median + 3 MAD on the log1p scale. If MAD is degenerate, empirical 0.5/99.5 percentiles are used as a fallback. These are distribution-derived per-sample bounds, not fixed universal cutoffs.",
    "",
    "| sample | nFeature lower-upper | nCount lower-upper | mt upper (%) | initial | QC remove | source doublet | union remove | proposed retain (%) |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(thresholds))) {
    x <- thresholds[i, ]
    md <- c(md, sprintf("| %s | %.0f-%.0f | %.0f-%.0f | %.3f | %d | %d (%.2f%%) | %d (%.2f%%) | %d (%.2f%%) | %d (%.2f%%) |",
                        x$sample, x$nFeature_lower, x$nFeature_upper,
                        x$nCount_lower, x$nCount_upper, x$percent_mt_upper,
                        x$n_cells, x$qc_candidate_remove, x$qc_candidate_remove_pct,
                        x$source_doublet_predicted, x$source_doublet_predicted_pct,
                        x$union_candidate_remove, x$union_candidate_remove_pct,
                        x$proposed_retain, x$proposed_retain_pct))
  }
  md <- c(md,
    "",
    paste0("Overall: initial ", format(total_cells, big.mark = ","),
           "; QC-metric candidates ", format(total_qc, big.mark = ","), " (", sprintf("%.2f", safe_pct(total_qc, total_cells)), "%)",
           "; source-predicted doublets ", format(total_doublet, big.mark = ","), " (", sprintf("%.2f", safe_pct(total_doublet, total_cells)), "%)",
           "; union candidates ", format(total_union, big.mark = ","), " (", sprintf("%.2f", safe_pct(total_union, total_cells)), "%)",
           "; proposed retained ", format(total_keep, big.mark = ","), " (", sprintf("%.2f", safe_pct(total_keep, total_cells)), "%)."),
    "",
    "## 5. Doublet risk",
    "",
    "- A score and Boolean prediction are supplied for every barcode. They were imported as source annotations; the file itself does not identify the algorithm or threshold provenance, so the calls should not be treated as independently validated.",
    "- The report compares observed source calls with the approximate 10x loading heuristic (~0.8% per 1,000 recovered cells).",
    "- Recommended independent validation after approval: scDblFinder per sample on raw counts. Installed dependencies are available (scDblFinder and SingleCellExperiment); no extra installation is currently required.",
    "- Doublet calls were not used to subset any object in this stage.",
    "",
    "## 6. Ambient RNA risk",
    "",
    "- These files contain only cell-called nuclei (5,072-9,818 barcodes per sample) and no unfiltered droplet matrix. Therefore EmptyDrops, SoupX, and CellBender cannot be run correctly from the available inputs.",
    "- A descriptive proxy was computed only: cardiac parenchymal-marker share (TTN/MYH6/MYH7/TNNT2/TNNI3/ACTC1/MYL2/MYL3/NPPA/NPPB/MB) in each sample's lowest-count decile. This is not a contamination estimate and was not used for filtering.",
    "- Samples with the highest relative proxy values are:",
    ""
  )
  for (i in seq_len(nrow(top_ambient))) {
    x <- top_ambient[i, ]
    md <- c(md, sprintf("  - %s: low-count cardiac-panel %.3f%%; other cells %.3f%%; ratio %.2f.",
                        x$sample, x$cardiac_panel_pct_low_count_proxy,
                        x$cardiac_panel_pct_other_cells, x$cardiac_panel_low_vs_other_ratio))
  }
  md <- c(md,
    "",
    "- If formal ambient correction is requested later: obtain raw/unfiltered droplets and use SoupX (R, requires raw + filtered matrices and cluster labels) or CellBender (Python/GPU, raw matrix). Do not stack ambient-correction tools. With filtered cells only, DecontX/celda is possible but adds Bioconductor dependencies and relies on cluster structure; it is not run in Stage 1.",
    "",
    "## 7. Methods and software checks",
    "",
    "- Local bioSkills checked: single-cell/data-io, single-cell/preprocessing, single-cell/doublet-detection, workflows/scrnaseq-pipeline.",
    "- External SkillHub search returned no relevant single-cell skill above the registry relevance threshold; no skill was installed.",
    paste0("- R: ", R.version.string, "; Seurat ", as.character(packageVersion("Seurat")),
           "; SeuratObject ", as.character(packageVersion("SeuratObject")),
           "; Matrix ", as.character(packageVersion("Matrix")), "."),
    "- Official references used:",
    "  - https://satijalab.org/seurat/reference/readmtx",
    "  - https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html",
    "  - https://satijalab.org/seurat/reference/percentagefeatureset",
    "  - https://satijalab.org/seurat/articles/pbmc3k_tutorial",
    "",
    "## 8. Decision required",
    "",
    "Please review the per-sample plots and threshold table, confirm whether to accept the proposed union rule (QC outliers OR source doublet prediction), and decide whether an independent scDblFinder validation should be run before filtering.",
    "",
    "**Filtering status: NOT EXECUTED.**"
  )
  con <- file(report_path, open = "wb")
  writeLines(enc2utf8(md), con, useBytes = TRUE)
  close(con)

  writeLines(c(
    paste0("finalized_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("samples=", nrow(thresholds)),
    paste0("initial_cells=", total_cells),
    paste0("qc_candidates=", total_qc),
    paste0("source_doublets=", total_doublet),
    paste0("union_candidates=", total_union),
    paste0("proposed_retain=", total_keep),
    "filtering_executed=FALSE"
  ), file.path(log_dir, "stage1_finalize.log"), useBytes = TRUE)
  message("Stage-1 QC finalized. No filtering executed.")
}

if (!mode %in% c("validate", "sample", "finalize")) stop("Unknown mode: ", mode)
