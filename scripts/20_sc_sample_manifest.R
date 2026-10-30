# =============================================================================
# 20_sc_sample_manifest.R  --  GSE292067 single-nucleus sample manifest builder
#
# Purpose
#   GSE292067 ships NO separate phenotype table in the raw_data folder. All
#   grouping information is encoded in the supplementary file names, e.g.
#       GSM8848306_AJAI-Zhen-Ctrl1_1666_barcodes.tsv.gz
#                          ^^^^^  ^^^^
#                      group token  subject id
#   This script parses that structure, maps each token onto the experimental
#   arm reported in the source paper (Circulation 2025, PMID 40808662), and
#   emits a standard single-cell sample manifest (sample / group / covariates)
#   for every downstream step (Read10X, integration, pseudobulk, DEG).
#
# Input  (read-only):  raw_data/GSE292067/GSM*_*.{barcodes,features,matrix}...gz
# Output (clean_data): GSE292067_sample_manifest.csv
#                      GSE292067_group_counts.csv
#                      GSE292067_sample_manifest_qc.txt
#
# NOTE: this source file is intentionally kept in pure ASCII. All Chinese
# strings are written with \u escapes so the script survives any locale.
# =============================================================================

rm(list = ls())

# ---------- 0. paths ----------
proj_dir  <- "E:/another data file/GEObulk analysis project"
raw_dir   <- file.path(proj_dir, "raw_data", "GSE292067")
clean_dir <- file.path(proj_dir, "clean_data")
if (!dir.exists(clean_dir)) dir.create(clean_dir, recursive = TRUE)

gse      <- "GSE292067"
platform <- "GPL24676"

t_start <- Sys.time()
hr  <- function() cat(strrep("-", 72), "\n", sep = "")
msg <- function(...) cat(..., "\n", sep = "")
LOG <- character(0)
log_add <- function(...) LOG <<- c(LOG, paste0(...))

if (!dir.exists(raw_dir)) stop("raw_data/GSE292067 not found")
msg(">>> GSE292067 sample manifest builder")
msg("    raw_dir : ", raw_dir)

# ---------- 1. group dictionary (token -> canonical level) ----------
# Source: GSE292067 = human LV snRNA-seq, 25 samples, 3 arms.
#   Ctrl1 = non-failing donor           (reference / control)
#   DoxT2 = doxorubicin cardiomyopathy  (DoxCM)
#   DCM4  = non-ischemic cardiomyopathy (NICM)
GROUP_MAP <- list(
  Ctrl = list(group     = "NF",
              full      = "non-failing donor",
              cn        = "\u975e\u8870\u7aed\u4f9b\u4f53\u5bf9\u7167",
              disease   = "non-failing donor",
              control   = TRUE),
  DoxT = list(group     = "DoxCM",
              full      = "doxorubicin-induced cardiomyopathy",
              cn        = "\u963f\u9709\u7d20\u6027\u5fc3\u808c\u75c5",
              disease   = "doxorubicin cardiomyopathy (DoxCM)",
              control   = FALSE),
  DCM  = list(group     = "NICM",
              full      = "non-ischemic cardiomyopathy",
              cn        = "\u975e\u7f3a\u8840\u6027\u6269\u5f20\u578b\u5fc3\u808c\u75c5",
              disease   = "non-ischemic cardiomyopathy (NICM)",
              control   = FALSE)
)

# down-stream contrast order: control first, matches GSE57338 NF/DCM/ICM habit
GROUP_LEVELS <- c("NF", "NICM", "DoxCM")

# ---------- 2. enumerate files & parse names ----------
files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)
files <- files[grepl("\\.gz$", files, ignore.case = TRUE)]
msg("    gzip files found: ", length(files))
if (length(files) == 0) stop("no .gz files under raw_data/GSE292067")

FTYPES <- c(barcodes = "barcodes.tsv.gz",
            features = "features.tsv.gz",
            matrix   = "matrix.mtx.gz",
            doublet  = "doublet_prediction.txt.gz")

parse_name <- function(bname) {
  # strip every known file-type suffix to get the sample stem
  stem <- bname
  stem <- sub("_barcodes\\.tsv\\.gz$|_features\\.tsv\\.gz$|_matrix\\.mtx\\.gz$|_doublet_prediction\\.txt\\.gz$",
              "", stem)
  tok  <- strsplit(stem, "_")[[1]]
  if (length(tok) < 3 || !grepl("^GSM[0-9]+$", tok[1])) return(NULL)

  gsm    <- tok[1]
  sid    <- tok[length(tok)]                       # trailing numeric = subject id
  if (!grepl("^[0-9]+$", sid)) return(NULL)
  title  <- paste(tok[-c(1, length(tok))], collapse = "_")  # e.g. AJAI-Zhen-Ctrl1
  pieces <- strsplit(title, "-")[[1]]
  token  <- pieces[length(pieces)]                 # group token, e.g. Ctrl1
  list(gsm = gsm, title = stem, subject_id = sid, group_token = token,
       arm_core = paste(pieces[-length(pieces)], collapse = "-"))
}

rec <- list()
for (f in files) {
  p <- parse_name(basename(f))
  if (is.null(p)) { warning("cannot parse: ", basename(f)); next }
  ftype <- names(Filter(function(x) grepl(x, basename(f), fixed = TRUE), FTYPES))
  if (length(ftype) != 1) { warning("ambiguous file type: ", basename(f)); next }
  rec[[length(rec) + 1]] <- data.frame(sample      = p$gsm,
                                       sample_title= p$title,
                                       subject_id  = p$subject_id,
                                       group_token = p$group_token,
                                       ftype       = ftype,
                                       file        = basename(f),
                                       stringsAsFactors = FALSE)
}
dat <- do.call(rbind, rec)
stopifnot(nrow(dat) > 0)

# wide: one row per sample, one column per file type
samples <- unique(dat$sample)
manifest <- data.frame(sample = sort(samples), stringsAsFactors = FALSE)
for (ft in names(FTYPES)) {
  sub <- dat[dat$ftype == ft, c("sample", "file")]
  manifest[[paste0("file_", ft)]] <- sub$file[match(manifest$sample, sub$sample)]
}
meta <- unique(dat[, c("sample", "sample_title", "subject_id", "group_token")])
manifest <- merge(manifest, meta, by = "sample", all.x = TRUE)
manifest$has_all_files <- rowSums(!is.na(manifest[, paste0("file_", names(FTYPES))])) == length(FTYPES)

msg("    samples parsed: ", nrow(manifest), " (files per sample: ",
    paste(table(table(dat$sample)), collapse = "/"), ")")
log_add("Parsed ", nrow(manifest), " samples from ", nrow(dat), " gz files")

# ---------- 3. intelligent grouping ----------
norm_token <- toupper(manifest$group_token)
# token = letters + optional trailing index, e.g. Ctrl1 / DoxT2 / DCM4
base <- sub("([A-Za-z]+).*", "\\1", manifest$group_token)
match_key <- character(nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  hits <- names(GROUP_MAP)[match(toupper(base[i]), toupper(names(GROUP_MAP)))]
  # prefix match so Ctrl1 / Control / Ctrl all land on the same arm
  if (is.na(hits)) {
    cand <- names(GROUP_MAP)[startsWith(toupper(base[i]), toupper(names(GROUP_MAP)))]
    hits <- if (length(cand) == 1) cand else NA
  }
  match_key[i] <- if (is.na(hits)) "" else hits
}
unresolved <- match_key == ""
if (any(unresolved))
  msg("    !! unresolved group tokens: ",
      paste(unique(manifest$group_token[unresolved]), collapse = ", "))

manifest$group        <- ifelse(unresolved, NA, sapply(match_key, function(k) if (nzchar(k)) GROUP_MAP[[k]]$group else NA))
manifest$group_full   <- ifelse(unresolved, NA, sapply(match_key, function(k) if (nzchar(k)) GROUP_MAP[[k]]$full  else NA))
manifest$group_cn     <- ifelse(unresolved, NA, sapply(match_key, function(k) if (nzchar(k)) GROUP_MAP[[k]]$cn    else NA))
manifest$disease_raw  <- ifelse(unresolved, NA, sapply(match_key, function(k) if (nzchar(k)) GROUP_MAP[[k]]$disease else NA))
manifest$group_source <- "sample title (AJAI-Zhen-<arm token>)"
manifest$condition    <- ifelse(is.na(manifest$group), NA,
                                ifelse(sapply(match_key, function(k) if (nzchar(k)) GROUP_MAP[[k]]$control else NA),
                                       "Control", "Disease"))
manifest$group        <- factor(manifest$group, levels = GROUP_LEVELS)
manifest$condition    <- factor(manifest$condition, levels = c("Control", "Disease"))

# static covariates (identical for every library in this series)
manifest$organism  <- "Homo sapiens"
manifest$tissue    <- "heart left ventricle"
manifest$cell_type <- "DRAQ5+ single nucleus"
manifest$assay     <- "snRNA-seq (10x Genomics Chromium Single Cell 5' v2)"
manifest$platform  <- platform
manifest$molecule  <- "nuclear RNA"
manifest$gse       <- gse

# ---------- 4. read matrix dimensions + feature checksums ----------
msg(">>> reading matrix headers & barcode counts ...")
n_cells <- n_genes <- n_entries <- integer(nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  fm <- file.path(raw_dir, manifest$file_matrix[i])
  hd <- suppressWarnings(readLines(gzfile(fm), n = 10L))
  # 10x .mtx starts with "#" / "%" comment lines, then <genes cells> <nnz>
  hd <- hd[!grepl("^[#%]", hd)]
  hd <- trimws(hd[nzchar(hd)])
  d  <- suppressWarnings(as.numeric(strsplit(hd[1], "\\s+")[[1]]))
  n_genes[i] <- d[1]; n_cells[i] <- d[2]
  # MatrixMarket dims line carries three values: <nrow ncol nnz>
  if (length(d) >= 3 && !is.na(d[3])) n_entries[i] <- d[3]
  fb <- file.path(raw_dir, manifest$file_barcodes[i])
  nb <- length(suppressWarnings(readLines(gzfile(fb))))
  if (!length(nb)) nb <- NA_integer_
  if (!is.na(nb) && nb != n_cells[i])
    warning("barcode lines (", nb, ") != matrix cols (", n_cells[i], ") for ", manifest$sample[i])
}
manifest$n_cells   <- n_cells
manifest$n_genes   <- n_genes
manifest$n_entries <- n_entries

# NB: md5sum() on the .gz itself is useless here -- every gzip header stores the
# original file name, so identical content yields 25 different checksums. Hash the
# decompressed gene list instead (also verifies feature order across libraries).
fhash <- function(p) {
  ln <- suppressWarnings(readLines(gzfile(p)))
  digest::digest(paste(ln, collapse = "\n"), algo = "md5", serialize = FALSE)
}
manifest$features_md5 <- vapply(file.path(raw_dir, manifest$file_features), fhash, "")
feature_identical <- length(unique(as.character(manifest$features_md5))) == 1

# one subject per library?
subject_unique <- !any(duplicated(manifest$subject_id))

# ---------- 5. write outputs ----------
# CSV is kept ASCII-only: every down-stream tool (pandas, Seurat, Excel) reads it
# without guessing an encoding. The bilingual group dictionary lives in the QC file.
out_cols <- c("sample", "sample_title", "subject_id", "patient_id",
              "group", "group_full", "condition",
              "group_token", "group_source", "disease_raw",
              "organism", "tissue", "cell_type", "assay", "platform",
              "molecule", "gse", "n_cells", "n_genes", "n_entries",
              "has_all_files", "features_md5",
              "file_barcodes", "file_features", "file_matrix", "file_doublet")
manifest$patient_id <- manifest$subject_id
out <- manifest[, out_cols]
out <- out[order(out$group, out$sample), ]
rownames(out) <- NULL

# manual serialiser: write.csv() would escape UTF-8 strings as <U+XXXX>
write_ascii_csv <- function(df, path) {
  # write factors as labels, not as their internal integer codes
  for (j in seq_along(df)) if (is.factor(df[[j]])) df[[j]] <- as.character(df[[j]])
  con <- file(path, open = "wb")
  # quote only when needed; keep numerics bare so pandas/read.csv type them as ints
  q <- function(x) {
    x <- as.character(x)
    needs <- grepl('[",\n\r]', x)
    ifelse(needs, paste0('"', gsub('"', '""', x, fixed = TRUE), '"'), x)
  }
  lines <- c(paste(q(names(df)), collapse = ","),
             vapply(seq_len(nrow(df)),
                    function(i) paste(q(unlist(df[i, ], use.names = FALSE)), collapse = ","),
                    ""))
  writeLines(lines, con, useBytes = TRUE)
  close(con)
  invisible(path)
}

f_manifest <- file.path(clean_dir, paste0(gse, "_sample_manifest.csv"))
write_ascii_csv(out, f_manifest)

gc <- data.frame(group = levels(droplevels(out$group)))
gc$group_full <- out$group_full[match(gc$group, as.character(out$group))]
gc$condition  <- as.character(out$condition)[match(gc$group, as.character(out$group))]
gc$n_samples  <- as.integer(table(as.character(out$group))[gc$group])
gc$n_subjects <- as.integer(tapply(out$subject_id, as.character(out$group),
                                   function(x) length(unique(x)))[gc$group])
gc$total_nuclei   <- as.integer(tapply(out$n_cells, as.character(out$group), sum)[gc$group])
gc$median_nuclei  <- as.integer(tapply(out$n_cells, as.character(out$group), median)[gc$group])
rownames(gc) <- NULL
f_counts <- file.path(clean_dir, paste0(gse, "_group_counts.csv"))
write_ascii_csv(gc, f_counts)

qc <- c(
  paste0("dataset            : ", gse),
  paste0("title              : Human single nuclear RNA sequencing identifies CD47 as a therapeutic target for doxorubicin-induced cardiomyopathy"),
  paste0("source             : Circulation 2025, PMID 40808662"),
  paste0("samples            : ", nrow(out)),
  paste0("files per sample   : ", sum(!is.na(out$file_barcodes)) + sum(!is.na(out$file_features)) +
           sum(!is.na(out$file_matrix)) + sum(!is.na(out$file_doublet)),
         " / ", nrow(out) * 4),
  paste0("all files present  : ", all(out$has_all_files)),
  paste0("groups             : ", paste(levels(droplevels(out$group)), collapse = ", ")),
  paste0("group n            : ", paste(paste0(levels(droplevels(out$group)), "=",
           as.integer(table(as.character(out$group))[levels(droplevels(out$group))])), collapse = ", ")),
  paste0("unresolved tokens  : ", sum(is.na(as.character(out$group)))),
  paste0("unique subjects    : ", length(unique(out$subject_id)),
         " (1 library per subject: ", subject_unique, ")"),
  paste0("features identical : ", feature_identical, " (unique md5 n=",
         length(unique(as.character(out$features_md5))), ")"),
  paste0("genes per sample   : ", unique(out$n_genes)[1]),
  paste0("total nuclei       : ", sum(out$n_cells)),
  paste0("nuclei / sample    : min=", min(out$n_cells), " median=", median(out$n_cells),
         " max=", max(out$n_cells)),
  "",
  "group mapping (parsed from file names):",
  "  Ctrl1 -> NF    non-failing donor",
  "  DoxT2 -> DoxCM doxorubicin-induced cardiomyopathy",
  "  DCM4  -> NICM  non-ischemic cardiomyopathy",
  "",
  "log:", LOG
)
f_qc <- file.path(clean_dir, paste0(gse, "_sample_manifest_qc.txt"))
writeLines(qc, f_qc, useBytes = TRUE)

hr()
msg("OUTPUT  ", f_manifest)
msg("        ", f_counts)
msg("        ", f_qc)
hr()
print(gc)
msg("\nelapsed ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s")
