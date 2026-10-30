# =============================================================================
# 80_ext_parse.R  —— 外部队列解析：GSE5406（HG-U133A 芯片）、GSE116250（RNA-seq RPKM）
# 只做"读进来 + 探针/转录本映射到基因符号 + 样本表"，不做任何建模
# 日志一律 ASCII，避免 C locale 下中文被写坏
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

root <- "E:/another data file/GEObulk analysis project"
raw  <- file.path(root, "raw_data", "external")
outd <- file.path(root, "results", "ExternalValidation")
dirs <- file.path(outd, c("01_data", "01_tables", "02_figures", "03_logs"))
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

logcon <- file(file.path(outd, "03_logs", "80_parse_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, logcon) }

genes4 <- c("AQP3", "GGT5", "ISLR", "LAD1")
extra  <- c("NPPB", "NPPA", "CD63", "CD163", "STAT3")

# ---------------------------------------------------------------- series matrix
read_matrix <- function(gz) {
  con <- gzfile(gz, "r"); ln <- readLines(con, warn = FALSE); close(con)
  bidx <- grep("^!series_matrix_table_begin", ln); eidx <- grep("^!series_matrix_table_end", ln)
  if (length(bidx) == 0) stop("no table in ", gz)
  if (length(eidx) == 0) eidx <- length(ln) + 1
  body <- ln[(bidx + 1):(eidx - 1)]; body <- body[nzchar(body)]
  read.delim(text = body, header = TRUE, check.names = FALSE, quote = "\"",
             stringsAsFactors = FALSE)
}
# 重复的 !Sample_ 字段拆成 _1 _2 _3 ...
sample_meta <- function(gz) {
  con <- gzfile(gz, "r"); ln <- readLines(con, warn = FALSE); close(con)
  ln <- ln[grepl("^!Sample_", ln)]
  key <- sub("^!Sample_([^\t]+)\t.*$", "\\1", ln)
  vals <- lapply(ln, function(x) {
    x <- sub("^!Sample_[^\t]+\t", "", x)
    gsub('^"|"$', "", strsplit(x, "\t", fixed = TRUE)[[1]])
  })
  n <- length(vals[[1]])
  res <- list()
  for (i in seq_along(key)) {
    k <- key[i]; v <- vals[[i]]
    if (length(v) != n) next
    if (is.null(res[[k]])) res[[k]] <- list(v)
    else {
      dup <- any(sapply(res[[k]], function(z) identical(z, v)))
      if (!dup) res[[k]] <- c(res[[k]], list(v))
    }
  }
  out <- list()
  for (k in names(res)) {
    L <- res[[k]]
    if (length(L) == 1) out[[k]] <- L[[1]]
    else for (j in seq_along(L)) out[[paste0(k, "_", j)]] <- L[[j]]
  }
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}
grp_from_text <- function(txt) {
  x <- tolower(txt); out <- rep(NA_character_, length(x))
  out[grepl("nonfail|non-fail|non failing|non-failing|normal|unused donor", x)] <- "NF"
  out[grepl("ischemic", x)] <- "ICM"
  out[grepl("idiopathic|dilated", x)] <- "DCM"
  out
}

# ============================================================ 1. GSE5406
say("================ GSE5406 (GPL96 / HG-U133A / RMA log2) ================")
m5406 <- read_matrix(file.path(raw, "GSE5406_series_matrix.txt.gz"))
say("raw matrix dim: ", nrow(m5406), " x ", ncol(m5406))
met5406 <- sample_meta(file.path(raw, "GSE5406_series_matrix.txt.gz"))
chr_col <- grep("characteristics", colnames(met5406), value = TRUE)
say("characteristics columns: ", paste(chr_col, collapse = ", "))
for (cc in chr_col) { tb <- table(met5406[[cc]]); for (nm in names(tb))
  say(sprintf("  [%s] %-60s n=%d", cc, substr(nm, 1, 60), tb[[nm]])) }

txt5406 <- apply(met5406[, chr_col, drop = FALSE], 1, paste, collapse = " ")
grp <- grp_from_text(txt5406)
gt  <- grp_from_text(met5406$title)
say("group table (characteristics):"); print(table(grp, useNA = "ifany"))
say("NA in characteristics-based group: ", sum(is.na(grp)),
    " ; rescued by title: ", sum(is.na(grp) & !is.na(gt)))
grp[is.na(grp)] <- gt[is.na(grp)]
say("final group table:"); print(table(grp, useNA = "ifany"))
say("age/sex available in GSE5406: ",
    any(grepl("age|sex|gender", tolower(colnames(met5406)))))

info5406 <- data.frame(sample = met5406$geo_accession, title = met5406$title,
                       group = grp, source = met5406$source_name_ch1,
                       stringsAsFactors = FALSE)

ann <- fread(file.path(raw, "GPL96.annot.gz"), sep = "\t", header = TRUE,
             comment.char = "", quote = "", skip = "ID\t", data.table = FALSE)
say("GPL96.annot dim: ", nrow(ann), " x ", ncol(ann))
sym_col <- grep("^Gene symbol$", colnames(ann), value = TRUE)[1]
expr5406 <- as.matrix(m5406[, -1, drop = FALSE]); rownames(expr5406) <- m5406[[1]]
say("value range: ", paste(round(range(expr5406, na.rm = TRUE), 2), collapse = " ~ "),
    " ; median = ", round(median(expr5406, na.rm = TRUE), 2))

ps <- setNames(ann[[sym_col]], ann$ID)[rownames(expr5406)]
ps[is.na(ps) | ps == "" | ps == "---"] <- NA
ps <- sub(" ?///.*$", "", ps)

target <- unique(c(genes4, extra))
collapse_gene <- function(expr, ps, symbols) {
  res <- matrix(NA_real_, length(symbols), ncol(expr), dimnames = list(symbols, colnames(expr)))
  nprobe <- setNames(integer(length(symbols)), symbols)
  for (g in symbols) {
    idx <- which(ps == g); if (!length(idx)) next
    sub <- expr[idx, , drop = FALSE]; k <- which.max(rowMeans(sub, na.rm = TRUE))
    res[g, ] <- sub[k, ]; nprobe[[g]] <- length(idx)
  }
  list(expr = res, nprobe = nprobe)
}
cg <- collapse_gene(expr5406, ps, target)
say("--- probes per gene on GPL96 ---")
for (g in target) say(sprintf("  %-6s np=%d  detected=%s", g, cg$nprobe[[g]],
                              ifelse(is.na(cg$expr[g, 1]), "NO", "YES")))
say("missing genes: ", paste(names(cg$nprobe)[cg$nprobe == 0], collapse = ","))
say("four-gene value range: ",
    paste(round(range(cg$expr[genes4, ], na.rm = TRUE), 2), collapse = " ~ "))

write.csv(cg$expr, file.path(outd, "01_data", "GSE5406_gene_expr.csv"))
write.csv(info5406, file.path(outd, "01_data", "GSE5406_sample_info.csv"), row.names = FALSE)

# ============================================================ 2. GSE116250
say("")
say("================ GSE116250 (GPL16791 / HiSeq2500 / RPKM) ================")
met11 <- sample_meta(file.path(raw, "GSE116250_series_matrix.txt.gz"))
chr11 <- grep("characteristics", colnames(met11), value = TRUE)
say("characteristics columns: ", paste(chr11, collapse = ", "))
for (cc in chr11) say(sprintf("  [%s] e.g.: %s", cc, met11[[cc]][1]))

rpkm <- fread(file.path(raw, "GSE116250_rpkm.txt.gz"), sep = "\t", header = TRUE,
              data.table = FALSE, check.names = FALSE)
say("RPKM dim: ", nrow(rpkm), " x ", ncol(rpkm),
    " ; first cols: ", paste(colnames(rpkm)[1:2], collapse = " | "))
smp <- setdiff(colnames(rpkm), c("Gene", "Common_name"))
say("n samples in RPKM: ", length(smp), " ; series matrix: ", nrow(met11),
    " ; names match: ", all(sort(smp) == sort(met11$title)))
rp <- as.matrix(rpkm[, smp, drop = FALSE])
say("RPKM max = ", max(rp, na.rm = TRUE), " ; n cells > 1e6: ", sum(rp > 1e6, na.rm = TRUE))

sym11 <- rpkm$Common_name
say("--- transcript matching by Common_name ---")
for (g in target) say(sprintf("  %-6s n=%d  %s", g, sum(sym11 == g),
                              paste(rpkm$Gene[sym11 == g], collapse = ",")))
pick <- function(g) {
  idx <- which(sym11 == g); if (!length(idx)) return(rep(NA_real_, length(smp)))
  if (length(idx) == 1) return(rp[idx, ])
  rp[idx[which.max(rowMeans(rp[idx, , drop = FALSE], na.rm = TRUE))], ]
}
expr11 <- t(sapply(target, pick)); colnames(expr11) <- smp
say("four-gene raw RPKM range: ",
    paste(round(range(expr11[genes4, ], na.rm = TRUE), 2), collapse = " ~ "))
say("four-gene any value > 1e6: ", any(expr11[genes4, ] > 1e6, na.rm = TRUE))
expr11_log <- log2(expr11 + 1)

dis_col <- chr11[grepl("disease", apply(met11[, chr11], 2, function(z) z[1]))][1]
grp11 <- sub("^disease:\\s*", "", met11[[dis_col]])
grp11[grepl("non-failing|nonfailing", grp11)] <- "NF"
grp11[grepl("dilated|idiopathic", grp11)] <- "DCM"
grp11[grepl("ischemic", grp11)] <- "ICM"
say("group table (GSE116250):"); print(table(grp11, useNA = "ifany"))

age_col <- chr11[grepl("age", apply(met11[, chr11], 2, function(z) z[1]))][1]
sex_col <- chr11[grepl("sex", apply(met11[, chr11], 2, function(z) z[1]), ignore.case = TRUE)][1]
age11 <- as.numeric(sub("^age:\\s*", "", met11[[age_col]]))
sex11 <- sub("^Sex:\\s*|^sex:\\s*", "", met11[[sex_col]])
say("age range: ", paste(range(age11, na.rm = TRUE), collapse = "-"),
    " ; sex table: ", paste(names(table(sex11)), table(sex11), collapse = " "))

info11 <- data.frame(sample = met11$geo_accession, title = met11$title, group = grp11,
                     age = age11, sex = sex11, stringsAsFactors = FALSE)
write.csv(expr11_log, file.path(outd, "01_data", "GSE116250_gene_expr_log2rpkm.csv"))
write.csv(expr11, file.path(outd, "01_data", "GSE116250_gene_rpkm_raw.csv"))
write.csv(info11, file.path(outd, "01_data", "GSE116250_sample_info.csv"), row.names = FALSE)

# ============================================================ 3. GSE55296（备用第三队列）
say("")
say("================ GSE55296 (RNA-seq, 10 Control / 13 Dilated / 13 Ischemic) ================")
rd <- fread(file.path(raw, "GSE55296_processed_data_readme.txt"), sep = "\t",
            header = TRUE, data.table = FALSE, check.names = FALSE)
say("readme dim: ", nrow(rd), " x ", ncol(rd), " ; cols: ", paste(colnames(rd), collapse = " | "))
cnt <- fread(file.path(raw, "GSE55296_count_data.txt.gz"), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE, fill = TRUE)
say("count table dim: ", nrow(cnt), " x ", ncol(cnt))
say("first two columns: ", paste(colnames(cnt)[1:2], collapse = " | "),
    " ; example ids: ", paste(cnt[1:2, 1], collapse = ","), " / ", paste(cnt[1:2, 2], collapse = ","))

col_map <- setNames(rd[[1]], rd[[3]])              # G-number -> GSM（注意方向）
grp_map <- ifelse(grepl("Control", rd[[2]]), "NF",
                  ifelse(grepl("Isch", rd[[2]]), "ICM", "DCM"))
names(grp_map) <- rd[[3]]                           # G-number -> group
keep <- intersect(rd[[3]], colnames(cnt))
say("matched count columns: ", length(keep), " / ", nrow(rd))

sym5 <- cnt[[2]]
say("--- target genes in GSE55296 (by symbol, 2nd column) ---")
for (g in unique(c(genes4, extra))) say(sprintf("  %-6s n=%d", g, sum(sym5 == g)))

e5 <- as.matrix(cnt[, keep, drop = FALSE])
e5 <- apply(e5, 2, function(z) suppressWarnings(as.numeric(z)))
rownames(e5) <- cnt[[1]]; colnames(e5) <- keep
say("NAs after numeric coercion: ", sum(is.na(e5)), " ; per-column max NA: ",
    max(colSums(is.na(e5))))
say("columns with any NA: ", paste(names(which(colSums(is.na(e5)) > 0)), collapse = ","))
say("unmapped count columns (not in readme): ",
    paste(setdiff(grep("^G[0-9]+$", colnames(cnt), value = TRUE), rd[[3]]), collapse = ", "))
e5 <- e5[sym5 != "" & !is.na(sym5), , drop = FALSE]
sym5b <- sym5[sym5 != "" & !is.na(sym5)]
collapse5 <- function(g) {
  idx <- which(sym5b == g); if (!length(idx)) return(rep(NA_real_, ncol(e5)))
  if (length(idx) == 1) return(e5[idx, ])
  e5[idx[which.max(rowMeans(e5[idx, , drop = FALSE], na.rm = TRUE))], ]
}
expr5 <- t(sapply(unique(c(genes4, extra)), collapse5))
colnames(expr5) <- keep
say("value range: ", paste(round(range(expr5, na.rm = TRUE), 2), collapse = " ~ "))
expr5_log <- log2(expr5 + 1)
info5 <- data.frame(sample = unname(col_map[keep]), title = keep, group = unname(grp_map[keep]),
                    stringsAsFactors = FALSE)
say("group table:"); for (g in c("NF", "DCM", "ICM")) say(sprintf("  %s = %d", g, sum(info5$group == g)))
write.csv(expr5_log, file.path(outd, "01_data", "GSE55296_gene_expr_log2.csv"))
write.csv(info5, file.path(outd, "01_data", "GSE55296_sample_info.csv"), row.names = FALSE)

# ============================================================ 4. 训练集
say("")
say("================ training cohort GSE57338 ================")
tr <- fread(file.path(root, "clean_data", "GSE57338_clean_expression_matrix.csv"),
            data.table = FALSE, check.names = FALSE)
say("train dim: ", nrow(tr), " x ", ncol(tr))
if (is.character(tr[[1]])) { rn <- tr[[1]]; trm <- as.matrix(tr[, -1, drop = FALSE]); rownames(trm) <- rn
} else trm <- as.matrix(tr)
say("all four genes in train: ", all(genes4 %in% rownames(trm)))
tr4 <- trm[genes4, , drop = FALSE]
say("train mean/sd:")
for (g in genes4) say(sprintf("  %-6s mean=%.3f sd=%.3f", g, mean(tr4[g, ]), sd(tr4[g, ])))
write.csv(tr4, file.path(outd, "01_data", "train_fourgene_expr.csv"))

say(""); say("DONE"); close(logcon)
