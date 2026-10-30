# =============================================================================
# 64_sc_fourgene_extract.R
# Stage 5 (A3) - slim extraction of the four ML genes + metadata.
#
# Input : results/Singlecel/objects/03_merged/05_final_annotated.rds
# Output: results/Singlecel/fourgene/objects/fourgene_slim.rds
#         results/Singlecel/fourgene/tables/64_meta_probe.csv
#         results/Singlecel/fourgene/tables/64_sample_group.csv
#         results/Singlecel/fourgene/64_extract.log
#
# Why a slim object: the full object is 1.6 GB on disk. Every downstream plot
# only needs the four features + meta.data + the harmony UMAP, so the heavy
# read happens exactly once and everything else runs on ~100 MB.
#
# Note: the RNA assay is keyed by Ensembl ID; gene symbols live in the
#       `gene_symbol` feature metadata (same convention as 60/61).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
obj_f <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
out_d <- file.path(proj, "results", "Singlecel", "fourgene")
tab_d <- file.path(out_d, "tables")
obj_d <- file.path(out_d, "objects")
for (d in c(tab_d, obj_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file.path(out_d, "64_extract.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")

lg("=== 64 four-gene extraction | seed ", SEED, " ===")
t0 <- Sys.time()
so <- readRDS(obj_f)
lg(sprintf("loaded %d genes x %d nuclei in %.1f min",
           nrow(so), ncol(so), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
lg("assays: ", paste(names(so@assays), collapse = ","),
   " | default: ", DefaultAssay(so))
lg("reductions: ", paste(names(so@reductions), collapse = ","))
lg("meta columns: ", paste(colnames(so@meta.data), collapse = ","))

# ------------------------------------------------------- feature lookup -----
gs_map <- so[["RNA"]][[]][["gene_symbol"]]
idx <- match(GENES, gs_map)
names(idx) <- GENES
if (any(is.na(idx))) {
  lg("!! missing symbols: ", paste(GENES[is.na(idx)], collapse = ","))
}
ens <- rownames(so)[idx]
names(ens) <- GENES
ens <- ens[!is.na(ens)]
lg("resolved: ", paste(sprintf("%s=%s", names(ens), ens), collapse = " | "))

# duplicate Ensembl sanity check (a symbol may map to several features)
dup <- table(gs_map[gs_map %in% GENES])
lg("feature count per symbol: ", paste(sprintf("%s:%d", names(dup), as.integer(dup)), collapse = " | "))

# ------------------------------------------------------------- matrices -----
cnt <- GetAssayData(so, assay = "RNA", layer = "counts")[unname(ens), , drop = FALSE]
dat <- GetAssayData(so, assay = "RNA", layer = "data")[unname(ens), , drop = FALSE]
rownames(cnt) <- names(ens)[match(rownames(cnt), unname(ens))]
rownames(dat) <- rownames(cnt)
lg("counts layer: ", nrow(cnt), " x ", ncol(cnt),
   " | max count ", max(cnt), " | data max ", round(max(dat), 3))

md <- so@meta.data
um <- Embeddings(so, "umap.harmony")
stopifnot(nrow(um) == ncol(so), identical(rownames(um), colnames(so)))
colnames(um) <- c("UMAP_1", "UMAP_2")

# ----------------------------------------------------------- meta probe -----
probe <- data.frame(
  column = colnames(md),
  class  = sapply(md, function(x) class(x)[1]),
  n_unique = sapply(md, function(x) length(unique(x))),
  example = sapply(md, function(x) paste(head(unique(as.character(x)), 4), collapse = "|")),
  stringsAsFactors = FALSE)
write.csv(probe, file.path(tab_d, "64_meta_probe.csv"), row.names = FALSE)

stopifnot(all(c("sample", "group", "celltype") %in% colnames(md)))
sg <- unique(md[, c("sample", "group", "condition")])
sg <- sg[order(sg$group, sg$sample), ]
write.csv(sg, file.path(tab_d, "64_sample_group.csv"), row.names = FALSE)
lg("sample x group:")
print(table(md$sample, md$group))

# --------------------------------------------------------- detection --------
det <- data.frame(
  gene = rownames(cnt),
  n_cells_total = ncol(cnt),
  n_detect_total = as.integer(Matrix::rowSums(cnt > 0)),
  pct_total = round(Matrix::rowMeans(cnt > 0) * 100, 2),
  mean_count = round(Matrix::rowMeans(cnt), 4),
  stringsAsFactors = FALSE)
write.csv(det, file.path(tab_d, "64_detection_overall.csv"), row.names = FALSE)
lg("overall detection: ", paste(sprintf("%s %.2f%%", det$gene, det$pct_total), collapse = " | "))

# ------------------------------------------------------------- save ---------
slim <- list(
  counts = cnt, data = dat, meta = md, umap = um,
  genes = names(ens), ensembl = ens,
  provenance = list(source = obj_f, n_cells = ncol(so), n_genes = nrow(so),
                    seed = SEED, built = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
saveRDS(slim, file.path(obj_d, "fourgene_slim.rds"))
lg("slim object saved: ", round(file.size(file.path(obj_d, "fourgene_slim.rds")) / 1e6, 1), " MB")

writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_64.txt"))
lg("=== done ===")
