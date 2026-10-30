# =============================================================================
# 54_sc_annotation_codet_check.R
# Recompute cross-lineage co-detection at the single-nucleus level.
#
# Why: the deposited doublet_like_index left TTN x PTPRC and LYZ x TNNT2 empty,
#      because PTPRC and TNNT2 are not in the lineage anchor set used there. The
#      cardiomyocyte x immune pair is the most likely residual doublet in heart
#      snRNA-seq, so it is retested here directly on the count matrix.
#
# Notes: rownames(so) are Ensembl IDs; gene symbols live in the `gene_symbol`
#        feature-metadata column, and Seurat reorders features on subsetting, so
#        rows are realigned by Ensembl ID before any arithmetic (same fix as 51).
#
# Output: results/annotation/tables/codet_recheck.csv
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

proj <- "E:/another data file/GEObulk analysis project"
obj  <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "04_integrated.rds")
outd <- file.path(proj, "results", "annotation", "tables")
CL_COL <- "harmony_snn_res.0.6"

pairs <- list(
  TTN_PTPRC    = c("TTN",    "PTPRC"),
  LYZ_TNNT2    = c("LYZ",    "TNNT2"),
  COL1A1_CD3D  = c("COL1A1", "CD3D"),
  PECAM1_MYH11 = c("PECAM1", "MYH11"))

so <- readRDS(obj)
cat("loaded:", nrow(so), "features x", ncol(so), "cells\n")

md <- so@meta.data
stopifnot(CL_COL %in% colnames(md))
cl_chr <- as.character(md[[CL_COL]])
lev    <- as.character(sort(unique(as.integer(cl_chr))))

gs <- so[["RNA"]][[]][["gene_symbol"]]
if (is.null(gs)) stop("gene_symbol feature metadata not found")

genes <- unique(unlist(pairs))
ens   <- rownames(so)[match(genes, gs)]
names(ens) <- genes
cat("missing:", if (anyNA(ens)) paste(genes[is.na(ens)], collapse = ",") else "none", "\n")
key <- ens[!is.na(ens)]

sub <- so[as.character(key), ]
rm(so); invisible(gc(verbose = FALSE))
if (length(Layers(sub, assay = "RNA")) > 1) sub <- JoinLayers(sub, assay = "RNA")
cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")

# order-safe alignment
ord <- match(as.character(key), rownames(cnt))
if (anyNA(ord)) stop("Ensembl IDs missing after subset")
cnt <- cnt[ord, , drop = FALSE]
rownames(cnt) <- names(key)
cat("anchors:", paste(rownames(cnt), collapse = ", "), "\n")

det <- as.matrix(cnt > 0)
rm(cnt, sub); invisible(gc(verbose = FALSE))

res <- do.call(rbind, lapply(names(pairs), function(p) {
  g <- pairs[[p]]
  if (!all(g %in% rownames(det))) return(NULL)
  a <- det[g[1], ]; b <- det[g[2], ]
  do.call(rbind, lapply(lev, function(k) {
    idx <- which(cl_chr == k)
    pa <- mean(a[idx]); pb <- mean(b[idx]); po <- mean(a[idx] & b[idx])
    data.frame(cluster = k, pair = p, n_cells = length(idx),
               pct_gene1 = round(pa * 100, 3), pct_gene2 = round(pb * 100, 3),
               observed_pct = round(po * 100, 3),
               expected_pct = round(pa * pb * 100, 3),
               ratio = round(po / (pa * pb), 3), stringsAsFactors = FALSE)
  }))
}))

res$verdict <- ifelse(res$observed_pct >= 5 & res$ratio >= 2,
                      "enriched beyond independence - doublet candidate",
                      "consistent with independent co-detection (ambient)")
res <- res[order(res$pair, -res$observed_pct), ]
write.csv(res, file.path(outd, "codet_recheck.csv"), row.names = FALSE)
cat("rows:", nrow(res), "| max observed:", max(res$observed_pct),
    "| max ratio:", max(res$ratio, na.rm = TRUE), "\n")
print(res[res$observed_pct >= 1, ], row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(proj, "results", "annotation", "sessionInfo_54.txt"))
