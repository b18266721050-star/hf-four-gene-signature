# =============================================================================
# 56_sc_muc3a_origin.R
# Where does the C3 MUC3A signal come from?
#
# Script 55 showed it is NOT diffuse (66.4% of C3 vs <8% elsewhere), NOT single
# read (only 32% of positives have exactly 1 count), NOT depth driven
# (rho 0.36 vs RGS5 0.46) and NOT explained by mucin paralogue multi-mapping
# (no co-detection excess; MUC17/MUC5B absent). Ensembl confirms ENSG00000169894
# is the canonical MUC3A (7:100,949,414-100,968,348, protein_coding).
#
# Remaining candidates tested here
#   A  sample / batch driven (a few donors or one run)
#   B  marks a coherent pericyte sub-state (real biology)
#   C  spatially random within C3 (low-level noise, no sub-structure)
#
# Output: results/annotation/tables/muc3a_by_sample.csv
#         results/annotation/tables/muc3a_vs_neg.csv
#         results/annotation/tables/muc3a_spatial.csv
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

proj <- "E:/another data file/GEObulk analysis project"
obj  <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "04_integrated.rds")
outd <- file.path(proj, "results", "annotation", "tables")
CL_COL <- "harmony_snn_res.0.6"

so <- readRDS(obj)
md <- so@meta.data
keep <- as.character(md[[CL_COL]]) == "3"
cat("C3 nuclei:", sum(keep), "\n")

genes <- c("MUC3A","KRT8","RGS5","PDGFRB","HIGD1B","KCNJ8","EGFLAM","NDUFA4L2",
           "NOTCH3","ABCC9","GUCY1A2","ACTA2","MYH11","DCN")
gs <- so[["RNA"]][[]][["gene_symbol"]]
ens <- rownames(so)[match(genes, gs)]
names(ens) <- genes
key <- ens[!is.na(ens)]
cat("missing:", paste(genes[is.na(ens)], collapse = ","), "\n")

um <- Embeddings(so, "umap.harmony")[keep, , drop = FALSE]
samp <- as.character(md[["sample"]][keep])
cond <- as.character(md[["condition"]][keep])
nC   <- md[["nCount_RNA"]][keep]
nF   <- md[["nFeature_RNA"]][keep]
pmt  <- md[["percent.mt"]][keep]

sub <- so[as.character(key), keep]
rm(so); invisible(gc(verbose = FALSE))
if (length(Layers(sub, assay = "RNA")) > 1) sub <- JoinLayers(sub, assay = "RNA")
cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")
ord <- match(as.character(key), rownames(cnt))
cnt <- cnt[ord, , drop = FALSE]
rownames(cnt) <- names(key)
det <- cnt > 0
m3  <- as.numeric(cnt["MUC3A", ])
pos <- m3 > 0
cat("MUC3A positive in C3:", sum(pos), sprintf("(%.1f%%)\n", mean(pos) * 100))

# ---- A: sample / condition -------------------------------------------------
sp <- unique(samp)
by_s <- data.frame(sample = sp,
                   n_cells = sapply(sp, function(s) sum(samp == s)),
                   pct_MUC3A = sapply(sp, function(s) round(mean(m3[samp == s] > 0) * 100, 2)),
                   mean_count_pos = sapply(sp, function(s) {
                     v <- m3[samp == s]; v <- v[v > 0]
                     if (!length(v)) NA_real_ else round(mean(v), 2) }),
                   stringsAsFactors = FALSE)
cd <- unique(cond)
by_c <- data.frame(sample = paste0("BY_", cd),
                   n_cells = sapply(cd, function(k) sum(cond == k)),
                   pct_MUC3A = sapply(cd, function(k) round(mean(m3[cond == k] > 0) * 100, 2)),
                   mean_count_pos = sapply(cd, function(k) {
                     v <- m3[cond == k]; v <- v[v > 0]
                     if (!length(v)) NA_real_ else round(mean(v), 2) }),
                   stringsAsFactors = FALSE)
by_s <- by_s[order(-by_s$pct_MUC3A), ]
cat("\n--- A: per sample ---\n"); print(by_s, row.names = FALSE)
cat("\n--- A: per condition ---\n"); print(by_c, row.names = FALSE)
write.csv(rbind(by_s, by_c), file.path(outd, "muc3a_by_sample.csv"), row.names = FALSE)
cat(sprintf("samples with any MUC3A+ cell: %d / %d\n",
            sum(by_s$pct_MUC3A > 0), nrow(by_s)))

# ---- B: does MUC3A split a pericyte sub-state? -----------------------------
mk <- setdiff(rownames(cnt), "MUC3A")
cmp <- data.frame(
  gene = mk,
  pct_in_MUC3A_pos = round(as.numeric(Matrix::rowMeans(det[mk, pos, drop = FALSE])) * 100, 2),
  pct_in_MUC3A_neg = round(as.numeric(Matrix::rowMeans(det[mk, !pos, drop = FALSE])) * 100, 2),
  stringsAsFactors = FALSE)
cmp$delta <- round(cmp$pct_in_MUC3A_pos - cmp$pct_in_MUC3A_neg, 2)
qc <- data.frame(gene = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                 pct_in_MUC3A_pos = c(round(mean(nC[pos]), 1), round(mean(nF[pos]), 1),
                                      round(mean(pmt[pos]), 3)),
                 pct_in_MUC3A_neg = c(round(mean(nC[!pos]), 1), round(mean(nF[!pos]), 1),
                                      round(mean(pmt[!pos]), 3)),
                 stringsAsFactors = FALSE)
qc$delta <- round(qc$pct_in_MUC3A_pos - qc$pct_in_MUC3A_neg, 2)
out_b <- rbind(cmp, qc)
cat("\n--- B: MUC3A+ vs MUC3A- inside C3 ---\n"); print(out_b, row.names = FALSE)
write.csv(out_b, file.path(outd, "muc3a_vs_neg.csv"), row.names = FALSE)

# ---- C: spatial coherence on the harmony UMAP ------------------------------
grid_n <- 30
gx <- cut(um[, 1], breaks = grid_n)
gy <- cut(um[, 2], breaks = grid_n)
bin <- paste(gx, gy)
tb  <- table(bin)
ok  <- names(tb)[tb >= 20]
p_hat <- mean(pos)
k <- sapply(ok, function(b) sum(pos[bin == b]))
n <- sapply(ok, function(b) sum(bin == b))
X2 <- sum((k - n * p_hat)^2 / (n * p_hat * (1 - p_hat)))
df <- length(ok)
spat <- data.frame(n_bins_used = df, n_cells = sum(n), p_global = round(p_hat, 4),
                   pearson_X2 = round(X2, 1), dispersion = round(X2 / df, 2),
                   stringsAsFactors = FALSE)
cat("\n--- C: spatial dispersion on UMAP (binomial null) ---\n"); print(spat, row.names = FALSE)
cat("dispersion  ~1 random | >1.5 spatially clustered\n")
write.csv(spat, file.path(outd, "muc3a_spatial.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(proj, "results", "annotation", "sessionInfo_56.txt"))
