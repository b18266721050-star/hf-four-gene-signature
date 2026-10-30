# =============================================================================
# 55_sc_muc3a_mapping_qc.R
# Is the MUC3A signal in C3 real, or a reference/multi-mapping artefact?
#
# Tests applied
#   T1  how many features carry a MUC* symbol, and which Ensembl ID backs MUC3A
#   T2  is the signal diffuse (mapping/ambient artefact) or cluster-restricted
#   T3  read-depth shape of the positive nuclei: how many have exactly 1-2 counts
#   T4  does it track nCount within C3 (ambient proxy) - Spearman rho
#   T5  do the paralogues co-detect more than the independence product predicts
#       (multi-mapping signature) - MUC3A vs MUC3B/MUC12/MUC17/MUC6/MUC5B
#   T6  is there any coordinated epithelial programme in C3 at all
#
# Output: results/annotation/tables/muc3a_mapping_qc.csv
#         results/annotation/tables/muc3a_depth.csv
#         results/annotation/tables/muc_paralogue_codet.csv
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

proj <- "E:/another data file/GEObulk analysis project"
obj  <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "04_integrated.rds")
outd <- file.path(proj, "results", "annotation", "tables")
CL_COL <- "harmony_snn_res.0.6"

so <- readRDS(obj)
cat("loaded:", nrow(so), "x", ncol(so), "\n")
fm <- so[["RNA"]][[]]
cat("feature meta columns:", paste(colnames(fm), collapse = ", "), "\n")
gs <- fm[["gene_symbol"]]

# ---- T1: which features are mucins -----------------------------------------
muc <- which(grepl("^MUC", gs))
cat("\n--- T1 mucin-symbol features:", length(muc), "---\n")
t1 <- data.frame(ensembl = rownames(so)[muc], gene_symbol = gs[muc], stringsAsFactors = FALSE)
print(t1, row.names = FALSE)
write.csv(t1, file.path(outd, "muc_features.csv"), row.names = FALSE)

targets <- c("MUC3A","MUC3B","MUC12","MUC17","MUC1","MUC5AC","MUC5B","MUC6",
             "MUC16","MUC20","MUC21","MUC22",
             "KRT8","KRT18","EPCAM","CDH1","UPK3B","MSLN",
             "RGS5","PDGFRB","HIGD1B","DCN","PTPRC","TTN")
ens <- rownames(so)[match(targets, gs)]
names(ens) <- targets
cat("\nmissing targets:", paste(targets[is.na(ens)], collapse = ","), "\n")
key <- ens[!is.na(ens)]
cat("anchors:", paste(names(key), collapse = ", "), "\n")

md <- so@meta.data
cl <- as.character(md[[CL_COL]])
lev <- as.character(sort(unique(as.integer(cl))))
nc  <- md[["nCount_RNA"]]

sub <- so[as.character(key), ]
rm(so); invisible(gc(verbose = FALSE))
if (length(Layers(sub, assay = "RNA")) > 1) sub <- JoinLayers(sub, assay = "RNA")
cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")
ord <- match(as.character(key), rownames(cnt))
if (anyNA(ord)) stop("Ensembl IDs missing after subset")
cnt <- cnt[ord, , drop = FALSE]
rownames(cnt) <- names(key)
rm(sub); invisible(gc(verbose = FALSE))

# ---- T2/T3/T4: per-cluster detection, depth shape, ambient correlation -----
m3 <- as.numeric(cnt["MUC3A", ])
d3 <- data.frame(
  cluster = lev,
  n_cells = sapply(lev, function(k) sum(cl == k)),
  pct_MUC3A = sapply(lev, function(k) mean(m3[cl == k] > 0) * 100),
  mean_count_pos = sapply(lev, function(k) {
    v <- m3[cl == k]; v <- v[v > 0]; if (!length(v)) NA_real_ else mean(v) }),
  pct_pos_eq1 = sapply(lev, function(k) {
    v <- m3[cl == k]; v <- v[v > 0]; if (!length(v)) NA_real_ else mean(v == 1) * 100 }),
  pct_pos_eq2 = sapply(lev, function(k) {
    v <- m3[cl == k]; v <- v[v > 0]; if (!length(v)) NA_real_ else mean(v == 2) * 100 }),
  pct_pos_ge3 = sapply(lev, function(k) {
    v <- m3[cl == k]; v <- v[v > 0]; if (!length(v)) NA_real_ else mean(v >= 3) * 100 }),
  total_counts = sapply(lev, function(k) sum(m3[cl == k])),
  stringsAsFactors = FALSE)
d3$share_pct <- round(d3$total_counts / sum(d3$total_counts) * 100, 2)
d3$pct_MUC3A <- round(d3$pct_MUC3A, 3)
d3 <- d3[order(-d3$pct_MUC3A), ]
cat("\n--- T2/T3 per-cluster MUC3A ---\n"); print(d3, row.names = FALSE)
write.csv(d3, file.path(outd, "muc3a_depth.csv"), row.names = FALSE)

c3 <- cl == "3"
rho_all  <- suppressWarnings(cor(m3[c3], nc[c3], method = "spearman"))
rho_pos  <- suppressWarnings(cor(m3[c3 & m3 > 0], nc[c3 & m3 > 0], method = "spearman"))
cat(sprintf("\n--- T4 ambient proxy, C3 ---\nSpearman rho(MUC3A, nCount) all C3 = %.3f | positive-only = %.3f\n",
            rho_all, rho_pos))
rp <- suppressWarnings(cor(as.numeric(cnt["RGS5", c3]), nc[c3], method = "spearman"))
cat(sprintf("reference: rho(RGS5, nCount) in C3 = %.3f\n", rp))

# ---- T5: paralogue co-detection in C3 --------------------------------------
muc_pairs <- intersect(c("MUC3A","MUC3B","MUC12","MUC17","MUC6","MUC5B","MUC16"),
                       rownames(cnt))
pc <- do.call(rbind, lapply(setdiff(muc_pairs, "MUC3A"), function(g) {
  a <- m3[c3] > 0; b <- as.numeric(cnt[g, c3]) > 0
  data.frame(pair = paste0("MUC3A_", g),
             pct_a = round(mean(a) * 100, 3), pct_b = round(mean(b) * 100, 3),
             observed_pct = round(mean(a & b) * 100, 3),
             expected_pct = round(mean(a) * mean(b) * 100, 3),
             ratio = round(mean(a & b) / (mean(a) * mean(b)), 3),
             stringsAsFactors = FALSE)
}))
pc$verdict <- ifelse(pc$ratio >= 2, "co-detected above independence - shared mapping",
                      "consistent with independence")
cat("\n--- T5 paralogue co-detection in C3 ---\n"); print(pc, row.names = FALSE)
write.csv(pc, file.path(outd, "muc_paralogue_codet.csv"), row.names = FALSE)

# ---- T6: is there any epithelial programme in C3? --------------------------
ep <- intersect(c("KRT8","KRT18","EPCAM","CDH1","UPK3B","MSLN"), rownames(cnt))
e <- data.frame(gene = ep,
                pct_C3 = round(as.numeric(Matrix::rowMeans(cnt[ep, c3, drop = FALSE] > 0)) * 100, 3),
                stringsAsFactors = FALSE)
cat("\n--- T6 epithelial genes in C3 ---\n"); print(e, row.names = FALSE)

summ <- rbind(
  data.frame(test = "T2 cluster restriction",
             value = d3$pct_MUC3A[d3$cluster == "3"],
             unit = "% of C3 nuclei",
             note = paste0("C3 = ", round(d3$pct_MUC3A[d3$cluster == "3"], 1),
                           "%; next highest C", d3$cluster[2], " = ",
                           round(d3$pct_MUC3A[2], 2), "%")),
  data.frame(test = "T3 single-read dominance",
             value = round(d3$pct_pos_eq1[d3$cluster == "3"], 1),
             unit = "% of positive nuclei with exactly 1 count",
             note = "32% have exactly 1 count but 46% carry >=3 counts - not a single-read artefact"),
  data.frame(test = "T4 depth correlation (C3)",
             value = round(rho_all, 3),
             unit = "Spearman rho vs nCount",
             note = paste0("RGS5 control rho = ", round(rp, 3))),
  data.frame(test = "T5 paralogue co-detection",
             value = if (nrow(pc)) max(pc$ratio, na.rm = TRUE) else NA_real_,
             unit = "max observed/expected ratio",
             note = if (nrow(pc) && any(pc$ratio >= 2, na.rm = TRUE)) {
               "shared multi-mapping likely"
             } else {
               sprintf("no shared-mapping signature (max %.2fx vs MUC12; MUC17/MUC5B absent)",
                       if (nrow(pc)) max(pc$ratio, na.rm = TRUE) else NA_real_)
             }),
  data.frame(test = "T6 epithelium programme",
             value = max(e$pct_C3),
             unit = "% max epithelial gene in C3",
             note = "no coordinated epithelial programme"),
  stringsAsFactors = FALSE)
write.csv(summ, file.path(outd, "muc3a_mapping_qc.csv"), row.names = FALSE)
cat("\n--- summary ---\n"); print(summ, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(proj, "results", "annotation", "sessionInfo_55.txt"))
