# =============================================================================
# 50_sc_annotation_panel.R
# Stage 3 - Step 1: quantitative marker-panel evidence for cell-type annotation.
#
# Reads the integrated object ONCE, extracts a curated canonical marker panel
# (average log-normalised expression + detection rate per cluster), per-cluster
# QC / doublet metrics, and lineage co-detection (doublet-like proxy).
# A small RDS with embedding + metadata is saved so that figures can be redrawn
# later without reloading the 7.8 GB object.
#
# Outputs -> results/Singlecel/annotation/
#   tables/panel_avg_expr.csv, panel_pct_expr.csv, panel_zscore.csv,
#   cluster_qc.csv, lineage_module_score.csv, doublet_like_index.csv
#   objects/annotation_small.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

SEED <- 20261028
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel")
out_d <- file.path(sdir, "annotation")
tab_d <- file.path(out_d, "tables")
obj_d <- file.path(out_d, "objects")
fig_d <- file.path(out_d, "figures")
for (d in c(out_d, tab_d, obj_d, fig_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file(file.path(out_d, "50_annotation_panel.log"), "wt")
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf); flush(lgf)
}
options(error = function() {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ERROR: ", geterrmessage(), "\n"),
      file = lgf, append = TRUE); flush(lgf); q(status = 1)
})

CL_COL <- "harmony_snn_res.0.6"
RED    <- "umap.harmony"
f_in   <- file.path(sdir, "objects", "03_merged", "04_integrated.rds")

# ------------------------------------------------------------- panel --------
panel <- list(
  CM     = c("TTN","MYH6","MYH7","TNNT2","TNNI3","ACTC1","MYL2","MYL3","NPPA","NPPB",
             "RYR2","PLN","MYBPC3","ACTN2","TRDN","CHRM2","ALPK2","ESRRG","RBM20","MYO18B"),
  FB     = c("DCN","LUM","COL1A1","COL1A2","COL3A1","PDGFRA","FBN1","MFAP5","SERPINE2",
             "POSTN","ABCA9","CFD","FBLN1","COL6A3"),
  EC     = c("PECAM1","CDH5","VWF","ERG","FLT1","KDR","EMCN","AQP1","CA4","BTNL9",
             "ADGRL4","F8","ACKR1","NR2F2"),
  LEC    = c("CCL21","PROX1","LYVE1","MMRN1","FLT4","TFF3","PKHD1L1","STAB2"),
  ENDO   = c("NPR3","HAPLN1","NRG1","NFATC1"),
  PC     = c("RGS5","KCNJ8","ABCC9","PDGFRB","HIGD1B","NDUFA4L2","NOTCH3","EGFLAM","GUCY1A2"),
  SMC    = c("MYH11","ACTA2","TAGLN","CNN1","LMOD1","CARMN","MYOCD","DES","NTRK3","SLIT3"),
  MYE    = c("PTPRC","LYZ","CD68","CD14","CSF1R","C1QA","MRC1","CD163","F13A1","MS4A4A",
             "LILRB5","ITGAX","FCN1","S100A8","STAB1","MS4A7","KYNU"),
  TNK    = c("CD3D","CD3E","CD2","IL7R","CD8A","NKG7","GNLY","KLRD1","SKAP1","ITK",
             "BCL11B","THEMIS","CD247","TRBC2"),
  BPL    = c("CD79A","MS4A1","CD19","JCHAIN","MZB1","DERL3","IGHG1"),
  MAST   = c("CPA3","TPSAB1","TPSB2","KIT","MS4A2","HDC","CMA1"),
  NEUR   = c("NRXN1","PLP1","S100B","SOX10","MPZ","PTPRZ1","XKR4","GPM6B","CDH19"),
  ADIPO  = c("PLIN1","ADIPOQ","GPD1","CIDEA","GPAM","LEP"),
  PROL   = c("MKI67","TOP2A","RRM2","BUB1B","CENPF"),
  MESO   = c("UPK3B","MSLN","WT1","KRT19"),
  AMBI   = c("MUC3A","KRT8","KRT18","EPCAM","HBB","HBA1","ALAS2")
)
panel_genes <- unique(unlist(panel, use.names = FALSE))
mod_of <- setNames(rep(names(panel), lengths(panel)), unlist(panel, use.names = FALSE))

anchors <- c(CM = "TTN", FB = "DCN", EC = "PECAM1", PC = "RGS5", SMC = "MYH11",
             MYE = "LYZ", TNK = "CD3D", MAST = "CPA3", ADIPO = "PLIN1", NEUR = "PLP1")

lg("=== 50 annotation panel (", length(panel_genes), " genes, ",
   length(panel), " modules) ===")
lg("reading ", f_in)
t0 <- Sys.time()
so <- readRDS(f_in)
lg(sprintf("loaded in %.1f min | %d genes x %d cells", 
           as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(so), ncol(so)))

md <- so@meta.data
stopifnot(CL_COL %in% colnames(md))
cl  <- as.character(md[[CL_COL]])
clu <- sort(unique(cl))
lg("clusters: ", paste(clu, collapse = ","))

emb <- Embeddings(so, RED)
gs  <- so[["RNA"]][[]][["gene_symbol"]]
if (is.null(gs)) stop("gene_symbol feature metadata not found")

ens <- rownames(so)[match(panel_genes, gs)]
names(ens) <- panel_genes
missing <- panel_genes[is.na(ens)]
lg("panel genes not found: ", if (length(missing)) paste(missing, collapse = ",") else "none")

# small per-cell annotation frame (kept for later figure scripts)
meta_small <- data.frame(
  cell = colnames(so),
  cluster = cl,
  sample = as.character(md[["sample"]]),
  donor = as.character(md[["donor"]]),
  condition = as.character(md[["condition"]]),
  nCount = md[["nCount_RNA"]],
  nFeature = md[["nFeature_RNA"]],
  pct_mt = md[["percent.mt"]],
  dbl_score = if (!is.null(md[["doublet_score_source"]])) md[["doublet_score_source"]] else NA_real_,
  stringsAsFactors = FALSE
)
meta_small$UMAP_1 <- emb[, 1]
meta_small$UMAP_2 <- emb[, 2]

keep <- ens[!is.na(ens)]
sub  <- so[as.character(keep), ]
rm(so); invisible(gc(verbose = FALSE))
lg("subset kept ", nrow(sub), " genes; full object released")

# layers were split per sample at merge time -> join on the small subset only
if (length(Layers(sub, assay = "RNA")) > 1) {
  sub <- JoinLayers(sub, assay = "RNA")
  lg("layers joined on panel subset: ", paste(Layers(sub, assay = "RNA"), collapse = ","))
}
cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")
dat <- GetAssayData(sub, assay = "RNA", layer = "data")
sym <- as.character(sub[["RNA"]][[]][["gene_symbol"]])
rownames(cnt) <- sym; rownames(dat) <- sym
rm(sub); invisible(gc(verbose = FALSE))
lg("matrices: counts ", nrow(cnt), "x", ncol(cnt), " | data nnz = ", length(dat@x))

# --------------------------------------------- per-cluster avg / pct ---------
avg <- vapply(clu, function(k) Matrix::rowMeans(dat[, cl == k, drop = FALSE]), numeric(nrow(dat)))
pct <- vapply(clu, function(k) Matrix::rowMeans(cnt[, cl == k, drop = FALSE] > 0), numeric(nrow(cnt)))
colnames(avg) <- clu; colnames(pct) <- clu
z <- t(scale(t(avg)))
z[is.na(z)] <- 0
colnames(z) <- clu

write.csv(data.frame(gene = rownames(avg), module = mod_of[rownames(avg)], avg, check.names = FALSE),
          file.path(tab_d, "panel_avg_expr.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(pct), module = mod_of[rownames(pct)], pct, check.names = FALSE),
          file.path(tab_d, "panel_pct_expr.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(z), module = mod_of[rownames(z)], z, check.names = FALSE),
          file.path(tab_d, "panel_zscore.csv"), row.names = FALSE)
lg("panel tables written")

# ------------------------------------------- module score per cluster --------
mods <- names(panel)
mod_score <- sapply(clu, function(k) sapply(mods, function(m) {
  g <- intersect(panel[[m]], rownames(z))
  if (!length(g)) return(NA_real_)
  mean(z[g, k])
}))
mod_score <- as.data.frame(mod_score)
mod_score <- data.frame(module = rownames(mod_score), mod_score, check.names = FALSE)
write.csv(mod_score, file.path(tab_d, "lineage_module_score.csv"), row.names = FALSE)

# top modules per cluster (z-based)
top_mod <- sapply(clu, function(k) {
  v <- mod_score[[k]]; names(v) <- mod_score$module
  v <- sort(v, decreasing = TRUE)
  paste0(names(v)[1:3], "(", sprintf("%.2f", v[1:3]), ")", collapse = " | ")
})
write.csv(data.frame(cluster = clu, top3_modules = top_mod, check.names = FALSE),
          file.path(tab_d, "cluster_top_modules.csv"), row.names = FALSE)
lg("module scores written")

# ------------------------------------------------- per-cluster QC ------------
qc <- data.frame(
  cluster = clu,
  n_cells = as.integer(table(cl)[clu]),
  mean_nCount = round(tapply(meta_small$nCount, cl, mean)[clu], 1),
  median_nCount = round(tapply(meta_small$nCount, cl, stats::median)[clu], 1),
  mean_nFeature = round(tapply(meta_small$nFeature, cl, mean)[clu], 1),
  mean_pct_mt = round(tapply(meta_small$pct_mt, cl, mean)[clu], 4),
  mean_dbl_score = if (all(is.na(meta_small$dbl_score))) NA_real_ else
    round(tapply(meta_small$dbl_score, cl, mean, na.rm = TRUE)[clu], 4),
  pct_disease = round(tapply(meta_small$condition == "Disease", cl, mean)[clu] * 100, 2),
  n_samples = as.integer(tapply(meta_small$sample, cl, function(x) length(unique(x)))[clu]),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(tab_d, "cluster_qc.csv"), row.names = FALSE)
lg("cluster QC written")

# ------------------------------------- lineage co-detection (doublet proxy) --
det <- cnt[intersect(anchors, rownames(cnt)), , drop = FALSE] > 0
n_anchor <- Matrix::colSums(det)
lg("cells detecting >=2 lineage anchors: ", sprintf("%.2f%%", mean(n_anchor >= 2) * 100))

pairs <- list(c("TTN","PTPRC"), c("COL1A1","CD3D"), c("PECAM1","MYH11"),
              c("LYZ","TNNT2"), c("DCN","CD3D"))
pair_lab <- vapply(pairs, function(p) paste(p, collapse = "_"), character(1))
pair_mat <- sapply(clu, function(k) {
  idx <- which(cl == k)
  vapply(pairs, function(p) {
    if (!all(p %in% rownames(det))) return(NA_real_)
    round(mean(det[p[1], idx] & det[p[2], idx]) * 100, 3)
  }, numeric(1))
})
rownames(pair_mat) <- pair_lab
pair_df <- as.data.frame(t(pair_mat))
colnames(pair_df) <- paste0("codet_", pair_lab)
dl <- data.frame(
  cluster = clu,
  pct_ge2_anchor = round(tapply(n_anchor >= 2, cl, mean)[clu] * 100, 2),
  mean_n_anchor = round(tapply(n_anchor, cl, mean)[clu], 3),
  pair_df,
  check.names = FALSE, row.names = NULL
)
write.csv(dl, file.path(tab_d, "doublet_like_index.csv"), row.names = FALSE)
lg("doublet-like index written")

saveRDS(list(meta = meta_small, avg = avg, pct = pct, z = z,
             mod_score = mod_score, qc = qc, dl = dl,
             panel = panel, cl_col = CL_COL, red = RED),
        file.path(obj_d, "annotation_small.rds"))
lg("saved annotation_small.rds (", sprintf("%.1f MB",
   file.size(file.path(obj_d, "annotation_small.rds")) / 1024^2), ")")

writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_50.txt"))
lg("=== 50 END ===")
close(lgf)
