# =============================================================================
# 51_sc_annotation_panel.R  (v2, order-safe)
# Stage 3 - Step 1: quantitative marker-panel evidence for cell-type annotation.
#
# v1 bug: after `so[features, ]` Seurat may return features in genome order, and
# gene symbols taken from feature metadata were assigned positionally -> gene
# rows were mislabelled. v2 aligns everything by Ensembl ID and asserts it.
#
# Reads the integrated object ONCE, extracts a curated canonical marker panel
# (mean log-normalised expression + detection rate per cluster), per-cluster
# QC / doublet metrics, lineage co-detection, and a small RDS for later figures.
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
for (d in c(out_d, tab_d, obj_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file(file.path(out_d, "51_annotation_panel.log"), "wt")
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
mod_of <- mod_of[!duplicated(names(mod_of))]

anchors <- c(CM = "TTN", FB = "DCN", EC = "PECAM1", PC = "RGS5", SMC = "MYH11",
             MYE = "LYZ", TNK = "CD3D", MAST = "CPA3", ADIPO = "PLIN1", NEUR = "PLP1")

lg("=== 51 annotation panel v2 (", length(panel_genes), " genes, ",
   length(panel), " modules) ===")
lg("reading ", f_in)
t0 <- Sys.time()
so <- readRDS(f_in)
lg(sprintf("loaded in %.1f min | %d genes x %d cells",
           as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(so), ncol(so)))

md <- so@meta.data
stopifnot(CL_COL %in% colnames(md))
cl_chr <- as.character(md[[CL_COL]])
clu    <- as.character(sort(unique(as.integer(cl_chr))))   # numeric order, explicit
lg("clusters (numeric order): ", paste(clu, collapse = ","))

emb <- Embeddings(so, RED)
gs  <- so[["RNA"]][[]][["gene_symbol"]]
if (is.null(gs)) stop("gene_symbol feature metadata not found")

ens <- rownames(so)[match(panel_genes, gs)]
names(ens) <- panel_genes
missing <- panel_genes[is.na(ens)]
lg("panel genes not found: ", if (length(missing)) paste(missing, collapse = ",") else "none")
key <- ens[!is.na(ens)]

meta_small <- data.frame(
  cell = colnames(so),
  cluster = cl_chr,
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

sub <- so[as.character(key), ]
rm(so); invisible(gc(verbose = FALSE))
lg("subset kept ", nrow(sub), " features")

if (length(Layers(sub, assay = "RNA")) > 1) {
  sub <- JoinLayers(sub, assay = "RNA")
  lg("layers joined: ", paste(Layers(sub, assay = "RNA"), collapse = ","))
}
cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")
dat <- GetAssayData(sub, assay = "RNA", layer = "data")

# ---- order-safe alignment: reorder every matrix to the panel key order ----
ord_c <- match(as.character(key), rownames(cnt))
ord_d <- match(as.character(key), rownames(dat))
if (anyNA(ord_c) || anyNA(ord_d)) stop("panel Ensembl IDs missing from the assay after subset")
if (!identical(ord_c, seq_along(ord_c))) lg("NOTE: counts rows were reordered by Seurat -> realigned")
cnt <- cnt[ord_c, , drop = FALSE]
dat <- dat[ord_d, , drop = FALSE]
symv <- names(key)
rownames(cnt) <- symv
rownames(dat) <- symv
rm(sub); invisible(gc(verbose = FALSE))
lg("matrices aligned: ", nrow(cnt), " genes x ", ncol(cnt), " cells")

# --------------------------------------------- per-cluster avg / pct ---------
avg <- vapply(clu, function(k) Matrix::rowMeans(dat[, cl_chr == k, drop = FALSE]), numeric(nrow(dat)))
pct <- vapply(clu, function(k) Matrix::rowMeans(cnt[, cl_chr == k, drop = FALSE] > 0), numeric(nrow(cnt)))
colnames(avg) <- paste0("C", clu)
colnames(pct) <- paste0("C", clu)
rownames(avg) <- symv
rownames(pct) <- symv
z <- t(scale(t(avg)))
z[is.na(z)] <- 0
colnames(z) <- colnames(avg)

# ---- self-check against the independent marker table ----------------------
chk <- read.csv(file.path(sdir, "tables", "04_markers_all.csv"), stringsAsFactors = FALSE)
chk <- chk[!is.na(chk$gene_symbol), ]
chk_pct <- tapply(chk$pct.1, list(chk$gene_symbol, chk$cluster), max)
common_g <- intersect(rownames(pct), rownames(chk_pct))
common_c <- intersect(clu, colnames(chk_pct))
cmp <- data.frame(panel = as.vector(pct[common_g, paste0("C", common_c)]),
                  marker = as.vector(chk_pct[common_g, common_c]))
cmp <- cmp[stats::complete.cases(cmp), ]
r <- suppressWarnings(stats::cor(cmp$panel, cmp$marker))
lg(sprintf("self-check vs 04_markers_all pct.1: n=%d  pearson r=%.3f", nrow(cmp), r))
probe <- c(PKHD1L1 = "8", MKI67 = "15", MYH11 = "11", PLIN1 = "16", CCL21 = "13",
           CPA3 = "14", CD3D = "10", BTNL9 = "2", DCN = "1")
for (g in names(probe)) {
  if (!g %in% rownames(pct)) next
  k <- probe[[g]]
  lg(sprintf("  %-8s cluster %-2s pct=%.3f (marker table expected ~%.3f)",
             g, k, pct[g, paste0("C", k)],
             if (g %in% rownames(chk_pct) && k %in% colnames(chk_pct)) chk_pct[g, k] else NA_real_))
}
if (!is.na(r) && r < 0.7) stop("self-check failed: panel/marker pct correlation too low")

write.csv(data.frame(gene = symv, module = mod_of[symv], avg, check.names = FALSE),
          file.path(tab_d, "panel_avg_expr.csv"), row.names = FALSE)
write.csv(data.frame(gene = symv, module = mod_of[symv], pct, check.names = FALSE),
          file.path(tab_d, "panel_pct_expr.csv"), row.names = FALSE)
write.csv(data.frame(gene = symv, module = mod_of[symv], z, check.names = FALSE),
          file.path(tab_d, "panel_zscore.csv"), row.names = FALSE)
write.csv(data.frame(cluster = clu, column = paste0("C", clu),
                     n_cells = as.integer(table(cl_chr)[clu]), check.names = FALSE),
          file.path(tab_d, "cluster_column_map.csv"), row.names = FALSE)
lg("panel tables written")

# ------------------------------------------- module score per cluster --------
mods <- names(panel)
ms <- sapply(clu, function(k) sapply(mods, function(m) {
  g <- intersect(panel[[m]], rownames(z))
  if (!length(g)) return(NA_real_)
  mean(z[g, paste0("C", k)])
}))
colnames(ms) <- paste0("C", clu)
ms <- data.frame(module = rownames(ms), ms, check.names = FALSE)
write.csv(ms, file.path(tab_d, "lineage_module_score.csv"), row.names = FALSE)

top_mod <- sapply(clu, function(k) {
  v <- ms[[paste0("C", k)]]; names(v) <- ms$module
  v <- sort(v, decreasing = TRUE)
  paste0(names(v)[1:3], "(", sprintf("%.2f", v[1:3]), ")", collapse = " | ")
})
write.csv(data.frame(cluster = clu, top3_modules = top_mod, check.names = FALSE),
          file.path(tab_d, "cluster_top_modules.csv"), row.names = FALSE)

# ------------------------------------------------- per-cluster QC ------------
qc <- data.frame(
  cluster = clu,
  n_cells = as.integer(table(cl_chr)[clu]),
  mean_nCount = round(tapply(meta_small$nCount, cl_chr, mean)[clu], 1),
  mean_nFeature = round(tapply(meta_small$nFeature, cl_chr, mean)[clu], 1),
  mean_pct_mt = round(tapply(meta_small$pct_mt, cl_chr, mean)[clu], 4),
  mean_dbl_score = if (all(is.na(meta_small$dbl_score))) NA_real_ else
    round(tapply(meta_small$dbl_score, cl_chr, mean, na.rm = TRUE)[clu], 4),
  pct_disease = round(tapply(meta_small$condition == "Disease", cl_chr, mean)[clu] * 100, 2),
  n_samples = as.integer(tapply(meta_small$sample, cl_chr,
                                function(x) length(unique(x)))[clu]),
  stringsAsFactors = FALSE, row.names = NULL
)
write.csv(qc, file.path(tab_d, "cluster_qc.csv"), row.names = FALSE)

# ------------------------------------- lineage co-detection (doublet proxy) --
det <- cnt[intersect(anchors, rownames(cnt)), , drop = FALSE] > 0
n_anchor <- Matrix::colSums(det)
lg("cells detecting >=2 lineage anchors: ", sprintf("%.2f%%", mean(n_anchor >= 2) * 100))

pairs <- list(c("TTN","PTPRC"), c("COL1A1","CD3D"), c("PECAM1","MYH11"),
              c("LYZ","TNNT2"), c("DCN","CD3D"))
pair_lab <- vapply(pairs, function(p) paste(p, collapse = "_"), character(1))
pair_mat <- sapply(clu, function(k) {
  idx <- which(cl_chr == k)
  vapply(pairs, function(p) {
    if (!all(p %in% rownames(det))) return(NA_real_)
    round(mean(det[p[1], idx] & det[p[2], idx]) * 100, 3)
  }, numeric(1))
})
rownames(pair_mat) <- pair_lab
pair_df <- as.data.frame(t(pair_mat))
colnames(pair_df) <- paste0("codet_", pair_lab)
dl <- data.frame(cluster = clu,
                 pct_ge2_anchor = round(tapply(n_anchor >= 2, cl_chr, mean)[clu] * 100, 2),
                 mean_n_anchor = round(tapply(n_anchor, cl_chr, mean)[clu], 3),
                 pair_df, check.names = FALSE, row.names = NULL)
write.csv(dl, file.path(tab_d, "doublet_like_index.csv"), row.names = FALSE)

saveRDS(list(meta = meta_small, avg = avg, pct = pct, z = z, mod_score = ms,
             qc = qc, dl = dl, panel = panel, cl_col = CL_COL, red = RED,
             check_r = r),
        file.path(obj_d, "annotation_small.rds"))
lg("saved annotation_small.rds (",
   sprintf("%.1f MB", file.size(file.path(obj_d, "annotation_small.rds")) / 1024^2), ")")
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_51.txt"))
lg("=== 51 END ===")
close(lgf)
