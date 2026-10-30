# =============================================================================
# 62_sc_C8_resolution.R
# Resolve C8: is it endothelial at the major-type level, and is the
# "endocardial-like" programme real?
#
# Tests
#   A  core endothelial identity: fraction of C8 nuclei positive for PECAM1 /
#      CDH5 / VWF / EMCN / ERG vs the blood-EC clusters C2 and C7
#   B  endocardial programme: NPR3 / NRG1 / NFATC1 / HAPLN1 detection and the
#      fraction of nuclei co-expressing NPR3 and NRG1
#   C  differential expression C8 vs (C2 + C7) - the direct question "what
#      makes C8 different from the other endothelial clusters"
#   D  alternative identities: lymphatic, mesenchymal/valve, cardiomyocyte
#
# Output: results/annotation/tables/c8_panel_detection.csv
#         results/annotation/tables/c8_cell_level.csv
#         results/annotation/tables/c8_vs_bloodEC_markers.csv
#         results/annotation/62_C8_resolution.txt
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

proj <- "E:/another data file/GEObulk analysis project"
obj  <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
outd <- file.path(proj, "results", "annotation", "tables")
so <- readRDS(obj)
cat("loaded", ncol(so), "nuclei\n")

cl  <- as.character(so$cluster)
gs  <- so[["RNA"]][[]][["gene_symbol"]]

panel <- list(
  EC_core   = c("PECAM1","CDH5","VWF","EMCN","ERG","FLT1","KDR","CD34"),
  EC_endo   = c("NPR3","NRG1","NFATC1","HAPLN1","PLVAP","ESAM","SMOC1","LGR5",
                "TBX3","KLF2","KLF4","SOX17","NOTCH1","DLL4"),
  EC_blood  = c("BTNL9","CA4","ACKR1","ADGRL4","F8","NR2F2","AQP1"),
  LEC       = c("PROX1","FLT4","CCL21","MMRN1","PKHD1L1","STAB2"),
  Mesench   = c("POSTN","COL1A1","COL3A1","FN1","DCN","VIM","SPARC"),
  CM        = c("TTN","MYH7","TNNT2","NPPA","MYL2"))

genes <- unique(unlist(panel))
have  <- genes %in% gs
cat("panel genes missing from object:", paste(genes[!have], collapse = ", "), "\n")
ens <- rownames(so)[match(genes[have], gs)]
names(ens) <- genes[have]

dat <- GetAssayData(so, assay = "RNA", layer = "data")[unname(ens), , drop = FALSE]
rownames(dat) <- names(ens)
cnt <- GetAssayData(so, assay = "RNA", layer = "counts")[unname(ens), , drop = FALSE]
rownames(cnt) <- names(ens)

# ---- A/B/D: per-cluster detection ------------------------------------------
lev <- sort(unique(as.integer(cl)))
det <- sapply(lev, function(k) round(Matrix::rowMeans(cnt[, cl == as.character(k), drop = FALSE] > 0) * 100, 2))
avg <- sapply(lev, function(k) round(Matrix::rowMeans(dat[, cl == as.character(k), drop = FALSE]), 3))
colnames(det) <- paste0("C", lev); colnames(avg) <- paste0("C", lev)

ord <- unlist(panel); ord <- ord[ord %in% rownames(det)]
tab <- data.frame(gene = ord,
                  module = rep(names(panel), lengths(panel))[match(ord, unlist(panel))],
                  det_C2 = det[ord, "C2"], det_C7 = det[ord, "C7"], det_C8 = det[ord, "C8"],
                  avg_C2 = avg[ord, "C2"], avg_C7 = avg[ord, "C7"], avg_C8 = avg[ord, "C8"],
                  stringsAsFactors = FALSE)
write.csv(tab, file.path(outd, "c8_panel_detection.csv"), row.names = FALSE)
cat("\n=== A/B/D panel (detection %, data-layer mean) ===\n")
print(tab, row.names = FALSE)

# ---- cell-level: how many C8 nuclei look endothelial -----------------------
c8  <- cl == "8"; c2 <- cl == "2"; c7 <- cl == "7"
pos <- function(g, idx) as.numeric(cnt[g, idx] > 0)
ec_any <- function(idx) {
  m <- sapply(intersect(c("PECAM1","CDH5","VWF","EMCN","ERG"), rownames(cnt)),
              function(g) as.numeric(cnt[g, idx] > 0))
  rowSums(m) >= 2
}
summ <- data.frame(
  cluster = c("C2", "C7", "C8"),
  n = c(sum(c2), sum(c7), sum(c8)),
  pct_EC_2of5 = round(c(mean(ec_any(c2)), mean(ec_any(c7)), mean(ec_any(c8))) * 100, 2),
  pct_NPR3 = round(c(mean(pos("NPR3", c2)), mean(pos("NPR3", c7)), mean(pos("NPR3", c8))) * 100, 2),
  pct_NRG1 = round(c(mean(pos("NRG1", c2)), mean(pos("NRG1", c7)), mean(pos("NRG1", c8))) * 100, 2),
  pct_NPR3_and_NRG1 = round(c(mean(pos("NPR3", c2) & pos("NRG1", c2)),
                              mean(pos("NPR3", c7) & pos("NRG1", c7)),
                              mean(pos("NPR3", c8) & pos("NRG1", c8))) * 100, 2),
  pct_POSTN = round(c(mean(pos("POSTN", c2)), mean(pos("POSTN", c7)), mean(pos("POSTN", c8))) * 100, 2),
  pct_PECAM1_and_POSTN = round(c(mean(pos("PECAM1", c2) & pos("POSTN", c2)),
                                 mean(pos("PECAM1", c7) & pos("POSTN", c7)),
                                 mean(pos("PECAM1", c8) & pos("POSTN", c8))) * 100, 2),
  stringsAsFactors = FALSE)
write.csv(summ, file.path(outd, "c8_cell_level.csv"), row.names = FALSE)
cat("\n=== cell-level summary ===\n"); print(summ, row.names = FALSE)

# ---- C: C8 vs blood EC -----------------------------------------------------
mk <- FindMarkers(so, ident.1 = colnames(so)[c8], ident.2 = colnames(so)[c2 | c7],
                  only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
mk$gene <- gs[match(rownames(mk), rownames(so))]
mk_out <- mk[order(-mk$avg_log2FC), c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
write.csv(head(mk_out, 60), file.path(outd, "c8_vs_bloodEC_markers.csv"), row.names = FALSE)
cat("\n=== C: top 30 up in C8 vs C2+C7 (by log2FC) ===\n")
print(head(mk_out, 30), row.names = FALSE)

cat("\n=== C: canonical genes in the DE table (if present) ===\n")
probe <- c("NPR3","NRG1","NFATC1","HAPLN1","POSTN","COL3A1","SMOC1","PLVAP","ESAM",
           "ACKR1","BTNL9","CA4","F8","MMRN1","PROX1","TTN","MYH7")
print(mk_out[mk_out$gene %in% probe, ], row.names = FALSE)

sink(file.path(proj, "results", "annotation", "62_C8_resolution.txt"))
cat("C8 resolution evidence\n\n")
print(summ, row.names = FALSE)
cat("\nTop 40 up in C8 vs C2+C7:\n")
print(head(mk_out, 40), row.names = FALSE)
sink()
cat("\nwritten tables + 62_C8_resolution.txt\n")
