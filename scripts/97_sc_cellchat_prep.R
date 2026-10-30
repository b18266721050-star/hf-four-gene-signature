# =============================================================================
# 97_sc_cellchat_prep.R
# Stage 6 (A4) - build the ligand-receptor input for cell-cell communication.
#
# Input : results/Singlecel/objects/03_merged/05_final_annotated.rds
#         raw_data/CellChatDB.human.rda   (official CellChatDB, jinworks/CellChat)
# Output: results/Singlecel/cellchat/objects/cc_input.rds
#         results/Singlecel/cellchat/tables/97_{db_summary,gene_coverage,meta_check}.csv
#         results/Singlecel/cellchat/logs/97_prep.log
#
# Why a slim object again: the whole-transcriptome object is 1.6 GB and the
# communication analysis only needs the 1.4k genes that appear in the database.
# Reading it once here keeps every later step under a minute.
#
# Deviations from a stock CellChat run (documented in the report):
#   * CellChat itself is not installable here (its src/ needs a C++ toolchain,
#     no Rtools on this machine), so the mass-action model is re-implemented
#     on top of the official CellChatDB interaction table. The formulas follow
#     CellChat v1: P = Hill(meanL * meanR, Kh = 0.5, n = 1).
#   * complexes are collapsed with min() over subunits (CellChat convention:
#     every subunit must be present).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

SEED <- 20251029
set.seed(SEED)

proj   <- "E:/another data file/GEObulk analysis project"
obj_f  <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
db_f   <- file.path(proj, "raw_data", "CellChatDB.human.rda")
out_d  <- file.path(proj, "results", "Singlecel", "cellchat")
tab_d  <- file.path(out_d, "tables"); obj_d <- file.path(out_d, "objects")
log_d  <- file.path(out_d, "logs")
for (d in c(tab_d, obj_d, log_d)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

lgf <- file.path(log_d, "97_prep.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

CT_ORDER <- c("Cardiomyocyte", "Fibroblast", "Endothelial cell",
              "Lymphatic endothelial cell", "Pericyte", "Smooth muscle cell",
              "Macrophage", "Monocyte/DC", "T/NK cell", "Mast cell",
              "Neural/glial cell", "Proliferating cell")
GRP <- c("NF", "NICM", "DoxCM")

lg("=== 97 cellchat prep | seed ", SEED, " ===")

# --------------------------------------------------------------- database ----
e <- new.env(); load(db_f, envir = e)
db <- get(ls(e)[1], envir = e)
inter <- db$interaction
cmplx <- db$complex
lg("CellChatDB interactions: ", nrow(inter), " | pathways: ", length(unique(inter$pathway_name)),
   " | annotation: ", paste(names(table(inter$annotation)), as.integer(table(inter$annotation)), sep = "=", collapse = " "))

# complex -> subunits (min() later on)
cmplx_list <- lapply(seq_len(nrow(cmplx)), function(i) {
  v <- as.character(unlist(cmplx[i, ])); v <- v[!is.na(v) & v != ""]; v
})
names(cmplx_list) <- rownames(cmplx)

split_mol <- function(x) {
  # "TGFbR1_R2" -> c("TGFBR1","TGFBR2") via the complex table; "TGFB1" -> itself
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  out <- unlist(lapply(parts, function(p) if (p %in% names(cmplx_list)) cmplx_list[[p]] else p))
  unique(toupper(out))
}
lig_list <- lapply(inter$ligand,   split_mol)
rec_list <- lapply(inter$receptor, split_mol)

# keep one row per unique (ligand-set, receptor-set, pathway)
key <- vapply(seq_len(nrow(inter)), function(i)
  paste(inter$pathway_name[i], paste(lig_list[[i]], collapse = "+"),
        paste(rec_list[[i]], collapse = "+"), sep = "|"), character(1))
keep <- !duplicated(key)
lg("unique LR records after complex collapse: ", sum(keep), " (from ", nrow(inter), ")")

pairs <- data.frame(
  interaction_name = inter$interaction_name_2,
  pathway          = inter$pathway_name,
  annotation       = inter$annotation,
  evidence         = inter$evidence,
  ligand_str       = inter$ligand,
  receptor_str     = inter$receptor,
  n_lig_subunit    = vapply(lig_list, length, integer(1)),
  n_rec_subunit    = vapply(rec_list, length, integer(1)),
  stringsAsFactors = FALSE)[keep, ]
pairs$ligand   <- I(lig_list[keep])
pairs$receptor <- I(rec_list[keep])
rownames(pairs) <- NULL

genes <- sort(unique(unlist(c(pairs$ligand, pairs$receptor))))
lg("unique signalling genes: ", length(genes))

# neurotransmitter-only pathways are irrelevant to ventricular tissue and
# dominate the pair count (Glutamate = 401 pairs); they are kept in the object
# but flagged so downstream steps can exclude them if needed.
NEURO <- c("Glutamate", "GABA-A", "GABA-B", "5-HT", "Ach", "Dopamine", "Serotonin",
           "Noradrenaline", "Adrenaline", "OPIOID", "Adenosine", "Prostaglandin", "NPY")
pairs$neuro <- pairs$pathway %in% NEURO
lg("neurotransmitter pairs flagged: ", sum(pairs$neuro))

# ------------------------------------------------------------------ object ---
t0 <- Sys.time()
so <- readRDS(obj_f)
lg(sprintf("loaded object %d genes x %d nuclei in %.1f min",
           nrow(so), ncol(so), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

gs <- so[["RNA"]][[]][["gene_symbol"]]
cnt_all <- GetAssayData(so, assay = "RNA", layer = "counts")
md <- so@meta.data

# symbol -> row(s); duplicated symbols are summed (pseudogene/duplicate probes)
present <- genes[genes %in% gs]
missing <- setdiff(genes, gs)
lg("gene coverage: ", length(present), "/", length(genes), " present (",
   round(length(present) / length(genes) * 100, 1), "%)")
lg("missing (first 40): ", paste(head(missing, 40), collapse = ","))

idx_by_gene <- split(seq_len(nrow(so)), gs)[present]
n_sub <- vapply(idx_by_gene, length, integer(1))
rows  <- unlist(idx_by_gene, use.names = FALSE)
cnt <- cnt_all[rows, , drop = FALSE]
if (any(n_sub > 1)) {
  # several Ensembl features per symbol -> sum them (keeps the symbol as unit)
  fac <- factor(rep(present, n_sub), levels = present)
  M <- Matrix::sparseMatrix(i = as.integer(fac), j = seq_along(fac), x = 1,
                            dims = c(length(present), length(fac)))
  cnt <- as(M %*% cnt, "CsparseMatrix")
}
rownames(cnt) <- present
rm(cnt_all); gc(verbose = FALSE)
lg("LR count matrix: ", nrow(cnt), " genes x ", ncol(cnt), " nuclei | nnz ", length(cnt@x),
   " | sparsity ", round((1 - length(cnt@x) / (nrow(cnt) * ncol(cnt))) * 100, 2), "%")

meta <- data.frame(
  cell   = colnames(so),
  sample = as.character(md$sample),
  donor  = as.character(md$donor),
  group  = factor(as.character(md$group), levels = GRP),
  condition = as.character(md$condition),
  celltype  = factor(as.character(md$celltype), levels = CT_ORDER),
  celltype_confidence = as.character(md$celltype_confidence),
  nCount_RNA = md$nCount_RNA,
  nFeature_RNA = md$nFeature_RNA,
  percent_mt = md$percent.mt,
  stringsAsFactors = FALSE)
stopifnot(identical(meta$cell, colnames(cnt)), !any(is.na(meta$celltype)))

# ----------------------------------------------------------------- checks ----
tbl <- table(meta$celltype, meta$group)
write.csv(as.data.frame.matrix(tbl), file.path(tab_d, "97_meta_check.csv"))
lg("nuclei per celltype x group:")
print(tbl)

depth <- aggregate(nCount_RNA ~ celltype + group, data = meta, FUN = median)
write.csv(depth, file.path(tab_d, "97_depth_by_celltype_group.csv"), row.names = FALSE)
lg("median nCount per celltype x group written (depth confound audit)")

db_sum <- data.frame(
  metric = c("interactions_raw", "interactions_unique", "pathways", "signalling_genes",
             "genes_present", "genes_missing", "nuclei", "celltypes",
             "secreted", "ecm_receptor", "cell_contact", "non_protein"),
  value = c(nrow(inter), nrow(pairs), length(unique(pairs$pathway)), length(genes),
            length(present), length(missing), ncol(cnt), length(unique(meta$celltype)),
            sum(pairs$annotation == "Secreted Signaling"),
            sum(pairs$annotation == "ECM-Receptor"),
            sum(pairs$annotation == "Cell-Cell Contact"),
            sum(pairs$annotation == "Non-protein Signaling")))
write.csv(db_sum, file.path(tab_d, "97_db_summary.csv"), row.names = FALSE)

cov <- data.frame(gene = c("AQP3", "GGT5", "ISLR", "LAD1",
                           "STAT3", "CD163", "FPR1", "IL1R2", "ANPEP", "S100A9", "VSIG4"),
                  in_db = NA, n_records = NA, records = NA, stringsAsFactors = FALSE)
for (i in seq_len(nrow(cov))) {
  g <- cov$gene[i]
  hit <- which(vapply(pairs$ligand, function(x) g %in% x, logical(1)) |
               vapply(pairs$receptor, function(x) g %in% x, logical(1)))
  cov$in_db[i] <- length(hit) > 0
  cov$n_records[i] <- length(hit)
  cov$records[i] <- paste(head(unique(pairs$interaction_name[hit]), 6), collapse = "; ")
}
write.csv(cov, file.path(tab_d, "97_gene_coverage.csv"), row.names = FALSE)
lg("four model genes / PPI hubs in CellChatDB:")
print(cov[, 1:3])

# ------------------------------------------------------------------ save -----
cc <- list(counts = cnt, meta = meta, pairs = pairs, pathways = sort(unique(pairs$pathway)),
           neuro = NEURO, ct_order = CT_ORDER, groups = GRP,
           params = list(Kh = 0.5, hill = 1, complex_rule = "min", trim = "mean",
                         normalization = "log1p(count/nCount*1e4)",
                         db = basename(db_f), db_version = unique(inter$version),
                         seed = SEED, built = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
saveRDS(cc, file.path(obj_d, "cc_input.rds"))
lg("saved cc_input.rds: ", round(file.size(file.path(obj_d, "cc_input.rds")) / 1e6, 1), " MB")

writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_97.txt"))
lg("=== done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min ===")
