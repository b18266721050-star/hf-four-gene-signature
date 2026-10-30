# =============================================================================
# 63_sc_C8_apply.R
# Stage 4b - apply the C8 resolution and refresh everything that depends on it.
#
# Evidence (script 62) showed C8 is endothelial beyond doubt and carries a
# coherent endocardial programme:
#   * 98.9% of C8 nuclei are positive for >=2 of PECAM1/CDH5/VWF/EMCN/ERG
#     (C2 96.2%, C7 97.3% - so C8 is as endothelial as the other EC clusters)
#   * C8 vs C2+C7: NRG1 log2FC 8.90 (46.0% vs 1.9%), SMOC1 8.65 (77.6% vs
#     0.5%), SALL1 8.23, SALL3 8.09, GATA5 7.78, NPR3 5.92, POSTN 5.56
#   * blood-EC markers are absent (BTNL9 1.0%, CA4 0.3%, ACKR1 1.7%)
#   * lymphatic markers are absent (PROX1 6.1%, FLT4 0.6%, CCL21 2.9%)
#
# Decision: major type = Endothelial cell (High). The remaining uncertainty was
# only ever about the subtype, and the subtype is now supported by a multi-gene
# endocardial programme, so it is recorded in `celltype_note` instead of being
# carried as `Uncertain` in the analysis label.
#
# Pipeline position: 60 (initial annotation) -> 63 (this file) -> 61 (figures)
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

proj <- "E:/another data file/GEObulk analysis project"
obj_f <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
adr   <- file.path(proj, "results", "annotation")
tab_d <- file.path(adr, "tables")

so <- readRDS(obj_f)
md <- so@meta.data
cat("loaded", ncol(so), "nuclei; celltypes:", paste(levels(so$celltype), collapse = " | "), "\n")

# ------------------------------------------------------- apply C8 -----------
c8 <- as.character(so$cluster) == "8"
if (all(as.character(so$celltype[c8]) == "Uncertain")) {
  so$celltype[c8] <- "Endothelial cell"
  so$celltype_confidence[c8] <- "High"
  cat("C8 relabelled: Uncertain -> Endothelial cell\n")
} else if (all(as.character(so$celltype[c8]) == "Endothelial cell")) {
  cat("C8 already resolved; refreshing notes and tables only (idempotent rerun)\n")
} else {
  stop("unexpected C8 label: ", paste(unique(so$celltype[c8]), collapse = ","))
}
lv <- c("Cardiomyocyte","Fibroblast","Endothelial cell","Pericyte",
        "Smooth muscle cell","Macrophage","Monocyte/DC","T/NK cell",
        "Lymphatic endothelial cell","Mast cell","Neural/glial cell",
        "Proliferating cell")
so$celltype <- factor(as.character(so$celltype), levels = lv)
so$celltype_confidence[c8] <- "High"
cat("celltypes now:", paste(levels(so$celltype), collapse = " | "), "\n")
print(table(so$celltype))

# per-cluster note (informational only - never used as a grouping variable)
note <- c("0" = "CM, atrial-like (MYH6/CHRM2/ALPK2)",
          "1" = "FB, quiescent (PDGFRA/ABCA9, MFAP5 low)",
          "2" = "EC, capillary/arterial (BTNL9/CA4, ACKR1-)",
          "3" = "Pericyte (RGS5/HIGD1B/EGFLAM, MYH11-)",
          "4" = "FB, activated (MFAP5/SERPINE2)",
          "5" = "Macrophage, resident (MRC1/CD163/F13A1/LYVE1)",
          "6" = "CM, ventricular (MYBPC3/RBM20/MYL2)",
          "7" = "EC, venous (ACKR1/NR2F2)",
          "8" = "EC, endocardial (NRG1/SMOC1/SALL1/GATA5/NPR3)",
          "9" = "Monocyte/DC, non-resident (LYZ/FCN1/ITGAX)",
          "10" = "T/NK (THEMIS/BCL11B/ITK/SKAP1)",
          "11" = "SMC (MYH11/MYOCD/LMOD1)",
          "12" = "Neural/glial (SOX10/PLP1/PTPRZ1)",
          "13" = "LEC (PROX1/CCL21/FLT4/MMRN1)",
          "14" = "Mast (KIT/CPA3/TPSB2/MS4A2)",
          "15" = "Proliferating, myeloid-biased (MKI67 + resident macrophage genes)")
so$celltype_note <- unname(note[as.character(so$cluster)])
stopifnot(!any(is.na(so$celltype_note)))

# --------------------------------------------- refresh composition tables ---
md <- so@meta.data            # refresh: the earlier copy still holds `Uncertain`
grp <- setNames(read.csv(file.path(proj, "clean_data", "GSE292067_sample_manifest.csv"),
                         stringsAsFactors = FALSE)$group,
                read.csv(file.path(proj, "clean_data", "GSE292067_sample_manifest.csv"),
                         stringsAsFactors = FALSE)$sample)

cnt_by_sample <- as.data.frame.matrix(table(md$sample, so$celltype))
cnt_by_sample <- data.frame(sample = rownames(cnt_by_sample),
                            group = unname(grp[rownames(cnt_by_sample)]),
                            cnt_by_sample, check.names = FALSE, row.names = NULL)
prop_by_sample <- cnt_by_sample
prop_by_sample[, lv] <- round(sweep(as.matrix(cnt_by_sample[, lv]), 1,
                                    rowSums(as.matrix(cnt_by_sample[, lv])), "/") * 100, 3)

cnt_by_group <- as.data.frame.matrix(table(so$group, so$celltype))
cnt_by_group <- data.frame(group = rownames(cnt_by_group), cnt_by_group,
                           check.names = FALSE, row.names = NULL)
prop_by_group <- cnt_by_group
prop_by_group[, lv] <- round(sweep(as.matrix(cnt_by_group[, lv]), 1,
                                   rowSums(as.matrix(cnt_by_group[, lv])), "/") * 100, 3)

overall <- data.frame(celltype = lv,
                      n_cells = as.integer(table(so$celltype)[lv]),
                      stringsAsFactors = FALSE)
overall$pct_of_all <- round(overall$n_cells / ncol(so) * 100, 3)
overall$confidence <- c("High","High","High","High","High","High","Moderate",
                        "High","High","High","High","Moderate")
overall$n_samples <- sapply(lv, function(k) length(unique(so$sample[so$celltype == k])))
overall$clusters <- c("C0, C6","C1, C4","C2, C7, C8","C3","C11","C5","C9","C10",
                      "C13","C14","C12","C15")
overall$note <- sapply(lv, function(k) paste(unique(so$celltype_note[so$celltype == k]),
                                             collapse = "; "))

write.csv(cnt_by_sample,  file.path(tab_d, "celltype_counts_by_sample.csv"), row.names = FALSE)
write.csv(prop_by_sample, file.path(tab_d, "celltype_proportions_by_sample.csv"), row.names = FALSE)
write.csv(cnt_by_group,   file.path(tab_d, "celltype_counts_by_group.csv"), row.names = FALSE)
write.csv(prop_by_group,  file.path(tab_d, "celltype_proportions_by_group.csv"), row.names = FALSE)
write.csv(overall,        file.path(tab_d, "celltype_overall.csv"), row.names = FALSE)

samp_lev <- unique(so$sample[order(match(so$group, c("NF","NICM","DoxCM")), so$sample)])
long_s <- do.call(rbind, lapply(lv, function(k) {
  data.frame(sample = cnt_by_sample$sample, group = cnt_by_sample$group, celltype = k,
             n = cnt_by_sample[[k]], pct = prop_by_sample[[k]], stringsAsFactors = FALSE)
}))
write.csv(long_s, file.path(tab_d, "celltype_composition_long.csv"), row.names = FALSE)

print(overall, row.names = FALSE)
print(prop_by_group, row.names = FALSE)

# --------------------------------------------- refresh mapping + decisions --
map <- read.csv(file.path(tab_d, "final_celltype_mapping.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
map$celltype[map$cluster == "8"] <- "Endothelial cell"
map$confidence[map$cluster == "8"] <- "High"
map$reason[map$cluster == "8"] <-
  "resolved in script 62: 98.9% EC-marker positive; NRG1 +8.90, SMOC1 +8.65, SALL1 +8.23, GATA5 +7.78, NPR3 +5.92 vs C2+C7; blood-EC and lymphatic markers absent"
write.csv(map, file.path(tab_d, "final_celltype_mapping.csv"), row.names = FALSE)

dec <- read.csv(file.path(tab_d, "annotation_decisions.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
dec$decision[dec$cluster == "8"] <- "include"
dec$decided_on[dec$cluster == "8"] <- "2026-10-29"
dec$note[dec$cluster == "8"] <-
  "Endothelial cell - resolved by differential expression vs C2+C7 (NRG1/SMOC1/SALL1/GATA5/NPR3); subtype recorded as endocardial in celltype_note"
write.csv(dec, file.path(tab_d, "annotation_decisions.csv"), row.names = FALSE)

# ------------------------------------------------------------- save ---------
saveRDS(so, obj_f)
cat("saved", obj_f, "|", round(file.size(obj_f) / 1024^2, 1), "MB\n")
cat("celltype_note levels:\n"); print(table(so$celltype_note))
