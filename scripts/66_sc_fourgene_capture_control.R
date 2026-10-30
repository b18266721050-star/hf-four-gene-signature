# =============================================================================
# 66_sc_fourgene_capture_control.R
# Stage 5 (A3) - is the disease-related drop gene-specific or a global
#                capture/library shift?
#
# The project already has a documented "everything goes down" problem at the
# bulk level (see 72_global_shift_check.R). Before reporting that AQP3/GGT5
# fall in disease, the snRNA-seq side needs the same control: a panel of
# stably expressed genes measured in exactly the same nuclei.
#
# Input : results/Singlecel/objects/03_merged/05_final_annotated.rds
# Output: results/Singlecel/fourgene/tables/66_capture_control_*.csv
#         results/Singlecel/fourgene/figures/66_capture_control.{png,pdf}
#         results/Singlecel/fourgene/66_capture_control.log
#
# For every donor x cell type:
#   det_t   = % nuclei with > 0 counts for the target gene
#   det_hk  = mean % nuclei with > 0 counts over the housekeeping panel
#   ratio   = det_t / det_hk          <- capture-corrected detection
# The ratio is what gets tested across groups; a uniform loss of RNA in one
# group cancels out, a genuine loss of the gene does not.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork)
})

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
obj_f <- file.path(proj, "results", "Singlecel", "objects", "03_merged", "05_final_annotated.rds")
out_d <- file.path(proj, "results", "Singlecel", "fourgene")
tab_d <- file.path(out_d, "tables"); fig_d <- file.path(out_d, "figures")
if (!dir.exists(fig_d)) dir.create(fig_d, recursive = TRUE)

lgf <- file.path(out_d, "66_capture_control.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

TARGET <- c("AQP3", "GGT5", "ISLR", "LAD1")
HK <- c("MALAT1", "GAPDH", "ACTB", "B2M", "TPT1", "EEF1A1",
        "RPLP0", "PPIA", "UBC", "FTL")
GRP <- c("NF", "NICM", "DoxCM")
GCOL <- c(NF = "#2C7FB8", NICM = "#E9A03B", DoxCM = "#D1495B")
MIN_CELLS <- 50
MIN_DONOR <- 4

lg("=== 66 capture control | seed ", SEED, " ===")
t0 <- Sys.time()
so <- readRDS(obj_f)
lg(sprintf("loaded %d genes x %d nuclei in %.1f min", nrow(so), ncol(so),
           as.numeric(difftime(Sys.time(), t0, units = "mins"))))

gs <- so[["RNA"]][[]][["gene_symbol"]]
want <- c(TARGET, HK)
idx <- match(want, gs)
names(idx) <- want
lg("missing symbols: ", paste(want[is.na(idx)], collapse = ","))
want <- want[!is.na(idx)]
idx  <- idx[!is.na(idx)]
cnt0 <- GetAssayData(so, assay = "RNA", layer = "counts")[rownames(so)[idx], , drop = FALSE]
sel <- setNames(rownames(so)[idx], want)
dup <- names(sel)[duplicated(names(sel))]
if (length(dup)) {
  drop <- unlist(lapply(dup, function(g) {
    rows <- which(names(sel) == g)
    rows[-which.max(Matrix::rowSums(cnt0[rows, , drop = FALSE]))]
  }))
  sel <- sel[-drop]
  lg("dropped duplicated-symbol features: ", paste(dup, collapse = ","))
}
cnt <- cnt0[unname(sel), , drop = FALSE]
rownames(cnt) <- names(sel)
md <- so@meta.data
md$group <- factor(as.character(md$group), levels = GRP)
md$celltype <- factor(as.character(md$celltype))
md$sample <- factor(as.character(md$sample))
CT <- levels(md$celltype)
rm(so); invisible(gc())

# ------------------------------------------------ donor x cell type rates ---
key <- factor(paste(md$sample, md$celltype, sep = "|"))
D <- sparse.model.matrix(~0 + key); colnames(D) <- levels(key)
n   <- setNames(as.numeric(colSums(D)), colnames(D))
det <- as.matrix(((cnt > 0) * 1) %*% D)
pct <- sweep(det, 2, n, "/") * 100
tot <- setNames(as.numeric(crossprod(D, md$nCount_RNA)), colnames(D))
cpm <- sweep(as.matrix(cnt %*% D), 2, tot, "/") * 1e6

cell <- data.frame(key = colnames(D), stringsAsFactors = FALSE)
cell$sample <- sub("\\|.*$", "", cell$key)
cell$celltype <- sub("^.*\\|", "", cell$key)
cell$group <- as.character(md$group[match(cell$sample, as.character(md$sample))])
cell$n_nuclei <- as.numeric(n[cell$key])
hkv <- intersect(HK, rownames(pct))
cell$det_hk <- colMeans(pct[hkv, , drop = FALSE])
cell$cpm_hk <- exp(colMeans(log1p(cpm[hkv, , drop = FALSE])))
for (g in TARGET) {
  cell[[paste0("det_", g)]] <- as.numeric(pct[g, cell$key])
  cell[[paste0("ratio_", g)]] <- as.numeric(pct[g, cell$key]) / cell$det_hk
}
cell$usable <- cell$n_nuclei >= MIN_CELLS
cell$group <- factor(cell$group, levels = GRP)
write.csv(cell, file.path(tab_d, "66_capture_control_sample_x_celltype.csv"), row.names = FALSE)
lg("housekeeping panel used: ", paste(hkv, collapse = ","))

# per-donor whole tissue (all nuclei)
key2 <- factor(md$sample)
D2 <- sparse.model.matrix(~0 + key2); colnames(D2) <- levels(key2)
n2 <- setNames(as.numeric(colSums(D2)), colnames(D2))
pct2 <- sweep(as.matrix(((cnt > 0) * 1) %*% D2), 2, n2, "/") * 100
DNA <- data.frame(sample = colnames(D2), stringsAsFactors = FALSE)
DNA$group <- factor(as.character(md$group[match(DNA$sample, as.character(md$sample))]), levels = GRP)
DNA$n_nuclei <- as.numeric(n2[DNA$sample])
DNA$det_hk <- colMeans(pct2[hkv, , drop = FALSE])
for (g in TARGET) {
  DNA[[paste0("det_", g)]] <- as.numeric(pct2[g, DNA$sample])
  DNA[[paste0("ratio_", g)]] <- DNA[[paste0("det_", g)]] / DNA$det_hk
}
write.csv(DNA, file.path(tab_d, "66_capture_control_donor.csv"), row.names = FALSE)

# per-donor QC (median per nucleus)
qc <- aggregate(cbind(nFeature_RNA, nCount_RNA, percent.mt) ~ sample + group,
                data = md, FUN = median)
write.csv(qc, file.path(tab_d, "66_donor_qc_medians.csv"), row.names = FALSE)
fmt <- function(v, d = 0) paste(sprintf(paste0("%s=%." , d, "f"), names(v), as.numeric(v)), collapse = " ")
lg("median nFeature by group: ", fmt(tapply(qc$nFeature_RNA, qc$group, median)))
lg("median nCount by group:   ", fmt(tapply(qc$nCount_RNA, qc$group, median)))
lg("median det_hk by group:   ", fmt(tapply(DNA$det_hk, DNA$group, median), 2))

# ---------------------------------------------------------------- tests -----
wt <- function(x, g) {
  a <- x[g == "NF"]; b <- x[g != "NF"]
  if (length(a) < MIN_DONOR || length(b) < MIN_DONOR) return(NA_real_)
  tryCatch(wilcox.test(b, a, exact = FALSE)$p.value, error = function(e) NA_real_)
}
rows <- do.call(rbind, lapply(TARGET, function(g) {
  r <- data.frame(gene = g,
                  det_NF = median(DNA[[paste0("det_", g)]][DNA$group == "NF"]),
                  det_Disease = median(DNA[[paste0("det_", g)]][DNA$group != "NF"]),
                  ratio_NF = median(DNA[[paste0("ratio_", g)]][DNA$group == "NF"]),
                  ratio_Disease = median(DNA[[paste0("ratio_", g)]][DNA$group != "NF"]),
                  p_raw_det = wt(DNA[[paste0("det_", g)]], DNA$group),
                  p_ratio_det = wt(DNA[[paste0("ratio_", g)]], DNA$group),
                  stringsAsFactors = FALSE)
  r$fold_ratio <- r$ratio_Disease / r$ratio_NF
  r
}))
rows$fdr_ratio <- p.adjust(rows$p_ratio_det, "BH")
write.csv(rows, file.path(tab_d, "66_capture_control_stats.csv"), row.names = FALSE)
lg("--- capture-corrected detection, donor level ---")
for (i in seq_len(nrow(rows)))
  lg(sprintf("  %-5s ratio %.4f -> %.4f (x%.2f) p=%.3g FDR=%.3g | raw det %.2f%% -> %.2f%% p=%.3g",
             rows$gene[i], rows$ratio_NF[i], rows$ratio_Disease[i], rows$fold_ratio[i],
             rows$p_ratio_det[i], rows$fdr_ratio[i], rows$det_NF[i], rows$det_Disease[i],
             rows$p_raw_det[i]))

# housekeeping panel stability (is the panel itself flat across groups?)
hk_long <- do.call(rbind, lapply(hkv, function(g) data.frame(
  gene = g, sample = DNA$sample, group = DNA$group, det = as.numeric(pct2[g, DNA$sample]))))
hk_stat <- do.call(rbind, lapply(hkv, function(g) {
  d <- hk_long[hk_long$gene == g, ]
  data.frame(gene = g, det_NF = median(d$det[d$group == "NF"]),
             det_Disease = median(d$det[d$group != "NF"]), p = wt(d$det, d$group))
}))
write.csv(hk_stat, file.path(tab_d, "66_housekeeping_stability.csv"), row.names = FALSE)
lg("HK genes with p<0.05 across groups: ", sum(hk_stat$p < 0.05, na.rm = TRUE), " / ", nrow(hk_stat))

# ---------------------------------------------------------------- figure ----
qc_long <- rbind(
  data.frame(metric = "median nFeature per nucleus", qc[, c("sample", "group")], value = qc$nFeature_RNA),
  data.frame(metric = "median nCount per nucleus", qc[, c("sample", "group")], value = qc$nCount_RNA),
  data.frame(metric = "housekeeping detection (%)", DNA[, c("sample", "group")], value = DNA$det_hk),
  data.frame(metric = "median mitochondrial %", qc[, c("sample", "group")], value = qc$percent.mt))
qc_long$metric <- factor(qc_long$metric, levels = unique(qc_long$metric))

p_qc <- ggplot(qc_long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, linewidth = 0.3) +
  geom_point(aes(colour = group), size = 1.2, stroke = 0.15,
             position = position_jitter(width = 0.12, seed = SEED)) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  labs(title = "Per-donor capture / RNA-quality metrics",
       subtitle = "if one group has systematically less RNA per nucleus, every gene looks down there",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(size = 7.4),
        strip.background = element_rect(fill = "grey94", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ratio_long <- do.call(rbind, lapply(TARGET, function(g) data.frame(
  gene = g, sample = DNA$sample, group = DNA$group,
  value = DNA[[paste0("ratio_", g)]])))
ratio_long$gene <- factor(ratio_long$gene, levels = TARGET)
lab <- data.frame(gene = factor(rows$gene, levels = TARGET),
                  value = Inf, txt = sprintf("FDR %s", ifelse(rows$fdr_ratio < 0.001,
                    formatC(rows$fdr_ratio, format = "e", digits = 1),
                    sprintf("%.3f", rows$fdr_ratio))))
p_ratio <- ggplot(ratio_long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, linewidth = 0.3) +
  geom_point(aes(colour = group), size = 1.3, stroke = 0.15,
             position = position_jitter(width = 0.12, seed = SEED)) +
  geom_text(data = lab, aes(x = 2, y = value, label = txt), vjust = 1.5,
            size = 2.7, inherit.aes = FALSE) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GCOL) + scale_colour_manual(values = GCOL) +
  labs(title = "Capture-corrected detection: target gene / housekeeping panel",
       subtitle = sprintf("housekeeping panel: %s\na flat ratio means the gene itself changed; a falling baseline would mean the library changed",
                          paste(hkv, collapse = ", ")),
       x = NULL, y = "detection ratio") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 8),
        strip.background = element_rect(fill = "grey94", colour = NA),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

p66 <- p_qc / p_ratio + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(title = "Is the disease-related loss gene-specific?",
                  subtitle = "donor-level; one point = one donor (25 donors)")
ggsave(file.path(fig_d, "66_capture_control.png"), p66, width = 10.6, height = 6.6, dpi = 300)
ggsave(file.path(fig_d, "66_capture_control.pdf"), p66, width = 10.6, height = 6.6)

writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_66.txt"))
lg("=== done ===")
