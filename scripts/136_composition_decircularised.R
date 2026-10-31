###############################################################################
## 136_composition_decircularised.R
##
## R1 of the pre-submission audit: remove the circularity of the
## composition-adjusted analysis.
##
## Problem
##   92_composition_adjust.R adjusted the differential expression of the 70
##   candidate genes for xCell compartment scores. Seven of those 70 genes
##   (CD163, VSIG4, ANPEP, S100A9, FPR1, IL1R2, STAT3) are canonical myeloid
##   markers that also contribute to the xCell macrophage signature. The
##   attenuation produced by that model is therefore partly mechanical: the
##   covariate is built from the outcome.
##
## What this script does
##   1. Quantifies the overlap between each xCell compartment signature and the
##      set of 70 candidate genes (counts, Jaccard, gene names).
##   2. Rebuilds the four compartment scores by ssGSEA, twice: once with the
##      complete signature (method control) and once with all 70 candidate
##      genes removed from the signature (de-circularised).
##   3. Fits four limma models on the same samples:
##        M0  group + sex + age                      (no composition)
##        M1  M0 + z(xCell scores)                   (published analysis)
##        M2  M0 + z(ssGSEA, complete signature)     (method control)
##        M3  M0 + z(ssGSEA, 70 genes removed)       (de-circularised)
##      and reports how much of the 70-gene effect survives each.
##
## Interpretation
##   If attenuation persists in M3, the compositional contribution is real.
##   If it collapses in M3 relative to M2, the published attenuation was
##   carried by the shared marker genes and must be reported as an upper bound.
##
## Inputs
##   clean_data/GSE57338_clean_expression_matrix.csv
##   results/Composition/01_tables/T04_composition_merged.csv
##   results/HubGene/GSE57338_candidate_genes.csv
##
## Outputs (results/Composition/)
##   01_tables/T21_decircularised_models.csv
##   01_tables/T22_decircularised_key_genes.csv
##   01_tables/T23_signature_overlap.csv
##   02_figures/F10_attenuation_decirculated.{pdf,png}
##   03_logs/136_run_log.txt , sessionInfo_136.txt
##
## Run:  Rscript scripts/136_composition_decircularised.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GSVA)
  library(GSEABase)
  library(xCell)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20261030)

proj <- "E:/another data file/GEObulk analysis project"
out_root <- file.path(proj, "results", "Composition")
d_tab <- file.path(out_root, "01_tables")
d_fig <- file.path(out_root, "02_figures")
d_log <- file.path(out_root, "03_logs")
for (d in c(d_tab, d_fig, d_log)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_log <- character(0)
log_msg <- function(...) {
  txt <- paste0(...); run_log <<- c(run_log, txt); cat(txt, "\n", sep = ""); flush.console()
}
hr <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")
write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

LFC_CUT <- 0.585
PADJ_CUT <- 0.05
GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")
HUBS <- c("CD163", "VSIG4", "S100A9", "FPR1", "STAT3", "ANPEP", "IL1R2")

## ===========================================================================
## 1. Data
## ===========================================================================
hr("=")
log_msg("R1  Composition adjustment without circularity (GSE57338)")
hr("=")

comp <- data.table::fread(file.path(d_tab, "T04_composition_merged.csv"))
expr_dt <- data.table::fread(file.path(proj, "clean_data", "GSE57338_clean_expression_matrix.csv"),
                             data.table = FALSE, showProgress = FALSE)
expr <- as.matrix(expr_dt[, -1, drop = FALSE]); rownames(expr) <- expr_dt[[1]]
rm(expr_dt); invisible(gc())
log_msg("[1] expression: ", nrow(expr), " genes x ", ncol(expr), " samples")

stopifnot(setequal(comp$sample, colnames(expr)))
comp <- comp[match(colnames(expr), sample)]
comp[, group2 := factor(ifelse(group == "NF", "Control", "Disease"),
                        levels = c("Control", "Disease"))]
comp[, sex := factor(sex, levels = c("Female", "Male"))]
if (anyNA(comp$age_num)) comp[is.na(age_num), age_num := median(comp$age_num, na.rm = TRUE)]

cand <- data.table::fread(file.path(proj, "results", "HubGene", "GSE57338_candidate_genes.csv"))
cand70 <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]
log_msg("    candidate genes: ", length(cand70))

## ===========================================================================
## 2. Signature overlap (the circularity, quantified)
## ===========================================================================
hr()
log_msg("[2] Overlap between xCell signatures and the 70 candidate genes")

gsc <- GSEABase::geneIds(xCell.data$signatures)
ct <- sub("%.*$", "", names(gsc))

comp_sets <- list(
  Macrophages  = grep("^Macrophages", ct),
  Fibroblasts  = grep("^Fibroblasts", ct),
  Endothelial  = grep("Endothelial cells$", ct),
  Myocytes     = grep("^Myocytes", ct)
)
comp_sets <- comp_sets[lengths(comp_sets) > 0]

overlap_rows <- list()
gene_sets_full <- list()
for (nm in names(comp_sets)) {
  idx <- comp_sets[[nm]]
  g <- unique(unlist(gsc[idx]))
  g <- intersect(g, rownames(expr))
  ov <- intersect(g, cand70)
  gene_sets_full[[nm]] <- g
  overlap_rows[[nm]] <- data.table(
    compartment = nm,
    n_sets_in_xcell = length(idx),
    n_genes_in_signature = length(g),
    n_overlap_with_70 = length(ov),
    jaccard = round(length(ov) / length(union(g, cand70)), 4),
    pct_of_signature = round(100 * length(ov) / length(g), 2),
    overlap_genes = paste(ov, collapse = ";")
  )
  log_msg(sprintf("    %-12s signature=%4d genes  overlap=%2d  jaccard=%.4f  [%s]",
                  nm, length(g), length(ov),
                  length(ov) / length(union(g, cand70)), paste(ov, collapse = ",")))
}
ov_tab <- rbindlist(overlap_rows)
write_csv(ov_tab, file.path(d_tab, "T23_signature_overlap.csv"))

## de-circularised gene sets: drop every one of the 70 candidate genes
gene_sets_dec <- lapply(gene_sets_full, function(g) setdiff(g, cand70))
for (nm in names(gene_sets_dec)) {
  log_msg(sprintf("    %-12s after removal: %d genes (-%d)", nm,
                  length(gene_sets_dec[[nm]]),
                  length(gene_sets_full[[nm]]) - length(gene_sets_dec[[nm]])))
}

## ===========================================================================
## 3. Rebuild compartment scores by ssGSEA (complete and de-circularised)
## ===========================================================================
hr()
log_msg("[3] ssGSEA compartment scores (twice)")

run_ssgsea <- function(gs_list, tag) {
  par <- ssgseaParam(expr, gs_list, normalize = TRUE, alpha = 0.25)
  s <- gsva(par, verbose = FALSE)
  s <- t(s)                        # samples x compartments
  colnames(s) <- paste0("ss_", tag, "_", colnames(s))
  log_msg("    ", tag, ": ", ncol(s), " scores, range ",
          paste(round(range(s), 2), collapse = " to "))
  s
}
ss_full <- run_ssgsea(gene_sets_full, "full")
ss_dec  <- run_ssgsea(gene_sets_dec,  "dec")

comp <- cbind(comp, as.data.table(ss_full), as.data.table(ss_dec))
log_msg("    correlation xCell vs ssGSEA(full), Macrophages: r = ",
        round(cor(comp$xcell_Macrophages, comp$ss_full_Macrophages), 3))
log_msg("    correlation ssGSEA(full) vs ssGSEA(dec), Macrophages: r = ",
        round(cor(comp$ss_full_Macrophages, comp$ss_dec_Macrophages), 3))

z <- function(v) as.numeric(scale(v))
xcell_map <- c(Macrophages = "xcell_Macrophages",
               Fibroblasts = "xcell_Fibroblasts",
               Endothelial = "xcell_Endothelial_cells",
               Myocytes    = "xcell_Myocytes")
stopifnot(all(xcell_map %in% names(comp)))
for (nm in names(xcell_map)) {
  comp[, (paste0("z_xc_", nm)) := z(get(xcell_map[[nm]]))]
  comp[, (paste0("z_sf_", nm)) := z(get(paste0("ss_full_", nm)))]
  comp[, (paste0("z_sd_", nm)) := z(get(paste0("ss_dec_", nm)))]
}

## ===========================================================================
## 4. Four limma models
## ===========================================================================
hr()
log_msg("[4] limma models")

fit_one <- function(rhs, tag) {
  design <- model.matrix(as.formula(paste("~", rhs)), data = comp)
  f <- eBayes(lmFit(expr, design), trend = TRUE, robust = TRUE)
  tt <- topTable(f, coef = "group2Disease", number = Inf, sort.by = "none",
                 adjust.method = "BH")
  d <- data.table(GeneSymbol = rownames(tt), log2FC = round(tt$logFC, 4),
                  padj = tt$adj.P.Val)
  d[, sig := !is.na(padj) & padj < PADJ_CUT & abs(log2FC) > LFC_CUT]
  log_msg(sprintf("    %-10s significant = %3d", tag, sum(d$sig)))
  d
}

rhs0 <- "group2 + sex + age_num"
cov_xc <- "z_xc_Macrophages + z_xc_Fibroblasts + z_xc_Endothelial + z_xc_Myocytes"
cov_sf <- "z_sf_Macrophages + z_sf_Fibroblasts + z_sf_Endothelial + z_sf_Myocytes"
cov_sd <- "z_sd_Macrophages + z_sd_Fibroblasts + z_sd_Endothelial + z_sd_Myocytes"

m0 <- fit_one(rhs0, "M0")
m1 <- fit_one(paste(rhs0, "+", cov_xc), "M1_xCell")
m2 <- fit_one(paste(rhs0, "+", cov_sf), "M2_ssGSEA")
m3 <- fit_one(paste(rhs0, "+", cov_sd), "M3_dec")

## ===========================================================================
## 5. Attenuation per model
## ===========================================================================
hr()
log_msg("[5] Attenuation of the candidate genes")

summarise_model <- function(d, tag) {
  d <- copy(d); setnames(d, c("log2FC", "padj", "sig"),
                         c("log2FC_m", "padj_m", "sig_m"))
  mm <- merge(m0[, .(GeneSymbol, log2FC_base = log2FC)],
              d[, .(GeneSymbol, log2FC_m, padj_m, sig_m)], by = "GeneSymbol")
  mm[, att := 100 * (1 - abs(log2FC_m) / abs(log2FC_base))]
  cand_m <- mm[GeneSymbol %in% cand70]
  hub_m <- mm[GeneSymbol %in% HUBS]
  gene_m <- mm[GeneSymbol %in% GENES]
  data.table(
    model = tag,
    n_sig = sum(d$sig),
    n_sig_70 = sum(cand_m$sig_m),
    n_sig_hub = sum(hub_m$sig_m),
    median_att_70 = round(median(cand_m$att, na.rm = TRUE), 1),
    median_att_hub = round(median(hub_m$att, na.rm = TRUE), 1),
    median_att_4genes = round(median(gene_m$att, na.rm = TRUE), 1),
    median_abs_log2fc_70 = round(median(abs(cand_m$log2FC_m), na.rm = TRUE), 4)
  )
}

models <- rbindlist(list(
  summarise_model(m0, "M0 no composition"),
  summarise_model(m1, "M1 xCell (published)"),
  summarise_model(m2, "M2 ssGSEA complete"),
  summarise_model(m3, "M3 ssGSEA 70 genes removed")
))
write_csv(models, file.path(d_tab, "T21_decircularised_models.csv"))
log_msg("")
for (i in seq_len(nrow(models))) {
  log_msg(sprintf("    %-26s sig=%3d  70 genes sig=%2d  hub sig=%d  median att 70=%5.1f%%  hub=%5.1f%%  model genes=%5.1f%%",
                  models$model[i], models$n_sig[i], models$n_sig_70[i], models$n_sig_hub[i],
                  models$median_att_70[i], models$median_att_hub[i], models$median_att_4genes[i]))
}

## key gene table: three models side by side
d0 <- m0[, list(GeneSymbol, log2FC_base = log2FC, padj_base = padj)]
d1 <- m1[, list(GeneSymbol, log2FC_xc = log2FC, padj_xc = padj, sig_xc = sig)]
d2 <- m2[, list(GeneSymbol, log2FC_sf = log2FC, padj_sf = padj, sig_sf = sig)]
d3 <- m3[, list(GeneSymbol, log2FC_sd = log2FC, padj_sd = padj, sig_sd = sig)]
key <- merge(d0, d1, by = "GeneSymbol", sort = FALSE)
key <- merge(key, d2, by = "GeneSymbol", sort = FALSE)
key <- merge(key, d3, by = "GeneSymbol", sort = FALSE)
key <- key[GeneSymbol %in% c(GENES, HUBS)]
key[, role := ifelse(GeneSymbol %in% GENES, "model gene", "PPI hub")]
key[, att_xc := round(100 * (1 - abs(log2FC_xc) / abs(log2FC_base)), 1)]
key[, att_sf := round(100 * (1 - abs(log2FC_sf) / abs(log2FC_base)), 1)]
key[, att_sd := round(100 * (1 - abs(log2FC_sd) / abs(log2FC_base)), 1)]
setorder(key, role, GeneSymbol)
write_csv(key, file.path(d_tab, "T22_decircularised_key_genes.csv"))
log_msg("")
for (i in seq_len(nrow(key))) {
  log_msg(sprintf("    %-7s %-10s raw=%+.3f  xCell=%+.3f (%+5.1f%%)  ssGSEA=%+.3f (%+5.1f%%)  dec=%+.3f (%+5.1f%%)",
                  key$GeneSymbol[i], key$role[i], key$log2FC_base[i],
                  key$log2FC_xc[i], key$att_xc[i], key$log2FC_sf[i], key$att_sf[i],
                  key$log2FC_sd[i], key$att_sd[i]))
}

## ===========================================================================
## 6. Figure: attenuation before and after de-circularisation
## ===========================================================================
hr()
log_msg("[6] Figure")

plot_dt <- key[, .(GeneSymbol, role, att_xc, att_sf, att_sd)]
long <- melt(plot_dt, id.vars = c("GeneSymbol", "role"),
             variable.name = "model", value.name = "attenuation")
long[, model := factor(model,
                       levels = c("att_xc", "att_sf", "att_sd"),
                       labels = c("xCell (published)", "ssGSEA complete", "ssGSEA 70 genes removed"))]
long[, GeneSymbol := factor(GeneSymbol, levels = rev(plot_dt[order(role, att_sd)]$GeneSymbol))]

p <- ggplot(long, aes(attenuation, GeneSymbol, colour = model, shape = role)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_line(aes(group = GeneSymbol), colour = "grey80", linewidth = 0.4) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_colour_manual(values = c("#9A9A9A", "#4C78A8", "#C53030"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = "Attenuation of log2 fold change after composition adjustment (%)",
       y = NULL,
       title = "Composition adjustment before and after removing the candidate genes from the signature",
       subtitle = "Negative values indicate that the gene changed more after adjustment") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(d_fig, "F10_attenuation_decirculated.pdf"), p, width = 7.4, height = 4.4)
ggsave(file.path(d_fig, "F10_attenuation_decirculated.png"), p, width = 7.4, height = 4.4, dpi = 600)

hr("=")
log_msg("Done.")
hr("=")
writeLines(run_log, file.path(d_log, "136_run_log.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(d_log, "sessionInfo_136.txt"))
