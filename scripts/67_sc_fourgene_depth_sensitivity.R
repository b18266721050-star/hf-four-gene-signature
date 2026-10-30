# =============================================================================
# 67_sc_fourgene_depth_sensitivity.R
# Stage 5 (A3) - triangulation of the donor-level group effect.
#
# 66 showed that NF nuclei carry distinctly more RNA than disease nuclei
# (median nFeature 1944 vs 1584 / 1909; housekeeping detection 55% vs 31% / 39%).
# Three independent ways of asking "does the gene still change once depth is
# accounted for":
#   (1) raw pseudobulk CPM, donor level            (65)
#   (2) detection ratio vs a housekeeping panel    (66)
#   (3) depth-residualised pseudobulk              (67, this script)
# plus a mixed model with the per-donor nFeature median as a covariate, and the
# attenuation of the disease coefficient when the covariate enters.
#
# Input : 65_pseudobulk_sample_x_celltype.csv, 65_stats_*.csv,
#         66_donor_qc_medians.csv, 66_capture_control_donor.csv,
#         66_capture_control_stats.csv
# Output: tables/67_depth_models.csv, tables/67_depth_residual_tests.csv,
#         tables/67_triangulation.csv, 67_depth_sensitivity.log
# =============================================================================

suppressPackageStartupMessages({ library(lme4); library(car) })

SEED <- 20251029
proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "fourgene")
tab_d <- file.path(out_d, "tables")
lgf <- file.path(out_d, "67_depth_sensitivity.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

pb  <- read.csv(file.path(tab_d, "65_pseudobulk_sample_x_celltype.csv"), stringsAsFactors = FALSE)
qc  <- read.csv(file.path(tab_d, "66_donor_qc_medians.csv"), stringsAsFactors = FALSE)
cap <- read.csv(file.path(tab_d, "66_capture_control_donor.csv"), stringsAsFactors = FALSE)
capstat <- read.csv(file.path(tab_d, "66_capture_control_stats.csv"), stringsAsFactors = FALSE)
st_ct <- read.csv(file.path(tab_d, "65_stats_celltype.csv"), stringsAsFactors = FALSE)
st_dn <- read.csv(file.path(tab_d, "65_stats_donor_whole_tissue.csv"), stringsAsFactors = FALSE)
GENES <- c("AQP3", "GGT5", "ISLR", "LAD1")

pb <- pb[pb$usable, ]
pb <- merge(pb, qc[, c("sample", "nFeature_RNA", "nCount_RNA", "percent.mt")],
            by = "sample", all.x = TRUE)
pb$group <- factor(pb$group, levels = c("NF", "NICM", "DoxCM"))
pb$grp2 <- factor(ifelse(pb$group == "NF", "NF", "Disease"), levels = c("NF", "Disease"))
pb$depth <- as.numeric(scale(pb$nFeature_RNA))
lg("usable rows: ", nrow(pb), "; donors: ", length(unique(pb$sample)))

# ------------------------------------------- 1. mixed models, +/- covariate --
mods <- do.call(rbind, lapply(GENES, function(g) {
  d <- pb[pb$gene == g, ]
  m0 <- lmer(log1p_cpm ~ grp2 + celltype + (1 | sample), data = d, REML = FALSE)
  m1 <- lmer(log1p_cpm ~ grp2 + celltype + depth + (1 | sample), data = d, REML = FALSE)
  b0 <- fixef(m0)["grp2Disease"]; b1 <- fixef(m1)["grp2Disease"]
  s0 <- summary(m0)$coefficients; s1 <- summary(m1)$coefficients
  tcol0 <- grep("(z|t) value", colnames(s0))[1]
  tcol1 <- grep("(z|t) value", colnames(s1))[1]
  p0 <- 2 * pnorm(-abs(s0["grp2Disease", tcol0]))
  p1 <- 2 * pnorm(-abs(s1["grp2Disease", tcol1]))
  pd <- 2 * pnorm(-abs(s1["depth", tcol1]))
  a0 <- car::Anova(m0, type = 2, test.statistic = "Chisq")
  a1 <- car::Anova(m1, type = 2, test.statistic = "Chisq")
  data.frame(gene = g,
             beta_unadj = b0, p_unadj = p0,
             chisq_group_unadj = a0$Chisq[rownames(a0) == "grp2"],
             p_group_unadj = a0$`Pr(>Chisq)`[rownames(a0) == "grp2"],
             beta_adj = b1, p_adj = p1,
             beta_depth = fixef(m1)["depth"], p_depth = pd,
             chisq_group_adj = a1$Chisq[rownames(a1) == "grp2"],
             p_group_adj = a1$`Pr(>Chisq)`[rownames(a1) == "grp2"],
             attenuation = 1 - b1 / b0,
             n_obs = nrow(d), n_donor = length(unique(d$sample)),
             stringsAsFactors = FALSE)
}))
write.csv(mods, file.path(tab_d, "67_depth_models.csv"), row.names = FALSE)
lg("--- mixed model: disease coefficient with / without depth covariate ---")
for (i in seq_len(nrow(mods)))
  lg(sprintf("  %-5s beta %+.3f (p=%.3g) -> %+.3f (p=%.3g) | attenuation %+.0f%% | depth p=%.3g",
             mods$gene[i], mods$beta_unadj[i], mods$p_unadj[i], mods$beta_adj[i],
             mods$p_adj[i], mods$attenuation[i] * 100, mods$p_depth[i]))

# --------------------------------- 2. depth-residualised donor-level tests ---
res <- do.call(rbind, lapply(GENES, function(g) {
  do.call(rbind, lapply(unique(pb$celltype), function(ct) {
    d <- pb[pb$gene == g & pb$celltype == ct, ]
    if (nrow(d) < 10 || length(unique(d$group)) < 3) return(NULL)
    fit <- lm(log1p_cpm ~ nFeature_RNA, data = d)
    d$resid <- resid(fit)
    nf <- d$resid[d$group == "NF"]; dis <- d$resid[d$group != "NF"]
    data.frame(gene = g, celltype = ct,
               rho_depth = suppressWarnings(cor(d$log1p_cpm, d$nFeature_RNA, method = "spearman")),
               delta_resid = median(dis) - median(nf),
               p_resid = tryCatch(wilcox.test(dis, nf, exact = FALSE)$p.value,
                                  error = function(e) NA_real_),
               n_NF = length(nf), n_Disease = length(dis),
               stringsAsFactors = FALSE)
  }))
}))
res$fdr_resid <- p.adjust(res$p_resid, "BH")
res <- res[order(res$gene, res$fdr_resid), ]
write.csv(res, file.path(tab_d, "67_depth_residual_tests.csv"), row.names = FALSE)
lg("depth-residualised celltype tests at FDR<0.05: ", sum(res$fdr_resid < 0.05, na.rm = TRUE),
   " / ", sum(!is.na(res$fdr_resid)))
lg("median |rho| between depth and expression: ", round(median(abs(res$rho_depth), na.rm = TRUE), 3))

# ------------------------------------------------ 3. triangulation table -----
tri <- do.call(rbind, lapply(GENES, function(g) {
  d <- st_dn[st_dn$gene == g, ]
  r <- capstat[capstat$gene == g, ]
  rsel <- res[res$celltype == "Endothelial cell" | TRUE, ]
  main_ct <- st_ct[st_ct$gene == g, ]
  main_ct <- main_ct[order(main_ct$fdr), ][1, ]
  rr <- res[res$gene == g & res$celltype == main_ct$celltype, ]
  data.frame(
    gene = g,
    bulk_log2FC = d$bulk_log2FC,
    pseudobulk_log2FC = d$log2FC_median,
    pseudobulk_FDR = d$fdr,
    main_cell_type = main_ct$celltype,
    main_ct_log2FC = main_ct$log2FC_median,
    main_ct_FDR = main_ct$fdr,
    hk_ratio_fold = r$fold_ratio,
    hk_ratio_FDR = r$fdr_ratio,
    depth_resid_main_delta = rr$delta_resid[1],
    depth_resid_main_FDR = rr$fdr_resid[1],
    stringsAsFactors = FALSE)
}))
tri$direction_bulk <- sign(tri$bulk_log2FC)
tri$concordant_pseudobulk <- sign(tri$pseudobulk_log2FC) == tri$direction_bulk
tri$concordant_hk <- sign(tri$hk_ratio_fold - 1) == tri$direction_bulk
tri$concordant_resid <- sign(tri$depth_resid_main_delta) == tri$direction_bulk
tri$robust <- tri$concordant_pseudobulk & tri$concordant_hk & tri$concordant_resid
write.csv(tri, file.path(tab_d, "67_triangulation.csv"), row.names = FALSE)

lg("--- triangulation (bulk direction vs three single-cell readouts) ---")
for (i in seq_len(nrow(tri)))
  lg(sprintf("  %-5s bulk %+.2f | pseudo %+.2f (FDR %.3g) | HK-ratio x%.2f (FDR %.3g) | depth-resid %+.2f (FDR %.3g) | robust=%s",
             tri$gene[i], tri$bulk_log2FC[i], tri$pseudobulk_log2FC[i], tri$pseudobulk_FDR[i],
             tri$hk_ratio_fold[i], tri$hk_ratio_FDR[i], tri$depth_resid_main_delta[i],
             tri$depth_resid_main_FDR[i], tri$robust[i]))
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_67.txt"))
lg("sessionInfo written: sessionInfo_67.txt")
lg("=== done ===")
