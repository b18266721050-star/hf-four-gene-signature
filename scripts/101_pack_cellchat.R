# =============================================================================
# 101_pack_cellchat.R
# Stage 6 (A4) - summary numbers for the report / manuscript text.
#
# Input : cellchat/tables/98_*.csv, 99_*.csv, 100_*.csv
#         cellchat/objects/cc_donor.rds
# Output: cellchat/tables/101_summary_numbers.csv, 101_focal_nrg_detail.csv
#         cellchat/logs/101_pack.log
# =============================================================================

suppressPackageStartupMessages({ library(openxlsx) })

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "cellchat")
tab_d <- file.path(out_d, "tables"); obj_d <- file.path(out_d, "objects")
log_d <- file.path(out_d, "logs")
lgf <- file.path(log_d, "101_pack.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}
lg("=== 101 pack ===")

t99 <- read.csv(file.path(tab_d, "99_pathway_tests.csv"), stringsAsFactors = FALSE)
e99 <- read.csv(file.path(tab_d, "99_focal_edge_tests.csv"), stringsAsFactors = FALSE)
p100 <- read.csv(file.path(tab_d, "100_focal_pair_detail.csv"), stringsAsFactors = FALSE)
gsum <- read.csv(file.path(tab_d, "98_group_summary.csv"), stringsAsFactors = FALSE)
don <- readRDS(file.path(obj_d, "cc_donor.rds"))
score <- don$score; grp <- as.character(don$donor_grp$group); dep <- don$depth$median_nCount

nn <- t99[!t99$neuro, ]
nn <- nn[!is.na(nn$p), ]                     # pathways with no detectable flow have no test
sig <- nn[nn$fdr < 0.05, ]
lg("non-neuro pathways tested: ", nrow(nn), " | FDR<0.05: ", nrow(sig),
   " (up ", sum(sig$log2_ratio > 0, na.rm = TRUE), " / down ",
   sum(sig$log2_ratio < 0, na.rm = TRUE), ")")
lg("direction concordance abs vs share (all non-neuro): ",
   round(mean(sign(nn$log2_ratio) == sign(nn$log2_ratio_share), na.rm = TRUE) * 100, 1), "%")
lg("depth-adjusted p<0.05: ", sum(nn$p_adj_depth < 0.05, na.rm = TRUE),
   " | of those also abs FDR<0.05: ", sum(nn$p_adj_depth < 0.05 & nn$fdr < 0.05, na.rm = TRUE))

# donor depth is not systematically lower in the disease group -> report the numbers
lg("median nCount: NF ", median(dep[grp == "NF"]), " | NICM ", median(dep[grp == "NICM"]),
   " | DoxCM ", median(dep[grp == "DoxCM"]))

# focal pathway, per group
FOCAL <- "NRG"
v <- score[, FOCAL]
q <- tapply(v, grp, function(x) sprintf("%.2f (%.2f-%.2f)", median(x), quantile(x, .25), quantile(x, .75)))
p_nicm <- wilcox.test(v[grp == "NICM"], v[grp == "NF"], exact = FALSE)$p.value
p_dox  <- wilcox.test(v[grp == "DoxCM"], v[grp == "NF"], exact = FALSE)$p.value
lg(FOCAL, " score by group: ", paste(names(q), q, collapse = " | "))
lg("NICM vs NF p=", signif(p_nicm, 3), " | DoxCM vs NF p=", signif(p_dox, 3))
lg("cor(NRG score, donor depth) = ", round(cor(v, dep), 3),
   " | rho = ", round(cor(v, dep, method = "spearman"), 3))

summ <- data.frame(
  item = c("pathways_tested_non_neuro", "pathways_fdr_lt_0.05", "pathways_up", "pathways_down",
           "direction_concordance_pct", "depth_adjusted_p_lt_0.05",
           "NRG_median_NF", "NRG_median_NICM", "NRG_median_DoxCM",
           "NRG_p_NICM_vs_NF", "NRG_p_DoxCM_vs_NF",
           "NRG_log2FC", "NRG_p", "NRG_fdr", "NRG_p_share", "NRG_p_adj_depth",
           "NRG_edges_fdr_lt_0.05", "NRG_top_edge",
           "detected_LR_NF", "detected_LR_NICM", "detected_LR_DoxCM",
           "total_strength_NF", "total_strength_NICM", "total_strength_DoxCM"),
  value = c(nrow(nn), nrow(sig), sum(sig$log2_ratio > 0, na.rm = TRUE),
            sum(sig$log2_ratio < 0, na.rm = TRUE),
            round(mean(sign(nn$log2_ratio) == sign(nn$log2_ratio_share), na.rm = TRUE) * 100, 1),
            sum(nn$p_adj_depth < 0.05, na.rm = TRUE),
            median(v[grp == "NF"]), median(v[grp == "NICM"]), median(v[grp == "DoxCM"]),
            signif(p_nicm, 3), signif(p_dox, 3),
            signif(t99$log2_ratio[t99$pathway == FOCAL], 3), signif(t99$p[t99$pathway == FOCAL], 3),
            signif(t99$fdr[t99$pathway == FOCAL], 3), signif(t99$p_share[t99$pathway == FOCAL], 3),
            signif(t99$p_adj_depth[t99$pathway == FOCAL], 3),
            sum(e99$fdr[e99$pathway == FOCAL] < 0.05),
            paste0(e99$sender[e99$pathway == FOCAL][which.min(e99$p[e99$pathway == FOCAL])], " -> ",
                   e99$receiver[e99$pathway == FOCAL][which.min(e99$p[e99$pathway == FOCAL])]),
            gsum$n_interactions_detected[gsum$group == "NF"],
            gsum$n_interactions_detected[gsum$group == "NICM"],
            gsum$n_interactions_detected[gsum$group == "DoxCM"],
            round(gsum$total_strength[gsum$group == "NF"]),
            round(gsum$total_strength[gsum$group == "NICM"]),
            round(gsum$total_strength[gsum$group == "DoxCM"])))
write.csv(summ, file.path(tab_d, "101_summary_numbers.csv"), row.names = FALSE)
print(summ)

nrg_pair <- p100[p100$pathway == FOCAL, ]
nrg_pair <- nrg_pair[order(nrg_pair$p), ]
write.csv(nrg_pair, file.path(tab_d, "101_focal_nrg_detail.csv"), row.names = FALSE)
lg("NRG pairs with FDR<0.05: ", sum(nrg_pair$fdr < 0.05), "/", nrow(nrg_pair))

# one workbook for the supplementary tables of this module
wb <- createWorkbook()
addWorksheet(wb, "pathway_tests"); writeData(wb, "pathway_tests", t99)
addWorksheet(wb, "NRG_edges");     writeData(wb, "NRG_edges", e99[e99$pathway == FOCAL, ])
addWorksheet(wb, "NRG_pairs");     writeData(wb, "NRG_pairs", nrg_pair)
addWorksheet(wb, "group_summary"); writeData(wb, "group_summary", gsum)
addWorksheet(wb, "summary");       writeData(wb, "summary", summ)
saveWorkbook(wb, file.path(out_d, "cellchat_supplementary_tables.xlsx"), overwrite = TRUE)
lg("workbook written")
lg("=== done ===")
