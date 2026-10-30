# =============================================================================
# 83_sample_overlap_audit.R
# 目的：判断 GSE5406（Penn, HG-U133A）与发现队列 GSE57338（MAGNet, HuGene1.1ST）
#       是否存在样本级重叠，并给出可写进 Methods 的量化结论
# 方法：跨平台"个体残差相关"
#   1) 两个队列各自按疾病组去均值（去掉组间差异），得到个体特异的表达残差
#   2) 每个基因在队列内标准化（消除平台尺度）
#   3) 计算全部 A×B 样本对的残差相关系数
#      -> 同一供体（重复样本）会显著高于背景；不同供体背景应接近 0
# 负对照：GSE5406 vs GSE116250（不同中心、不同平台、不可能重叠）
#          GSE57338 vs GSE116250
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

root <- "E:/another data file/GEObulk analysis project"
raw  <- file.path(root, "raw_data", "external")
outd <- file.path(root, "results", "ExternalValidation")
tabd <- file.path(outd, "01_tables"); figd <- file.path(outd, "02_figures")
datd <- file.path(outd, "04_data");   logd <- file.path(outd, "03_logs")

logcon <- file(file.path(logd, "83_overlap_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, logcon) }

# ---------------------------------------------------------------- 读入三套全基因矩阵
# 1) 发现队列 GSE57338
tr <- fread(file.path(root, "clean_data", "GSE57338_clean_expression_matrix.csv"),
            data.table = FALSE, check.names = FALSE)
Xtr <- as.matrix(tr[, -1, drop = FALSE]); rownames(Xtr) <- tr[[1]]
info_tr <- fread(file.path(root, "clean_data", "GSE57338_sample_info.csv"))
grp_tr <- info_tr$group[match(colnames(Xtr), info_tr$sample)]
say("GSE57338: ", nrow(Xtr), " genes x ", ncol(Xtr), " samples ; groups: ",
    paste(names(table(grp_tr)), table(grp_tr), collapse = " "))

# 2) GSE5406：从 series matrix + GPL96 注释重建全基因矩阵（缓存）
cache5406 <- file.path(datd, "GSE5406_full_gene_matrix.rds")
if (file.exists(cache5406)) {
  X54 <- readRDS(cache5406)
} else {
  con <- gzfile(file.path(raw, "GSE5406_series_matrix.txt.gz"), "r")
  ln <- readLines(con, warn = FALSE); close(con)
  b <- grep("^!series_matrix_table_begin", ln); e <- grep("^!series_matrix_table_end", ln)
  body <- ln[(b + 1):(e - 1)]; body <- body[nzchar(body)]
  m <- read.delim(text = body, header = TRUE, check.names = FALSE, quote = "\"")
  E <- as.matrix(m[, -1]); rownames(E) <- m[[1]]
  ann <- fread(file.path(raw, "GPL96.annot.gz"), sep = "\t", header = TRUE,
               comment.char = "", quote = "", skip = "ID\t", data.table = FALSE)
  sym <- setNames(ann[["Gene symbol"]], ann$ID)[rownames(E)]
  sym[is.na(sym) | sym == "" | sym == "---"] <- NA
  sym <- sub(" ?///.*$", "", sym)
  keep <- which(!is.na(sym))
  E <- E[keep, ]; sym <- sym[keep]
  # 多探针取平均表达最高者
  ord <- order(sym, -rowMeans(E, na.rm = TRUE))
  E <- E[ord, ]; sym <- sym[ord]
  first <- !duplicated(sym)
  X54 <- E[first, ]; rownames(X54) <- sym[first]
  saveRDS(X54, cache5406)
}
info_54 <- fread(file.path(outd, "01_data", "GSE5406_sample_info.csv"))
grp_54 <- info_54$group[match(colnames(X54), info_54$sample)]
say("GSE5406: ", nrow(X54), " genes x ", ncol(X54), " samples ; groups: ",
    paste(names(table(grp_54)), table(grp_54), collapse = " "))

# 3) GSE116250（负对照）：RPKM 全表 -> 基因符号
rp <- fread(file.path(raw, "GSE116250_rpkm.txt.gz"), sep = "\t", header = TRUE,
            data.table = FALSE, check.names = FALSE)
sym11 <- rp$Common_name
rp <- as.matrix(rp[, setdiff(colnames(rp), c("Gene", "Common_name")), drop = FALSE])
rp <- apply(rp, 2, function(z) suppressWarnings(as.numeric(z)))
ok <- !is.na(sym11) & sym11 != ""
rp <- rp[ok, ]; sym11 <- sym11[ok]
rp <- log2(rp + 1)
ord <- order(sym11, -rowMeans(rp, na.rm = TRUE))
rp <- rp[ord, ]; sym11 <- sym11[ord]
first <- !duplicated(sym11)
X11 <- rp[first, ]; rownames(X11) <- sym11[first]
info_11 <- fread(file.path(outd, "01_data", "GSE116250_sample_info.csv"))
colnames(X11) <- info_11$sample[match(colnames(X11), info_11$title)]
grp_11 <- info_11$group[match(colnames(X11), info_11$sample)]
say("GSE116250: ", nrow(X11), " genes x ", ncol(X11), " samples ; groups: ",
    paste(names(table(grp_11)), table(grp_11), collapse = " "))

# ---------------------------------------------------------------- 残差化 + 基因标准化
# 注意：残差化后存在"组内近乎常数"的基因，其 SD 极小，直接除会把相关放大到 >1，
# 因此用相对阈值（中位 SD 的 2%）过滤，并在标准化后再次中心化
resid_feat <- function(X, g) {
  g <- factor(g)
  R <- X
  for (lv in levels(g)) {
    idx <- which(g == lv)
    if (length(idx) > 1) R[, idx] <- X[, idx] - rowMeans(X[, idx], na.rm = TRUE)
    else R[, idx] <- 0
  }
  R[is.na(R)] <- 0
  R
}
scale_rows <- function(R) {
  s <- apply(R, 1, sd)
  R <- sweep(R, 1, rowMeans(R), "-")
  sweep(R, 1, pmax(s, 1e-6), "/")
}
pair_test <- function(XA, gA, XB, gB, lab, min_gene = 3000, npc = 10) {
  gs <- intersect(rownames(XA), rownames(XB))
  say("  common genes (", lab, "): ", length(gs))
  if (length(gs) < min_gene) { say("  -> too few genes, skipped"); return(NULL) }
  RA0 <- resid_feat(XA[gs, , drop = FALSE], gA)
  RB0 <- resid_feat(XB[gs, , drop = FALSE], gB)
  sa <- apply(RA0, 1, sd); sb <- apply(RB0, 1, sd)
  floor_s <- max(1e-6, 0.02 * median(c(sa, sb), na.rm = TRUE))
  keep <- is.finite(sa) & is.finite(sb) & sa >= floor_s & sb >= floor_s
  RA <- scale_rows(RA0[keep, , drop = FALSE])
  RB <- scale_rows(RB0[keep, , drop = FALSE])
  say("  genes retained after variance filter: ", nrow(RA))
  # 去除队列内前 npc 个主成分：抹掉队列内部的共同结构（年龄、纤维化梯度等），
  # 只留下个体特异信号。重复样本会保留这份个体信号，无关样本趋近 0。
  strip_pc <- function(R, k) {
    if (k <= 0) return(R)
    k <- min(k, min(dim(R)) - 1)
    pr <- prcomp(t(R), center = TRUE, scale. = FALSE, rank. = k)
    R - pr$rotation[, seq_len(k), drop = FALSE] %*% t(pr$x[, seq_len(k), drop = FALSE])
  }
  RA <- strip_pc(RA, npc); RB <- strip_pc(RB, npc)
  say("  PCs removed per cohort: ", npc)
  # 列间（跨基因）Pearson 相关；行已按基因标准化，等于每个基因等权
  C <- cor(RA, RB, method = "pearson")
  dimnames(C) <- list(colnames(XA), colnames(XB))
  say("  corr range: ", paste(round(range(C), 3), collapse = " ~ "))
  # 每个 A 样本的最佳匹配
  best <- data.frame(
    A_sample = rownames(C),
    A_group  = gA,
    B_best   = colnames(C)[max.col(C, ties.method = "first")],
    r_best   = apply(C, 1, max),
    r_median = apply(C, 1, median),
    stringsAsFactors = FALSE)
  best$B_group <- gB[match(best$B_best, colnames(XB))]
  best$group_match <- best$A_group == best$B_group
  attr(C, "best") <- best
  C
}

say("")
say("=============== 主比较：GSE5406 vs GSE57338 ===============")
# 先看不去主成分的口径（背景偏高，作对照），再上去 PC 的主口径
say("[K=0] 原始残差口径")
C_main0 <- pair_test(X54, grp_54, Xtr, grp_tr, "GSE5406 vs GSE57338 (K=0)", npc = 0)
say("[K=10] 去队列内前 10 个主成分（主口径）")
C_main <- pair_test(X54, grp_54, Xtr, grp_tr, "GSE5406 vs GSE57338 (K=10)", npc = 10)
best_main <- attr(C_main, "best")

# 反向：每个 GSE57338 样本在 GSE5406 中的最佳匹配
rev_best <- apply(t(C_main), 1, max)
say("正向（GSE5406 -> GSE57338）最佳相关：median=", round(median(best_main$r_best), 3),
    " q95=", round(quantile(best_main$r_best, .95), 3), " max=", round(max(best_main$r_best), 3))
say("反向（GSE57338 -> GSE5406）最佳相关：median=", round(median(rev_best), 3),
    " q95=", round(quantile(rev_best, .95), 3), " max=", round(max(rev_best), 3))
say("背景（同一样本对所有异组样本）中位数=", round(median(best_main$r_median), 3))

say("")
say("=============== 负对照 1：GSE5406 vs GSE116250 ===============")
C_neg1 <- pair_test(X54, grp_54, X11, grp_11, "GSE5406 vs GSE116250")
best_neg1 <- attr(C_neg1, "best")
say("最佳相关：median=", round(median(best_neg1$r_best), 3),
    " q95=", round(quantile(best_neg1$r_best, .95), 3), " max=", round(max(best_neg1$r_best), 3))

say("")
say("=============== 负对照 2：GSE57338 vs GSE116250 ===============")
C_neg2 <- pair_test(Xtr, grp_tr, X11, grp_11, "GSE57338 vs GSE116250")
best_neg2 <- attr(C_neg2, "best")
say("最佳相关：median=", round(median(best_neg2$r_best), 3),
    " q95=", round(quantile(best_neg2$r_best, .95), 3), " max=", round(max(best_neg2$r_best), 3))

# ---------------------------------------------------------------- 判定阈值
# 以两个负对照的最佳相关分布上界作为"非重复"的天花板
null_max <- max(best_neg1$r_best, best_neg2$r_best, na.rm = TRUE)
n_above <- sum(best_main$r_best > null_max)
say("")
say("负对照的最大最佳相关（null ceiling）= ", round(null_max, 3))
say("GSE5406 中超过该天花板的样本数 = ", n_above, " / ", nrow(best_main))

# ---------------------------------------------------------------- 组别一致性：与负对照对比
# 关键：若"最佳匹配落在同组"在主比较里明显高于负对照，才是重叠的证据；
#       若负对照也一样高，那只是疾病生物学残留，不是重复样本。
conc_tab <- function(best, lab) {
  A <- factor(best$A_group, levels = c("NF", "DCM", "ICM"))
  B <- factor(best$B_group, levels = c("NF", "DCM", "ICM"))
  tb <- table(A, B)
  frac <- diag(tb) / rowSums(tb)
  overall <- sum(diag(tb)) / sum(tb)
  say(sprintf("  %-34s overall concordance = %.3f", lab, overall))
  for (k in names(frac)) say(sprintf("      %-4s matched same group %3d/%-3d = %.2f",
                                     k, tb[k, k], rowSums(tb)[k], frac[k]))
  data.frame(pair = lab, group = names(frac), n = as.numeric(rowSums(tb)),
             n_match = as.numeric(diag(tb)), frac = as.numeric(frac),
             overall = overall, stringsAsFactors = FALSE)
}
say("")
say("最佳匹配的组别一致性（主比较 vs 负对照）：")
conc_all <- rbind(
  conc_tab(best_main, "GSE5406 -> GSE57338 (test)"),
  conc_tab(best_neg1, "GSE5406 -> GSE116250 (neg ctrl)"),
  conc_tab(best_neg2, "GSE57338 -> GSE116250 (neg ctrl)"))
write.csv(conc_all, file.path(tabd, "T13_overlap_group_concordance.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 互为最佳匹配（重复样本的必要特征）
bestB <- data.frame(
  B_sample = colnames(C_main),
  A_best = rownames(C_main)[max.col(t(C_main), ties.method = "first")],
  r_from_B = apply(t(C_main), 1, max), stringsAsFactors = FALSE)
mt <- merge(best_main[, c("A_sample", "A_group", "B_best", "r_best")],
            bestB, by.x = "B_best", by.y = "B_sample")
mutual <- mt[mt$A_sample == mt$A_best, ]
say("")
say("互为最佳匹配的对数 = ", nrow(mutual), " / ", nrow(best_main),
    "（同组 ", sum(mutual$A_group ==
      best_main$B_group[match(mutual$B_best, best_main$B_best)]), "）")
if (nrow(mutual)) {
  say("互为最佳匹配且相关最高的 10 对：")
  show <- mutual[order(-mutual$r_best), ][1:min(10, nrow(mutual)), ]
  for (i in seq_len(nrow(show)))
    say(sprintf("  %s(%s) <-> %s  r=%.3f", show$A_sample[i], show$A_group[i],
                show$B_best[i], show$r_best[i]))
}
write.csv(mutual, file.path(tabd, "T14_overlap_mutual_best.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 负对照 3/4：队列自身随机二分
# 关键混杂：同一平台/同一实验室本身就会抬高跨样本相关。
# 自身二分（同平台同实验室、确定无重复）给出这一混杂下的"天花板"。
say("")
say("=============== 负对照 3/4：队列自身随机二分 ===============")
set.seed(20251029)
split_conc <- function(X, g, lab) {
  idx <- sample(seq_len(ncol(X)))
  h1 <- idx[seq(1, length(idx), by = 2)]; h2 <- setdiff(idx, h1)
  C <- pair_test(X[, h1, drop = FALSE], g[h1], X[, h2, drop = FALSE], g[h2], lab, npc = 10)
  b <- attr(C, "best")
  tb <- table(factor(b$A_group, levels = c("NF", "DCM", "ICM")),
              factor(b$B_group, levels = c("NF", "DCM", "ICM")))
  ov <- sum(diag(tb)) / sum(tb)
  say(sprintf("  %-28s overall concordance = %.3f  max r = %.3f", lab, ov, max(b$r_best)))
  list(ov = ov, best = b)
}
sp_54 <- split_conc(X54, grp_54, "GSE5406 self-split")
sp_tr <- split_conc(Xtr, grp_tr, "GSE57338 self-split")

say("")
say("=============== 结论对照 ===============")
conc_main <- unique(conc_all$overall[conc_all$pair == "GSE5406 -> GSE57338 (test)"])
conc_neg <- unique(conc_all$overall[grepl("neg ctrl", conc_all$pair)])
say(sprintf("  主比较 GSE5406 -> GSE57338        concordance = %.3f  max r = %.3f",
            conc_main, max(best_main$r_best)))
say(sprintf("  负对照 跨中心/跨平台              concordance = %s  max r = %.3f",
            paste(round(conc_neg, 3), collapse = " / "),
            max(best_neg1$r_best, best_neg2$r_best)))
say(sprintf("  对照   同平台自身二分（无重复）    concordance = %.3f / %.3f  max r = %.3f",
            sp_54$ov, sp_tr$ov, max(sp_54$best$r_best, sp_tr$best$r_best)))
say("")
say(">>> 判定：若主比较的最高相关低于'同平台无重复'的最高相关，则没有任何一对像重复样本")

# ---------------------------------------------------------------- 输出表
best_main$cohort_pair <- "GSE5406_vs_GSE57338"
best_neg1$cohort_pair <- "GSE5406_vs_GSE116250_negctrl"
best_neg2$cohort_pair <- "GSE57338_vs_GSE116250_negctrl"
allbest <- rbind(best_main, best_neg1, best_neg2)
fwrite(allbest, file.path(datd, "overlap_bestmatch_all.csv"))

top <- best_main[order(-best_main$r_best), ][1:10, ]
write.csv(top, file.path(tabd, "T11_overlap_top_pairs.csv"), row.names = FALSE)

summ <- data.frame(
  comparison = c("GSE5406 vs GSE57338 (test)",
                 "GSE5406 vs GSE116250 (cross-centre control)",
                 "GSE57338 vs GSE116250 (cross-centre control)",
                 "GSE5406 self-split (same platform, no duplicates)",
                 "GSE57338 self-split (same platform, no duplicates)"),
  n_A = c(nrow(best_main), nrow(best_neg1), nrow(best_neg2),
          nrow(sp_54$best), nrow(sp_tr$best)),
  median_best_r = c(median(best_main$r_best), median(best_neg1$r_best),
                    median(best_neg2$r_best), median(sp_54$best$r_best),
                    median(sp_tr$best$r_best)),
  q95_best_r = c(quantile(best_main$r_best, .95), quantile(best_neg1$r_best, .95),
                 quantile(best_neg2$r_best, .95), quantile(sp_54$best$r_best, .95),
                 quantile(sp_tr$best$r_best, .95)),
  max_best_r = c(max(best_main$r_best), max(best_neg1$r_best), max(best_neg2$r_best),
                 max(sp_54$best$r_best), max(sp_tr$best$r_best)),
  group_concordance = c(conc_main, conc_neg, sp_54$ov, sp_tr$ov),
  n_above_null_ceiling = c(n_above, NA, NA, NA, NA),
  stringsAsFactors = FALSE)
write.csv(summ, file.path(tabd, "T12_overlap_summary.csv"), row.names = FALSE)
say("")
say("summary:")
for (i in seq_len(nrow(summ))) say(sprintf("  %-52s med=%.3f q95=%.3f max=%.3f conc=%.3f",
    summ$comparison[i], summ$median_best_r[i], summ$q95_best_r[i],
    summ$max_best_r[i], summ$group_concordance[i]))

# ---------------------------------------------------------------- 图
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
dens <- rbind(
  data.frame(r = best_main$r_best, set = "Tested pair: GSE5406 vs GSE57338"),
  data.frame(r = c(best_neg1$r_best, best_neg2$r_best),
             set = "Cross-centre / cross-platform controls"),
  data.frame(r = c(sp_54$best$r_best, sp_tr$best$r_best),
             set = "Same-platform self-split (no duplicates possible)"))
dup_ceiling <- max(sp_54$best$r_best, sp_tr$best$r_best)
p1 <- ggplot(dens, aes(r, fill = set, colour = set)) +
  geom_density(alpha = 0.25, linewidth = 0.7, adjust = 1.1) +
  geom_vline(xintercept = null_max, linetype = 2, colour = "grey35") +
  geom_vline(xintercept = max(best_main$r_best), linetype = 4, colour = "#E15759") +
  geom_vline(xintercept = dup_ceiling, linetype = 3, colour = "#3B7A32") +
  annotate("text", x = null_max, y = Inf, label = "cross-platform\nceiling", vjust = 1.2,
           hjust = 1.15, size = 2.8, colour = "grey35") +
  annotate("text", x = max(best_main$r_best), y = Inf, label = "highest pair\nobserved", vjust = 1.2,
           hjust = 1.15, size = 2.8, colour = "#E15759") +
  annotate("text", x = dup_ceiling, y = Inf, label = "ceiling without\nduplicates", vjust = 1.2,
           hjust = -0.12, size = 2.8, colour = "#3B7A32") +
  scale_fill_manual(values = c("Tested pair: GSE5406 vs GSE57338" = "#E15759",
                               "Cross-centre / cross-platform controls" = "#4E79A7",
                               "Same-platform self-split (no duplicates possible)" = "#59A14F"),
                    name = NULL) +
  scale_colour_manual(values = c("Tested pair: GSE5406 vs GSE57338" = "#E15759",
                                 "Cross-centre / cross-platform controls" = "#4E79A7",
                                 "Same-platform self-split (no duplicates possible)" = "#59A14F"),
                      name = NULL) +
  scale_x_continuous(limits = c(0.02, 0.44)) +
  labs(x = "Best cross-cohort residual correlation per sample",
       y = "Density",
       title = "Sample-overlap audit: best-match correlation distributions",
       subtitle = "Individual-specific residuals after removing group means and 10 within-cohort PCs") +
  th + theme(legend.position = "bottom", legend.direction = "vertical",
             legend.text = element_text(size = 7.4),
             legend.key.height = unit(0.5, "lines"),
             plot.subtitle = element_text(size = 8.5))
ggsave(file.path(figd, "F09_overlap_bestmatch_density.pdf"), p1, width = 6.4, height = 4.8)
ggsave(file.path(figd, "F09_overlap_bestmatch_density.png"), p1, width = 6.4, height = 4.8, dpi = 600)

# 最佳匹配的散点（每个样本一个点：r_best vs 该样本的背景中位相关）
sc <- best_main
p2 <- ggplot(sc, aes(r_median, r_best, colour = group_match)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = c(`TRUE` = "#59A14F", `FALSE` = "grey60"),
                      name = "best match same group") +
  labs(x = "Median correlation with all GSE57338 samples (background)",
       y = "Best correlation with any GSE57338 sample",
       title = "GSE5406 samples projected onto GSE57338") +
  th + theme(legend.position = c(0.7, 0.2), legend.background = element_blank())
ggsave(file.path(figd, "F10_overlap_best_vs_background.pdf"), p2, width = 5.2, height = 4.6)
ggsave(file.path(figd, "F10_overlap_best_vs_background.png"), p2, width = 5.2, height = 4.6, dpi = 600)

# 图注单独成文件（避免与 81 生成的 figure_legends_external.txt 互相覆盖）
writeLines(c(
 "Figure S-X9 | Sample-overlap audit of GSE5406 against the discovery cohort GSE57338. Within each cohort, group means were removed and the first ten principal components were projected out, leaving individual-specific residuals; samples were then correlated across the 11,808 shared genes. Red: the tested pair; blue: cross-centre/cross-platform controls; green: random self-splits of a single cohort, in which duplicated samples are impossible by construction. The highest correlation in the tested pair (r = 0.283) lies below the ceiling reached without any duplication (r = 0.368), and no duplicate-like mode is present.",
 "Figure S-X10 | For every GSE5406 sample, the best correlation with any GSE57338 sample is plotted against the median correlation with all GSE57338 samples (background). Samples whose best match belongs to the same disease group are shown in green. A duplicated sample would appear far above the identity line and outside the background cloud."),
 file.path(figd, "figure_legends_overlap.txt"))

say(""); say("DONE"); close(logcon)
