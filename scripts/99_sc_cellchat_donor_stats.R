# =============================================================================
# 99_sc_cellchat_donor_stats.R
# Stage 6 (A4) - donor-level communication statistics (the inferential step).
#
# Input : cellchat/objects/{cc_input.rds, cc_networks.rds}
# Output: cellchat/objects/cc_donor.rds
#         cellchat/tables/99_{donor_pathway_scores,pathway_tests,donor_depth}.csv
#         cellchat/logs/99_donor_stats.log
#
# Why donor level: 25 donors (NF 8 / NICM 6 / DoxCM 11). Pooling 150k nuclei and
# running CellChat's cell-level permutation test would inflate the sample size by
# four orders of magnitude. Here each donor contributes one pseudobulk value per
# (celltype, gene), the communication score is computed per donor, and the test
# is a rank test over 25 values. Nuclei are never replicates.
#
# Three read-outs per pathway, because the disease group has a different depth
# and cell-composition profile:
#   abs      raw summed probability (confounded by depth / composition)
#   share    pathway score / donor's non-neurotransmitter total (compositional,
#            removes a global expression shift shared by all pathways)
#   depth-adjusted linear model on log10(median nCount)
# =============================================================================

suppressPackageStartupMessages({ library(Matrix) })

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "cellchat")
tab_d <- file.path(out_d, "tables"); obj_d <- file.path(out_d, "objects")
log_d <- file.path(out_d, "logs")
lgf <- file.path(log_d, "99_donor_stats.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

KH <- 0.5; EXPR_FLOOR <- 0.05; MIN_CELLS <- 20
GRP <- c("NF", "NICM", "DoxCM")
lg("=== 99 donor-level communication stats | Kh=", KH, " min_cells=", MIN_CELLS, " ===")

cc <- readRDS(file.path(obj_d, "cc_input.rds"))
net <- readRDS(file.path(obj_d, "cc_networks.rds"))
cnt <- cc$counts; meta <- cc$meta; pairs <- cc$pairs; CT <- cc$ct_order
common <- net$common
lg("LR records used: ", length(common), " | donors: ", length(unique(meta$sample)))

# ------------------------------------------------------------ pseudobulk -----
lib <- meta$nCount_RNA / 1e4
dat <- cnt %*% Diagonal(x = 1 / lib); dat@x <- log1p(dat@x)

samples <- sort(unique(meta$sample))
donor_grp <- unique(meta[, c("sample", "group", "condition")])
donor_grp <- donor_grp[match(samples, donor_grp$sample), ]
stopifnot(!any(is.na(donor_grp$group)))

strata <- interaction(meta$sample, meta$celltype, drop = TRUE, lex.order = TRUE)
S2 <- sparseMatrix(i = seq_len(nrow(meta)), j = as.integer(strata), x = 1,
                   dims = c(nrow(meta), nlevels(strata)))
n_str <- as.numeric(colSums(S2))
sum_expr <- as.matrix(dat %*% S2)
mean_expr <- sweep(sum_expr, 2, pmax(n_str, 1), "/")
lg("pseudobulk strata: ", ncol(mean_expr), " (sample x celltype) | cells median ",
   median(n_str), " | min ", min(n_str))

# split strata label back into sample / celltype
lv <- levels(strata)
sp <- sub("\\..*$", "", lv)
ct <- sub("^.*\\.", "", lv)
sp <- factor(sp, levels = samples); ct <- factor(ct, levels = CT)
usable <- n_str >= MIN_CELLS
lg("usable strata (>= ", MIN_CELLS, " cells): ", sum(usable), "/", length(usable))

# genes x celltype x donor array
E <- array(0, dim = c(nrow(mean_expr), length(CT), length(samples)),
           dimnames = list(rownames(mean_expr), CT, samples))
for (k in which(usable)) E[, as.integer(ct[k]), as.integer(sp[k])] <- mean_expr[, k]
lg("donor expression array: ", paste(dim(E), collapse = " x "))

# ------------------------------------------------------------ probabilities --
pair_ids <- match(common, pairs$interaction_name)
lig_list <- pairs$ligand[pair_ids]; rec_list <- pairs$receptor[pair_ids]
pw <- pairs$pathway[pair_ids]
np <- length(common)

sub_expr <- function(subs) {
  subs <- subs[subs %in% dimnames(E)[[1]]]
  if (length(subs) == 0) return(NULL)
  if (length(subs) == 1) return(E[subs, , ])
  apply(E[subs, , , drop = FALSE], c(2, 3), min)
}

L <- array(0, dim = c(np, length(CT), length(samples)))   # ligand per sender
R <- array(0, dim = c(np, length(CT), length(samples)))   # receptor per receiver
for (i in seq_len(np)) {
  a <- sub_expr(lig_list[[i]]); b <- sub_expr(rec_list[[i]])
  if (!is.null(a)) L[i, , ] <- a
  if (!is.null(b)) R[i, , ] <- b
  if (i %% 500 == 0) lg("  ... subunit collapse ", i, "/", np)
}

prob_donor <- array(0, dim = c(np, length(CT), length(CT), length(samples)))
for (d in seq_len(length(samples))) {
  for (i in seq_len(length(CT))) for (j in seq_len(length(CT))) {
    x <- L[, i, d] * R[, j, d]
    p <- x / (x + KH)
    p[x < EXPR_FLOOR] <- 0
    prob_donor[, i, j, d] <- p
  }
}
lg("donor probability array: ", paste(dim(prob_donor), collapse = " x "),
   " (", round(object.size(prob_donor) / 1e6, 1), " MB )")

# -------------------------------------------------------- pathway scores -----
pw_lev <- sort(unique(pw))
pw_idx <- match(pw, pw_lev)
score <- matrix(0, nrow = length(samples), ncol = length(pw_lev),
                dimnames = list(samples, pw_lev))
for (d in seq_len(length(samples))) {
  s <- apply(prob_donor[, , , d], 1, sum)                  # per LR record
  score[d, ] <- as.numeric(rowsum(s, pw_idx, reorder = TRUE)[, 1])
}
NEURO <- cc$neuro
non_neuro <- setdiff(pw_lev, NEURO)
tot_non_neuro <- rowSums(score[, non_neuro])
share <- score / tot_non_neuro

donor_depth <- data.frame(
  sample = samples,
  group  = as.character(donor_grp$group),
  n_nuclei = as.integer(table(meta$sample)[samples]),
  median_nCount = as.numeric(tapply(meta$nCount_RNA, meta$sample, median)[samples]),
  total_non_neuro_score = as.numeric(tot_non_neuro[samples]),
  stringsAsFactors = FALSE)
write.csv(donor_depth, file.path(tab_d, "99_donor_depth.csv"), row.names = FALSE)
lg("donor depth table:"); print(donor_depth)

# --------------------------------------------------------------- tests -------
cliff <- function(x, y) {
  o <- outer(x, y, "-"); mean(sign(o))
}
pfmt <- function(p) ifelse(is.na(p), NA, signif(p, 3))

test_block <- function(mat, label) {
  out <- lapply(colnames(mat), function(p) {
    v <- mat[, p]
    if (length(unique(v)) < 2) {
      return(data.frame(pathway = p, n_NF = sum(donor_grp$group == "NF"),
                        n_disease = sum(donor_grp$group != "NF"),
                        median_NF = median(v[donor_grp$group == "NF"]),
                        median_disease = median(v[donor_grp$group != "NF"]),
                        log2_ratio = NA, delta = NA, p = NA, stringsAsFactors = FALSE))
    }
    x <- v[donor_grp$group == "NF"]; y <- v[donor_grp$group != "NF"]
    data.frame(pathway = p, n_NF = length(x), n_disease = length(y),
               median_NF = median(x), median_disease = median(y),
               log2_ratio = log2((median(y) + 1e-9) / (median(x) + 1e-9)),
               delta = cliff(x, y),
               p = suppressWarnings(wilcox.test(y, x, exact = FALSE)$p.value),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$fdr <- p.adjust(out$p, method = "BH")
  out$readout <- label
  out[, c("readout", "pathway", "n_NF", "n_disease", "median_NF", "median_disease",
          "log2_ratio", "delta", "p", "fdr")]
}

t_abs   <- test_block(score, "abs")
t_share <- test_block(share, "share")

# depth-adjusted linear model (disease indicator + log10 median depth)
grp2 <- as.integer(donor_grp$group != "NF")
dep  <- log10(donor_depth$median_nCount)
adj <- do.call(rbind, lapply(pw_lev, function(p) {
  v <- score[, p]
  if (length(unique(v)) < 2) return(data.frame(pathway = p, beta_disease = NA, p_adj_depth = NA))
  fit <- try(lm(v ~ grp2 + dep), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.frame(pathway = p, beta_disease = NA, p_adj_depth = NA))
  cf <- summary(fit)$coefficients
  data.frame(pathway = p, beta_disease = cf["grp2", "Estimate"],
             p_adj_depth = cf["grp2", "Pr(>|t|)"], stringsAsFactors = FALSE)
}))

tests <- merge(t_abs, t_share[, c("pathway", "log2_ratio", "delta", "p", "fdr")],
               by = "pathway", suffixes = c("", "_share"))
tests <- merge(tests, adj, by = "pathway", all.x = TRUE)
tests$fdr_adj_depth <- p.adjust(tests$p_adj_depth, method = "BH")
tests$neuro <- tests$pathway %in% NEURO
tests <- tests[order(tests$p), ]
write.csv(tests, file.path(tab_d, "99_pathway_tests.csv"), row.names = FALSE)

write.csv(data.frame(sample = rownames(score), donor_depth[, c("group", "median_nCount")],
                     round(score, 4), check.names = FALSE),
          file.path(tab_d, "99_donor_pathway_scores.csv"), row.names = FALSE)

lg("tests written | pathways with p<0.05 (abs): ", sum(tests$p < 0.05, na.rm = TRUE),
   " | FDR<0.05: ", sum(tests$fdr < 0.05, na.rm = TRUE))
lg("top 15 by abs p:")
print(head(tests[, c("pathway", "median_NF", "median_disease", "log2_ratio", "delta",
                     "p", "fdr", "p_share", "p_adj_depth", "neuro")], 15))

# pre-specified interest set (from the A4 plan; POSTN has no own pathway in
# CellChatDB v1 - it is annotated inside COLLAGEN / FN1, so both are kept)
FOCAL_SET <- c("NRG", "SPP1", "TGFb", "FN1", "COLLAGEN", "THBS", "LAMININ", "VEGF", "PDGF", "BMP")
focal_tab <- tests[tests$pathway %in% FOCAL_SET, ]
write.csv(focal_tab, file.path(tab_d, "99_focal_pathways.csv"), row.names = FALSE)
lg("pre-specified pathways:"); print(focal_tab[, c("pathway", "log2_ratio", "delta", "p", "fdr", "p_share", "p_adj_depth")])

# ------------------------------------------------- per celltype-pair detail --
# For the pre-specified pathways, test every sender -> receiver edge at donor
# level so the figure can show *who* talks to *whom*, not just a pathway total.
edge_tests <- do.call(rbind, lapply(FOCAL_SET, function(p) {
  sel <- which(pw == p)
  if (length(sel) == 0) return(NULL)
  em <- array(0, dim = c(length(CT), length(CT), length(samples)))
  for (d in seq_len(length(samples))) {
    em[, , d] <- apply(prob_donor[sel, , , d, drop = FALSE], c(2, 3), sum)
  }
  do.call(rbind, lapply(seq_along(CT), function(i) do.call(rbind, lapply(seq_along(CT), function(j) {
    v <- em[i, j, ]
    if (length(unique(v)) < 2) return(NULL)
    x <- v[donor_grp$group == "NF"]; y <- v[donor_grp$group != "NF"]
    data.frame(pathway = p, sender = CT[i], receiver = CT[j],
               median_NF = median(x), median_disease = median(y),
               log2_ratio = log2((median(y) + 1e-9) / (median(x) + 1e-9)),
               delta = cliff(x, y),
               p = suppressWarnings(wilcox.test(y, x, exact = FALSE)$p.value),
               stringsAsFactors = FALSE)
  }))))
}))
edge_tests$fdr <- p.adjust(edge_tests$p, method = "BH")
edge_tests <- edge_tests[order(edge_tests$p), ]
write.csv(edge_tests, file.path(tab_d, "99_focal_edge_tests.csv"), row.names = FALSE)
lg("edge tests: ", nrow(edge_tests), " | FDR<0.1: ", sum(edge_tests$fdr < 0.1, na.rm = TRUE))

saveRDS(list(score = score, share = share, donor_grp = donor_grp, depth = donor_depth,
             tests = tests, edge_tests = edge_tests, prob_donor = NULL,
             params = list(Kh = KH, floor = EXPR_FLOOR, min_cells = MIN_CELLS,
                           seed = SEED, neuro = NEURO, ct = CT, focal_set = FOCAL_SET)),
        file.path(obj_d, "cc_donor.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_99.txt"))
lg("=== done ===")
