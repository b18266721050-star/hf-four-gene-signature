# =============================================================================
# 98_sc_cellchat_network.R
# Stage 6 (A4) - group-level communication networks (CellChat mass-action model).
#
# Input : cellchat/objects/cc_input.rds (97)
# Output: cellchat/objects/cc_networks.rds
#         cellchat/tables/98_{pathway_flow,network_by_celltype,group_summary,top_lr_pairs}.csv
#         cellchat/logs/98_network.log
#
# Model (CellChat v1 formula):
#   x_ij = mean_expr(ligand | celltype i, group g) * mean_expr(receptor | celltype j, group g)
#   P_ij = x_ij / (x_ij + Kh),  Kh = 0.5
#   complex subunits are combined with min()
#   pairs whose mean expression < EXPR_FLOOR on either side are set to 0
#
# Group-level networks are descriptive (all nuclei of a group are pooled).
# Every p-value in this project comes from the donor-level step (99), never
# from pooling nuclei, which would treat 150k nuclei as independent replicates.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix) })

SEED <- 20251029
set.seed(SEED)

proj  <- "E:/another data file/GEObulk analysis project"
out_d <- file.path(proj, "results", "Singlecel", "cellchat")
tab_d <- file.path(out_d, "tables"); obj_d <- file.path(out_d, "objects")
log_d <- file.path(out_d, "logs")
lgf <- file.path(log_d, "98_network.log")
if (file.exists(lgf)) invisible(file.remove(lgf))
lg <- function(...) {
  m <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = lgf, append = TRUE)
}

KH <- 0.5
EXPR_FLOOR <- 0.05
GRP <- c("NF", "NICM", "DoxCM")
lg("=== 98 group-level communication networks | Kh=", KH, " floor=", EXPR_FLOOR, " ===")

cc <- readRDS(file.path(obj_d, "cc_input.rds"))
cnt <- cc$counts; meta <- cc$meta; pairs <- cc$pairs; CT <- cc$ct_order
ct_lv <- CT[CT %in% unique(as.character(meta$celltype))]
lg("input: ", nrow(cnt), " genes x ", ncol(cnt), " nuclei | ", nrow(pairs), " LR records")

# ---------------------------------------------------------- normalisation ----
lib <- meta$nCount_RNA / 1e4
dat <- cnt %*% Diagonal(x = 1 / lib)
dat@x <- log1p(dat@x)
lg("normalised: log1p(count/nCount*1e4) | max ", round(max(dat@x), 2))

# --------------------------------------------- mean expression per celltype --
S <- sparseMatrix(i = seq_len(nrow(meta)), j = as.integer(factor(meta$celltype, levels = ct_lv)),
                  x = 1, dims = c(nrow(meta), length(ct_lv)))
colnames(S) <- ct_lv
n_cells <- as.numeric(colSums(S))
mean_ct <- as.matrix(dat %*% S) / matrix(n_cells, nrow = nrow(dat), ncol = length(ct_lv), byrow = TRUE)
rownames(mean_ct) <- rownames(dat); colnames(mean_ct) <- ct_lv
lg("mean expression matrix: ", nrow(mean_ct), " genes x ", ncol(mean_ct), " celltypes | range ",
   paste(round(range(mean_ct), 2), collapse = "-"))

write.csv(data.frame(gene = rownames(mean_ct), round(mean_ct, 4), check.names = FALSE),
          file.path(tab_d, "98_mean_expr_by_celltype.csv"), row.names = FALSE)

# --------------------------------------------------------- pair evaluation --
# subunits -> vector over celltypes (min across subunits; NA if any missing)
subunit_expr <- function(subs, mat) {
  subs <- subs[subs %in% rownames(mat)]
  if (length(subs) == 0) return(NULL)
  if (length(subs) == 1) return(mat[subs, ])
  apply(mat[subs, , drop = FALSE], 2, min)
}

eval_pair <- function(i, mat) {
  L <- subunit_expr(pairs$ligand[[i]], mat)
  R <- subunit_expr(pairs$receptor[[i]], mat)
  if (is.null(L) || is.null(R)) return(NULL)
  x <- outer(L, R)                      # rows = sender, cols = receiver
  p <- x / (x + KH)
  p[outer(L, R, function(a, b) a < EXPR_FLOOR | b < EXPR_FLOOR)] <- 0
  p
}

build_prob <- function(mat) {
  out  <- lapply(seq_len(nrow(pairs)), function(i) eval_pair(i, mat))
  keep <- which(!vapply(out, is.null, logical(1)))
  a <- array(0, dim = c(length(keep), length(ct_lv), length(ct_lv)))
  for (m in seq_along(keep)) a[m, , ] <- out[[keep[m]]]
  dimnames(a) <- list(pairs$interaction_name[keep], ct_lv, ct_lv)
  list(prob = a, keep_ids = keep)
}

compute_group <- function(g) {
  idx <- which(as.character(meta$group) == g)
  if (length(idx) == 0) return(NULL)
  Sg <- S[idx, , drop = FALSE]
  ng <- as.numeric(colSums(Sg))
  mg <- as.matrix(dat[, idx, drop = FALSE] %*% Sg) /
        matrix(ng, nrow = nrow(dat), ncol = length(ct_lv), byrow = TRUE)
  rownames(mg) <- rownames(dat); colnames(mg) <- ct_lv
  b <- build_prob(mg)
  list(prob = b$prob, keep = rep(FALSE, nrow(pairs)) | seq_len(nrow(pairs)) %in% b$keep_ids,
       n = sum(ng), expr = mg)
}

res <- list()
res[["All"]] <- local({
  b <- build_prob(mean_ct)
  list(prob = b$prob, keep = seq_len(nrow(pairs)) %in% b$keep_ids, n = nrow(meta), expr = mean_ct)
})
for (g in GRP) res[[g]] <- compute_group(g)
lg("networks built for: ", paste(names(res), collapse = ","), " | pairs kept ",
   paste(sprintf("%s=%d", names(res), vapply(res, function(x) dim(x$prob)[1], integer(1))), collapse = " "))

# --------------------------------------------------------- aggregation ------
common <- Reduce(intersect, lapply(res, function(x) dimnames(x$prob)[[1]]))
lg("LR records present in all groups: ", length(common), " / ", nrow(pairs))
pw_of <- setNames(pairs$pathway, pairs$interaction_name)

pw_sum <- function(p) {
  # sum probability over all (sender, receiver) pairs for every pathway
  s <- apply(p, 1, sum)                              # per LR record total
  tapply(s[common], pw_of[common], sum, na.rm = TRUE)
}
flow <- sapply(res, function(x) pw_sum(x$prob[common, , , drop = FALSE]))
flow <- as.data.frame(flow)
flow$pathway <- rownames(flow)
flow <- flow[, c("pathway", names(res))]
flow <- flow[order(-flow[["All"]]), ]
write.csv(flow, file.path(tab_d, "98_pathway_flow.csv"), row.names = FALSE)

net_ct <- do.call(rbind, lapply(names(res), function(g) {
  p <- res[[g]]$prob[common, , , drop = FALSE]
  m <- apply(p, c(2, 3), sum)
  data.frame(group = g, sender = rep(rownames(m), times = ncol(m)),
             receiver = rep(colnames(m), each = nrow(m)),
             strength = as.vector(m), stringsAsFactors = FALSE)
}))
write.csv(net_ct, file.path(tab_d, "98_network_by_celltype.csv"), row.names = FALSE)

gsum <- do.call(rbind, lapply(names(res), function(g) {
  p <- res[[g]]$prob[common, , , drop = FALSE]
  data.frame(group = g, nuclei = res[[g]]$n,
             n_LR_tested = dim(p)[1],
             n_interactions_detected = sum(apply(p, 1, max) > 0),
             total_strength = sum(p),
             mean_strength_per_pair = sum(p) / dim(p)[1],
             n_pathways_active = sum(tapply(apply(p, 1, max), pw_of[common], max) > 0),
             stringsAsFactors = FALSE)
}))
write.csv(gsum, file.path(tab_d, "98_group_summary.csv"), row.names = FALSE)
lg("group summary:"); print(gsum)

# top pairs per group (for the bubble / table)
top_lr <- do.call(rbind, lapply(GRP, function(g) {
  p <- res[[g]]$prob[common, , , drop = FALSE]
  s <- apply(p, 1, max)
  o <- order(-s)[seq_len(min(40, length(s)))]
  data.frame(group = g, interaction = common[o], pathway = pw_of[common[o]],
             max_prob = round(s[o], 4), stringsAsFactors = FALSE)
}))
write.csv(top_lr, file.path(tab_d, "98_top_lr_pairs.csv"), row.names = FALSE)

saveRDS(list(prob = lapply(res, `[[`, "prob"), expr = lapply(res, `[[`, "expr"),
             common = common, flow = flow, ct = ct_lv, params = cc$params),
        file.path(obj_d, "cc_networks.rds"))
lg("saved cc_networks.rds: ", round(file.size(file.path(obj_d, "cc_networks.rds")) / 1e6, 1), " MB")

writeLines(capture.output(sessionInfo()), file.path(out_d, "sessionInfo_98.txt"))
lg("=== done ===")
