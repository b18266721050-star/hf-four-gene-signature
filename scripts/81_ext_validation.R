# =============================================================================
# 81_ext_validation.R —— 外部队列验证（GSE5406 芯片 / GSE116250 RNA-seq）
# 冻结训练集的 z-score 参数与 M2 系数，不重新拟合
#   主分析 V2：分位数映射到训练分布后代入冻结参数（跨平台可迁移）
#   敏感性 V3：外部队列内 z 标准化（检验判别力是否依赖尺度对齐）
#   诊断   V1：直接套用训练参数（不做任何尺度对齐，仅用于展示平台尺度差异）
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(pROC); library(ggplot2)
})

root <- "E:/another data file/GEObulk analysis project"
outd <- file.path(root, "results", "ExternalValidation")
tabd <- file.path(outd, "01_tables"); figd <- file.path(outd, "02_figures")
datd <- file.path(outd, "04_data");   logd <- file.path(outd, "03_logs")
for (d in c(tabd, figd, datd, logd)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

logcon <- file(file.path(logd, "81_validation_log.txt"), open = "wt")
say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, logcon) }

genes4 <- c("AQP3", "GGT5", "ISLR", "LAD1")

# ---------------------------------------------------------------- 冻结参数
zpar <- fread(file.path(root, "results/Nomogram/01_tables/T06_zscore_parameters_train_only.csv"))
coef_tab <- fread(file.path(root, "results/Nomogram/01_tables/T14_coefficients_M2_gene.csv"))
beta <- setNames(coef_tab$beta, coef_tab$term)
b0 <- beta[["Intercept"]]; w <- beta[paste0(genes4, "_z")]
names(w) <- genes4          # 便于按基因名取子集（留一基因敏感性分析）
say("frozen intercept: ", b0)
say("frozen weights: ", paste(sprintf("%s=%.4f", names(w), w), collapse = ", "))
tr_mean <- setNames(zpar$train_mean, zpar$gene)[genes4]
tr_sd   <- setNames(zpar$train_sd,   zpar$gene)[genes4]
say("frozen train mean: ", paste(sprintf("%s=%.4f", genes4, tr_mean), collapse = ", "))
say("frozen train sd:   ", paste(sprintf("%s=%.4f", genes4, tr_sd), collapse = ", "))

tr4 <- as.matrix(read.csv(file.path(outd, "01_data", "train_fourgene_expr.csv"),
                          row.names = 1, check.names = FALSE))
say("train 4-gene matrix: ", nrow(tr4), " x ", ncol(tr4))

# ---------------------------------------------------------------- 训练集内部判别力（复算，作对照）
mp <- fread(file.path(root, "results/Nomogram/04_data/modelpredictions.csv"))
say("train predictions cols ok: ", all(c("lp_M2", "prob_M2_gene", "HF", "set") %in% names(mp)),
    " ; n=", nrow(mp))
auc_int <- function(d, lab) {
  r <- roc(d$HF, d$prob_M2_gene, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(cohort = lab, n = nrow(d), n_case = sum(d$HF), n_ctrl = sum(d$HF == 0),
             auc = as.numeric(auc(r)), low = ci[1], high = ci[3], stringsAsFactors = FALSE)
}
int_tab <- rbind(
  auc_int(mp[set == "train"], "GSE57338 training (internal, 219)"),
  auc_int(mp[set == "validation"], "GSE57338 internal random split (94)")
)
say("internal AUC recomputed:")
for (i in seq_len(nrow(int_tab))) say(sprintf("  %-42s n=%3d AUC=%.4f (%.4f-%.4f)",
    int_tab$cohort[i], int_tab$n[i], int_tab$auc[i], int_tab$low[i], int_tab$high[i]))

# ---------------------------------------------------------------- 工具函数
qp_map <- function(x, ref) {              # 经验分位数映射到参考分布
  ok <- !is.na(x); n <- sum(ok)
  q <- (rank(x[ok], ties.method = "average") - 0.5) / n
  out <- x
  out[ok] <- as.numeric(quantile(ref, probs = q, type = 7, names = FALSE))
  out
}
lp_from_z <- function(zmat) as.numeric(b0 + t(as.matrix(zmat)) %*% as.numeric(w))
to_prob <- function(lp) 1 / (1 + exp(-lp))
# 固定宽度分箱的校准曲线（预测概率高度集中时，十分位分箱会抖到 0/1 两端）
cal_bins <- function(p, y, width = 0.1, min_n = 5) {
  brk <- seq(0, 1, by = width)
  b <- cut(p, breaks = brk, include.lowest = TRUE, labels = FALSE)
  d <- data.frame(b = b, obs = y, pred = p)
  ag <- aggregate(cbind(obs, pred) ~ b, data = d, FUN = mean)
  ag$n <- as.numeric(table(b)[as.character(ag$b)])
  ag <- ag[!is.na(ag$b) & ag$n >= min_n, ]
  if (nrow(ag) < 2) {                       # 退化时改用三分位
    q <- unique(quantile(p, c(0, 1/3, 2/3, 1)))
    if (length(q) >= 3) {
      b2 <- cut(p, breaks = q, include.lowest = TRUE, labels = FALSE)
      d2 <- data.frame(b = b2, obs = y, pred = p)
      ag <- aggregate(cbind(obs, pred) ~ b, data = d2, FUN = mean)
      ag$n <- as.numeric(table(b2)[as.character(ag$b)])
      ag <- ag[!is.na(ag$b), ]
    }
  }
  ag
}

# 三种口径的 z 矩阵
z_v1 <- function(x) sweep(sweep(x, 1, tr_mean, "-"), 1, tr_sd, "/")
z_v2 <- function(x) {                      # 分位数映射后再用冻结参数
  xm <- t(sapply(genes4, function(g) qp_map(x[g, ], tr4[g, ])))
  rownames(xm) <- genes4; colnames(xm) <- colnames(x)
  sweep(sweep(xm, 1, tr_mean, "-"), 1, tr_sd, "/")
}
z_v3 <- function(x) {                      # 队列内标准化
  m <- rowMeans(x, na.rm = TRUE); s <- apply(x, 1, sd, na.rm = TRUE)
  sweep(sweep(x, 1, m, "-"), 1, s, "/")
}

# ---------------------------------------------------------------- 队列装载
coh <- list()
# GSE5406
e1 <- as.matrix(read.csv(file.path(outd, "01_data", "GSE5406_gene_expr.csv"), row.names = 1, check.names = FALSE))
i1 <- fread(file.path(outd, "01_data", "GSE5406_sample_info.csv"))
coh[["GSE5406"]] <- list(expr = e1[genes4, ], info = i1, platform = "HG-U133A (RMA log2)",
                         type = "microarray", age_sex = FALSE)
# GSE116250
e2 <- as.matrix(read.csv(file.path(outd, "01_data", "GSE116250_gene_expr_log2rpkm.csv"), row.names = 1, check.names = FALSE))
i2 <- fread(file.path(outd, "01_data", "GSE116250_sample_info.csv"))
i2[, title := as.character(title)]
# 统一用 GSM 号做列名，便于与 sample_info 对齐
colnames(e2) <- i2$sample[match(colnames(e2), i2$title)]
coh[["GSE116250"]] <- list(expr = e2[genes4, ], info = i2, platform = "Illumina HiSeq2500 (log2 RPKM)",
                           type = "RNA-seq", age_sex = TRUE)

# GSE55296（第三队列，备用）
f5 <- file.path(outd, "01_data", "GSE55296_gene_expr_log2.csv")
if (file.exists(f5)) {
  e5 <- as.matrix(read.csv(f5, row.names = 1, check.names = FALSE))
  i5 <- fread(file.path(outd, "01_data", "GSE55296_sample_info.csv"))
  colnames(e5) <- i5$sample[match(colnames(e5), i5$title)]
  coh[["GSE55296"]] <- list(expr = e5[genes4, ], info = i5,
                            platform = "RNA-seq (normalized counts, log2)",
                            type = "RNA-seq", age_sex = FALSE)
}
bad <- which(coh[["GSE116250"]]$expr["LAD1", ] > 20)
if (length(bad)) {
  say("GSE116250 LAD1 placeholder detected in: ",
      paste(colnames(coh[["GSE116250"]]$expr)[bad], collapse = ", "), " -> set NA")
  coh[["GSE116250"]]$expr["LAD1", bad] <- NA
}
# 通用异常值巡检：log2 尺度下 >25 视为占位符/伪值
for (cn in names(coh)) {
  M <- coh[[cn]]$expr
  big <- which(M > 25 | M < -10, arr.ind = TRUE)
  if (nrow(big)) {
    say("WARN ", cn, " : ", nrow(big), " extreme cells -> set NA: ",
        paste(sprintf("%s/%s=%.1f", rownames(M)[big[, 1]], colnames(M)[big[, 2]],
                      M[big]), collapse = ", "))
    M[big] <- NA
    coh[[cn]]$expr <- M
  } else say(cn, " : no extreme value (|log2| sane)")
}
say("GSE5406 any |value|>20: ", any(abs(coh[["GSE5406"]]$expr) > 20, na.rm = TRUE))

# ---------------------------------------------------------------- 主循环
auc_rows <- list(); cal_rows <- list(); dir_rows <- list(); scale_rows <- list()
pred_store <- list()

for (cn in names(coh)) {
  X <- coh[[cn]]$expr
  info <- coh[[cn]]$info
  nm <- colnames(X)
  info <- info[match(nm, info$sample)]
  y <- ifelse(info$group %in% c("DCM", "ICM"), 1L, 0L)
  ok <- complete.cases(t(X))
  say("")
  say("============== ", cn, " (", coh[[cn]]$platform, ") ==============")
  say("samples: ", length(nm), " ; complete 4-gene: ", sum(ok),
      " ; cases: ", sum(y[ok]), " controls: ", sum(y[ok] == 0))
  say("group x complete:")
  for (g in c("NF", "DCM", "ICM")) say(sprintf("  %-4s n=%d complete=%d", g,
      sum(info$group == g), sum(info$group == g & ok)))

  # --- 尺度诊断
  for (g in genes4) scale_rows[[paste(cn, g)]] <- data.frame(
    cohort = cn, gene = g, mean_external = mean(X[g, ], na.rm = TRUE),
    sd_external = sd(X[g, ], na.rm = TRUE), mean_train = tr_mean[g], sd_train = tr_sd[g],
    sd_ratio = sd(X[g, ], na.rm = TRUE) / tr_sd[g])

  Z <- list(V1 = z_v1(X), V2 = z_v2(X), V3 = z_v3(X))
  for (v in names(Z)) {
    lp <- lp_from_z(Z[[v]]); p <- to_prob(lp)
    pred_store[[paste(cn, v)]] <- data.frame(
      cohort = cn, version = v, sample = nm, title = info$title, group = info$group,
      y = y, lp = lp, prob = p,
      z_AQP3 = Z[[v]]["AQP3", ], z_GGT5 = Z[[v]]["GGT5", ],
      z_ISLR = Z[[v]]["ISLR", ], z_LAD1 = Z[[v]]["LAD1", ])
  }
  P <- pred_store[[paste(cn, "V2")]]

  # --- 判别力（主口径 V2 + 敏感性 V1/V3）
  sub_auc <- function(d, lab, sel) {
    dd <- d[sel & ok, ]
    if (length(unique(dd$y)) < 2) return(NULL)
    r <- roc(dd$y, dd$prob, quiet = TRUE, direction = "<")
    ci <- as.numeric(ci.auc(r, method = "delong"))
    data.frame(cohort = cn, subset = lab, version = unique(d$version), n = nrow(dd),
               n_case = sum(dd$y), n_ctrl = sum(dd$y == 0),
               auc = as.numeric(auc(r)), low = ci[1], high = ci[3], stringsAsFactors = FALSE)
  }
  for (v in names(pred_store)) {
    if (!grepl(paste0("^", cn, " "), v)) next
    d <- pred_store[[v]]
    if (!identical(d$version[1], c(V1 = "V1", V2 = "V2", V3 = "V3")[[d$version[1]]])) next
    auc_rows[[paste(cn, d$version[1])]] <- rbind(
      sub_auc(d, "all disease vs NF", rep(TRUE, nrow(d))),
      sub_auc(d, "DCM vs NF", d$group %in% c("DCM", "NF")),
      sub_auc(d, "ICM vs NF", d$group %in% c("ICM", "NF")))
  }

  # --- 单基因 AUC（V2 口径，仅用于说明哪个基因在扛）
  for (g in genes4) {
    dd <- P[!is.na(P[[paste0("z_", g)]]) & ok, ]
    if (length(unique(dd$y)) < 2) next
    r <- roc(dd$y, dd[[paste0("z_", g)]], quiet = TRUE)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    auc_rows[[paste(cn, "single", g)]] <- data.frame(
      cohort = cn, subset = paste0("single gene ", g), version = "V2",
      n = nrow(dd), n_case = sum(dd$y), n_ctrl = sum(dd$y == 0),
      auc = as.numeric(auc(r)), low = ci[1], high = ci[3], stringsAsFactors = FALSE)
  }

  # --- 留一基因敏感性：系数冻结，只去掉某一项
  mods <- list(
    "score: all 4 genes"            = genes4,
    "score: 3 genes (drop LAD1)"    = c("AQP3", "GGT5", "ISLR"),
    "score: 3 genes (drop ISLR)"    = c("AQP3", "GGT5", "LAD1"),
    "score: 2 genes (AQP3+GGT5)"    = c("AQP3", "GGT5"))
  for (mn in names(mods)) {
    mm <- mods[[mn]]
    z <- Z[["V2"]][mm, , drop = FALSE]
    lp <- as.numeric(b0 + t(z) %*% as.numeric(w[mm]))
    pr <- to_prob(lp)
    okk <- complete.cases(t(z))
    n_ctrl2 <- sum(okk & y == 0); n_case2 <- sum(okk & y == 1)
    say(sprintf("  [%s] %-28s n=%d case=%d ctrl=%d", cn, mn, sum(okk), n_case2, n_ctrl2))
    if (n_ctrl2 < 2 || n_case2 < 2) { say("    -> skipped (too few in one arm)"); next }
    r <- roc(y[okk], pr[okk], quiet = TRUE, direction = "<")
    ci <- as.numeric(ci.auc(r, method = "delong"))
    auc_rows[[paste(cn, mn)]] <- data.frame(
      cohort = cn, subset = mn, version = "V2", n = sum(okk),
      n_case = sum(y[okk]), n_ctrl = sum(y[okk] == 0),
      auc = as.numeric(auc(r)), low = ci[1], high = ci[3], stringsAsFactors = FALSE)
  }

  # --- 校准与整体表现（V2）
  dd <- P[ok, ]
  p <- dd$prob; yy <- dd$y
  brier <- mean((p - yy)^2)
  lp_obs <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  citl <- tryCatch(suppressWarnings(coef(glm(yy ~ offset(lp_obs), family = binomial))[1]),
                   error = function(e) NA)
  slope <- tryCatch(suppressWarnings(coef(glm(yy ~ lp_obs, family = binomial))[2]),
                    error = function(e) NA)
  if (!is.na(slope) && abs(slope) > 20) slope <- NA_real_   # 完全分离时斜率不可估
  # 校准分箱（固定 0.1 宽，n>=5）
  cald <- cal_bins(p, yy)
  names(cald)[names(cald) == "b"] <- "dec"
  cald$cohort <- cn
  cal_rows[[cn]] <- data.frame(cohort = cn, n = nrow(dd), auc = as.numeric(auc(roc(yy, p, quiet = TRUE))),
                               brier = brier, cal_intercept = as.numeric(citl),
                               cal_slope = as.numeric(slope),
                               mean_pred = mean(p), obs_rate = mean(yy))
  write.csv(cald, file.path(datd, paste0("calibration_deciles_", cn, ".csv")), row.names = FALSE)

  # --- 方向一致性（limma ~ group）
  grp <- factor(info$group, levels = c("NF", "DCM", "ICM"))
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(X, design)
  ctr <- makeContrasts(DCM - NF, ICM - NF,
                       levels = design)
  fit2 <- eBayes(contrasts.fit(fit, ctr))
  for (g in genes4) {
    for (k in 1:2) {
      tt <- topTable(fit2, coef = k, number = Inf)
      a <- X[g, grp == c("DCM", "ICM")[k]]; b <- X[g, grp == "NF"]
      wn <- suppressWarnings(wilcox.test(a, b))
      sp <- sqrt(((length(a) - 1) * var(a, na.rm = TRUE) +
                    (length(b) - 1) * var(b, na.rm = TRUE)) / (length(a) + length(b) - 2))
      dir_rows[[paste(cn, g, k)]] <- data.frame(
        cohort = cn, gene = g, contrast = colnames(ctr)[k],
        log2FC = tt[g, "logFC"], t = tt[g, "t"], p = tt[g, "P.Value"],
        cohens_d = (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sp,
        wilcox_p = wn$p.value,
        n_case = sum(grp == c("DCM", "ICM")[k]), n_ctrl = sum(grp == "NF"),
        stringsAsFactors = FALSE)
    }
  }
  # 阳性对照 NPPB：若队列无 NPPB 会跳过
  if ("NPPB" %in% rownames(X)) {
    tt <- topTable(fit2, coef = 1, number = Inf)
    say("positive control NPPB DCM vs NF: log2FC=", round(tt["NPPB", "logFC"], 3),
        " p=", signif(tt["NPPB", "P.Value"], 3))
  }
  # 队列内四基因组间中位数（描述）
  say("--- median by group (log2 scale) ---")
  for (g in genes4) say(sprintf("  %-5s NF=%6.2f DCM=%6.2f ICM=%6.2f", g,
      median(X[g, grp == "NF"], na.rm = TRUE), median(X[g, grp == "DCM"], na.rm = TRUE),
      median(X[g, grp == "ICM"], na.rm = TRUE)))
}

auc_tab <- rbindlist(auc_rows)
cal_tab <- rbindlist(cal_rows)
dir_tab <- rbindlist(dir_rows)
dir_tab[, fdr := p.adjust(p, method = "BH"), by = .(cohort, contrast)]
scale_tab <- rbindlist(scale_rows)

# ---------------------------------------------------------------- DeLong：外部 vs 内部验证
rv <- mp[set == "validation"]
r_int <- roc(rv$HF, rv$prob_M2_gene, quiet = TRUE, direction = "<")
dl <- list()
for (cn in names(coh)) {
  d <- pred_store[[paste(cn, "V2")]]
  okk <- !is.na(d$prob) & !is.na(d$z_LAD1)
  r_ext <- roc(d$y[okk], d$prob[okk], quiet = TRUE, direction = "<")
  tt <- roc.test(r_int, r_ext, method = "delong")
  dl[[cn]] <- data.frame(
    comparison = paste0(cn, " external vs GSE57338 internal split"),
    n_internal = length(r_int$response), n_external = sum(okk),
    auc_internal = as.numeric(auc(r_int)), auc_external = as.numeric(auc(r_ext)),
    delta = as.numeric(auc(r_int)) - as.numeric(auc(r_ext)),
    p_delong = tt$p.value, stringsAsFactors = FALSE)
}
dl_tab <- rbindlist(dl)
write.csv(dl_tab, file.path(tabd, "T07_delong_internal_vs_external.csv"), row.names = FALSE)
say("DeLong internal vs external:")
for (i in seq_len(nrow(dl_tab)))
  say(sprintf("  %-52s AUC %.3f vs %.3f, p=%.3g", dl_tab$comparison[i],
              dl_tab$auc_internal[i], dl_tab$auc_external[i], dl_tab$p_delong[i]))

# ---------------------------------------------------------------- 分层 bootstrap AUC CI（主口径 V2）
set.seed(20251029)
boot_ci <- function(y, p, B = 2000) {
  idx1 <- which(y == 1); idx0 <- which(y == 0)
  out <- numeric(B)
  for (b in seq_len(B)) {
    s <- c(sample(idx1, replace = TRUE), sample(idx0, replace = TRUE))
    out[b] <- as.numeric(auc(roc(y[s], p[s], quiet = TRUE, direction = "<")))
  }
  c(auc = as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<"))),
    low = as.numeric(quantile(out, 0.025, names = FALSE)),
    high = as.numeric(quantile(out, 0.975, names = FALSE)))
}
boot_rows <- list()
for (cn in names(coh)) {
  d <- pred_store[[paste(cn, "V2")]]
  okk <- !is.na(d$prob) & !is.na(d$z_LAD1)
  bc <- boot_ci(d$y[okk], d$prob[okk])
  boot_rows[[cn]] <- data.frame(cohort = cn, n = sum(okk), auc = bc[1],
                                boot_low = bc[2], boot_high = bc[3], B = 2000)
}
boot_tab <- rbindlist(boot_rows)
write.csv(boot_tab, file.path(tabd, "T09_bootstrap_auc.csv"), row.names = FALSE)
say("bootstrap AUC (2000 stratified resamples):")
for (i in seq_len(nrow(boot_tab))) say(sprintf("  %-10s %.3f (%.3f-%.3f)",
    boot_tab$cohort[i], boot_tab$auc[i], boot_tab$boot_low[i], boot_tab$boot_high[i]))

# ---------------------------------------------------------------- 外部队列 DCA
dca_all <- list()
for (cn in names(coh)) {
  d <- pred_store[[paste(cn, "V2")]]
  okk <- !is.na(d$prob) & !is.na(d$z_LAD1)
  yv <- d$y[okk]; pv <- d$prob[okk]; nv <- length(yv); prev <- mean(yv)
  th <- seq(0.01, 0.99, by = 0.01)
  nb_model <- sapply(th, function(t) {
    tp <- sum(pv >= t & yv == 1); fp <- sum(pv >= t & yv == 0)
    tp / nv - fp / nv * (t / (1 - t))
  })
  nb_all <- prev - (1 - prev) * (th / (1 - th))
  # 敏感性：队列内逻辑再校准（乐观，仅用于判断"排序信息是否可用"）
  rc <- tryCatch(suppressWarnings(glm(yv ~ qlogis(pmin(pmax(pv, 1e-6), 1 - 1e-6)),
                                      family = binomial)), error = function(e) NULL)
  p_re <- if (is.null(rc)) pv else predict(rc, type = "response")
  nb_re <- sapply(th, function(t) {
    tp <- sum(p_re >= t & yv == 1); fp <- sum(p_re >= t & yv == 0)
    tp / nv - fp / nv * (t / (1 - t))
  })
  dca_all[[cn]] <- rbind(
    data.frame(cohort = cn, threshold = th, strategy = "Frozen 4-gene model",
               net_benefit = nb_model),
    data.frame(cohort = cn, threshold = th, strategy = "Frozen + in-cohort recalibration (optimistic)",
               net_benefit = nb_re),
    data.frame(cohort = cn, threshold = th, strategy = "Treat all", net_benefit = nb_all),
    data.frame(cohort = cn, threshold = th, strategy = "Treat none", net_benefit = 0))
  write.csv(dca_all[[cn]], file.path(datd, paste0("dca_", cn, ".csv")), row.names = FALSE)
}
dca_tab <- rbindlist(dca_all)
say("DCA computed for: ", paste(names(dca_all), collapse = ", "))

# ---------------------------------------------------------------- 队列基线表
base_rows <- list()
for (cn in names(coh)) {
  info <- coh[[cn]]$info
  ag <- if ("age" %in% names(info)) suppressWarnings(as.numeric(info$age)) else rep(NA, nrow(info))
  sx <- if ("sex" %in% names(info)) info$sex else rep(NA_character_, nrow(info))
  base_rows[[cn]] <- data.frame(
    cohort = cn, platform = coh[[cn]]$platform, n = nrow(info),
    NF = sum(info$group == "NF"), DCM = sum(info$group == "DCM"),
    ICM = sum(info$group == "ICM"),
    age_median = if (all(is.na(ag))) NA else round(median(ag, na.rm = TRUE), 1),
    age_range = if (all(is.na(ag))) "n/a" else paste(range(ag, na.rm = TRUE), collapse = "-"),
    n_male = if (all(is.na(sx))) NA else sum(grepl("male", sx, ignore.case = TRUE) &
                                              !grepl("female", sx, ignore.case = TRUE)),
    n_female = if (all(is.na(sx))) NA else sum(grepl("female", sx, ignore.case = TRUE)),
    note = c(GSE5406 = "no age/sex in GEO record",
             GSE116250 = "age and sex available",
             GSE55296 = "no age/sex in GEO record")[[cn]],
    stringsAsFactors = FALSE)
}
base_tab <- rbindlist(base_rows)
write.csv(base_tab, file.path(tabd, "T10_cohort_baseline.csv"), row.names = FALSE)
say("baseline table:")
for (i in seq_len(nrow(base_tab))) say(sprintf("  %-10s n=%3d NF=%2d DCM=%2d ICM=%3d age=%s",
    base_tab$cohort[i], base_tab$n[i], base_tab$NF[i], base_tab$DCM[i], base_tab$ICM[i],
    base_tab$age_range[i]))

# ---------------------------------------------------------------- 发现队列的校准分箱（供同图对比）
dv <- mp[set == "validation"]
cald0 <- cal_bins(dv$prob_M2_gene, dv$HF)
names(cald0)[names(cald0) == "b"] <- "dec"
cald0$cohort <- "GSE57338 internal split"
write.csv(cald0, file.path(datd, "calibration_deciles_GSE57338_internal.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 方向一致性汇总
train_fc <- data.frame(gene = genes4,
                       train_log2FC = c(-0.991, -0.811, 0.821, -0.667),
                       stringsAsFactors = FALSE)   # 来自 results/DEG（DCM+ICM vs NF）
wide <- dcast(dir_tab, gene ~ paste(cohort, contrast, sep = "."), value.var = "log2FC")
dir_sum <- merge(train_fc, wide, by = "gene", all.x = TRUE)
val_cols <- grep("^GSE", names(dir_sum), value = TRUE)
dir_sum$sign_train <- sign(dir_sum$train_log2FC)
for (cc in val_cols) dir_sum[[paste0("sign_", cc)]] <- sign(dir_sum[[cc]])
sign_cols <- paste0("sign_", val_cols)
dir_sum$n_consistent <- rowSums(
  sapply(sign_cols, function(k) dir_sum[[k]] == dir_sum$sign_train), na.rm = TRUE)
dir_sum$n_tested <- rowSums(!is.na(dir_sum[, val_cols, drop = FALSE]))

# ---------------------------------------------------------------- 汇总表
write.csv(int_tab, file.path(tabd, "T01_internal_auc_recomputed.csv"), row.names = FALSE)
write.csv(auc_tab, file.path(tabd, "T02_external_auc.csv"), row.names = FALSE)
write.csv(dir_tab, file.path(tabd, "T03_direction_external.csv"), row.names = FALSE)
write.csv(dir_sum, file.path(tabd, "T04_direction_summary.csv"), row.names = FALSE)
write.csv(cal_tab, file.path(tabd, "T05_calibration_brier.csv"), row.names = FALSE)
write.csv(scale_tab, file.path(tabd, "T06_scale_diagnostics.csv"), row.names = FALSE)
allpred <- rbindlist(pred_store)
fwrite(allpred, file.path(datd, "external_predictions_all_versions.csv"))

say("")
say("================ AUC (V2 primary) ================")
print(auc_tab[version == "V2"])
say("================ direction summary ================")
print(dir_sum)
say("================ calibration ================")
print(cal_tab)

# ---------------------------------------------------------------- 图
pal <- c("GSE57338 internal split" = "#4E79A7", "GSE5406" = "#E15759",
         "GSE116250" = "#59A14F", "GSE55296" = "#B07AA1")
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

# F01 ROC 叠加
rocdf <- list()
rv <- mp[set == "validation"]
r0 <- roc(rv$HF, rv$prob_M2_gene, quiet = TRUE, direction = "<")
rocdf[[1]] <- data.frame(fpr = 1 - r0$specificities, tpr = r0$sensitivities,
                         cohort = "GSE57338 internal split")
for (cn in names(coh)) {
  d <- pred_store[[paste(cn, "V2")]]
  okc <- !is.na(d$prob) & !is.na(d$z_LAD1)
  rr <- roc(d$y[okc], d$prob[okc], quiet = TRUE, direction = "<")
  rocdf[[cn]] <- data.frame(fpr = 1 - rr$specificities, tpr = rr$sensitivities, cohort = cn)
}
rocdf <- rbindlist(rocdf)
lab <- auc_tab[version == "V2" & subset == "all disease vs NF"]
lab <- rbind(data.frame(cohort = "GSE57338 internal split",
                        auc = as.numeric(auc(r0))), lab[, .(cohort, auc)])
lab$txt <- sprintf("%s: AUC %.3f", lab$cohort, lab$auc)
p1 <- ggplot(rocdf, aes(fpr, tpr, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_step(linewidth = 0.9, direction = "hv") +
  scale_colour_manual(values = c(pal, "GSE57338 internal split" = pal[[1]]), name = NULL,
                      labels = setNames(lab$txt, lab$cohort)) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity") +
  annotate("text", x = 0.62, y = 0.08,
           label = "frozen 4-gene model (M2)\nno refitting in external cohorts",
           size = 3, colour = "grey30") +
  th + theme(legend.position = c(0.62, 0.32), legend.background = element_blank())
ggsave(file.path(figd, "F01_ROC_internal_vs_external.pdf"), p1, width = 5.2, height = 5)
ggsave(file.path(figd, "F01_ROC_internal_vs_external.png"), p1, width = 5.2, height = 5, dpi = 600)

# F02 校准曲线（含发现队列作参照）
calall <- rbindlist(c(
  list(fread(file.path(datd, "calibration_deciles_GSE57338_internal.csv"))),
  lapply(names(coh), function(cn)
    fread(file.path(datd, paste0("calibration_deciles_", cn, ".csv"))))))
p2 <- ggplot(calall, aes(pred, obs, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(aes(size = n), alpha = 0.9) +
  scale_size_continuous(range = c(1.4, 4.2), name = "n per bin") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  coord_equal() +
  labs(x = "Predicted probability (frozen model)", y = "Observed proportion",
       subtitle = "0.1-wide probability bins with n \u2265 5; dashed line = perfect calibration") +
  th + theme(legend.position = "bottom", legend.box = "horizontal",
             legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
             plot.subtitle = element_text(size = 8)) +
  guides(size = guide_legend(order = 1), colour = guide_legend(order = 2, ncol = 2))
ggsave(file.path(figd, "F02_calibration_external.pdf"), p2, width = 5.4, height = 5.4)
ggsave(file.path(figd, "F02_calibration_external.png"), p2, width = 5.4, height = 5.4, dpi = 600)

# F03 四基因箱线图
bx <- list()
for (cn in names(coh)) {
  X <- coh[[cn]]$expr; info <- coh[[cn]]$info
  info <- info[match(colnames(X), info$sample)]
  for (g in genes4) bx[[paste(cn, g)]] <- data.frame(
    cohort = cn, gene = g, group = factor(info$group, levels = c("NF", "DCM", "ICM")),
    value = as.numeric(X[g, ]))
}
bx <- rbindlist(bx)
bx <- bx[!is.na(bx$value), ]
p3 <- ggplot(bx, aes(group, value, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.65) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.55) +
  facet_grid(gene ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c(NF = "#8CD17D", DCM = "#F1CE63", ICM = "#E15759"), guide = "none") +
  labs(x = NULL, y = expression(log[2]~"expression")) + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(figd, "F03_fourgene_boxplot_external.pdf"), p3, width = 7.2, height = 6.4)
ggsave(file.path(figd, "F03_fourgene_boxplot_external.png"), p3, width = 7.2, height = 6.4, dpi = 600)

# F04 AUC 森林图（只放主模型：全病例 / DCM / ICM）
main_sub <- c("all disease vs NF", "DCM vs NF", "ICM vs NF")
fp <- rbind(
  data.frame(cohort = "GSE57338 internal split", subset = "all disease vs NF",
             auc = as.numeric(auc(r0)),
             low = as.numeric(ci.auc(r0, method = "delong"))[1],
             high = as.numeric(ci.auc(r0, method = "delong"))[3]),
  as.data.frame(auc_tab[version == "V2" & subset %in% main_sub, .(cohort, subset, auc, low, high)]))
fp$cohort <- factor(fp$cohort, levels = c("GSE57338 internal split", names(coh)))
fp$subset <- factor(fp$subset, levels = main_sub)
fp$lab <- paste0(fp$cohort, " | ", fp$subset)
fp <- fp[order(fp$cohort, fp$subset), ]
fp$lab <- factor(fp$lab, levels = rev(unique(fp$lab)))
p4 <- ggplot(fp, aes(auc, lab, colour = cohort)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.16, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(x = high + 0.02, label = sprintf("%.3f", auc)), hjust = 0, size = 3,
            show.legend = FALSE) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(0.40, 1.20), breaks = seq(0.4, 1.0, 0.2)) +
  labs(x = "AUC (95% CI, DeLong)", y = NULL) + th
ggsave(file.path(figd, "F04_AUC_forest.pdf"), p4, width = 7.2, height = 3.4)
ggsave(file.path(figd, "F04_AUC_forest.png"), p4, width = 7.2, height = 3.4, dpi = 600)

# F05 方向一致性热图（discovery 作参照列，队列顺序固定）
hm <- as.data.frame(copy(dir_tab))
hm$cell <- paste0(hm$cohort, "\n", hm$contrast)
hm$star <- ifelse(hm$p < 0.001, "***", ifelse(hm$p < 0.01, "**", ifelse(hm$p < 0.05, "*", "")))
hm$gene <- factor(hm$gene, levels = rev(genes4))
ref <- data.frame(cell = "GSE57338\ndiscovery", gene = factor(genes4, levels = rev(genes4)),
                  log2FC = c(-0.991, -0.811, 0.821, -0.667), star = "ref",
                  stringsAsFactors = FALSE)
hm2 <- rbind(hm[, c("cell", "gene", "log2FC", "star")], ref[, c("cell", "gene", "log2FC", "star")])
ord <- c("GSE57338\ndiscovery", "GSE5406\nDCM - NF", "GSE5406\nICM - NF",
         "GSE116250\nDCM - NF", "GSE116250\nICM - NF", "GSE55296\nDCM - NF", "GSE55296\nICM - NF")
hm2$cell <- factor(hm2$cell, levels = ord)
p5 <- ggplot(hm2, aes(cell, gene, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(star == "ref", sprintf("%.2f", log2FC),
                               paste0(sprintf("%.2f", log2FC), star))), size = 3.0) +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", midpoint = 0,
                       name = expression(log[2]~FC)) +
  labs(x = NULL, y = NULL,
       title = "Direction of the four genes (cases vs non-failing)",
       subtitle = "leftmost column = discovery cohort; * p<0.05, ** p<0.01, *** p<0.001 (limma)") +
  th + theme(plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 8.5),
             axis.text.x = element_text(size = 8.2))
ggsave(file.path(figd, "F05_direction_heatmap.pdf"), p5, width = 8, height = 3.6)
ggsave(file.path(figd, "F05_direction_heatmap.png"), p5, width = 8, height = 3.6, dpi = 600)

# F07 单基因与留一基因 AUC 热图（补充材料）
keep_sub <- c("single gene AQP3", "single gene GGT5", "single gene ISLR", "single gene LAD1",
              "score: 2 genes (AQP3+GGT5)", "score: 3 genes (drop LAD1)",
              "score: 3 genes (drop ISLR)", "score: all 4 genes")
tl <- as.data.frame(auc_tab[version == "V2" & subset %in% keep_sub,
                            .(cohort, subset, auc)])
tl$subset <- factor(tl$subset, levels = keep_sub)
tl$cohort <- factor(tl$cohort, levels = names(coh))
p7 <- ggplot(tl, aes(subset, cohort, fill = auc)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", auc)), size = 3.1) +
  scale_fill_gradient(low = "#f7f7f7", high = "#E15759", limits = c(0.45, 1),
                      name = "AUC") +
  labs(x = NULL, y = NULL, title = "Single-gene and leave-one-gene-out discrimination",
       subtitle = "frozen coefficients; only the terms shown are retained") +
  th + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8.3),
             plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 8.5))
ggsave(file.path(figd, "F07_single_and_leaveoneout_AUC.pdf"), p7, width = 7.4, height = 3.6)
ggsave(file.path(figd, "F07_single_and_leaveoneout_AUC.png"), p7, width = 7.4, height = 3.6, dpi = 600)

# F06 尺度诊断
sc <- as.data.frame(copy(scale_tab)); sc$gene <- factor(sc$gene, levels = genes4)
p6 <- ggplot(sc, aes(mean_train, mean_external, label = cohort, colour = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 2.6) +
  geom_text(aes(label = paste0(cohort, "\nSD ratio ", sprintf("%.2f", sd_ratio))),
            size = 2.7, vjust = -0.7, show.legend = FALSE) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = "Training cohort mean (log2)", y = "External cohort mean (log2)") + th
ggsave(file.path(figd, "F06_scale_diagnostics.pdf"), p6, width = 5.4, height = 4.2)
ggsave(file.path(figd, "F06_scale_diagnostics.png"), p6, width = 5.4, height = 4.2, dpi = 600)

# F08 决策曲线（外部队列，冻结模型 vs 全治 / 不治）
dca_tab[, strategy := factor(strategy, levels = c("Frozen 4-gene model", "Treat all", "Treat none"))]
dca_tab[, cohort := factor(cohort, levels = names(coh))]
dca_fig <- dca_tab[!is.na(strategy)]
p8 <- ggplot(dca_fig, aes(threshold, net_benefit, colour = strategy, linetype = strategy)) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = c("Frozen 4-gene model" = "#E15759",
                                 "Treat all" = "grey45", "Treat none" = "grey75"), name = NULL) +
  scale_linetype_manual(values = c("Frozen 4-gene model" = 1, "Treat all" = 2,
                                   "Treat none" = 3), name = NULL) +
  coord_cartesian(ylim = c(-0.05, max(dca_fig$net_benefit, na.rm = TRUE) * 1.05)) +
  labs(x = "Threshold probability", y = "Net benefit") + th +
  theme(legend.position = "bottom", legend.key.width = unit(1.1, "lines"))
ggsave(file.path(figd, "F08_decision_curve_external.pdf"), p8, width = 8.4, height = 3.2)
ggsave(file.path(figd, "F08_decision_curve_external.png"), p8, width = 8.4, height = 3.2, dpi = 600)

# ---------------------------------------------------------------- 图注（英文，可直接贴投稿系统）
leg <- c(
 "Figure S-X1 | ROC curves of the frozen four-gene model. Discovery cohort (internal random split, n = 94) and three independent external cohorts (GSE5406, n = 210; GSE116250, n = 62; GSE55296, n = 33). No coefficient or z-score parameter was re-estimated in the external cohorts. Numbers in the legend are AUCs.",
 "Figure S-X2 | Calibration of the frozen model. Decile-based observed proportion versus predicted probability. The dashed line is the identity line. Slopes < 1 in the external cohorts indicate overfitting of the discovery model; absolute predicted probabilities are not transportable because the external cohorts were case-enriched (70-92% cases versus 57% in the discovery set).",
 "Figure S-X3 | Expression of AQP3, GGT5, ISLR and LAD1 by disease group in the three external cohorts (log2 scale; boxes: median and IQR; points: individual samples).",
 "Figure S-X4 | Discrimination of the frozen model overall and by aetiology. Points are AUCs with DeLong 95% confidence intervals; the dashed line marks AUC 0.5.",
 "Figure S-X5 | Directional agreement of the four genes. The leftmost column reproduces the discovery cohort; the remaining columns show log2 fold changes (cases versus non-failing) in the external cohorts. Asterisks denote limma p values (* p<0.05, ** p<0.01, *** p<0.001); no covariate adjustment was possible in GSE5406 and GSE55296, and none was applied in GSE116250 for comparability.",
 "Figure S-X6 | Platform scale diagnostics. Mean expression of each gene in the external cohort versus the training cohort; the dashed line is identity. SD ratios far from 1 indicate that raw values cannot be plugged into the frozen model without scale alignment.",
 "Figure S-X7 | Single-gene and leave-one-gene-out discrimination with frozen coefficients. Only the terms shown were retained in the linear predictor; this panel quantifies how much of the panel's performance is carried by ISLR.",
 "Figure S-X8 | Decision curve analysis of the frozen model in the external cohorts. Net benefit of the model versus treating all and treating none across threshold probabilities.")
writeLines(leg, file.path(figd, "figure_legends_external.txt"))

say(""); say("DONE; figures written to ", figd); close(logcon)
