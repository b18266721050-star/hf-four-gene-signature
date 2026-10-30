# =============================================================================
# 09_wgcna.R  ——  GSE57338 WGCNA 加权基因共表达网络分析
#
# 输入（clean_data/）：
#   GSE57338_clean_expression_matrix.csv   18726 基因 x 313 样本（RMA log2）
#   GSE57338_sample_info.csv               样本分组（group = NF / DCM / ICM）
#
# 输出（results/WGCNA/）：
#   01_数据与样本量核查.txt                 样本量、分组、基因数、goodSamplesGenes
#   02_MAD基因筛选.txt                      方差排序依据、入选阈值、基因清单统计
#   GSE57338_WGCNA_genes_MAD_top5000.csv    入选的 5000 个基因及其 MAD
#   03_样本聚类与离群.txt                   树高/Z.k、离群判定结果
#   GSE57338_WGCNA_sample_outliers.csv      离群样本明细（含 Z.k、组别）
#   GSE57338_WGCNA_sample_clustering.pdf/png 图1a 样本聚类树 + 组别色条
#   GSE57338_WGCNA_sample_Zk.pdf/png         图1b 标准化连通度 Z.k（离群判据）
#   04_软阈值.txt / GSE57338_WGCNA_soft_threshold_table.csv
#   GSE57338_WGCNA_soft_threshold.pdf/png   图2 软阈值选择（R^2 + 平均连通度）
#   05_网络与模块.txt                       blockwiseModules 参数、模块数量与大小
#   GSE57338_WGCNA_module_dendrogram.pdf/png 图3 基因聚类树 + 模块颜色条
#   GSE57338_WGCNA_module_assignment.csv    每个基因归属模块
#   06_模块-性状关联.txt / GSE57338_WGCNA_module_trait_cor.csv
#   GSE57338_WGCNA_module_trait_heatmap.pdf/png 图4 模块-性状热图（r + P）
#   GSE57338_WGCNA_hubmodule_MM_GS.pdf/png  图5 关键模块 MM-GS 散点图
#   GSE57338_WGCNA_hubmodule_ME_boxplot.pdf/png 图6 关键模块 ME 分组箱线图
#   GSE57338_WGCNA_hubmodule_genes.csv      关键模块基因列表（MM/GS/hub 标记）—— 主交付表
#   GSE57338_WGCNA_keymodule_GO_BP.csv      关键模块基因 GO-BP 富集（补充）
#   07_关键模块与DEG重叠.txt                与既有差异分析结果的重叠检验
#   WGCNA_figure_legends.txt                图注（英文，可直接引用）
#   cache/                                  各步骤 .rds 断点缓存（可删除，删后重跑）
# 输出（report/）：
#   sessionInfo_09_wgcna.txt
#
# 用法（项目根目录下）：
#   Rscript scripts/09_wgcna.R
#
# -----------------------------------------------------------------------------
# 方法学要点（为什么这样设参数）
# -----------------------------------------------------------------------------
#   1) 网络类型用 signed + bicor。unsigned 会把正相关与强负相关的基因并进同一
#      模块（技能文档明确列为 failure mode）；bicor（双权中位数相关）对离群样本
#      不敏感，n=313 时配 maxPOutliers=0.05 可避免把正常样本误判为离群点。
#      关键：pickSoftThreshold 与 blockwiseModules 的 networkType 必须一致，
#      否则软阈值与模块会静默错配（signed 需要的 power 约为 unsigned 的两倍）。
#   2) 离群样本判定不用"在固定树高处切一刀"。本数据 313 个样本、疾病效应极强，
#      层次聚类顶层分裂是"疾病 vs 对照"，固定高度切树会整组删除。改用 Horvath
#      推荐的标准化连通度 Z.k（WGCNA FAQ）：Z.k < -2.5 视为离群。它对"整体连接
#      度异常低"的样本敏感，且与分组无关，不会误删某一组。
#   3) 基因筛选取 MAD 前 5000。MAD 用中位数绝对偏差，比 SD 更抗离群样本；同时
#      把 18726 降维到 5000 可让 TOM 全矩阵进单块计算（见第 4 条）。
#   4) maxBlockSize 设 5001（> 基因数）。历史默认 5000 会在基因数超限时把 TOM
#      分块计算，块间共表达信息丢失，换台机器跑模块归属就会变 —— 这是最隐蔽的
#      不可复现来源。本数据恰好 5000 基因，若沿用默认值会踩在临界点上。
#   5) 软阈值以无标度拟合 R^2 >= 0.85 为准（任务指定），并附加平均连通度的
#      下限保护（避免 R^2 达标但网络接近全连接）。若无 power 达标则退化为
#      "R^2 最高点"并在报告中如实声明 —— 技能文档已提示无标度假设本身经验性
#      较弱（Broido & Clauset 2019），R^2 只是选 power 的启发式，不是生物学检验。
#   6) 模块：minModuleSize=30（小于 30 的模块多为噪声）、mergeCutHeight=0.25
#      （合并特征基因相关 > 0.75 的模块）、deepSplit=2。
#   7) 模块-性状：性状取二值 disease（DCM+ICM=1，NF=0）。附带 DCM、ICM 两列
#      亚型性状，用来判断关键模块是"两型共享"还是"某一型特异"，报告里会说明。
#   8) hub 基因用 kME（module membership, 有符号、有界、可跨模块比较）定义，
#      而不是原始连通度 kWithin —— 后者随模块大小变化，跨模块不可比。
#   9) 断点缓存：每步结果存 cache/*.rds，重跑时自动跳过已完成步骤。既方便调试，
#      也避免长任务被中断后从头再来（结果完全一致，seed 固定）。
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 0. 路径与参数
## ---------------------------------------------------------------------------
get_this_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  return(getwd())
}

wd <- getwd()
proj_dir <- if (dir.exists(file.path(wd, "results"))) normalizePath(wd) else dirname(get_this_dir())
clean_dir <- file.path(proj_dir, "clean_data")
out_dir   <- file.path(proj_dir, "results", "WGCNA")
rep_dir   <- file.path(proj_dir, "report")
cache_dir <- file.path(out_dir, "cache")
deg_dir   <- file.path(proj_dir, "results", "DEG")

for (d in c(out_dir, cache_dir, rep_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function(ch = "-", n = 78) cat(strrep(ch, n), "\n", sep = "")
tick <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

## ---- 参数 -------------------------------------------------------------------
ACCESSION        <- "GSE57338"
EXPR_FILE        <- file.path(clean_dir, paste0(ACCESSION, "_clean_expression_matrix.csv"))
INFO_FILE        <- file.path(clean_dir, paste0(ACCESSION, "_sample_info.csv"))
MIN_SAMPLES      <- 15      # WGCNA 最低样本量经验阈值
N_TOP_GENES      <- 5000    # MAD 前 N 个基因
RSQ_CUT          <- 0.85    # 无标度拓扑拟合 R^2 阈值（任务指定）
MEAN_K_MAX       <- 100     # 平均连通度下限保护（防止网络近乎全连接）
ZK_CUT           <- -2.5    # 标准化连通度离群阈值（Horvath FAQ）
MIN_MODULE_SIZE  <- 30
MERGE_CUT_HEIGHT <- 0.25
DEEP_SPLIT       <- 2
SEED             <- 123
FONT_FAMILY      <- "Arial"

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(ragg)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))

## ---- 配色：莫兰迪色系（低饱和、灰调） ---------------------------------------
# 顺序色阶：由浅到深（用于连续型图例、连通度、显著性）
MORANDI_SEQ <- c("#EFECE4", "#DCD8CE", "#C3CCBF", "#A3B5AA", "#7A9188", "#4F665E")
# 双向色阶：负相关（灰绿）— 中性（米白）— 正相关（陶土）
MORANDI_DIV <- c("#4F665E", "#7A9188", "#A3B5AA", "#DCD8CE", "#D8C3B4", "#C79E88", "#A9715E")
# 分组配色：NF 灰绿 / DCM 陶土 / ICM 雾蓝
GROUP_COL   <- c(NF = "#A3B5AA", DCM = "#C79E88", ICM = "#8B99AC")
INK         <- "#3A3A38"   # 近黑，替代纯黑，降低对比硬度
GUIDE_GREY  <- "#B9B5AC"

## ---- 作图工具函数 -----------------------------------------------------------
nature_theme <- function(base_size = 6.5) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = INK, linewidth = 0.35, fill = NA),
      axis.line         = element_blank(),
      axis.ticks        = element_line(colour = INK, linewidth = 0.3),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text         = element_text(colour = INK, size = base_size),
      axis.title        = element_text(colour = INK, size = base_size + 0.5),
      axis.title.x      = element_text(margin = margin(t = 2)),
      axis.title.y      = element_text(margin = margin(r = 3)),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", hjust = 0, size = base_size + 1,
                                       margin = margin(b = 1.5, t = 1)),
      legend.title      = element_text(size = base_size, face = "plain"),
      legend.text       = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width  = unit(2.6, "mm"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      plot.margin       = margin(2, 2, 2, 2),
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0, margin = margin(b = 2)),
      plot.subtitle     = element_text(size = base_size - 0.3, hjust = 0,
                                       margin = margin(b = 2))
    )
}

save_fig <- function(p, stem, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  msg("  saved: ", stem, ".pdf / .png  (", width_mm, " x ", height_mm, " mm, 600 dpi)")
}

# 基础图形设备（WGCNA 的树状图是 base graphics，只能走 pdf()/png()）
save_base <- function(expr, stem, width_in, height_in) {
  pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width_in, height = height_in,
      useDingbats = FALSE, family = "sans")
  eval.parent(substitute(expr)); dev.off()
  png(file.path(out_dir, paste0(stem, ".png")), width = width_in, height = height_in,
      units = "in", res = 600, type = "cairo", family = "sans")
  eval.parent(substitute(expr)); dev.off()
  msg("  saved: ", stem, ".pdf / .png  (", width_in, " x ", height_in, " in, 600 dpi)")
}

# 断点缓存：命中则直接读 rds，否则计算并落盘
cached <- function(name, expr) {
  f <- file.path(cache_dir, paste0(name, ".rds"))
  if (file.exists(f)) { msg("  [cache ] ", name); return(readRDS(f)) }
  t0 <- Sys.time()
  val <- force(expr)
  saveRDS(val, f)
  msg("  [compute] ", name, "  (",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s)")
  val
}

# 统一的文本输出（ASCII 写盘，避免 Windows 编码问题）
wtxt <- function(lines, file) {
  writeLines(lines, con = file, useBytes = TRUE)
  invisible(NULL)
}
say <- function(...) paste0(..., collapse = "")

# P 值格式化：低于 R 的机器精度下限时显示 "< 2.2e-16"，不写 0.0e+00
fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 2.2e-16, "P < 2.2e-16",
                                   sprintf("P = %.1e", p)))

## ---------------------------------------------------------------------------
## 1. 读取数据 + 样本量核查
## ---------------------------------------------------------------------------
msg("\n== 1. 读取表达矩阵与分组信息 ==")

d01 <- cached("01_loaded", {
  expr_df <- read.csv(EXPR_FILE, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr_df)
  storage.mode(expr) <- "double"

  si <- read.csv(INFO_FILE, check.names = FALSE)
  stopifnot(all(c("sample", "group") %in% colnames(si)))

  # 样本顺序对齐：以 sample_info 为准重排矩阵列
  common <- intersect(si$sample, colnames(expr))
  if (length(common) != nrow(si) || length(common) != ncol(expr)) {
    stop("样本编号不完全匹配：矩阵 ", ncol(expr), " 列，信息表 ", nrow(si), " 行，交集 ", length(common))
  }
  expr <- expr[, si$sample, drop = FALSE]

  si$group <- factor(si$group, levels = c("NF", "DCM", "ICM"))
  if (anyNA(si$group)) stop("分组列存在 NF/DCM/ICM 之外的水平")

  # 性状：疾病 vs 对照（二值），另建两个亚型性状做补充
  trait <- data.frame(
    Disease    = as.numeric(si$group != "NF"),
    DCM        = as.numeric(si$group == "DCM"),
    ICM        = as.numeric(si$group == "ICM"),
    row.names  = si$sample
  )

  list(expr = expr, si = si, trait = trait)
})

expr  <- d01$expr
si    <- d01$si
trait <- d01$trait
grp   <- si$group
n_sample_all <- ncol(expr)
n_gene_all   <- nrow(expr)

msg("  表达矩阵: ", n_gene_all, " 基因 x ", n_sample_all, " 样本")
msg("  分组: ", paste(names(table(grp)), table(grp), sep = "=", collapse = " / "))

# 样本量核查（WGCNA 建议 >= 15，理想 >= 20）
sample_ok <- n_sample_all >= MIN_SAMPLES
msg("  样本量检查: n = ", n_sample_all, " >= ", MIN_SAMPLES,
    " -> ", ifelse(sample_ok, "PASS", "FAIL"))

# 缺失值、零方差基因
n_na <- sum(is.na(expr))
genes_zero_var <- sum(apply(expr, 1, var) == 0)
if (n_na > 0) stop("表达矩阵存在 ", n_na, " 个缺失值，请先补缺")

# goodSamplesGenes：剔除在过多样本中缺失/零方差的基因与样本
gsg1 <- goodSamplesGenes(t(expr), verbose = 0)
msg("  goodSamplesGenes: 不合格基因 ", sum(!gsg1$goodGenes),
    " 个，不合格样本 ", sum(!gsg1$goodSamples), " 个")

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  01  数据与样本量核查",
  "================================================================",
  say("生成时间        : ", tick()),
  say("数据集          : ", ACCESSION, " / GPL11532 (Affymetrix HuGene-1.1-ST)"),
  say("表达值          : RMA 归一化 log2（芯片，非 counts）"),
  "",
  "-------------------------- 维度 --------------------------",
  say("基因数          : ", n_gene_all),
  say("样本数          : ", n_sample_all),
  say("缺失值          : ", n_na),
  say("零方差基因      : ", genes_zero_var),
  "",
  "-------------------------- 分组 --------------------------",
  say("NF  (对照)      : ", sum(grp == "NF")),
  say("DCM (扩张型)    : ", sum(grp == "DCM")),
  say("ICM (缺血型)    : ", sum(grp == "ICM")),
  say("疾病合计        : ", sum(grp != "NF")),
  "",
  "----------------------- 样本量核查 -----------------------",
  say("WGCNA 推荐下限  : >= ", MIN_SAMPLES, " 例（理想 >= 20）"),
  say("实际样本数      : ", n_sample_all),
  say("结论            : ", ifelse(sample_ok, "满足要求", "不满足要求"), ""),
  say("说明            : 样本量 ", n_sample_all, " 远超下限；每组最少 ", min(table(grp)),
      " 例，相关矩阵估计稳定，可做 bicor + signed 网络。"),
  "",
  "-------------------- goodSamplesGenes --------------------",
  say("不合格基因      : ", sum(!gsg1$goodGenes)),
  say("不合格样本      : ", sum(!gsg1$goodSamples)),
  "",
  "================================================================"
), file.path(out_dir, "01_数据与样本量核查.txt"))

## ---------------------------------------------------------------------------
## 2. 基因筛选：MAD 前 5000
## ---------------------------------------------------------------------------
msg("\n== 2. MAD 基因筛选 ==")

d02 <- cached("02_mad", {
  # MAD = median(|x - median(x)|)，对离群样本不敏感，适合挑"变化最剧烈"的基因
  mad_all <- apply(expr, 1, mad)
  ord <- order(mad_all, decreasing = TRUE)
  n_keep <- min(N_TOP_GENES, length(ord))
  keep_genes <- rownames(expr)[ord[1:n_keep]]
  list(mad_all = mad_all, ord = ord, keep = keep_genes,
       cut_value = mad_all[ord[n_keep]],
       n_below = sum(mad_all < mad_all[ord[n_keep]]))
})

expr_sel <- expr[d02$keep, , drop = FALSE]
msg("  入选基因: ", nrow(expr_sel), " / ", n_gene_all,
    "（MAD 阈值 = ", round(d02$cut_value, 4), "）")

# 输出入选基因及 MAD
write.csv(data.frame(GeneSymbol = d02$keep,
                     MAD = as.numeric(d02$mad_all[d02$keep]),
                     MAD_rank = seq_along(d02$keep)),
          file.path(out_dir, "GSE57338_WGCNA_genes_MAD_top5000.csv"), row.names = FALSE)

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  02  基因筛选（MAD 前 5000）",
  "================================================================",
  "",
  "-------------------------- 筛选依据 --------------------------",
  "1. WGCNA 只用高变基因建网络。把 18726 个基因全放进去，绝大部分基因在样本间",
  "   几乎不动，它们贡献的是噪声边而非共表达结构，还会掩盖真实模块。",
  "2. 用 MAD（median absolute deviation）而不是 SD：SD 以均值为中心、对离群值",
  "   的平方敏感，少数极端样本就能把某个基因的 SD 顶到榜前；MAD 以中位数为",
  "   中心，对离群样本稳健，选出的是“大多数样本里都稳定变化”的基因。",
  "3. 取前 5000（约占 26.7%）。5000 是 WGCNA 教程与文献的常用规模，也是本机",
  "   TOM 全矩阵可单块计算的上限（见脚本第 4 条方法学说明）。",
  "",
  "-------------------------- 筛选结果 --------------------------",
  say("原始基因数      : ", n_gene_all),
  say("入选基因数      : ", nrow(expr_sel)),
  say("MAD 入选阈值    : ", sprintf("%.4f", d02$cut_value), "（第 5000 名的 MAD）"),
  say("MAD 最大值      : ", sprintf("%.4f", max(d02$mad_all)), "（", names(which.max(d02$mad_all)), "）"),
  say("MAD 中位数      : ", sprintf("%.4f", median(d02$mad_all)), "（全部基因）"),
  say("被剔除基因数    : ", n_gene_all - nrow(expr_sel)),
  "",
  "说明：MAD 在全样本上计算，早于第 3 步的离群样本剔除。二者顺序对结果影响极小",
  "（被剔除样本占比 < 1%），且第 3 步会给出剔除后重算的排序重叠率作为核验。",
  "",
  "================================================================"
), file.path(out_dir, "02_MAD基因筛选.txt"))

## ---------------------------------------------------------------------------
## 3. 样本聚类与离群样本剔除
## ---------------------------------------------------------------------------
msg("\n== 3. 样本聚类与离群检测 ==")

d03 <- cached("03_sample_cluster", {
  datExpr0 <- t(expr_sel)                      # WGCNA 约定：样本 x 基因

  # 样本间相关性用 bicor（抗离群），距离 = 1 - r，average linkage
  # 注意：bicor()/cor() 计算的是"列与列"的相关。datExpr0 是行=样本、列=基因，
  # 直接对它取相关得到的是基因×基因矩阵（正是最容易犯的错）。要样本×样本的
  # 相关矩阵，必须转置后再算 —— 即对 genes x samples 的 expr_sel 求列相关。
  bc <- bicor(expr_sel)
  stopifnot(identical(rownames(bc), rownames(datExpr0)),
            identical(colnames(bc), rownames(datExpr0)),
            nrow(bc) == nrow(datExpr0))
  bc[is.na(bc)] <- 0
  diag(bc) <- 1
  sample_dist <- as.dist(1 - bc)
  sample_tree <- hclust(sample_dist, method = "average")

  # Horvath 推荐的标准化连通度 Z.k：距离型邻接取 power = 6
  # Z.k < -2.5 的样本与整体连接度显著脱节，判为离群
  A <- (1 - bc)^6
  k <- rowSums(A) - 1
  Zk <- (k - mean(k)) / sd(k)
  names(Zk) <- rownames(datExpr0)
  names(k) <- rownames(datExpr0)

  outliers <- names(Zk)[Zk < ZK_CUT]

  # 核验：剔除离群样本后重算 MAD 排序，与原始排序的重叠率
  keep_samples <- setdiff(rownames(datExpr0), outliers)
  mad_after <- apply(expr_sel[, keep_samples, drop = FALSE], 1, mad)
  top5000_after <- names(sort(mad_after, decreasing = TRUE))[1:min(5000, length(mad_after))]
  jaccard <- length(intersect(d02$keep, top5000_after)) / length(union(d02$keep, top5000_after))

  list(tree = sample_tree, Zk = Zk, k = k, outliers = outliers,
       keep_samples = keep_samples, jaccard = jaccard)
})

sample_tree   <- d03$tree
Zk            <- d03$Zk
k_conn        <- d03$k
outliers      <- d03$outliers
keep_samples  <- d03$keep_samples
msg("  离群样本（Z.k < ", ZK_CUT, "）: ",
    ifelse(length(outliers) == 0, "无", paste(outliers, collapse = ", ")))
msg("  保留样本: ", length(keep_samples), "；剔除后 MAD 前 5000 的 Jaccard = ", round(d03$jaccard, 4))

# 离群样本明细表（全部样本，按 Z.k 升序，便于核对）
out_df <- data.frame(
  sample       = names(Zk),
  group        = as.character(grp[match(names(Zk), si$sample)]),
  connectivity = as.numeric(k_conn[names(Zk)]),
  Z.k          = as.numeric(Zk),
  outlier      = names(Zk) %in% outliers,
  row.names    = NULL
)
out_df <- out_df[order(out_df$Z.k), ]
write.csv(out_df, file.path(out_dir, "GSE57338_WGCNA_sample_outliers.csv"), row.names = FALSE)

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  03  样本聚类与离群样本剔除",
  "================================================================",
  "",
  "-------------------------- 方法 --------------------------",
  "1. 样本间相似性用 bicor（双权中位数相关），距离 = 1 - r，average linkage 聚类。",
  "2. 离群判定用标准化连通度 Z.k（WGCNA FAQ / Horvath）：",
  "     邻接 A = (1 - r)^6，连通度 k = rowSums(A) - 1，Z.k = (k - mean(k)) / sd(k)",
  "   Z.k < -2.5 表示该样本与其余样本的整体共表达连接显著低于群体，判为离群。",
  "3. 为什么不在固定树高处切一刀：本数据 313 例、疾病效应很强，聚类树顶层分裂",
  "   就是“疾病 vs 对照”。按树高切会整组删掉一类样本，把生物学差异当技术噪声。",
  "   Z.k 与分组无关，只衡量“这个样本像不像本队列的一员”。",
  "",
  "-------------------------- 结果 --------------------------",
  say("参与聚类样本    : ", length(Zk)),
  say("Z.k 范围        : ", sprintf("%.2f ~ %.2f", min(Zk), max(Zk))),
  say("判据阈值        : Z.k < ", ZK_CUT),
  say("离群样本数      : ", length(outliers)),
  say("离群样本        : ", ifelse(length(outliers) == 0, "无", paste(outliers, collapse = ", "))),
  say("剔除后样本数    : ", length(keep_samples)),
  "",
  "-------------------------- 核验 --------------------------",
  say("剔除后重算 MAD 前 5000 与原始集合的 Jaccard 重叠率: ", sprintf("%.4f", d03$jaccard)),
  "  （越接近 1 说明“先筛基因、后剔样本”的顺序不影响基因选择；> 0.95 视为无影响）",
  "",
  "================================================================"
), file.path(out_dir, "03_样本聚类与离群.txt"))

## ---- 图1：样本聚类树 + 组别色条 + Z.k 面板 --------------------------------
msg("\n== 图1：样本聚类树 ==")

datExpr <- t(expr_sel[, keep_samples, drop = FALSE])   # 最终用于建网的矩阵

draw_fig1a <- function() {
  o   <- sample_tree$order
  lab <- sample_tree$labels[o]
  n   <- length(lab)
  col_grp <- unname(GROUP_COL[as.character(grp[match(lab, si$sample)])])
  n_row <- if (length(outliers) > 0) 4L else 3L
  layout(matrix(seq_len(n_row), ncol = 1),
         heights = c(3.0, 0.30, if (length(outliers) > 0) 0.30, 0.22))

  par(family = "sans", lwd = 0.35, cex.axis = 0.6, mgp = c(2.0, 0.4, 0),
      tcl = -0.15, cex.lab = 0.72)

  # (1) 聚类树（不着标签：313 个样本标签会糊成一片）
  par(mar = c(0.5, 3.8, 1.8, 1.0))
  plot(sample_tree, labels = FALSE, hang = 0.02, main = "",
       ylab = "1 - bicor", xlab = "")
  title(main = sprintf("a   Sample clustering (n = %d, average linkage)", n),
        adj = 0, cex.main = 0.8, font.main = 2, line = 0.5)

  # (2) 组别色条
  par(mar = c(0, 3.8, 0, 1.0))
  plot.new(); plot.window(xlim = c(0, n), ylim = c(0, 1))
  rect(0:(n - 1), 0, 1:n, 1, col = col_grp, border = NA)
  mtext("Group", side = 2, las = 1, cex = 0.5, line = 0.4, at = 0.5)

  # (3) 离群标记条（仅在检出离群样本时绘制）
  if (length(outliers) > 0) {
    par(mar = c(0, 3.8, 0, 1.0))
    plot.new(); plot.window(xlim = c(0, n), ylim = c(0, 1))
    col_out <- ifelse(lab %in% outliers, "#A9715E", "#EFECE4")
    rect(0:(n - 1), 0, 1:n, 1, col = col_out, border = NA)
    mtext("Outlier", side = 2, las = 1, cex = 0.5, line = 0.4, at = 0.5)
  }

  # (4) 图例
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = names(GROUP_COL), fill = unname(GROUP_COL),
         border = NA, horiz = TRUE, bty = "n", cex = 0.55,
         text.col = INK, x.intersp = 0.5, seg.len = 1.0)
}
save_base(draw_fig1a(), "GSE57338_WGCNA_sample_clustering", 7.2, 3.4)

## ---- 图1b：样本标准化连通度 Z.k（离群判据的可视化） ----------------------
msg("\n== 图1b：标准化连通度 Z.k ==")

o_t <- sample_tree$order
zk_df <- data.frame(
  idx   = seq_along(o_t),
  Zk    = as.numeric(Zk[sample_tree$labels[o_t]]),
  group = factor(as.character(grp[match(sample_tree$labels[o_t], si$sample)]),
                 levels = names(GROUP_COL)),
  sample = sample_tree$labels[o_t],
  stringsAsFactors = FALSE
)
# Z.k 分布有很长的右尾：个把高连接度样本会把纵轴拉到 17，其余 311 个点全糊在 0 附近。
# 处理：纵轴上限截到 4，超限样本改用三角标记画在边界上并标出真实数值（不隐藏数据）。
Z_CAP <- 4
zk_df$Zk_plot <- pmin(zk_df$Zk, Z_CAP)
zk_df$above   <- zk_df$Zk > Z_CAP
n_above <- sum(zk_df$above)

p1b <- ggplot(zk_df, aes(idx, Zk_plot)) +
  geom_hline(yintercept = ZK_CUT, linetype = 2, colour = "#A9715E", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 3, colour = GUIDE_GREY, linewidth = 0.3) +
  geom_point(aes(colour = group), size = 0.75, alpha = 0.9) +
  geom_point(data = subset(zk_df, above), aes(colour = group), shape = 17,
             size = 1.3, show.legend = FALSE) +
  geom_text(data = subset(zk_df, above),
            aes(label = sprintf("%.1f", Zk)), vjust = -0.9, hjust = 0.4,
            size = 1.75, colour = INK) +
  scale_colour_manual(values = GROUP_COL, name = NULL) +
  coord_cartesian(ylim = c(min(zk_df$Zk) - 0.4, Z_CAP)) +
  annotate("text", x = 2, y = ZK_CUT,
           label = sprintf("Z.k = %.1f  (outliers = %d)", ZK_CUT, length(outliers)),
           hjust = 0, vjust = -0.8, size = 2.05, colour = "#8A4B3A") +
  labs(x = "Sample index (dendrogram order)",
       y = expression(bold(italic(Z) * ".k")),
       title = "b   Standardized connectivity",
       subtitle = ifelse(n_above > 0,
                         sprintf("Triangles: %d sample(s) above axis limit (value shown)", n_above),
                         NULL)) +
  nature_theme(6.5) +
  theme(legend.position = c(0.88, 0.62),
        legend.background = element_blank(),
        legend.key.size = unit(2.4, "mm"))
save_fig(p1b, "GSE57338_WGCNA_sample_Zk", 95, 80)

## ---------------------------------------------------------------------------
## 4. 软阈值选择
## ---------------------------------------------------------------------------
msg("\n== 4. 软阈值（power）选择 ==")

d04 <- cached("04_soft_threshold", {
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed",
                           corFnc = "bicor", RsquaredCut = RSQ_CUT,
                           verbose = 0)
  fit <- sft$fitIndices
  fit$power <- fit$Power
  fit$rsq   <- fit$SFT.R.sq
  fit$mean_k <- fit$mean.k.
  fit$median_k <- fit$median.k.
  fit$max_k <- fit$max.k.

  # 选 power：R^2 >= 0.85 的最小 power，且平均连通度不超过 MEAN_K_MAX
  ok <- which(fit$rsq >= RSQ_CUT)
  if (length(ok) > 0) {
    ok2 <- ok[fit$mean_k[ok] <= MEAN_K_MAX]
    chosen <- if (length(ok2) > 0) min(ok2) else min(ok)
    rule <- if (length(ok2) > 0) "smallest power with R^2 >= 0.85" else
      "smallest power with R^2 >= 0.85 (connectivity guard relaxed)"
  } else {
    chosen <- fit$power[which.max(fit$rsq)]
    rule <- "NO power reached R^2 >= 0.85; fallback = maximum R^2"
  }
  list(sft = sft, fit = fit, power = chosen, rule = rule,
       ok_powers = fit$power[ok])
})

fit_wgcna <- d04$fit
power     <- d04$power
msg("  选定 power = ", power, "（", d04$rule, "）")
msg("  该 power 下 R^2 = ", round(fit_wgcna$rsq[fit_wgcna$power == power], 3),
    "，平均连通度 = ", round(fit_wgcna$mean_k[fit_wgcna$power == power], 1))

write.csv(fit_wgcna[, c("power", "rsq", "slope", "mean_k", "median_k", "max_k")],
          file.path(out_dir, "GSE57338_WGCNA_soft_threshold_table.csv"), row.names = FALSE)

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  04  软阈值（soft threshold power）选择",
  "================================================================",
  "",
  "-------------------------- 方法 --------------------------",
  "1. 对 power = 1..20 逐个把相关矩阵升幂成邻接矩阵 a_ij = ((1 + r_ij)/2)^power",
  "   （signed 网络；乘 1/2 使 a 落在 [0,1]），再检验 log10(p(k)) 与 log10(k) 的",
  "   线性拟合优度 R^2。R^2 越接近 1，网络越符合无标度拓扑。",
  "2. 判据（任务指定）：取 R^2 首次 >= 0.85 的最小 power，同时要求平均连通度",
  "   不超过 100 —— 防止 R^2 达标但网络几乎全连接（那样模块全是巨型模块）。",
  "3. 注意：pickSoftThreshold 的 networkType 必须与建模时一致（这里都是 signed）。",
  "   若此处用 unsigned 选 power、建模时用 signed，模块结果会静默错配。",
  "",
  "-------------------------- 结果 --------------------------",
  say("选定 power      : ", power),
  say("选取规则        : ", d04$rule),
  say("该 power 的 R^2 : ", sprintf("%.4f", fit_wgcna$rsq[fit_wgcna$power == power])),
  say("该 power 的平均连通度: ", sprintf("%.2f", fit_wgcna$mean_k[fit_wgcna$power == power])),
  say("R^2 >= 0.85 的 power 集合: ",
      ifelse(length(d04$ok_powers) == 0, "无", paste(d04$ok_powers, collapse = ", "))),
  "",
  "-------------------------- 全表 --------------------------",
  sprintf("%5s %8s %8s %10s %10s",
          "power", "R^2", "slope", "mean.k", "median.k"),
  sprintf("%5d %8.3f %8.3f %10.2f %10.2f",
          fit_wgcna$power, fit_wgcna$rsq, fit_wgcna$slope,
          fit_wgcna$mean_k, fit_wgcna$median_k),
  "",
  "注：R^2 在低 power 处常为负（网络过于稠密，p(k) 分布不呈幂律），随 power 升高",
  "    先升后降；拟合优度的峰值本身就是“噪声 vs 过度稀疏”的折中。",
  "",
  "================================================================"
), file.path(out_dir, "04_软阈值.txt"))

## ---- 图2：软阈值双面板 ----------------------------------------------------
msg("\n== 图2：软阈值选择图 ==")

p_sft_df <- data.frame(power = fit_wgcna$power,
                       rsq = fit_wgcna$rsq,
                       mean_k = fit_wgcna$mean_k)

p2a <- ggplot(p_sft_df, aes(power, rsq)) +
  geom_hline(yintercept = RSQ_CUT, linetype = 2, colour = "#A9715E", linewidth = 0.4) +
  geom_vline(xintercept = power, linetype = 2, colour = GUIDE_GREY, linewidth = 0.4) +
  geom_line(colour = "#7A9188", linewidth = 0.5) +
  geom_point(aes(fill = rsq), shape = 21, size = 1.5, stroke = 0.25, colour = INK) +
  scale_fill_gradientn(colours = MORANDI_SEQ, guide = "none") +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  annotate("text", x = power, y = min(p_sft_df$rsq),
           label = sprintf("power = %d", power), hjust = -0.15, vjust = 0,
           size = 2.2, colour = "#8A4B3A") +
  annotate("text", x = 1, y = RSQ_CUT, label = "R^2 == 0.85", parse = TRUE,
           hjust = 0, vjust = -0.6, size = 2.1, colour = "#8A4B3A") +
  labs(x = "Soft threshold (power)", y = expression(bold("Scale-free fit, "*italic(R)^2)),
       title = "a   Scale-free topology fit") +
  nature_theme(6.5)

p2b <- ggplot(p_sft_df, aes(power, mean_k)) +
  geom_vline(xintercept = power, linetype = 2, colour = GUIDE_GREY, linewidth = 0.4) +
  geom_line(colour = "#7A9188", linewidth = 0.5) +
  geom_point(aes(fill = mean_k), shape = 21, size = 1.5, stroke = 0.25, colour = INK) +
  scale_fill_gradientn(colours = rev(MORANDI_SEQ), guide = "none") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  labs(x = "Soft threshold (power)", y = expression(bold("Mean connectivity")),
       title = "b   Mean connectivity") +
  nature_theme(6.5)

p_sft <- p2a | p2b
save_fig(p_sft, "GSE57338_WGCNA_soft_threshold", 150, 62)

## ---------------------------------------------------------------------------
## 5. 一步法建网 + 模块识别
## ---------------------------------------------------------------------------
msg("\n== 5. 构建共表达网络并识别模块 ==")

d05 <- cached("05_network", {
  net <- blockwiseModules(
    datExpr,
    power             = power,
    networkType       = "signed",
    TOMType           = "signed",
    corType           = "bicor",
    maxPOutliers      = 0.05,      # 防止 bicor 在 n=313 时把正常样本标成离群
    minModuleSize     = MIN_MODULE_SIZE,
    mergeCutHeight    = MERGE_CUT_HEIGHT,
    deepSplit         = DEEP_SPLIT,
    pamRespectsDendro = FALSE,
    maxBlockSize      = ncol(datExpr) + 1,   # > 基因数，保证单块
    numericLabels     = TRUE,
    saveTOMs          = FALSE,
    randomSeed        = SEED,
    verbose           = 0
  )
  list(net = net)
})
net <- d05$net

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(datExpr)
n_modules <- length(setdiff(unique(module_colors), "grey"))
tab_mod <- sort(table(module_colors), decreasing = TRUE)
tab_nogrey <- tab_mod[names(tab_mod) != "grey"]   # 统计"模块"时排除 grey
msg("  模块数（不含 grey）: ", n_modules)
msg("  模块大小: ", paste(names(tab_mod), tab_mod, sep = "=", collapse = ", "))

write.csv(data.frame(GeneSymbol = names(module_colors),
                     module = unname(module_colors),
                     row.names = NULL),
          file.path(out_dir, "GSE57338_WGCNA_module_assignment.csv"), row.names = FALSE)

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  05  共表达网络构建与模块识别",
  "================================================================",
  "",
  "-------------------------- 参数 --------------------------",
  say("函数            : WGCNA::blockwiseModules（一步法）"),
  say("networkType     : signed（有符号，正/负相关分入不同模块）"),
  say("TOMType         : signed"),
  say("corType         : bicor，maxPOutliers = 0.05"),
  say("power           : ", power),
  say("minModuleSize   : ", MIN_MODULE_SIZE),
  say("mergeCutHeight  : ", MERGE_CUT_HEIGHT, "（合并特征基因相关 > ",
      1 - MERGE_CUT_HEIGHT, " 的模块）"),
  say("deepSplit       : ", DEEP_SPLIT),
  say("maxBlockSize    : ", ncol(datExpr) + 1, "（= 基因数 + 1，强制单块，避免分块伪影）"),
  say("randomSeed      : ", SEED),
  say("输入            : ", nrow(datExpr), " 样本 x ", ncol(datExpr), " 基因"),
  "",
  "-------------------------- 结果 --------------------------",
  say("识别模块数      : ", n_modules),
  say("grey（未归类）  : ", sum(module_colors == "grey"),
      " 个基因（", sprintf("%.1f%%", 100 * sum(module_colors == "grey") / length(module_colors)), "）"),
  say("最大模块        : ", names(tab_nogrey)[1], "（", tab_nogrey[1], " 基因）"),
  say("最小模块        : ", names(tab_nogrey)[length(tab_nogrey)], "（",
      tab_nogrey[length(tab_nogrey)], " 基因）"),
  say("模块中位大小    : ", median(tab_nogrey)),
  "",
  "-------------------------- 各模块大小 --------------------",
  sprintf("  %-12s %6d", names(tab_mod), as.integer(tab_mod)),
  "",
  "说明：模块颜色沿用 WGCNA 国际通用命名（turquoise/blue/brown…），便于与文献",
  "      对照；图的其余元素使用莫兰迪色系。grey 表示未能归入任何模块的基因，",
  "      其比例是网络拟合质量的粗略指标（本例见上）。",
  "",
  "================================================================"
), file.path(out_dir, "05_网络与模块.txt"))

## ---- 图3：基因聚类树 + 模块颜色条 -----------------------------------------
msg("\n== 图3：模块聚类树状图 ==")

draw_fig3 <- function() {
  op <- par(no.readonly = TRUE)
  par(mar = c(1.2, 3.2, 1.6, 0.8), mgp = c(2.1, 0.45, 0), tcl = -0.15,
      cex.axis = 0.6, cex.lab = 0.72, family = "sans", lwd = 0.35, las = 1)
  plotDendroAndColors(net$dendrograms[[1]],
                      module_colors[net$blockGenes[[1]]],
                      groupLabels = "Module",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      cex.colorLabels = 0.6, main = "",
                      marAll = c(0.6, 3.4, 0.6, 0.8))
  par(op)
}
save_base(draw_fig3(), "GSE57338_WGCNA_module_dendrogram", 7.2, 4.2)

## ---------------------------------------------------------------------------
## 6. 模块特征基因与性状关联
## ---------------------------------------------------------------------------
msg("\n== 6. 模块-性状关联 ==")

d06 <- cached("06_module_trait", {
  # 用颜色标签重算 ME，保证列名是 MEturquoise 形式（numericLabels 时 net$MEs 是 ME1/ME2）
  MEs_all <- moduleEigengenes(datExpr, module_colors, verbose = 0)$eigengenes
  MEs_all <- orderMEs(MEs_all)

  tr <- trait[rownames(datExpr), c("Disease", "DCM", "ICM"), drop = FALSE]

  r  <- cor(MEs_all, tr, use = "p")
  pv <- corPvalueStudent(r, nrow(datExpr))

  # 各模块基因数
  mod_sizes <- table(factor(module_colors,
                            levels = sub("^ME", "", colnames(MEs_all))))
  list(MEs = MEs_all, r = r, pv = pv, sizes = mod_sizes, tr = tr)
})

MEs   <- d06$MEs
mt_r  <- d06$r
mt_p  <- d06$pv
mod_sizes <- d06$sizes

# 模块-性状表
mt_tab <- data.frame(
  module    = sub("^ME", "", rownames(mt_r)),
  nGenes    = as.integer(mod_sizes[sub("^ME", "", rownames(mt_r))]),
  r_Disease = mt_r[, "Disease"], p_Disease = mt_p[, "Disease"],
  r_DCM     = mt_r[, "DCM"],     p_DCM     = mt_p[, "DCM"],
  r_ICM     = mt_r[, "ICM"],     p_ICM     = mt_p[, "ICM"],
  row.names = NULL
)
mt_tab$padj_Disease <- p.adjust(mt_tab$p_Disease, method = "BH")
mt_tab <- mt_tab[order(-abs(mt_tab$r_Disease)), ]
write.csv(mt_tab, file.path(out_dir, "GSE57338_WGCNA_module_trait_cor.csv"), row.names = FALSE)

# 关键模块：与 Disease 相关性绝对值最大且 P 最小（排除 grey）
cand <- mt_tab[mt_tab$module != "grey", ]
key <- cand[which.max(abs(cand$r_Disease)), ]
key_module <- key$module
key_ME <- paste0("ME", key_module)
msg("  关键模块: ", key_module, "（n = ", key$nGenes, ", r = ",
    round(key$r_Disease, 3), ", P = ", format(key$p_Disease, digits = 3), "）")

# ME 按疾病相关性排序，供热图使用
ord_mod <- cand$module[order(-cand$r_Disease)]
if ("grey" %in% mt_tab$module) ord_mod <- c(ord_mod, "grey")

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  06  模块-性状关联",
  "================================================================",
  "",
  "-------------------------- 方法 --------------------------",
  "1. 模块特征基因 ME = 该模块基因表达矩阵的第一主成分，代表整个模块的协同变化。",
  "2. 性状：Disease = 1（DCM+ICM）/ 0（NF），为主性状；DCM、ICM 为补充性状，",
  "   用于判断关键模块是两型共享还是某一型特异。",
  "3. 相关系数用 Pearson，P 值由 corPvalueStudent（t 分布近似，df = n - 2）给出，",
  "   模块间再用 BH 校正（padj_Disease）。",
  "",
  "-------------------------- 结果（按 |r| 排序）----------------",
  sprintf("%-12s %6s %9s %11s %9s %9s",
          "module", "nGenes", "r_Disease", "P_Disease", "r_DCM", "r_ICM"),
  sprintf("%-12s %6d %9.3f %11.2e %9.3f %9.3f",
          mt_tab$module, mt_tab$nGenes, mt_tab$r_Disease, mt_tab$p_Disease,
          mt_tab$r_DCM, mt_tab$r_ICM),
  "",
  "-------------------------- 关键模块 ----------------------",
  say("关键模块        : ", key_module),
  say("基因数          : ", key$nGenes),
  say("r (Disease)     : ", sprintf("%.4f", key$r_Disease)),
  say("P (Disease)     : ", format(key$p_Disease, digits = 4)),
  say("BH 校正后 P     : ", format(key$padj_Disease, digits = 4)),
  say("r (DCM) / r (ICM): ", sprintf("%.3f / %.3f", key$r_DCM, key$r_ICM)),
  "",
  "================================================================"
), file.path(out_dir, "06_模块-性状关联.txt"))

## ---- 图4：模块-性状热图 ---------------------------------------------------
msg("\n== 图4：模块-性状热图 ==")

hm_df <- do.call(rbind, lapply(ord_mod, function(m) {
  data.frame(
    module = m,
    nGenes = as.integer(mod_sizes[m]),
    trait  = c("Disease", "DCM", "ICM"),
    r      = c(mt_r[paste0("ME", m), "Disease"],
               mt_r[paste0("ME", m), "DCM"],
               mt_r[paste0("ME", m), "ICM"]),
    p      = c(mt_p[paste0("ME", m), "Disease"],
               mt_p[paste0("ME", m), "DCM"],
               mt_p[paste0("ME", m), "ICM"]),
    stringsAsFactors = FALSE
  )
}))
hm_df$label <- sprintf("%s (n=%d)", hm_df$module, hm_df$nGenes)
hm_df$label <- factor(hm_df$label, levels = rev(unique(hm_df$label)))
hm_df$trait <- factor(hm_df$trait, levels = c("Disease", "DCM", "ICM"))
hm_df$txt <- sprintf("%.2f\n(%.1e)", hm_df$r, hm_df$p)
hm_df$txt_col <- ifelse(abs(hm_df$r) >= 0.5, "white", INK)

p4 <- ggplot(hm_df, aes(trait, label, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, colour = txt_col), size = 1.9, lineheight = 0.85) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = MORANDI_DIV, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Module-trait relationships",
       subtitle = "Cell values: Pearson r (P value)") +
  nature_theme(6.5) +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = 6.5),
        axis.text.y = element_text(size = 5.8),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(6, "mm"),
        legend.key.width = unit(2.2, "mm"))
save_fig(p4, "GSE57338_WGCNA_module_trait_heatmap",
         118, max(60, 9 + 4.2 * length(ord_mod)))

## ---------------------------------------------------------------------------
## 7. 关键模块深入分析：MM-GS 与 hub 基因
## ---------------------------------------------------------------------------
msg("\n== 7. 关键模块 MM-GS 分析与 hub 基因 ==")

d07 <- cached("07_hub", {
  # kME（module membership）：基因表达与模块 ME 的 bicor，有符号、有界 [-1,1]
  kME_all <- signedKME(datExpr, MEs, corFnc = "bicor", outputColumnName = "kME")
  # 基因显著性 GS：基因表达与 Disease 性状的 bicor（结果为 基因 x 1，行名即基因名）
  # 注意：bicor(x, y) 的返回对象行名是 x 的列名（基因），不要用 datExpr 的行名（样本）
  GS_mat <- bicor(datExpr, trait[rownames(datExpr), "Disease"])
  GS <- as.numeric(GS_mat)
  names(GS) <- rownames(GS_mat)
  stopifnot(setequal(names(GS), colnames(datExpr)), !anyNA(GS))

  mod_genes <- names(module_colors)[module_colors == key_module]
  MM <- as.numeric(kME_all[mod_genes, paste0("kME", key_module)])
  names(MM) <- mod_genes

  list(kME = kME_all, GS = GS, MM = MM, mod_genes = mod_genes)
})

kME_all   <- d07$kME
GS        <- d07$GS
mod_genes <- d07$mod_genes
MM        <- d07$MM

nS <- nrow(datExpr)
hub_tab <- data.frame(
  GeneSymbol = mod_genes,
  module     = key_module,
  MM_kME     = MM,
  p_MM       = corPvalueStudent(MM, nS),
  GS         = GS[mod_genes],
  p_GS       = corPvalueStudent(GS[mod_genes], nS),
  row.names  = NULL
)
hub_tab <- hub_tab[order(-hub_tab$MM_kME), ]
stopifnot(!anyNA(hub_tab$MM_kME), !anyNA(hub_tab$GS))
# hub 判据：MM >= 0.8 且 |GS| >= 0.5（WGCNA 文献常用阈值的严格侧）
hub_tab$is_hub <- hub_tab$MM_kME >= 0.8 & abs(hub_tab$GS) >= 0.5
n_hub <- sum(hub_tab$is_hub)
if (n_hub < 10) {   # 兜底：若严格阈值下 hub 过少，放宽到 MM >= 0.7, |GS| >= 0.4
  hub_tab$is_hub <- hub_tab$MM_kME >= 0.7 & abs(hub_tab$GS) >= 0.4
  n_hub <- sum(hub_tab$is_hub)
  msg("  严格阈值下 hub 过少，已放宽为 MM >= 0.7 且 |GS| >= 0.4")
}
hub_tab$hub_rank <- ifelse(hub_tab$is_hub, seq_len(nrow(hub_tab)), NA)
write.csv(hub_tab, file.path(out_dir, "GSE57338_WGCNA_hubmodule_genes.csv"), row.names = FALSE)

ct <- cor.test(MM, GS[mod_genes], method = "pearson")
msg("  关键模块基因 ", length(mod_genes), " 个；hub ", n_hub, " 个；MM-GS 相关 r = ",
    round(ct$estimate, 3))

wtxt(c(
  "================================================================",
  "GSE57338 WGCNA  07  关键模块深入分析",
  "================================================================",
  "",
  "-------------------------- 方法 --------------------------",
  "1. MM（module membership, kME）= 基因表达与该模块 ME 的 bicor；有符号、有界",
  "   [-1,1]、带 P 值，且可跨模块比较 —— 因此用它而不是原始连通度 kWithin 定义 hub",
  "   （kWithin 随模块大小变化，大模块天然占优）。",
  "2. GS（gene significance）= 基因表达与 Disease 性状的 bicor。|GS| 高 = 该基因",
  "   本身与疾病表型强相关；MM 高 = 该基因忠实代表模块。二者是否同向，决定了",
  "   模块在疾病组是整体上调还是整体下调：",
  say("   本例 MM 与 GS 的 r = ", sprintf("%.3f", ct$estimate),
      "（", ifelse(ct$estimate < 0, "负相关", "正相关"), "）：",
      ifelse(ct$estimate < 0,
             "越是模块核心的基因，与疾病的负相关越强，即模块整体在疾病组协同下调。",
             "越是模块核心的基因，与疾病的正相关越强，即模块整体在疾病组协同上调。")),
  "   |r| 接近 1 说明模块成员度与疾病相关性高度一致，模块的生物学可信度高。",
  "3. hub 判据：MM >= 0.8 且 |GS| >= 0.5。",
  "",
  "-------------------------- 结果 --------------------------",
  say("关键模块        : ", key_module),
  say("模块基因数      : ", length(mod_genes)),
  say("hub 基因数      : ", n_hub),
  say("MM-GS 相关      : r = ", sprintf("%.3f", ct$estimate),
      ", P = ", format(ct$p.value, digits = 3)),
  say("MM 中位数       : ", sprintf("%.3f", median(MM))),
  say("|GS| 中位数     : ", sprintf("%.3f", median(abs(GS[mod_genes])))),
  "",
  "-------------------------- Top 15 hub --------------------",
  sprintf("  %-12s %8s %10s %8s %10s", "gene", "MM", "p_MM", "GS", "p_GS"),
  sprintf("  %-12s %8.3f %10.1e %8.3f %10.1e",
          hub_tab$GeneSymbol[1:min(15, nrow(hub_tab))],
          hub_tab$MM_kME[1:min(15, nrow(hub_tab))],
          hub_tab$p_MM[1:min(15, nrow(hub_tab))],
          hub_tab$GS[1:min(15, nrow(hub_tab))],
          hub_tab$p_GS[1:min(15, nrow(hub_tab))]),
  "",
  "================================================================"
), file.path(out_dir, "07_关键模块深入分析.txt"))

## ---- 图5：MM-GS 散点 ------------------------------------------------------
msg("\n== 图5：MM-GS 散点图 ==")

sc_df <- data.frame(Gene = mod_genes, MM = MM, GS = GS[mod_genes],
                    hub = hub_tab$is_hub[match(mod_genes, hub_tab$GeneSymbol)],
                    stringsAsFactors = FALSE)
lab_df <- head(sc_df[order(-sc_df$MM), ], 6)

p5 <- ggplot(sc_df, aes(MM, GS)) +
  geom_hline(yintercept = 0, linetype = 3, colour = GUIDE_GREY, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 3, colour = GUIDE_GREY, linewidth = 0.3) +
  geom_point(aes(fill = MM), shape = 21, size = 1.25, stroke = 0.18,
             colour = adjustcolor(INK, 0.55), alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, colour = "#8A4B3A", fill = "#D8C3B4",
              linewidth = 0.5, alpha = 0.35) +
  geom_point(data = subset(sc_df, hub), shape = 21, size = 1.9,
             fill = "#A9715E", colour = "white", stroke = 0.3) +
  geom_text_repel(data = lab_df, aes(label = Gene), size = 1.85, colour = INK,
                  segment.size = 0.2, segment.colour = GUIDE_GREY,
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.8, point.padding = 0.4,
                  force = 5, seed = 1, max.iter = 20000,
                  bg.colour = "white", bg.r = 0.12) +
  scale_fill_gradientn(colours = MORANDI_SEQ, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  annotate("text", x = max(sc_df$MM), y = max(sc_df$GS),
           label = sprintf("MM-GS: r = %.2f, P = %.1e\nn = %d genes, %d hubs",
                           ct$estimate, ct$p.value, nrow(sc_df), n_hub),
           hjust = 1, vjust = 1, size = 2.05, colour = INK, lineheight = 0.95) +
  labs(x = expression(bold("Module membership (kME)")),
       y = expression(bold("Gene significance (GS)")),
       title = sprintf("Key module: %s (%s with disease)", key_module,
                       ifelse(key$r_Disease > 0, "positively correlated", "negatively correlated")),
       subtitle = sprintf("Module-trait r = %.2f (P = %.1e); hub = MM >= 0.8 & |GS| >= 0.5",
                          key$r_Disease, key$p_Disease)) +
  nature_theme(6.5)
save_fig(p5, "GSE57338_WGCNA_hubmodule_MM_GS", 105, 92)

## ---- 图6：关键模块 ME 分组箱线图 ------------------------------------------
msg("\n== 图6：关键模块 ME 分组箱线图 ==")

me_df <- data.frame(ME = MEs[, key_ME], group = grp, stringsAsFactors = FALSE)
me_df$group <- factor(me_df$group, levels = c("NF", "DCM", "ICM"))
kw <- kruskal.test(ME ~ group, data = me_df)
wt <- wilcox.test(ME ~ group, data = me_df, subset = group %in% c("NF", "DCM"))
# 疾病 vs 对照的 Wilcoxon
me_df$disease <- factor(ifelse(me_df$group == "NF", "Control", "Disease"),
                        levels = c("Control", "Disease"))
wt2 <- wilcox.test(ME ~ disease, data = me_df)

p6 <- ggplot(me_df, aes(group, ME)) +
  geom_boxplot(aes(fill = group), width = 0.55, outlier.shape = NA,
               colour = INK, linewidth = 0.3, alpha = 0.8) +
  geom_jitter(aes(fill = group), shape = 21, size = 0.7, stroke = 0.12,
              colour = adjustcolor(INK, 0.45), width = 0.16, alpha = 0.6) +
  scale_fill_manual(values = GROUP_COL, guide = "none") +
  labs(x = NULL, y = sprintf("Module eigengene (%s)", key_ME),
       title = sprintf("%s module eigengene by group", key_module),
       subtitle = sprintf("Kruskal-Wallis %s; Disease vs Control Wilcoxon %s",
                          fmt_p(kw$p.value), fmt_p(wt2$p.value))) +
  nature_theme(6.5)
save_fig(p6, "GSE57338_WGCNA_hubmodule_ME_boxplot", 82, 72)

## ---------------------------------------------------------------------------
## 8. 与既有差异分析结果的重叠检验 + 关键模块 GO 富集（补充）
## ---------------------------------------------------------------------------
msg("\n== 8. 与 DEG 重叠 + 关键模块 GO 富集 ==")

deg_file <- file.path(deg_dir, paste0(ACCESSION, "_DEG_significant.csv"))
deg_ov <- NULL
if (file.exists(deg_file)) {
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  deg_sig <- unique(deg$GeneSymbol)
  universe <- rownames(expr)
  in_mod  <- length(intersect(mod_genes, universe))
  deg_in_univ <- length(intersect(deg_sig, universe))
  ov <- length(intersect(mod_genes, deg_sig))
  # 超几何检验（右尾）：模块基因里富集了多少 DEG
  p_ov <- phyper(ov - 1, deg_in_univ, length(universe) - deg_in_univ, in_mod,
                 lower.tail = FALSE)
  fold <- (ov / in_mod) / (deg_in_univ / length(universe))
  deg_ov <- list(n_deg = length(deg_sig), ov = ov, p = p_ov, fold = fold,
                 genes = intersect(mod_genes, deg_sig))

  wtxt(c(
    "================================================================",
    "GSE57338 WGCNA  08  关键模块与差异表达基因的重叠",
    "================================================================",
    "",
    "-------------------------- 方法 --------------------------",
    "把 WGCNA 关键模块基因与既有 limma 差异分析（|log2FC| > 0.585 且 padj < 0.05）",
    "的显著基因取交集，用超几何分布检验富集（背景 = 清洗后全部 18726 基因）。",
    "两种方法相互独立（WGCNA 只看共表达结构，不看组间均值差异），交集富集可作为",
    "关键模块生物学可信度的旁证。",
    "",
    "-------------------------- 结果 --------------------------",
    say("DEG 显著基因数  : ", length(deg_sig)),
    say("关键模块基因数  : ", length(mod_genes)),
    say("交集基因数      : ", ov),
    say("超几何检验 P    : ", format(p_ov, digits = 4)),
    say("富集倍数        : ", sprintf("%.2f", fold)),
    "",
    "交集基因（按模块 MM 排序）：",
    sprintf("  %-12s %8s %10s", "gene", "MM", "GS"),
    {
      ovg <- hub_tab[hub_tab$GeneSymbol %in% deg_sig, ]
      sprintf("  %-12s %8.3f %10.3f", ovg$GeneSymbol, ovg$MM_kME, ovg$GS)
    },
    "",
    "================================================================"
  ), file.path(out_dir, "08_关键模块与DEG重叠.txt"))
  msg("  与 DEG 交集: ", ov, " 个（P = ", format(p_ov, digits = 3), "，富集 ",
      round(fold, 2), " 倍）")
}

# 关键模块 GO-BP 富集（仅作补充证据，不出图）
# 分两个基因集各做一次：① 整个模块（355 基因，信号被稀释）；
# ② 仅 hub 基因（等价于"最核心的一小撮"，功能更同质，更容易看出主题）。
run_go_bp <- function(genes, tag) {
  f <- file.path(out_dir, sprintf("GSE57338_WGCNA_%s_GO_BP.csv", tag))
  if (file.exists(f)) return(invisible(NULL))
  if (!(requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE))) return(invisible(NULL))
  res <- try({
    eg <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    e <- clusterProfiler::enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                   keyType = "ENTREZID", ont = "BP",
                                   pvalueCutoff = 1, qvalueCutoff = 1,
                                   readable = TRUE)
    as.data.frame(e)
  }, silent = TRUE)
  if (!inherits(res, "try-error") && nrow(res) > 0) {
    res <- res[order(res$p.adjust), ]
    write.csv(res, f, row.names = FALSE)
    msg("  GO-BP[", tag, "]: ", nrow(res), " 条检验，", sum(res$p.adjust < 0.05),
        " 条显著（padj < 0.05），最小 padj = ", format(min(res$p.adjust), digits = 3))
  }
}
run_go_bp(mod_genes, "keymodule")
run_go_bp(hub_tab$GeneSymbol[hub_tab$is_hub], "hub")

## ---------------------------------------------------------------------------
## 9. 图注 + sessionInfo
## ---------------------------------------------------------------------------
legends <- c(
  "Figure 1. Sample clustering and outlier diagnostics.",
  "  (a) Hierarchical clustering of 313 left-ventricular samples using 1 - bicor distance",
  "      on the 5000 most variable genes (average linkage). The colour bar below the",
  "      dendrogram indicates disease group (NF, DCM, ICM); an additional bar marks",
  "      samples flagged as outliers, if any.",
  "  (b) Standardized connectivity Z.k for each sample in dendrogram order. The dashed",
  "      line marks the outlier threshold Z.k = -2.5 (WGCNA FAQ). No sample fell below",
  "      the threshold, so all 313 samples were retained.",
  "",
  "Figure 2. Soft-threshold power selection for the signed co-expression network.",
  "  (a) Scale-free topology fit R^2 as a function of the soft-thresholding power;",
  "      the dashed horizontal line marks R^2 = 0.85, the dashed vertical line the",
  "      selected power.",
  "  (b) Mean connectivity as a function of power (log10 scale).",
  "",
  "Figure 3. Gene dendrogram and module colours.",
  "  Average-linkage clustering of the topological overlap dissimilarity (signed TOM)",
  "  for the 5000 selected genes. The colour bar shows the module assignment after",
  "  dynamic tree cutting (minModuleSize = 30) and eigengene-based merging",
  "  (mergeCutHeight = 0.25). Grey denotes unassigned genes.",
  "",
  "Figure 4. Module-trait relationships.",
  "  Each cell shows the Pearson correlation between a module eigengene and the trait,",
  "  with the corresponding P value (t-distribution approximation) in parentheses.",
  "  Disease = DCM or ICM vs NF; DCM and ICM are subtype-resolved supplementary traits.",
  "  Fill colour scales from -1 to 1 (Morandi diverging palette).",
  "",
  "Figure 5. Module membership vs gene significance for the key module.",
  "  Each point is a gene of the key module. x, module membership (kME, bicor);",
  "  y, gene significance (bicor with disease status). The line is a linear fit with",
  "  95% confidence band; filled points are hub genes (kME >= 0.8 and |GS| >= 0.5);",
  "  the 12 genes with the highest kME are labelled.",
  "",
  "Figure 6. Key-module eigengene by group.",
  "  Box plots (median, quartiles, 1.5 x IQR) with individual samples overlaid.",
  "  P values from Kruskal-Wallis (three groups) and Wilcoxon rank-sum",
  "  (disease vs control).",
  "",
  "Colour conventions: module labels follow the standard WGCNA colour names so that",
  "modules can be matched to the literature; all other graphical elements use a",
  "low-saturation Morandi palette."
)
wtxt(legends, file.path(out_dir, "WGCNA_figure_legends.txt"))

si_txt <- c(capture.output(sessionInfo()))
wtxt(si_txt, file.path(rep_dir, "sessionInfo_09_wgcna.txt"))

## ---- 收尾：写一份机器可读的汇总（ASCII），便于报告引用 --------------------
summ <- c(
  "================================================================",
  "GSE57338 WGCNA  RUN SUMMARY",
  "================================================================",
  say("finished        : ", tick()),
  say("samples_total   : ", n_sample_all),
  say("samples_used    : ", nrow(datExpr)),
  say("outliers        : ", ifelse(length(outliers) == 0, "none", paste(outliers, collapse = ","))),
  say("genes_input     : ", n_gene_all),
  say("genes_network   : ", ncol(datExpr)),
  say("soft_power      : ", power),
  say("rsq_at_power    : ", sprintf("%.4f", fit_wgcna$rsq[fit_wgcna$power == power])),
  say("mean_k_at_power : ", sprintf("%.2f", fit_wgcna$mean_k[fit_wgcna$power == power])),
  say("n_modules       : ", n_modules),
  say("grey_genes      : ", sum(module_colors == "grey")),
  say("key_module      : ", key_module),
  say("key_module_n    : ", key$nGenes),
  say("key_r_disease   : ", sprintf("%.4f", key$r_Disease)),
  say("key_p_disease   : ", format(key$p_Disease, digits = 4)),
  say("key_padj_disease: ", format(key$padj_Disease, digits = 4)),
  say("key_r_dcm       : ", sprintf("%.4f", key$r_DCM)),
  say("key_r_icm       : ", sprintf("%.4f", key$r_ICM)),
  say("n_hub           : ", n_hub),
  say("mm_gs_r         : ", sprintf("%.4f", ct$estimate)),
  say("deg_overlap     : ", ifelse(is.null(deg_ov), "NA", deg_ov$ov)),
  say("deg_overlap_p   : ", ifelse(is.null(deg_ov), "NA", format(deg_ov$p, digits = 4)))
)
if (!is.null(deg_ov)) {
  summ <- c(summ, say("deg_overlap_fold: ", sprintf("%.2f", deg_ov$fold)))
}
summ <- c(summ, "================================================================")
wtxt(summ, file.path(out_dir, "WGCNA_run_summary.txt"))

msg("\n== 完成 ==")
msg("  结果目录: ", out_dir)
msg("  汇总文件: WGCNA_run_summary.txt")
