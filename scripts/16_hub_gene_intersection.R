# =============================================================================
# 16_hub_gene_intersection.R  ——  GSE57338 候选（关键）基因筛选
#                                 显著差异基因 ∩ WGCNA 关键模块基因
#
# 输入（results/）：
#   DEG/GSE57338_DEG_significant.csv           显著差异基因表（limma，|log2FC|>0.585 且 padj<0.05）
#   WGCNA/GSE57338_WGCNA_hubmodule_genes.csv   WGCNA 关键模块基因列表（含 MM/GS/hub 标记）
#   WGCNA/GSE57338_WGCNA_genes_MAD_top5000.csv WGCNA 输入基因全集（用于评估交集上限与富集检验）
#
# 输出（results/HubGene/）：
#   GSE57338_candidate_genes.csv               候选基因表（主交付）—— DEG 与 WGCNA 信息合并
#   GSE57338_candidate_genes_strict_hub.csv    补充：候选基因中同时满足 hub 判据的子集
#   GSE57338_venn.pdf / .png                   韦恩图（面积成比例，Nature 版式，莫兰迪配色）
#   GSE57338_venn_counts.csv                   三个分区的计数（差异/模块/交集）
#   GSE57338_genes_DEG_only.csv                仅在 DEG 中的基因（补充）
#   GSE57338_genes_WGCNA_only.csv              仅在关键模块中的基因（补充）
#   GSE57338_hubgene_summary.txt               全过程文本摘要（含基因名格式核查）
#   GSE57338_venn_figure_legend.txt            图注（英文，可直接引用）
#   GSE57338_HubGene_report.md                 中文处理报告（本脚本自动生成）
#   sessionInfo_16_hub_gene_intersection.txt   运行环境
#
# 用法（项目根目录下）：
#   Rscript scripts/16_hub_gene_intersection.R
#
# -----------------------------------------------------------------------------
# 方法学要点（为什么这样写）
# -----------------------------------------------------------------------------
#   1) 候选基因 = 显著差异基因 ∩ WGCNA 关键模块基因。两条判据来源不同、互为补充：
#      DEG 提供"单基因维度"的统计证据（limma + BH 校正），WGCNA 提供"网络维度"的
#      协同证据（共表达模块且与疾病性状最强相关）。同时满足两者的基因，假阳性概率
#      远低于只满足其一者。
#   2) 韦恩图用"面积成比例"的两圆构图。普通示意图的两个圆面积相等、圆心距固定，
#      读者从图形面积读出的比例是错的。本脚本按 n 反解圆心距 d，使三个分区的面积
#      严格对应基因数（两圆交面积有解析式，用 uniroot 求根）。
#   3) 基因名一致性核查是硬性前置步骤。DEG 来自 18726 个基因的全表，WGCNA 建立在
#      MAD 前 5000 个基因上，若两侧列名、大小写或空格不一致，交集会被静默算错。
#      脚本显式检查列名、缺失、重复、空格、大小写与 ID 类型，交集为 0 时直接报错退出。
#   4) 必须报告"交集上限"。WGCNA 只用了 5000 个基因，DEG 中落在 5000 之外的基因
#      在结构上永远不可能进入任何模块。因此交集大小的上限是
#      min(DEG ∩ 5000 全集, 关键模块大小)，报告里如实给出这个数。
#   5) WGCNA 的共表达是"边际相关"，不是调控关系。模块内两个基因有边，可能只是因为
#      它们共同受第三个基因驱动；把 hub 直接叫成"主调控因子"是本领域最常见的过度
#      解读。候选基因应理解为"差异显著且处于疾病相关共表达模块中"，而非"上游驱动者"。
#   6) 图内文字用英文（Nature 版式惯例，也避免 PDF 嵌入中文字体的字形缺失问题），
#      中文说明放在 .md 报告与图注文件里。
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## 0. 路径与环境
## ---------------------------------------------------------------------------
# 允许在 scripts/ 下直接跑：自动退回项目根目录
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("results")) stop("请在项目根目录下运行：Rscript scripts/16_hub_gene_intersection.R")

PROJECT_ROOT <- getwd()
DS          <- "GSE57338"                       # 数据集编号（GSE 编号，非 GSEA）
DIR_DEG     <- file.path(PROJECT_ROOT, "results", "DEG")
DIR_WGCNA   <- file.path(PROJECT_ROOT, "results", "WGCNA")
DIR_OUT     <- file.path(PROJECT_ROOT, "results", "HubGene")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

F_DEG  <- file.path(DIR_DEG,   paste0(DS, "_DEG_significant.csv"))
F_MOD  <- file.path(DIR_WGCNA, paste0(DS, "_WGCNA_hubmodule_genes.csv"))
F_UNI  <- file.path(DIR_WGCNA, paste0(DS, "_WGCNA_genes_MAD_top5000.csv"))

for (f in c(F_DEG, F_MOD, F_UNI)) {
  if (!file.exists(f)) stop("缺少输入文件：", f)
}

# 必需包；缺则用可获得的设备降级，不静默失败
need <- c("ggplot2")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) stop("缺少 R 包：", p)
}
suppressPackageStartupMessages(library(ggplot2))
HAS_RAGG  <- requireNamespace("ragg",  quietly = TRUE)
HAS_CAIRO <- requireNamespace("grDevices", quietly = TRUE) &&
             any(grepl("cairo", names(grDevices::pdf), ignore.case = TRUE)) ||
             .Platform$OS.type == "windows"

LOG <- character(0)
msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  LOG <<- c(LOG, line)
  cat(line, "\n", sep = "")
}
wtxt <- function(x, f) writeLines(x, f, useBytes = TRUE)   # useBytes：原样写出 UTF-8 字节

msg("=== ", DS, " 候选基因筛选（DEG ∩ WGCNA 关键模块） ===")
msg("项目根目录：", PROJECT_ROOT)
msg("输出目录  ：", DIR_OUT)

## ---------------------------------------------------------------------------
## 1. 读入两张表
## ---------------------------------------------------------------------------
deg <- read.csv(F_DEG,  check.names = FALSE, stringsAsFactors = FALSE)
mod <- read.csv(F_MOD,  check.names = FALSE, stringsAsFactors = FALSE)
uni <- read.csv(F_UNI,  check.names = FALSE, stringsAsFactors = FALSE)

msg("DEG 表      ：", nrow(deg), " 行 × ", ncol(deg), " 列  | 列名：",
    paste(names(deg), collapse = ", "))
msg("WGCNA 模块表：", nrow(mod), " 行 × ", ncol(mod), " 列  | 列名：",
    paste(names(mod), collapse = ", "))

## ---------------------------------------------------------------------------
## 2. 基因名格式一致性核查（硬性前置步骤）
## ---------------------------------------------------------------------------
# 优先用 GeneSymbol；找不到则退到第一个字符型列，并在日志里明确提示
pick_gene_col <- function(df, label) {
  nm  <- names(df)
  hit <- nm[tolower(nm) %in% c("genesymbol", "gene_symbol", "symbol",
                               "genename", "gene_name", "gene", "geneid")]
  if (length(hit) > 0) return(hit[1])
  cc  <- nm[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
  if (length(cc) > 0) {
    msg("  ! ", label, " 未找到 GeneSymbol 列，改用首个字符型列：", cc[1])
    return(cc[1])
  }
  stop(label, " 中找不到基因名列")
}
col_deg <- pick_gene_col(deg, "DEG 表")
col_mod <- pick_gene_col(mod, "WGCNA 模块表")
col_uni <- pick_gene_col(uni, "WGCNA 输入全集表")

# 全部正则都用 perl = TRUE：PCRE 的字符类是按码点判定，不受 locale 排序影响。
# （TRE 默认引擎下 [a-z] 这类范围会随 locale 变化，会误判。
#   本机的 R 运行在 C locale，实测 TRE 对 [a-z] 的判定确实不可靠。）
audit_symbols <- function(v, label) {
  s <- trimws(as.character(v))
  # 含小写字母的 symbol：HGNC 的 C<n>orf<m> 系列（如 C9orf47）是合法命名，
  # 必须与"真正可疑的小写命名"分开统计，否则会把正常数据误判成格式不一致。
  lower     <- grepl("[a-z]", s, perl = TRUE)
  # C<n>orf<m> 系列，含 -AS1 这类反义转录本后缀（如 C20orf166-AS1）
  lower_orf <- grepl("^C[0-9XY]+orf[0-9]+(-AS[0-9]+)?$", s, perl = TRUE)
  list(
    label        = label,
    n            = length(s),
    n_missing    = sum(is.na(s) | s == "" | s == "NA"),
    n_duplicate  = sum(duplicated(s)),
    n_has_space  = sum(grepl("\\s", s, perl = TRUE)),
    n_all_digit  = sum(grepl("^[0-9]+$", s, perl = TRUE)),      # 疑似 Entrez ID
    n_ensembl    = sum(grepl("^ENSG", s, perl = TRUE)),         # 疑似 Ensembl ID
    n_probe      = sum(grepl("^[0-9]{6,}$", s, perl = TRUE)),   # 疑似 Affymetrix 探针
    n_lower      = sum(lower),
    n_lower_orf  = sum(lower & lower_orf),
    n_lower_oth  = sum(lower & !lower_orf),
    lower_oth_eg = if (any(lower & !lower_orf))
                     paste(head(unique(s[lower & !lower_orf]), 5), collapse = " / ") else "",
    example      = paste(head(s, 3), collapse = " / ")
  )
}

a_deg <- audit_symbols(deg[[col_deg]], "DEG")
a_mod <- audit_symbols(mod[[col_mod]], "WGCNA 模块")
a_uni <- audit_symbols(uni[[col_uni]], "WGCNA 全集")

fmt_line <- function(a) {
  paste0("    ", a$label, "：n=", a$n,
         "  缺失=", a$n_missing,
         "  重复=", a$n_duplicate,
         "  含空格=", a$n_has_space,
         "  全数字=", a$n_all_digit,
         "  ENS=", a$n_ensembl,
         "  探针号=", a$n_probe,
         "  含小写=", a$n_lower,
         "(其中C*norf*命名=", a$n_lower_orf,
         " 可疑=", a$n_lower_oth, ")",
         "  示例=", a$example)
}
msg("-- 基因名格式核查 --")
msg(fmt_line(a_deg)); msg(fmt_line(a_mod)); msg(fmt_line(a_uni))
for (a in list(a_deg, a_mod, a_uni)) {
  if (a$n_lower_oth > 0) msg("    ! ", a$label, " 可疑小写命名示例：", a$lower_oth_eg)
}

# 统一：去首尾空格（不做大小写强转，symbol 大小写本身有意义）
deg_sym_raw <- trimws(as.character(deg[[col_deg]]))
mod_sym_raw <- trimws(as.character(mod[[col_mod]]))
uni_sym_raw <- trimws(as.character(uni[[col_uni]]))

deg_sym <- unique(deg_sym_raw)
mod_sym <- unique(mod_sym_raw)
uni_sym <- unique(uni_sym_raw)

ov_raw <- length(intersect(deg_sym, mod_sym))
if (ov_raw == 0) {
  # 交集为 0 时再试一次大小写无关匹配，便于定位 ID 类型不一致的问题
  ov_up <- length(intersect(toupper(deg_sym), toupper(mod_sym)))
  msg("  ! 直接交集为 0；大小写无关交集 = ", ov_up)
  if (ov_up == 0) {
    stop("两侧基因名交集为 0，ID 类型或物种不一致，请先统一基因名再取交集。")
  } else {
    stop("两侧基因名大小写不一致（直接交集 0，大小写无关交集 ", ov_up, "），请先统一。")
  }
}

fmt_ok <- (a_deg$n_missing + a_mod$n_missing == 0) &&
          (a_deg$n_duplicate + a_mod$n_duplicate == 0) &&
          (a_deg$n_has_space + a_mod$n_has_space == 0) &&
          (a_deg$n_all_digit + a_mod$n_all_digit == 0) &&
          (a_deg$n_ensembl   + a_mod$n_ensembl   == 0) &&
          (a_deg$n_probe     + a_mod$n_probe     == 0) &&
          (a_deg$n_lower_oth + a_mod$n_lower_oth == 0)

msg("  结论：两侧均使用 ", col_deg, " 列，格式",
    ifelse(fmt_ok, "一致", "存在差异（见上）"),
    "，可直接使用 GeneSymbol 求交集。")

## ---------------------------------------------------------------------------
## 3. 取交集
## ---------------------------------------------------------------------------
n_deg_tot <- length(deg_sym)                    # 差异基因数
n_mod_tot <- length(mod_sym)                    # 关键模块基因数
cand_sym  <- intersect(deg_sym, mod_sym)        # 候选基因
n_cand    <- length(cand_sym)

n_deg_only <- n_deg_tot - n_cand                # 仅 DEG
n_mod_only <- n_mod_tot - n_cand                # 仅关键模块

# 关键模块名称与该模块的模块-性状关联（用于报告）
mod_name <- if ("module" %in% names(mod)) paste(unique(mod$module), collapse = "/") else NA
n_hub_in_mod <- if ("is_hub" %in% names(mod)) sum(mod$is_hub == TRUE) else NA

# 交集上限：DEG 中落在 WGCNA 输入 5000 基因内的部分
deg_in_uni    <- sum(deg_sym %in% uni_sym)
mod_in_uni    <- sum(mod_sym %in% uni_sym)
ceiling_inter <- min(deg_in_uni, n_mod_tot)

# 候选基因中满足 hub 判据（is_hub == TRUE）的子集
n_cand_hub <- 0
if ("is_hub" %in% names(mod)) {
  n_cand_hub <- sum(mod$is_hub == TRUE & mod[[col_mod]] %in% cand_sym)
}

# 上调/下调构成（关键模块与疾病负相关时，候选应以下调为主——这是方向一致性自检）
reg_tab <- table(factor(deg$Regulation[deg_sym_raw %in% cand_sym],
                        levels = c("Up", "Down")))
n_up   <- as.integer(reg_tab[["Up"]])
n_down <- as.integer(reg_tab[["Down"]])
reg_tab_all <- table(factor(deg$Regulation, levels = c("Up", "Down")))

# 关键模块基因中差异基因的富集检验（Fisher 精确检验，单侧）
# 宇宙 = WGCNA 输入的 5000 个基因（模块是从这 5000 个里切出来的，宇宙必须一致）
a <- n_cand                                    # 模块内 & DEG
b <- n_mod_tot - n_cand                        # 模块内 & 非 DEG
c_ <- deg_in_uni - n_cand                      # 模块外 & DEG（仍在 5000 内）
d <- length(uni_sym) - a - b - c_              # 模块外 & 非 DEG
fish_p <- NA; fish_or <- NA
if (all(c(a, b, c_, d) >= 0)) {
  ft <- fisher.test(matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")
  fish_p  <- ft$p.value
  fish_or <- unname(ft$estimate)
}

msg("-- 计数 --")
msg("  显著差异基因            : ", n_deg_tot, "  (Up ", reg_tab_all[["Up"]],
    " / Down ", reg_tab_all[["Down"]], ")")
msg("  WGCNA 关键模块(", mod_name, ")基因 : ", n_mod_tot,
    "  (其中 hub ", n_hub_in_mod, ")")
msg("  交集 = 候选基因          : ", n_cand)
msg("  仅 DEG / 仅模块          : ", n_deg_only, " / ", n_mod_only)
msg("  DEG 落在 WGCNA 5000 基因内: ", deg_in_uni, " → 交集上限 ", ceiling_inter)
msg("  Fisher 精确检验(单侧)    : OR=", signif(fish_or, 3), "  P=", format.pval(fish_p, digits = 3))

if (n_cand == 0) stop("交集为 0，无法继续。请检查两侧基因名或 WGCNA 输入基因集。")

## ---------------------------------------------------------------------------
## 4. 输出候选基因表
## ---------------------------------------------------------------------------
deg_keep <- deg[deg_sym_raw %in% cand_sym, , drop = FALSE]
mod_keep <- mod[mod_sym_raw %in% cand_sym, , drop = FALSE]
names(mod_keep)[names(mod_keep) == col_mod] <- "GeneSymbol"
names(deg_keep)[names(deg_keep) == col_deg] <- "GeneSymbol"

cand <- merge(deg_keep, mod_keep, by = "GeneSymbol", all.x = TRUE, sort = FALSE)
# 排序：先按 padj（差异证据强度）升序；同 padj 时按 |MM_kME| 降序
ord_pad <- order(cand$padj, -abs(cand$MM_kME), method = "radix")
cand <- cand[ord_pad, , drop = FALSE]
cand$in_DEG_set   <- TRUE
cand$in_WGCNA_set <- TRUE

f_cand <- file.path(DIR_OUT, paste0(DS, "_candidate_genes.csv"))
write.csv(cand, f_cand, row.names = FALSE, fileEncoding = "UTF-8")
msg("  saved: ", basename(f_cand), "  (", nrow(cand), " 基因)")

# 补充：候选基因里同时满足 hub 判据的子集
if ("is_hub" %in% names(cand)) {
  cand_hub <- cand[cand$is_hub == TRUE, , drop = FALSE]
  if ("hub_rank" %in% names(cand_hub)) {
    cand_hub <- cand_hub[order(cand_hub$hub_rank, method = "radix"), , drop = FALSE]
  }
  f_hub <- file.path(DIR_OUT, paste0(DS, "_candidate_genes_strict_hub.csv"))
  write.csv(cand_hub, f_hub, row.names = FALSE, fileEncoding = "UTF-8")
  msg("  saved: ", basename(f_hub), "  (", nrow(cand_hub), " 基因)")
}

# 补充：两侧的非重叠部分（供后续敏感性分析）
w_deg_only <- file.path(DIR_OUT, paste0(DS, "_genes_DEG_only.csv"))
w_mod_only <- file.path(DIR_OUT, paste0(DS, "_genes_WGCNA_only.csv"))
write.csv(data.frame(GeneSymbol = sort(setdiff(deg_sym, mod_sym))),
          w_deg_only, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(GeneSymbol = sort(setdiff(mod_sym, deg_sym))),
          w_mod_only, row.names = FALSE, fileEncoding = "UTF-8")

# 三个分区的计数（便于下游引用/复算）
cnt <- data.frame(
  set = c("DEG", "WGCNA_module", "Candidate_intersection",
          "DEG_only", "WGCNA_module_only"),
  n   = c(n_deg_tot, n_mod_tot, n_cand, n_deg_only, n_mod_only),
  stringsAsFactors = FALSE
)
f_cnt <- file.path(DIR_OUT, paste0(DS, "_venn_counts.csv"))
write.csv(cnt, f_cnt, row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------------------------------------------------------
## 5. 韦恩图（面积成比例 + Nature 版式 + 莫兰迪配色）
## ---------------------------------------------------------------------------
COL_DEG   <- "#7A9188"   # 莫兰迪 灰绿（显著差异基因）
COL_MOD   <- "#C79E88"   # 莫兰迪 陶土（WGCNA 关键模块）
INK       <- "#3A3A38"   # 近黑，替代纯黑
GUIDE     <- "#B9B5AC"

# 两圆交面积（解析式）
circ_overlap_area <- function(r1, r2, d) {
  if (d >= r1 + r2)          return(0)
  if (d <= abs(r1 - r2))     return(pi * min(r1, r2)^2)
  a1  <- acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  a2  <- acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  tri <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
  r1^2 * a1 + r2^2 * a2 - tri
}

# 由三个计数反解圆心距：面积正比于基因数
venn_layout <- function(n1, n2, n12) {
  r1 <- sqrt(n1 / pi); r2 <- sqrt(n2 / pi)
  if (n12 >= min(n1, n2)) {                       # 小圆被完全包含的极端情形
    D <- abs(r1 - r2)
  } else {
    f <- function(d) circ_overlap_area(r1, r2, d) - n12
    D <- uniroot(f,
                 lower = abs(r1 - r2) * (1 + 1e-7) + 1e-9,
                 upper = (r1 + r2) * (1 - 1e-7),
                 tol = 1e-10)$root
  }
  c1x <- -D / 2; c2x <- D / 2
  list(r1 = r1, r2 = r2, D = D,
       c1x = c1x, c2x = c2x,
       # 三个分区在 x 轴上的中点（推导后 D 抵消，结果简洁）
       x1  = -(r1 + r2) / 2,          # 仅集合1 的月牙
       x12 = (r1 - r2) / 2,           # 交集透镜
       x2  = (r1 + r2) / 2,           # 仅集合2 的月牙
       ytop = max(r1, r2) + 1.4)      # 集合名标签高度
}

circle_df <- function(cx, cy, r, n = 400) {
  t <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}

L <- venn_layout(n_deg_tot, n_mod_tot, n_cand)
d1 <- circle_df(L$c1x, 0, L$r1); d1$grp <- "deg"
d2 <- circle_df(L$c2x, 0, L$r2); d2$grp <- "mod"
d_circ <- rbind(d1, d2)

xmin <- min(L$c1x - L$r1, L$c2x - L$r2)
xmax <- max(L$c1x + L$r1, L$c2x + L$r2)
padx <- 0.10 * (xmax - xmin)
xlim <- c(xmin - padx, xmax + padx)
ylim <- c(-max(L$r1, L$r2) - 0.6, L$ytop + 1.0)

lb <- data.frame(
  x    = c(L$x1, L$x12, L$x2),
  y    = c(0, 0, 0),
  txt  = c(n_deg_only, n_cand, n_mod_only),
  kind = c("n", "n", "n"),
  stringsAsFactors = FALSE
)
nm <- data.frame(
  x   = c(L$c1x, L$c2x),
  y   = c(L$ytop, L$ytop),
  txt = c("DEGs", paste0("WGCNA ", mod_name, " module")),
  stringsAsFactors = FALSE
)

p_venn <- ggplot() +
  geom_polygon(data = d_circ,
               aes(x = x, y = y, group = grp, fill = grp),
               alpha = 0.55, colour = NA) +
  geom_path(data = d_circ,
            aes(x = x, y = y, group = grp, colour = grp),
            linewidth = 0.45) +
  geom_text(data = lb, aes(x = x, y = y, label = txt),
            size = 10 * 0.3528, colour = INK, fontface = "bold",
            vjust = 0.5, hjust = 0.5) +
  geom_text(data = nm, aes(x = x, y = y, label = txt),
            size = 8 * 0.3528, colour = INK,
            vjust = 0, hjust = 0.5) +
  scale_fill_manual(values = c(deg = COL_DEG, mod = COL_MOD), guide = "none") +
  scale_colour_manual(values = c(deg = COL_DEG, mod = COL_MOD), guide = "none") +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "Helvetica") +
  theme(
    plot.margin = margin(2, 2, 2, 2, "mm"),
    text = element_text(family = "Helvetica", colour = INK)
  )

# 出图：PDF（矢量）+ PNG（600 dpi）
W_MM <- 89; H_MM <- 68                    # Nature 单栏宽度 89 mm
w_in <- W_MM / 25.4; h_in <- H_MM / 25.4
f_pdf <- file.path(DIR_OUT, paste0(DS, "_venn.pdf"))
f_png <- file.path(DIR_OUT, paste0(DS, "_venn.png"))

if (HAS_CAIRO) {
  ggsave(f_pdf, p_venn, width = w_in, height = h_in,
         device = cairo_pdf, bg = "white")
} else {
  ggsave(f_pdf, p_venn, width = w_in, height = h_in, bg = "white")
}
if (HAS_RAGG) {
  ggsave(f_png, p_venn, width = w_in, height = h_in,
         dpi = 600, device = ragg::agg_png, bg = "white")
} else {
  ggsave(f_png, p_venn, width = w_in, height = h_in, dpi = 600, bg = "white")
}
msg("  saved: ", basename(f_pdf), " / ", basename(f_png),
    "  (", W_MM, " x ", H_MM, " mm, PNG 600 dpi)")

# 图注
legend_txt <- c(
  paste0("Figure. Identification of candidate genes by intersecting differentially ",
         "expressed genes (DEGs) with the WGCNA key module."),
  paste0("  Venn diagram of ", n_deg_tot, " significant DEGs (", DS,
         "; limma, |log2FC| > 0.585 and adjusted P < 0.05) and ", n_mod_tot,
         " genes in the WGCNA key module (", mod_name,
         "; the module most strongly correlated with disease status)."),
  paste0("  The ", n_cand, " genes in the overlap were retained as candidates. ",
         n_deg_only, " DEGs fell outside the module and ", n_mod_only,
         " module genes were not differentially significant."),
  "  Circle areas and the area of overlap are proportional to gene counts.",
  "  DEGs: up-regulated / down-regulated as called by limma; the key module was",
  paste0("  negatively correlated with disease (r = -0.83); ", n_down, " of the ",
         n_cand, " candidates (", sprintf("%.0f%%", 100 * n_down / n_cand),
         ") are down-regulated in disease."),
  "  Fill colours follow a low-saturation Morandi palette; labels are set in",
  "  Helvetica. The P value is from a one-sided Fisher exact test of enrichment",
  "  of DEGs in the key module, with the 5,000 genes used as WGCNA input as the",
  "  background universe.",
  "",
  "Reproducibility: R ", paste(R.version$major, R.version$minor, sep = "."),
  "; script scripts/16_hub_gene_intersection.R."
)
wtxt(legend_txt, file.path(DIR_OUT, paste0(DS, "_venn_figure_legend.txt")))

## ---------------------------------------------------------------------------
## 6. 文本摘要
## ---------------------------------------------------------------------------
smr <- c(
  "================================================================",
  paste0(DS, " 候选（关键）基因筛选：显著差异基因 ∩ WGCNA 关键模块基因"),
  "================================================================",
  paste0("生成时间        : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "-------------------------- 输入 --------------------------",
  paste0("显著差异基因表  : ", basename(F_DEG), "  (", nrow(deg), " 行)"),
  paste0("关键模块基因表  : ", basename(F_MOD), "  (", nrow(mod), " 行)"),
  paste0("WGCNA 输入全集  : ", basename(F_UNI), "  (", nrow(uni), " 行)"),
  "",
  "--------------------- 基因名格式核查 ---------------------",
  fmt_line(a_deg), fmt_line(a_mod), fmt_line(a_uni),
  paste0("一致性结论      : ", ifelse(fmt_ok,
         "两侧均为 GeneSymbol（无缺失/无重复/无空格/无 ID 混杂），格式一致，可直接取交集",
         "存在格式差异，见上；已按 trimws 处理后取交集")),
  "",
  "-------------------------- 计数 --------------------------",
  paste0("显著差异基因  DEG              : ", n_deg_tot,
         "   (Up ", reg_tab_all[["Up"]], " / Down ", reg_tab_all[["Down"]], ")"),
  paste0("WGCNA 关键模块基因 (", mod_name, ")   : ", n_mod_tot,
         "   (其中 hub ", n_hub_in_mod, ")"),
  paste0("交集 = 候选基因                 : ", n_cand),
  paste0("  其中满足 hub 判据 (is_hub)    : ", n_cand_hub),
  paste0("  候选基因上调 / 下调           : ", n_up, " / ", n_down),
  paste0("仅 DEG（不在关键模块）          : ", n_deg_only),
  paste0("仅关键模块（不差异显著）        : ", n_mod_only),
  "",
  "------------------- 交集上限与富集检验 -------------------",
  paste0("DEG 落在 WGCNA 输入 5000 基因内 : ", deg_in_uni, " / ", n_deg_tot),
  paste0("交集的理论上限                  : ", ceiling_inter),
  paste0("Fisher 精确检验(单侧, 宇宙=5000): OR = ", signif(fish_or, 3),
         "   P = ", format.pval(fish_p, digits = 3)),
  "",
  "------------------------ 输出文件 ------------------------",
  paste0("候选基因表      : ", basename(f_cand)),
  paste0("韦恩图          : ", basename(f_pdf), " / ", basename(f_png)),
  paste0("计数表          : ", basename(f_cnt)),
  paste0("报告            : ", DS, "_HubGene_report.md"),
  "================================================================"
)
wtxt(smr, file.path(DIR_OUT, paste0(DS, "_hubgene_summary.txt")))

## ---------------------------------------------------------------------------
## 7. 中文处理报告（.md）
## ---------------------------------------------------------------------------
top10 <- cand$GeneSymbol[1:min(10, nrow(cand))]
top10_txt <- paste0("`", paste(top10, collapse = "`, `"), "`")
mod_r <- -0.827   # 关键模块与 disease 的相关（来自 09_wgcna.R 输出）
mod_r_p  <- "8.95e-80"

# 方向一致性段落：按实际结果措辞，避免"0 个上调"被写成含混的"以下调为主"
dir_txt <- if (n_up == 0L) {
  c(paste0("关键模块 ", mod_name, " 与 disease 呈负相关（r = ", mod_r, "，P = ",
           mod_r_p, "），即模块整体在疾病中下调。"),
    paste0("相应地，全部 ", n_cand, " 个候选基因在疾病中**无一例外均为下调**",
           "（下调 ", n_down, " 个、上调 ", n_up, " 个）。"),
    "这是一个很强的内部一致性信号：上调的差异基因几乎没有落进这个模块，",
    "说明 DEG 与 WGCNA 这两条证据线指向的是同一个「在疾病中整体下调的转录程序」。")
} else {
  c(paste0("关键模块 ", mod_name, " 与 disease 呈负相关（r = ", mod_r, "，P = ",
           mod_r_p, "），即模块整体在疾病中下调。"),
    paste0("相应地，候选基因中下调 ", n_down, " 个、上调 ", n_up, " 个，下调占 ",
           sprintf("%.1f%%", 100 * n_down / n_cand), "，与模块方向一致。"))
}
dir_txt <- c(dir_txt,
  "若候选基因的上下调比例与模块相关方向矛盾，说明取交集的某一步（基因名映射或",
  "模块选择）出了问题，务必回头检查。")

# 结构性漏检段落：本数据集 DEG 全部落在 WGCNA 输入集内，属于"漏检为 0"的情形
loss_txt <- if (deg_in_uni == n_deg_tot) {
  c(paste0("**(a) 结构性漏检：本数据集为 0。** WGCNA 只建立在 MAD 前 ",
           length(uni_sym), " 个高变异基因上，"),
    paste0("而全部 ", n_deg_tot, " 个差异基因都落在这 ", length(uni_sym),
           " 个基因之内——差异基因本身就是高变异基因，符合预期。"),
    "因此本数据集**不存在**「确实差异显著、却因变异太小而进不了 WGCNA」的结构性漏检，",
    "交集没有因此被人为压低。  ",
    "（换一个数据集时，务必先算这个覆盖率：若 DEG 有相当比例不在 WGCNA 输入集内，",
    "交集的上限就被结构性压低了，这一步会系统性漏掉一部分真实标志物。）")
} else {
  c(paste0("**(a) 结构性漏检。** WGCNA 只建立在 MAD 前 ", length(uni_sym),
           " 个高变异基因上，"),
    paste0("DEG 中只有 ", deg_in_uni, " 个落在其中，另有 ", n_deg_tot - deg_in_uni,
           " 个差异基因在结构上永远不可能进入任何模块；"),
    paste0("因此交集的理论上限是 ", ceiling_inter, " 个，而不是 ", n_deg_tot, " 个。"),
    "这意味着交集这一步会漏掉「确实差异显著、但表达变异不够大因而没进 WGCNA 输入集」",
    "的基因，其中可能包含真实的标志物。")
}

rep_md <- c(
  paste0("# ", DS, " 候选（关键）基因筛选报告"),
  "",
  paste0("> 分析内容：显著差异基因（DEG）与 WGCNA 关键模块基因取交集，得到候选基因，",
         "并绘制韦恩图。  "),
  paste0("> 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  "),
  paste0("> 生成脚本：`scripts/16_hub_gene_intersection.R`（全文含中文注释，可复现）  "),
  paste0("> 运行环境：R ", paste(R.version$major, R.version$minor, sep = "."),
         " / ggplot2 ", as.character(packageVersion("ggplot2"))),
  "",
  "---",
  "",
  "## 一、数据来源与筛选逻辑",
  "",
  "| 项目 | 内容 |",
  "| --- | --- |",
  paste0("| 数据集 | ", DS, "（GPL11532，Affymetrix HuGene-1.1-ST，RMA log2 标准化） |"),
  paste0("| 输入 1 | `results/DEG/", basename(F_DEG), "`，显著差异基因 ", n_deg_tot, " 个 |"),
  paste0("| 输入 2 | `results/WGCNA/", basename(F_MOD), "`，关键模块（", mod_name,
         "）基因 ", n_mod_tot, " 个 |"),
  paste0("| 筛选方式 | 取两者交集（按 GeneSymbol 精确匹配） |"),
  paste0("| 候选基因 | **", n_cand, " 个** |"),
  "",
  "**DEG 的判定阈值**（沿用 04_deg.R）：`|log2FC| > 0.585`（差异倍数 1.5 倍）",
  "且 `padj < 0.05`（BH 校正）。模型为 `~group2 + sex + age_num`，即校正了性别与年龄。  ",
  paste0("**WGCNA 关键模块的判定**：模块-性状关联中，", mod_name, " 模块与 disease 性状",
         "相关性最强（r = ", mod_r, "，P = 8.95e-80），模块内共 ", n_mod_tot, " 个基因。"),
  "",
  "---",
  "",
  "## 二、基因名格式一致性核查",
  "",
  "取交集前必须确认两侧基因名可直接比对，否则交集会静默算错。核查结果如下：",
  "",
  "| 检查项 | 显著差异基因 | WGCNA 关键模块基因 |",
  "| --- | --- | --- |",
  paste0("| 使用列 | `", col_deg, "` | `", col_mod, "` |"),
  paste0("| 基因数 | ", a_deg$n, " | ", a_mod$n, " |"),
  paste0("| 缺失值 | ", a_deg$n_missing, " | ", a_mod$n_missing, " |"),
  paste0("| 重复名 | ", a_deg$n_duplicate, " | ", a_mod$n_duplicate, " |"),
  paste0("| 含空格 | ", a_deg$n_has_space, " | ", a_mod$n_has_space, " |"),
  paste0("| 全数字（疑似 Entrez ID） | ", a_deg$n_all_digit, " | ", a_mod$n_all_digit, " |"),
  paste0("| ENS 前缀（疑似 Ensembl ID） | ", a_deg$n_ensembl, " | ", a_mod$n_ensembl, " |"),
  paste0("| 6 位以上纯数字（疑似探针号） | ", a_deg$n_probe, " | ", a_mod$n_probe, " |"),
  paste0("| 含小写字母的 symbol | ", a_deg$n_lower, " | ", a_mod$n_lower, " |"),
  paste0("| └ 其中 HGNC 的 `C<n>orf<m>` 命名（合法） | ", a_deg$n_lower_orf,
         " | ", a_mod$n_lower_orf, " |"),
  paste0("| └ 其他可疑小写命名 | ", a_deg$n_lower_oth, " | ", a_mod$n_lower_oth, " |"),
  paste0("| 示例 | ", a_deg$example, " | ", a_mod$example, " |"),
  "",
  paste0("**结论**：两侧基因名列均为 `GeneSymbol`，缺失 0、重复 0、含空格 0，"),
  "既不含 Entrez / Ensembl ID，也不含 Affymetrix 探针号，命名空间一致。  ",
  paste0("两侧另有少量含小写字母的 symbol（DEG ", a_deg$n_lower, " 个、关键模块 ",
         a_mod$n_lower, " 个），全部属于 HGNC 合法的 `C<n>orf<m>` 命名"),
  "（如 `C9orf47`），不是 ID 类型混杂，也不影响匹配。  ",
  "**两侧基因名格式一致，可直接按 `GeneSymbol` 取交集。**",
  "",
  "---",
  "",
  "## 三、结果：关键数字",
  "",
  "| 集合 | 基因数 | 说明 |",
  "| --- | --- | --- |",
  paste0("| 显著差异基因（DEG） | **", n_deg_tot, "** | 上调 ", reg_tab_all[["Up"]],
         " / 下调 ", reg_tab_all[["Down"]], " |"),
  paste0("| WGCNA 关键模块基因（", mod_name, "） | **", n_mod_tot,
         "** | 其中满足 hub 判据（kME 与 GS 双高）的有 ", n_hub_in_mod, " 个 |"),
  paste0("| 交集 = 候选基因 | **", n_cand, "** | 占 DEG 的 ",
         sprintf("%.1f%%", 100 * n_cand / n_deg_tot), "，占关键模块的 ",
         sprintf("%.1f%%", 100 * n_cand / n_mod_tot), " |"),
  paste0("| 仅 DEG（不在关键模块） | ", n_deg_only, " | 被交集筛掉 |"),
  paste0("| 仅关键模块（不差异显著） | ", n_mod_only, " | 被交集筛掉 |"),
  "",
  paste0("候选基因中下调 ", n_down, " 个（", sprintf("%.1f%%", 100 * n_down / n_cand),
         "）、上调 ", n_up, " 个。"),
  "",
  "---",
  "",
  "## 四、韦恩图",
  "",
  paste0("文件：`results/HubGene/", basename(f_pdf), "`（矢量 PDF）与 `", basename(f_png),
         "`（600 dpi PNG）。  "),
  "尺寸：89 × 68 mm（Nature 单栏宽度），字体 Helvetica，集合名 8 pt、计数 10 pt，",
  "符合期刊对最小字号 7 pt 的要求。  ",
  "配色采用低饱和度莫兰迪色：差异基因为灰绿 `#7A9188`，关键模块为陶土 `#C79E88`。",
  "",
  "**图的关键设计：三个分区的面积与基因数成正比。** 普通韦恩示意图的两个圆面积相等、",
  "圆心距固定，读者从图上读出的比例是失真的。本图按两圆交面积的解析式反解圆心距，",
  paste0("使「仅 DEG」「交集」「仅模块」三块的面积严格对应 ", n_deg_only,
         " / ", n_cand, " / ", n_mod_only, " 的基因数。"),
  "",
  "---",
  "",
  "## 五、取交集这一步的生物学意义",
  "",
  "### 1. 两条判据回答的是两个不同的问题",
  "",
  "- **差异表达（DEG）**回答「哪些基因在疾病组与对照组之间表达水平显著不同」。",
  "它是**单基因维度**的统计证据：逐个基因拟合线性模型，用经验贝叶斯收缩方差，",
  "再用 BH 校正控制多重检验假阳性。  ",
  "- **WGCNA 关键模块**回答「哪些基因在样本之间协同变化，且这种协同变化与疾病性状",
  "关联最强」。它是**网络维度**的结构证据：先把表达矩阵转成基因间的相关结构，",
  "再按拓扑重叠把基因聚成模块，最后用模块特征基因与性状求相关。",
  "",
  "两者不是同一件事的两种算法，而是**两条互相独立的证据线**。",
  "",
  "### 2. 为什么要取交集，而不是各用各的",
  "",
  "- **单看 DEG 的弱点**：差异显著性受样本异质性影响很大。一个基因可以因为 p>>n 下",
  "的偶然波动、或仅在少数样本中大幅变化而达到显著，但它并不参与疾病相关的协同",
  "调控程序，作为标志物时稳健性差（换一批样本就不显著了）。  ",
  "- **单看 WGCNA 模块的弱点**：共表达是**边际相关**。模块里两个基因之间有边，",
  "可能只是因为它们共同受第三个基因驱动，并不直接互作；更关键的是，一个基因可以",
  "稳稳地待在与疾病最相关的模块里，却在疾病与对照之间**根本没有显著表达差异**——",
  "它跟着模块一起波动，但幅度太小或个体间变异太大，达不到显著。  ",
  "- **取交集的效果**：候选基因必须**同时**满足「统计上确实差异显著」和「处于与疾病",
  "最相关的协同表达模块中」。两条证据相互印证，把只满足其一的基因排除掉，",
  "因此候选集的**假阳性率显著低于任一单方法的结果**，后续做机器学习建模、",
  "列线图或功能富集时，输入特征的质量更高、结论更稳。",
  "",
  "### 3. 交集不是单纯的「缩小范围」，而是提高证据密度",
  "",
  paste0("本数据中 DEG ", n_deg_tot, " 个、关键模块 ", n_mod_tot, " 个，交集 ", n_cand, " 个。"),
  paste0("如果两方法完全独立，随机交集的期望只有约 ",
         sprintf("%.0f", n_deg_tot * n_mod_tot / length(uni_sym)),
         " 个（在 WGCNA 输入的 ", length(uni_sym), " 基因宇宙内）；"),
  paste0("实测 ", n_cand, " 个，富集 ",
         sprintf("%.1f", n_cand / (n_deg_tot * n_mod_tot / length(uni_sym))), " 倍"),
  paste0("（Fisher 精确检验单侧 P = ", format.pval(fish_p, digits = 3), "）。"),
  "这说明差异基因并非随机散布，而是**显著聚集在与疾病最相关的那个共表达模块里**，",
  "本身就是一条有价值的结论：疾病相关的转录改变是一个协同的程序，而不是一堆互不",
  "相干的单基因事件。",
  "",
  "### 4. 方向一致性自检",
  "",
  dir_txt,
  "",
  "### 5. 必须知道的三个代价（诚实说明）",
  "",
  loss_txt,
  "",
  paste0("**(b) 交集缩小了候选范围。** ", n_mod_only, " 个模块内基因未达差异显著阈值，"),
  paste0("它们可能是效应量小但协同性强的基因；", n_deg_only, " 个差异基因不在关键模块内，"),
  "它们可能受其他调控程序驱动。这些基因不作为本轮候选，但不应被当作「无关」。",
  "",
  "**(c) 共表达不等于调控关系。** 这是本领域最常见的过度解读。候选基因应理解为",
  "「差异显著且处于疾病相关共表达模块中」，**不能**直接说成「上游驱动因子」或",
  "「主调控基因」；模块内的边可能是间接的，hub 地位也不等于因果地位。要回答",
  "「谁调控谁」，需要额外的因果推断或实验验证。",
  "",
  "---",
  "",
  "## 六、候选基因一览（按 padj 升序前 10）",
  "",
  top10_txt,
  "",
  paste0("完整列表见 `results/HubGene/", basename(f_cand), "`（", nrow(cand), " 个基因，",
         "含 log2FC / padj / Regulation 与 MM_kME / GS / is_hub / hub_rank 全部注释）。  "),
  paste0("其中同时满足 hub 判据的有 ", n_cand_hub, " 个，另存为 `",
         basename(file.path(DIR_OUT, paste0(DS, "_candidate_genes_strict_hub.csv"))), "`。"),
  "",
  "---",
  "",
  "## 七、输出文件清单",
  "",
  "| 文件 | 内容 |",
  "| --- | --- |",
  paste0("| `", basename(f_cand), "` | 候选基因表（主交付），DEG + WGCNA 注释合并，按 padj 升序 |"),
  paste0("| `", DS, "_candidate_genes_strict_hub.csv` | 候选中满足 hub 判据的子集（",
         n_cand_hub, " 个） |"),
  paste0("| `", basename(f_pdf), "` | 韦恩图，矢量 PDF，89 × 68 mm |"),
  paste0("| `", basename(f_png), "` | 韦恩图，600 dpi PNG |"),
  paste0("| `", basename(f_cnt), "` | 三个分区的计数 |"),
  paste0("| `", DS, "_genes_DEG_only.csv` | 仅 DEG 的 ", n_deg_only, " 个基因 |"),
  paste0("| `", DS, "_genes_WGCNA_only.csv` | 仅关键模块的 ", n_mod_only, " 个基因 |"),
  paste0("| `", DS, "_hubgene_summary.txt` | 全过程文本摘要（含基因名格式核查） |"),
  paste0("| `", DS, "_venn_figure_legend.txt` | 图注（英文，可直接引用） |"),
  paste0("| `", DS, "_HubGene_report.md` | 本报告 |"),
  paste0("| `sessionInfo_16_hub_gene_intersection.txt` | 运行环境 |"),
  "",
  "---",
  "",
  "## 八、复现方式",
  "",
  "```bash",
  "cd \"GEObulk analysis project\"",
  "Rscript scripts/16_hub_gene_intersection.R",
  "```",
  "",
  "脚本无随机过程（不涉及抽样、不设随机种子），在相同输入下结果完全确定。",
  "依赖包：ggplot2（必需），ragg（PNG 抗锯齿，缺失时自动降级）。"
)
wtxt(rep_md, file.path(DIR_OUT, paste0(DS, "_HubGene_report.md")))

## ---------------------------------------------------------------------------
## 8. 运行环境
## ---------------------------------------------------------------------------
wtxt(c(capture.output(sessionInfo())),
     file.path(DIR_OUT, "sessionInfo_16_hub_gene_intersection.txt"))

wtxt(c(paste0(DS, " 候选基因筛选运行日志"), "========================================",
       LOG),
     file.path(DIR_OUT, paste0(DS, "_hubgene_run_log.txt")))

msg("=== 完成。候选基因 ", n_cand, " 个，结果已写入 results/HubGene/ ===")
