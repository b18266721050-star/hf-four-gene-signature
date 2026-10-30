## 128_locus_coloc_figure.R -- Figure S72: 9p13.3 colocalisation and instrument specificity
## A: cis-eQTL -log10 p for AQP3 and two neighbours
## B: heart-failure -log10 p in the same window
## C: PP.H4 for every gene x outcome pair (primary prior)
## D: cis-MR estimates with all, block-shared and AQP3-specific instruments
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

ROOT <- "E:/another data file/GEObulk analysis project"
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
TAB  <- file.path(ROOT, "results", "Coloc", "01_tables")
FIG  <- file.path(ROOT, "results", "Coloc", "02_figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

res <- read.csv(file.path(TAB, "T01_coloc_results.csv"), stringsAsFactors = FALSE)
prim <- res[res$sdY == 0.15 & res$p12 == 1e-5, ]

rd <- function(g) {
  f <- file.path(DAT, sprintf("coloc_%s_HERMES.csv", gsub("\\.", "_", g)))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$mlog10p_eqtl <- -pnorm(-abs(d$z1), log.p = TRUE) / log(10)
  d$mlog10p_hf   <- -pnorm(-abs(d$z2), log.p = TRUE) / log(10)
  d$mb <- d$pos / 1e6
  d
}
rd2 <- function(g, o) {
  f <- file.path(DAT, sprintf("coloc_%s_%s.csv", gsub("[.]", "_", g), o))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$mlog10p_hf <- -pnorm(-abs(d$z2), log.p = TRUE) / log(10)
  d$mb <- d$pos / 1e6
  d
}
aq <- rd("AQP3"); ln <- rd("RP11-176F3.7"); ga <- rd("GALT")
fq <- rd2("AQP3", "FinnGen")

pa <- ggplot() +
  geom_point(data = aq, aes(mb, mlog10p_eqtl, colour = "AQP3"), size = .9, alpha = .7) +
  { if (!is.null(ln)) geom_point(data = ln, aes(mb, mlog10p_eqtl, colour = "RP11-176F3.7"),
                                 size = .9, alpha = .7) } +
  { if (!is.null(ga)) geom_point(data = ga, aes(mb, mlog10p_eqtl, colour = "GALT"),
                                 size = .9, alpha = .7) } +
  scale_colour_manual(values = c(AQP3 = "#2E5E8C", "RP11-176F3.7" = "#C4622D", GALT = "#4E8A5A"),
                      name = "cis-eQTL") +
  coord_cartesian(xlim = c(32.9, 34.9), ylim = c(0, 120)) +
  labs(x = NULL, y = expression(-log[10]~italic(p)[eQTL]),
       title = "A  Expression (eQTLGen, blood)") +
  theme_bw(base_size = 8) +
  theme(legend.position = "right", legend.key.size = unit(3, "mm"),
        plot.title = element_text(face = "bold", size = 8))

pb <- ggplot() +
  geom_point(data = aq, aes(mb, mlog10p_hf, colour = "HERMES"), size = .9, alpha = .7) +
  { if (!is.null(fq)) geom_point(data = fq, aes(mb, mlog10p_hf, colour = "FinnGen"),
                                 size = .9, alpha = .7) } +
  geom_hline(yintercept = -log10(5e-8), linetype = 2, linewidth = .25, colour = "grey20") +
  scale_colour_manual(values = c(HERMES = "#2E5E8C", FinnGen = "#B03A2E"), name = "HF GWAS") +
  coord_cartesian(xlim = c(32.9, 34.9)) +
  labs(x = "Position on chromosome 9 (Mb)", y = expression(-log[10]~italic(p)[HF]),
       title = "B  Heart failure in the AQP3 window") +
  theme_bw(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 8))

prim$lab <- paste0(prim$gene, " | ", prim$outcome)
prim$lab <- factor(prim$lab, levels = prim$lab[order(prim$gene, prim$outcome)])
pc <- ggplot(prim, aes(x = PP.H4, y = lab, fill = outcome)) +
  geom_col(width = .68) +
  geom_vline(xintercept = .8, linetype = 2, linewidth = .25) +
  scale_fill_manual(values = c(HERMES = "#2E5E8C", FinnGen = "#B03A2E"), name = NULL) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = expression(PP.H[4]~"(shared causal variant)"), y = NULL,
       title = "C  Colocalisation with heart failure") +
  theme_bw(base_size = 8) +
  theme(legend.position = "right", legend.key.size = unit(3, "mm"),
        plot.title = element_text(face = "bold", size = 8))

fp <- read.csv(file.path(TAB, "T04_instrument_specificity.csv"), stringsAsFactors = FALSE)
fp <- fp[fp$n_snp > 0, ]
fp$row <- paste0(fp$outcome, ": ", fp$set, " (", fp$n_snp, " SNP", ifelse(fp$n_snp > 1, "s", ""), ")")
fp$row <- factor(fp$row, levels = rev(fp$row))
fp$OR <- as.numeric(fp$OR); fp$lo <- as.numeric(fp$lo); fp$hi <- as.numeric(fp$hi)
pd <- ggplot(fp, aes(x = OR, y = row, colour = outcome)) +
  geom_vline(xintercept = 1, linewidth = .25) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = .28, linewidth = .5) +
  scale_x_log10() +
  scale_colour_manual(values = c(HERMES = "#2E5E8C", FinnGen = "#B03A2E"), guide = "none") +
  labs(x = "OR of heart failure per 1-SD genetically predicted AQP3", y = NULL,
       title = "D  Which instruments carry the estimate") +
  theme_bw(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 8))

p <- (pa / pb) | (pc / pd) + plot_layout(widths = c(1.45, 1))
ggsave(file.path(FIG, "FigureS72_locus_coloc.pdf"), p,
       width = 180, height = 150, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG, "FigureS72_locus_coloc.png"), p,
       width = 180, height = 150, units = "mm", dpi = 600)
cat("written FigureS72_locus_coloc.{pdf,png}\n")
