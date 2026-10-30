# 123_composite_fig2_3_4.R
# Build single-file composites for Figure 2 (DEG), Figure 3 (enrichment) and
# Figure 4 (co-expression + feature selection) from existing 600 dpi PNG panels.
# Layout: 2 x 2 grid, panel labels A-D, final width 180 mm at 600 dpi.
# ASCII only.

suppressPackageStartupMessages(library(magick))

root <- Sys.getenv("GEOBULK_ROOT", unset = "E:/another data file/GEObulk analysis project")
r <- function(...) file.path(root, "results", ...)
outd <- r("Manuscript", "figures")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

panel <- function(f, h = 1600) {
  img <- image_read(f)
  img <- image_trim(img)
  img <- image_flatten(img)
  img <- image_resize(img, paste0("x", h))
  image_border(img, "white", "60x60")
}
label <- function(img, lab) {
  tryCatch(image_annotate(img, lab, size = 110, gravity = "northwest",
                          color = "black", location = "+20+14"),
           error = function(e) img)
}
blank <- function(h, w) image_blank(w, h, color = "white")

build <- function(files, labs, out_stem) {
  ps <- lapply(files, panel)
  ps <- Map(label, ps, labs)
  # pad to 4 panels
  while (length(ps) < 4) {
    w <- image_info(ps[[1]])$width
    h <- image_info(ps[[1]])$height
    ps[[length(ps) + 1]] <- blank(h, w)
  }
  row1 <- image_append(image_join(ps[1:2]), stack = FALSE)
  row2 <- image_append(image_join(ps[3:4]), stack = FALSE)
  w <- min(image_info(row1)$width, image_info(row2)$width)
  row1 <- image_resize(row1, paste0(w, "x"))
  row2 <- image_resize(row2, paste0(w, "x"))
  fig <- image_append(c(row1, row2), stack = TRUE)
  fig <- image_background(fig, "white", flatten = TRUE)
  fig <- image_resize(fig, paste0(round(180 / 25.4 * 600), "x"))
  png <- file.path(outd, paste0(out_stem, ".png"))
  pdf <- file.path(outd, paste0(out_stem, ".pdf"))
  image_write(fig, png, format = "png", density = 600)
  image_write(fig, pdf, format = "pdf", density = 600)
  info <- image_info(fig)
  cat(out_stem, ":", info$width, "x", info$height, "px\n")
}

# Figure 2: DEG
build(
  c(r("DEG", "GSE57338_volcano.png"),
    r("DEG", "GSE57338_heatmap.png"),
    r("DEG", "GSE57338_QC_diagnostics.png")),
  c("A", "B", "C"),
  "Figure2_differential_expression"
)

# Figure 3: enrichment
build(
  c(r("Enrich", "GO_dotplot.png"),
    r("Enrich", "KEGG_dotplot.png"),
    r("GSEA", "GSEA_bubble_overview.png"),
    r("GSEA", "GSEA_runningscore_top5.png")),
  c("A", "B", "C", "D"),
  "Figure3_functional_enrichment"
)

# Figure 4: WGCNA + machine learning
build(
  c(r("WGCNA", "GSE57338_WGCNA_soft_threshold.png"),
    r("WGCNA", "GSE57338_WGCNA_module_trait_heatmap.png"),
    r("ML", "05_comparison", "GSE57338_ML_venn.png"),
    r("ML", "02_lasso", "GSE57338_LASSO_coef_barplot.png")),
  c("A", "B", "C", "D"),
  "Figure4_module_feature_selection"
)

cat("done\n")
