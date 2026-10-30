# 122_figure5_composite.R
# Assemble Figure 5 (model + external validation) with M2 as the primary model.
# Panels: A M2 nomogram, B M2 bootstrap calibration, C DCA, D external ROC,
#         E external calibration, F AUC forest.
# Output: results/Manuscript/figures/Figure5_model_validation.{pdf,png}
# ASCII only (R locale on this machine is not UTF-8).

suppressPackageStartupMessages(library(magick))

root <- Sys.getenv("GEOBULK_ROOT", unset = "E:/another data file/GEObulk analysis project")
nom  <- file.path(root, "results/Nomogram/02_figures")
ext  <- file.path(root, "results/ExternalValidation/02_figures")
outd <- file.path(root, "results/Manuscript/figures")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

files <- c(
  A = file.path(nom, "F02_nomogram_M2_gene.png"),
  B = file.path(nom, "F07_calibration_train_bootstrap_M2.png"),
  C = file.path(nom, "F10_DCA_validation.png"),
  D = file.path(ext, "F01_ROC_internal_vs_external.png"),
  E = file.path(ext, "F02_calibration_external.png"),
  F = file.path(ext, "F04_AUC_forest.png")
)
stopifnot(all(file.exists(files)))

# common height for each row's panels, keep aspect ratio
panel <- function(f, h) {
  img <- image_read(f)
  img <- image_trim(img)
  img <- image_resize(img, paste0("x", h))
  image_border(img, "white", "50x50")
}

row1 <- lapply(files[c("A", "B", "C")], panel, h = 1400)
row2 <- lapply(files[c("D", "E", "F")], panel, h = 1400)

add_label <- function(img, lab) {
  tryCatch(
    image_annotate(img, lab, size = 90, gravity = "northwest",
                   color = "black", location = "+18+12"),
    error = function(e) img
  )
}
row1 <- Map(add_label, row1, c("A", "B", "C"))
row2 <- Map(add_label, row2, c("D", "E", "F"))

r1 <- image_append(image_join(row1), stack = FALSE)
r2 <- image_append(image_join(row2), stack = FALSE)
w  <- min(image_info(r1)$width, image_info(r2)$width)
r1 <- image_resize(r1, paste0(w, "x"))
r2 <- image_resize(r2, paste0(w, "x"))
fig <- image_append(c(r1, r2), stack = TRUE)
fig <- image_background(fig, "white", flatten = TRUE)

# fix final physical width at the journal text width: 180 mm at 600 dpi = 4252 px
target_w <- round(180 / 25.4 * 600)
fig <- image_resize(fig, paste0(target_w, "x"))

png_out <- file.path(outd, "Figure5_model_validation.png")
pdf_out <- file.path(outd, "Figure5_model_validation.pdf")
image_write(fig, png_out, format = "png", density = 600)
image_write(fig, pdf_out, format = "pdf", density = 600)

info <- image_info(fig)
cat("panel sizes (px):\n")
print(image_info(fig)[, c("width", "height")])
cat("written:", png_out, "\n")
cat("written:", pdf_out, "\n")
