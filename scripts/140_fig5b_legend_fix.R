## Redraw the Figure 5B calibration panel with spelled-out legend labels.
## The published panel came from rms::plot.calibrate, whose built-in legend
## abbreviates confidence limits as "C.L."; reviewers read that as an
## unexplained abbreviation. Everything else is kept identical.

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(rms) })

root <- "E:/another data file/GEObulk analysis project"
rds  <- file.path(root, "results/Nomogram/03_models/nomogram_models.rds")
outd <- file.path(root, "results/Nomogram/02_figures")

x <- readRDS(rds)
cal <- x$calibration_models
cat("calibration_models:", paste(names(cal), collapse = ", "), "\n")

get_cal <- function() {
  for (nm in names(cal)) if (grepl("M2|gene", nm, ignore.case = TRUE)) return(cal[[nm]])
  NULL
}
cal_m2 <- get_cal()
if (is.null(cal_m2)) { cat("M2 calibration object not found\n"); quit(status = 1) }
cat("using:", class(cal_m2)[1], "\n")

lab <- c("Apparent", "Bias-corrected", "Ideal", "95% confidence limits")
lty <- c(3, 1, 2, 1)
col <- c("black", "black", "grey40", "grey75")

draw <- function() {
  plot(cal_m2, legend = FALSE, subtitles = FALSE,
       xlab = "Predicted probability", ylab = "Observed proportion",
       main = "Calibration (key gene model, bootstrap bias-corrected)")
  abline(0, 1, lty = 2, col = "grey40")
  mtext("B = 1000 bootstrap repetitions; mean absolute error = 0.022; n = 219",
        side = 1, line = 3.7, cex = 0.72, col = "grey25")
  legend("bottomright", legend = lab, lty = lty, col = col, lwd = 1,
         bty = "n", cex = 0.92, seg.len = 2.4)
}

png(file.path(outd, "F07_calibration_train_bootstrap_M2_v2.png"),
    width = 6.5, height = 6, units = "in", res = 600)
draw(); dev.off()

pdf(file.path(outd, "F07_calibration_train_bootstrap_M2_v2.pdf"),
    width = 6.5, height = 6)
draw(); dev.off()

cat("written v2 panel\n")
