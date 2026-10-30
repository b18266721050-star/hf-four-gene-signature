# =============================================================================
# 40_qc_plot_check.R  (v2 - memory-safe, incremental)
# =============================================================================

rm(list = ls())
proj_dir <- "E:/another data file/GEObulk analysis project"
plot_dir <- file.path(proj_dir, "results", "Singlecel", "qc_plots")
rep_dir <- file.path(proj_dir, "results", "Singlecel", "logs")

one <- function(f) {
  img <- tryCatch(png::readPNG(f), error = function(e) NULL)
  if (is.null(img)) {
    return(data.frame(file = basename(f), w = NA, h = NA, mean_lum = NA,
                      white_frac = NA, dark_frac = NA, ok = FALSE))
  }
  h <- dim(img)[1]; w <- dim(img)[2]
  # subsample for statistics only (keeps memory low)
  rs <- seq(1, h, by = max(1, floor(h / 300)))
  cs <- seq(1, w, by = max(1, floor(w / 500)))
  sub <- img[rs, cs, , drop = FALSE]
  lum <- if (length(dim(sub)) == 3) {
    0.299 * sub[, , 1] + 0.587 * sub[, , 2] + 0.114 * sub[, , 3]
  } else sub
  wf <- mean(lum > 0.97)
  df_ <- mean(lum < 0.7)
  rm(img, sub, lum); invisible(gc(FALSE))
  data.frame(file = basename(f), w = w, h = h,
             mean_lum = round(mean(c(wf, 1 - df_)), 3),
             white_frac = round(wf, 3), dark_frac = round(df_, 4),
             ok = (wf < 0.995 && df_ > 0.002), stringsAsFactors = FALSE)
}

for (g in c("violin", "scatter")) {
  fs <- sort(list.files(file.path(plot_dir, g), pattern = "[.]png$", full.names = TRUE))
  st <- data.frame()
  for (f in fs) {
    st <- rbind(st, one(f))
    cat(basename(f), "ok\n", file = file.path(rep_dir, "plot_check_progress.log"), append = TRUE)
  }
  st$type <- g
  assign(paste0("st_", g), st)
  cat("\n==", g, ":", nrow(st), "files\n")
  cat("   sizes:", paste(unique(paste0(st$w, "x", st$h)), collapse = ", "), "\n")
  cat("   white_frac:", sprintf("%.3f-%.3f", min(st$white_frac), max(st$white_frac)),
      "| dark_frac:", sprintf("%.4f-%.4f", min(st$dark_frac), max(st$dark_frac)), "\n")
  cat("   all ok:", all(st$ok), "\n")
  if (any(!st$ok)) print(st[!st$ok, ])
}

top <- sort(list.files(plot_dir, pattern = "^all_samples.*[.]png$", full.names = TRUE))
if (length(top)) {
  st_top <- do.call(rbind, lapply(top, one)); st_top$type <- "overview"
  cat("\n== overview\n"); print(st_top)
} else st_top <- NULL

all_st <- rbind(st_violin, st_scatter, st_top)
write.csv(all_st, file.path(rep_dir, "qc_plot_check.csv"), row.names = FALSE)
cat("\nwritten", file.path(rep_dir, "qc_plot_check.csv"), "\n")
