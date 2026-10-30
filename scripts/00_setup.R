## 00_setup.R -- single entry point for paths, packages and plotting defaults
## ---------------------------------------------------------------------------
## Usage (top of every analysis script):
##   source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
##
## This project was developed on a Windows machine whose R library sits in
## D:/R/library/4.6 (517 packages) instead of the default R library. Override
## with the environment variables GEOBULK_LIB / GEOBULK_ROOT when cloning
## elsewhere, e.g.  Sys.setenv(GEOBULK_ROOT = "/data/hf_project").
## ---------------------------------------------------------------------------

R_LIB <- Sys.getenv("GEOBULK_LIB", "D:/R/library/4.6")
if (dir.exists(R_LIB)) .libPaths(c(R_LIB, .libPaths()))

PROJ <- Sys.getenv("GEOBULK_ROOT", "E:/another data file/GEObulk analysis project")

DIR <- list(
  root    = PROJ,
  raw     = file.path(PROJ, "raw_data"),
  clean   = file.path(PROJ, "clean_data"),
  scripts = file.path(PROJ, "scripts"),
  results = file.path(PROJ, "results"),
  report  = file.path(PROJ, "report"),
  tmp     = file.path(PROJ, "_tmp")
)

## ---- packages -------------------------------------------------------------
need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("missing packages: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  invisible(TRUE)
}

use <- function(...) {
  pkgs <- unlist(list(...))
  need(pkgs)
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(pkgs)
}

## ---- directories ----------------------------------------------------------
new_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

## ---- plotting defaults (journal style: 8 pt, black-white, no minor grid) --
theme_pub <- function(base_size = 8, legend = "bottom") {
  need("ggplot2")
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", size = base_size, hjust = 0),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "grey20"),
      legend.key.size = grid::unit(8, "pt"),
      legend.position = legend
    )
}

## always write both a vector PDF and a 600 dpi PNG
save_pub <- function(plot, name, dir = DIR$results, width = 7.1, height = 5,
                     dpi = 600) {
  need(c("ggplot2"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot,
                  width = width, height = height, device = "pdf")
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = dpi)
  invisible(file.path(dir, paste0(name, c(".pdf", ".png"))))
}

## ---- provenance -----------------------------------------------------------
write_session <- function(dir, tag, note = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0("sessionInfo_", tag, ".txt"))
  con <- file(f, open = "wt")
  writeLines(c(paste0("## ", tag, "  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               if (!is.null(note)) paste0("## ", note) else NULL,
               "## .libPaths():", paste0("##   ", .libPaths()),
               "## ----", capture.output(sessionInfo())), con)
  close(con)
  invisible(f)
}

## ---- small helpers used across modules ------------------------------------
`%||%` <- function(a, b) if (is.null(a) || !length(a) || is.na(a[1])) b else a

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = digits)))
}

message("00_setup loaded | root = ", PROJ, " | library = ", R_LIB)
