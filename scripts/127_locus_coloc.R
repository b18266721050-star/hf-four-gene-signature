## 127_locus_coloc.R -- single-causal-variant colocalisation at the 9p13.3 block
## Runs coloc.abf (Giambartolomei 2014) for every gene x outcome pair under a grid
## of priors, so the conclusion does not hinge on one prior choice.
## Requires: coloc (installed in the project-local library below).
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(coloc))

ROOT <- "E:/another data file/GEObulk analysis project"
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
OUT  <- file.path(ROOT, "results", "Coloc")
TAB  <- file.path(OUT, "01_tables"); dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
LOG  <- file.path(OUT, "03_logs");   dir.create(LOG,  recursive = TRUE, showWarnings = FALSE)
RDS  <- file.path(OUT, "04_objects");dir.create(RDS,  recursive = TRUE, showWarnings = FALSE)

GENES <- c("AQP3", "RP11-176F3.7", "DCTN3", "GALT", "NUDT2", "RP11-384P7.7", "BAG1", "NOL6")
OUTS  <- list(HERMES  = list(n = 977323, s = 47309/977323),
              FinnGen = list(n = 453733, s = 33250/453733))
SDY   <- c(0.15, 1.0)
P12   <- c(1e-4, 1e-5, 1e-6)

sink(file.path(LOG, "127_locus_coloc.log"), split = TRUE)

res <- list(); k <- 0
for (g in GENES) {
  for (o in names(OUTS)) {
    fn <- file.path(DAT, sprintf("coloc_%s_%s.csv", gsub("\\.", "_", g), o))
    if (!file.exists(fn)) { cat("MISSING", fn, "\n"); next }
    d <- read.csv(fn, stringsAsFactors = FALSE)
    d <- d[!duplicated(d$snp) & is.finite(d$beta1) & is.finite(d$beta2) &
             d$varbeta1 > 0 & d$varbeta2 > 0, ]
    if (nrow(d) < 10) { cat("TOO FEW SNPs:", g, o, nrow(d), "\n"); next }
    N1 <- 31684
    for (sdY in SDY) {
      d1 <- list(snp = d$snp, beta = d$beta1, varbeta = d$varbeta1,
                 type = "quant", N = N1, sdY = sdY)
      d2 <- list(snp = d$snp, beta = d$beta2, varbeta = d$varbeta2,
                 type = "cc", N = OUTS[[o]]$n, s = OUTS[[o]]$s)
      for (p12 in P12) {
        r <- try(coloc.abf(d1, d2, p12 = p12), silent = TRUE)
        if (inherits(r, "try-error")) { cat("FAIL", g, o, sdY, p12, "\n"); next }
        s <- r$summary
        k <- k + 1
        res[[k]] <- data.frame(
          gene = g, outcome = o, nsnp = as.integer(s[["nsnps"]]),
          sdY = sdY, p12 = p12,
          PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]], PP.H2 = s[["PP.H2.abf"]],
          PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
          lead_eqtl = d$snp[which.max(abs(d$z1))],
          lead_eqtl_pos = d$pos[which.max(abs(d$z1))],
          lead_out = d$snp[which.max(abs(d$z2))],
          lead_out_pos = d$pos[which.max(abs(d$z2))],
          stringsAsFactors = FALSE)
        saveRDS(r, file.path(RDS, sprintf("coloc_%s_%s_sdY%s_p12%s.rds",
                                          gsub("\\.", "_", g), o, sdY, format(p12, scientific = TRUE))))
        cat(sprintf("%-14s %-8s sdY=%-4s p12=%-8s n=%5d  H3=%.3f H4=%.3f\n",
                    g, o, sdY, format(p12, scientific = TRUE), nrow(d),
                    s[["PP.H3.abf"]], s[["PP.H4.abf"]]))
      }
    }
  }
}
out <- do.call(rbind, res)
out <- out[order(out$gene, out$outcome, out$sdY, -out$p12), ]
write.csv(out, file.path(TAB, "T01_coloc_results.csv"), row.names = FALSE)

## region-wide outcome lead (independent of any gene)
for (o in names(OUTS)) {
  fn <- if (o == "HERMES") "hermes_region.csv" else "finngen_region.csv"
  f <- read.csv(file.path(DAT, fn), stringsAsFactors = FALSE)
  f$p <- suppressWarnings(as.numeric(f$p))
  f <- f[is.finite(f$p), ]
  b <- f[which.min(f$p), ]
  cat(sprintf("REGION LEAD %-8s %s pos=%s p=%.3g beta=%.3f\n", o, b$snp, b$pos, b$p, as.numeric(b$beta)))
}

cat("\nDone:", nrow(out), "coloc runs\n")
sink()
