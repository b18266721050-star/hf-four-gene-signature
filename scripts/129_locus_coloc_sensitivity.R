## 129_locus_coloc_sensitivity.R
## (1) leave-top-k-out sensitivity: is any PP.H4 carried by a single SNP?
## (2) eQTL x eQTL colocalisation inside the block: can these variants separate
##     AQP3 from its neighbours at all?
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(coloc))

ROOT <- "E:/another data file/GEObulk analysis project"
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
TAB  <- file.path(ROOT, "results", "Coloc", "01_tables")
LOG  <- file.path(ROOT, "results", "Coloc", "03_logs")
sink(file.path(LOG, "129_locus_coloc_sensitivity.log"), split = TRUE)

GENES <- c("AQP3", "RP11-176F3.7", "DCTN3", "GALT", "NUDT2", "RP11-384P7.7", "BAG1", "NOL6")
OUTS  <- list(HERMES  = list(n = 977323, s = 47309/977323),
              FinnGen = list(n = 453733, s = 33250/453733))
N1 <- 31684; SDY <- 0.15; P12 <- 1e-5

rd <- function(g, o) {
  f <- file.path(DAT, sprintf("coloc_%s_%s.csv", gsub("\\.", "_", g), o))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}

run <- function(d1, d2) {
  r <- try(coloc.abf(d1, d2, p12 = P12), silent = TRUE)
  if (inherits(r, "try-error")) return(c(NA, NA, NA, NA, NA))
  s <- r$summary
  c(s[["PP.H1.abf"]], s[["PP.H2.abf"]], s[["PP.H3.abf"]], s[["PP.H4.abf"]],
    s[["PP.H4.abf"]] / max(1e-12, s[["PP.H3.abf"]] + s[["PP.H4.abf"]]))
}

## ---------- (1) leave-top-k-out ----------
rows <- list(); k <- 0
for (g in GENES) for (o in names(OUTS)) {
  d <- rd(g, o); if (is.null(d)) next
  ord <- order(-abs(d$z1))
  for (drop in c(0, 1, 3, 5, 10)) {
    dd <- if (drop == 0) d else d[-ord[seq_len(drop)], ]
    if (nrow(dd) < 10) next
    v <- run(list(snp = dd$snp, beta = dd$beta1, varbeta = dd$varbeta1,
                  type = "quant", N = N1, sdY = SDY),
             list(snp = dd$snp, beta = dd$beta2, varbeta = dd$varbeta2,
                  type = "cc", N = OUTS[[o]]$n, s = OUTS[[o]]$s))
    k <- k + 1
    rows[[k]] <- data.frame(gene = g, outcome = o, dropped_top = drop, nsnp = nrow(dd),
                            PP.H1 = v[1], PP.H2 = v[2], PP.H3 = v[3], PP.H4 = v[4],
                            H4_cond = v[5],
                            min_p_outcome = 2 * pnorm(-max(abs(dd$z2))),
                            stringsAsFactors = FALSE)
  }
}
t2 <- do.call(rbind, rows)
write.csv(t2, file.path(TAB, "T02_coloc_leave_top_k_out.csv"), row.names = FALSE)
cat("\n== PP.H4 by number of top eQTL SNPs removed ==\n")
print(reshape(t2[, c("gene", "outcome", "dropped_top", "PP.H4")],
              idvar = c("gene", "outcome"), timevar = "dropped_top", direction = "wide"), row.names = FALSE)

## ---------- (2) eQTL x eQTL ----------
pairs <- list(); m <- 0
for (i in 1:(length(GENES) - 1)) for (j in (i + 1):length(GENES)) {
  a <- rd(GENES[i], "FinnGen"); b <- rd(GENES[j], "FinnGen")
  if (is.null(a) || is.null(b)) next
  dd <- merge(a[, c("snp", "pos", "beta1", "varbeta1", "z1")],
              b[, c("snp", "beta1", "varbeta1", "z1")], by = "snp", suffixes = c(".a", ".b"))
  if (nrow(dd) < 20) next
  v <- run(list(snp = dd$snp, beta = dd$beta1.a, varbeta = dd$varbeta1.a,
                type = "quant", N = N1, sdY = SDY),
           list(snp = dd$snp, beta = dd$beta1.b, varbeta = dd$varbeta1.b,
                type = "quant", N = N1, sdY = SDY))
  m <- m + 1
  pairs[[m]] <- data.frame(gene_a = GENES[i], gene_b = GENES[j], nsnp = nrow(dd),
                           PP.H3 = v[3], PP.H4 = v[4], H4_cond = v[5],
                           stringsAsFactors = FALSE)
  cat(sprintf("%-14s x %-14s n=%4d  H4=%.3f\n", GENES[i], GENES[j], nrow(dd), v[4]))
}
t3 <- do.call(rbind, pairs)
t3 <- t3[order(-t3$PP.H4), ]
write.csv(t3, file.path(TAB, "T03_coloc_eqtl_vs_eqtl.csv"), row.names = FALSE)
cat("\n== shared cis-eQTL signal (AQP3 vs neighbours) ==\n")
print(t3[t3$gene_a == "AQP3" | t3$gene_b == "AQP3", ], row.names = FALSE)
sink()
