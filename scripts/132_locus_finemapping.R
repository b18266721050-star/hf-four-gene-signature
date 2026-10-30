## 132_locus_finemapping.R -- SuSiE fine-mapping at 9p13.3 with 1000G EUR LD
## For each gene: fine-map the cis-eQTL and both heart-failure outcomes on the same
## SNP set and LD matrix, then run coloc.susie wherever both traits yield a
## credible set. Also reports the posterior inclusion probability of every MR
## instrument, which is the quantitative version of the "not gene-specific" claim.
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(susieR); library(coloc) })

ROOT <- "E:/another data file/GEObulk analysis project"
LD   <- file.path(ROOT, "results", "Coloc", "05_ld")
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
TAB  <- file.path(ROOT, "results", "Coloc", "01_tables")
LOG  <- file.path(ROOT, "results", "Coloc", "03_logs")
sink(file.path(LOG, "132_locus_finemapping.log"), split = TRUE)

GENES <- c("AQP3", "RP11-176F3.7", "GALT", "DCTN3", "NOL6", "RP11-384P7.7", "BAG1")
OUTS  <- list(HERMES = list(n = 977323, s = 47309/977323),
              FinnGen = list(n = 453733, s = 33250/453733))
N_EXPR <- 31684
L <- 10

align <- read.csv(file.path(LD, "align.csv"), stringsAsFactors = FALSE)
flipmap <- setNames(align$flip, align$snp)

fixR <- function(R) {
  R <- as.matrix(R)
  e <- eigen(R, symmetric = TRUE)
  if (min(e$values) < 1e-4) {
    v <- e$values; v[v < 1e-4] <- 1e-4
    R <- e$vectors %*% (v * t(e$vectors))
    R <- (R + t(R)) / 2
    d <- sqrt(diag(R)); R <- R / outer(d, d); diag(R) <- 1
  }
  R
}

cs_table <- function(s, gene, outcome, trait) {
  out <- list()
  if (is.null(s$sets$cs) || length(s$sets$cs) == 0) return(NULL)
  for (k in seq_along(s$sets$cs)) {
    ii <- s$sets$cs[[k]]
    out[[k]] <- data.frame(
      gene = gene, outcome = outcome, trait = trait, cs = k, cs_size = length(ii),
      lead_snp = s$snp[s$sets$cs[[k]][which.max(s$pip[ii])]],
      max_pip = max(s$pip[ii]), mean_pip = mean(s$pip[ii]),
      min_purity = min(s$sets$purity$min.abs.corr[k]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

csrows <- list(); instrows <- list(); susierows <- list()
INSTR <- c("rs10971815", "rs116972937", "rs12551364", "rs3814511", "rs573538", "rs626069")
pipmat <- data.frame(snp = INSTR)

for (g in GENES) {
  fn <- file.path(LD, sprintf("ld_%s.csv", gsub("[.]", "_", g)))
  if (!file.exists(fn)) { cat("no LD for", g, "\n"); next }
  snps <- read.csv(file.path(LD, sprintf("ld_%s_snps.csv", gsub("[.]", "_", g))),
                   stringsAsFactors = FALSE)$snp
  R <- fixR(read.csv(fn, header = FALSE))
  cat(sprintf("\n== %s: %d SNPs, LD rank %d\n", g, length(snps),
              sum(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 1e-3)))
  for (o in names(OUTS)) {
    d <- read.csv(file.path(DAT, sprintf("coloc_%s_%s.csv", gsub("[.]", "_", g), o)),
                  stringsAsFactors = FALSE)
    m <- match(snps, d$snp)
    ok <- !is.na(m)
    if (sum(ok) < 20) { cat("  too few overlapping SNPs:", o, sum(ok), "\n"); next }
    Rsub <- R[ok, ok]
    z1 <- d$z1[m[ok]] * flipmap[snps[ok]]
    z2 <- d$z2[m[ok]] * flipmap[snps[ok]]
    sn <- snps[ok]
    s1 <- try(susie_rss(z = z1, R = Rsub, n = N_EXPR, L = L), silent = TRUE)
    s2 <- try(susie_rss(z = z2, R = Rsub, n = OUTS[[o]]$n, L = L), silent = TRUE)
    if (inherits(s1, "try-error")) { cat("  susie eQTL failed\n"); next }
    if (inherits(s2, "try-error")) { cat("  susie outcome failed\n"); next }
    s1$snp <- sn; s2$snp <- sn
    t1 <- cs_table(s1, g, o, "eQTL"); t2 <- cs_table(s2, g, o, "outcome")
    if (!is.null(t1)) csrows[[length(csrows) + 1]] <- t1
    if (!is.null(t2)) csrows[[length(csrows) + 1]] <- t2
    cat(sprintf("  %-8s eQTL CS=%d (maxPIP %.3f)  outcome CS=%d (maxPIP %.3f)\n",
                o, ifelse(is.null(t1), 0, nrow(t1)), max(s1$pip),
                ifelse(is.null(t2), 0, nrow(t2)), max(s2$pip)))
    if (!is.null(t1) && !is.null(t2)) {
      cr <- try(coloc.susie(s1, s2), silent = TRUE)
      if (!inherits(cr, "try-error") && !is.null(cr$summary)) {
        sm <- cr$summary
        sm$gene <- g; sm$outcome <- o
        susierows[[length(susierows) + 1]] <- sm
        cat("    coloc.susie best PP.H4 =", round(max(sm$PP.H4.abf), 3), "\n")
      }
    }
    # instrument PIPs in the eQTL fine-mapping (the gene-specificity question)
    pip <- setNames(s1$pip, sn)
    for (r in INSTR) {
      if (!is.na(pip[r])) {
        instrows[[length(instrows) + 1]] <- data.frame(
          gene = g, outcome = o, snp = r, pip = as.numeric(pip[r]),
          in_95cs = any(sapply(s1$sets$cs, function(ii) r %in% sn[ii])),
          stringsAsFactors = FALSE)
      }
    }
  }
}

if (length(csrows)) {
  cs <- do.call(rbind, csrows)
  write.csv(cs, file.path(TAB, "T06_susie_credible_sets.csv"), row.names = FALSE)
  cat("\n== credible sets ==\n"); print(cs, row.names = FALSE)
}
if (length(susierows)) {
  ss <- do.call(rbind, susierows)
  write.csv(ss, file.path(TAB, "T07_coloc_susie.csv"), row.names = FALSE)
} else {
  cat("\ncoloc.susie: no gene x outcome pair produced credible sets in BOTH traits\n")
}
if (length(instrows)) {
  ia <- do.call(rbind, instrows)
  w <- reshape(ia[, c("gene", "snp", "pip")], idvar = "snp", timevar = "gene", direction = "wide")
  names(w) <- sub("pip[.]", "pip_", names(w))
  write.csv(w, file.path(TAB, "T08_instrument_pip.csv"), row.names = FALSE)
  cat("\n== MR instrument PIPs per gene ==\n"); print(w, row.names = FALSE)
  csflag <- ia[ia$in_95cs, c("gene", "snp")]
  cat("\ninstruments inside a 95% credible set:\n"); print(csflag, row.names = FALSE)
}
sink()
