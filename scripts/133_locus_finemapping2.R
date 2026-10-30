## 133_locus_finemapping2.R -- SuSiE fine-mapping, second pass.
## Fixes relative to 132: (i) fine-map a +/-250 kb window around the lead SNP
## instead of the whole cis window, so that the 503-sample LD matrix is close to
## full rank; (ii) LD-prune at r2 > 0.95 to remove near-duplicates; (iii) L = 5,
## max_iter = 2000, and convergence is checked before any credible set is reported;
## (iv) traits without a genome-wide significant signal in the window are NOT
## fine-mapped (a credible set cannot be defined without an association);
## (v) instead of PIPs, the instruments are related to each gene's lead SNP by r2,
## which is robust and directly answers the gene-specificity question.
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(susieR) })

ROOT <- "E:/another data file/GEObulk analysis project"
LD   <- file.path(ROOT, "results", "Coloc", "05_ld")
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
TAB  <- file.path(ROOT, "results", "Coloc", "01_tables")
LOG  <- file.path(ROOT, "results", "Coloc", "03_logs")
sink(file.path(LOG, "133_locus_finemapping2.log"), split = TRUE)

GENES <- c("AQP3", "RP11-176F3.7", "GALT", "DCTN3", "NOL6", "RP11-384P7.7", "BAG1")
OUTS  <- list(HERMES = list(n = 977323, s = 47309/977323),
              FinnGen = list(n = 453733, s = 33250/453733))
N_EXPR <- 31684
HALF <- 250000
R2_PRUNE <- 0.95
GW <- 5.45e-8
INSTR <- c("rs10971815", "rs116972937", "rs12551364", "rs3814511", "rs573538", "rs626069")

align <- read.csv(file.path(LD, "align.csv"), stringsAsFactors = FALSE)
flipmap <- setNames(align$flip, align$snp)

fixR <- function(R) {
  e <- eigen(R, symmetric = TRUE)
  if (min(e$values) < 1e-4) {
    v <- e$values; v[v < 1e-4] <- 1e-4
    R <- e$vectors %*% (v * t(e$vectors)); R <- (R + t(R)) / 2
    d <- sqrt(diag(R)); R <- R / outer(d, d); diag(R) <- 1
  }
  R
}

prune <- function(R, z, thr = R2_PRUNE) {
  ord <- order(-abs(z))
  keep <- integer(0)
  for (i in ord) {
    if (length(keep) && max(R[i, keep]^2) > thr) next
    keep <- c(keep, i)
  }
  sort(keep)
}

res_cs <- list(); res_inst <- list(); res_w <- list()

for (g in GENES) {
  fn <- file.path(LD, sprintf("ld_%s.csv", gsub("[.]", "_", g)))
  if (!file.exists(fn)) { cat("no LD for", g, "\n"); next }
  snps <- read.csv(file.path(LD, sprintf("ld_%s_snps.csv", gsub("[.]", "_", g))),
                   stringsAsFactors = FALSE)$snp
  Rfull <- fixR(as.matrix(read.csv(fn, header = FALSE)))
  if (nrow(Rfull) != length(snps)) { cat("LD/SNP mismatch for", g, "\n"); next }

  eqt <- read.csv(file.path(DAT, "eqtl_region_long.csv"), stringsAsFactors = FALSE)
  eqt <- eqt[!duplicated(eqt$snp), ]
  pos <- eqt$pos[match(snps, eqt$snp)]
  cat(sprintf("\n===== %s: %d SNPs, cis window %.2f-%.2f Mb =====\n",
              g, length(snps), min(pos, na.rm = TRUE)/1e6, max(pos, na.rm = TRUE)/1e6))

  for (o in names(OUTS)) {
    d <- read.csv(file.path(DAT, sprintf("coloc_%s_%s.csv", gsub("[.]", "_", g), o)),
                  stringsAsFactors = FALSE)
    m <- match(snps, d$snp)
    ok <- !is.na(m)
    z1 <- rep(NA_real_, length(snps)); z2 <- z1
    z1[ok] <- d$z1[m[ok]] * flipmap[snps[ok]]
    z2[ok] <- d$z2[m[ok]] * flipmap[snps[ok]]
    lead1 <- which.max(abs(z1))
    win <- which(!is.na(z1) & abs(pos - pos[lead1]) <= HALF)
    R <- Rfull[win, win]; s <- snps[win]
    z1w <- z1[win]; z2w <- z2[win]
    kp <- prune(R, z1w)
    R <- R[kp, kp]; s <- s[kp]; z1w <- z1w[kp]; z2w <- z2w[kp]
    dimnames(R) <- list(s, s)
    nrank <- sum(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 1e-3)
    p2 <- 2 * pnorm(-max(abs(z2w)))
    cat(sprintf("  %-8s window=%d  pruned=%d  LD rank=%d  max|z| eQTL=%.1f  HF=%.2f (min p=%.2g)\n",
                o, length(win), length(kp), nrank, max(abs(z1w)), max(abs(z2w)), p2))

    # ---- eQTL fine-mapping
    f1 <- try(susie_rss(z = z1w, R = R, n = N_EXPR, L = 5, max_iter = 2000), silent = TRUE)
    conv1 <- !inherits(f1, "try-error") && isTRUE(f1$converged)
    if (!inherits(f1, "try-error")) {
      cat(sprintf("    eQTL: converged=%s  CS=%d  maxPIP=%.3f\n", conv1,
                  length(f1$sets$cs), max(f1$pip)))
      if (conv1 && length(f1$sets$cs)) {
        for (k in seq_along(f1$sets$cs)) {
          ii <- f1$sets$cs[[k]]
          res_cs[[length(res_cs) + 1]] <- data.frame(
            gene = g, trait = "eQTL", outcome = o, cs = k, size = length(ii),
            lead = s[ii][which.max(f1$pip[ii])], max_pip = max(f1$pip[ii]),
            min_purity = min(f1$sets$purity$min.abs.corr[k]),
            stringsAsFactors = FALSE)
        }
      }
      # instruments vs lead SNP of this gene (robust alternative to PIP)
      lp <- s[which.max(abs(z1w))]
      for (r in INSTR) {
        if (!(r %in% s)) next
        res_inst[[length(res_inst) + 1]] <- data.frame(
          gene = g, lead_snp = lp, instrument = r,
          r2 = R[r, lp]^2, dist_kb = (pos[match(r, snps)] - pos[match(lp, snps)]) / 1000,
          stringsAsFactors = FALSE)
      }
    } else cat("    eQTL: susie failed\n")

    # ---- outcome fine-mapping (only if the window actually contains a signal)
    if (p2 < GW) {
      f2 <- try(susie_rss(z = z2w, R = R, n = OUTS[[o]]$n, L = 5, max_iter = 2000), silent = TRUE)
      conv2 <- !inherits(f2, "try-error") && isTRUE(f2$converged)
      cat(sprintf("    HF  : converged=%s  CS=%d\n", conv2, if (inherits(f2, "try-error")) 0 else length(f2$sets$cs)))
      if (conv2 && length(f2$sets$cs)) {
        for (k in seq_along(f2$sets$cs)) {
          ii <- f2$sets$cs[[k]]
          res_cs[[length(res_cs) + 1]] <- data.frame(
            gene = g, trait = "outcome", outcome = o, cs = k, size = length(ii),
            lead = s[ii][which.max(f2$pip[ii])], max_pip = max(f2$pip[ii]),
            min_purity = min(f2$sets$purity$min.abs.corr[k]),
            stringsAsFactors = FALSE)
        }
      }
    } else {
      cat(sprintf("    HF  : no genome-wide significant signal in window (min p=%.2g) -> not fine-mapped\n", p2))
      res_w[[length(res_w) + 1]] <- data.frame(gene = g, outcome = o, min_p_hf = p2,
                                               max_z = max(abs(z2w)), n_snp = length(kp),
                                               stringsAsFactors = FALSE)
    }
  }
}

if (length(res_cs)) {
  cs <- do.call(rbind, res_cs)
  write.csv(cs, file.path(TAB, "T06_susie_credible_sets.csv"), row.names = FALSE)
  cat("\n== credible sets (converged runs only) ==\n"); print(cs, row.names = FALSE)
} else cat("\nno credible sets reported\n")

if (length(res_inst)) {
  ia <- do.call(rbind, res_inst)
  write.csv(ia, file.path(TAB, "T08_instrument_vs_lead.csv"), row.names = FALSE)
  cat("\n== instruments vs each gene's lead eQTL SNP ==\n"); print(ia, row.names = FALSE, digits = 3)
}
if (length(res_w)) {
  w <- do.call(rbind, res_w)
  write.csv(w, file.path(TAB, "T09_outcome_window_power.csv"), row.names = FALSE)
}
sink()
