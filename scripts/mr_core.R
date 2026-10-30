## mr_core.R -- self-contained Mendelian randomisation core (pure R, no external MR package)
## Implements: allele harmonisation, Wald ratio, IVW (fixed + multiplicative random effects),
## MR-Egger (with residual-SE correction and intercept test), weighted median (bootstrap SE),
## Cochran Q, leave-one-out, F statistic, Steiger directionality.

comp_strand <- function(x) chartr("ACGT", "TGCA", toupper(x))

harmonise_mr <- function(exp, out, label = "") {
  ## exp: rsid, effect_allele, other_allele, beta, se, pval, gene, fstat
  ## out: rsid, effect_allele, other_allele, beta, se, pval
  e <- as.data.table(exp); o <- as.data.table(out)
  setnames(e, c("effect_allele","other_allele"), c("eA","oA"), skip_absent = TRUE)
  setnames(o, c("effect_allele","other_allele"), c("oeA","ooA"), skip_absent = TRUE)
  m <- merge(e, o, by = "rsid", suffixes = c("_exp","_out"))
  if (!nrow(m)) return(list(data = m, summary = data.table(label = label, n_in = nrow(e), n_matched = 0)))
  m[, `:=`(eA = toupper(eA), oA = toupper(oA), oeA = toupper(oeA), ooA = toupper(ooA))]
  m[, palin := paste0(eA, oA) %in% c("AT","TA","CG","GC")]
  m[, status := NA_character_]
  m[eA == oeA & oA == ooA, `:=`(status = "match", flip = 1L)]
  m[eA == ooA & oA == oeA, `:=`(status = "swap", flip = -1L)]
  m[is.na(status) & comp_strand(eA) == oeA & comp_strand(oA) == ooA, `:=`(status = "strand_match", flip = 1L)]
  m[is.na(status) & comp_strand(eA) == ooA & comp_strand(oA) == oeA, `:=`(status = "strand_swap", flip = -1L)]
  m[is.na(status), `:=`(status = "allele_mismatch", flip = NA_integer_)]
  ## palindromic SNPs: outcome allele frequency unavailable in HERMES -> drop conservatively
  m[status %in% c("match","swap","strand_match","strand_swap") & palin == TRUE,
    `:=`(status = "palindromic_dropped", flip = NA_integer_)]
  ok <- m[status %in% c("match","swap","strand_match","strand_swap")]
  ok[, `:=`(bx = beta_exp, by = beta_out * flip,
            sx = se_exp,   sy = se_out,
            effect_allele = eA, other_allele = oA)]
  summ <- data.table(label = label, n_in = nrow(e), n_matched = nrow(m), n_used = nrow(ok),
                     n_swap = sum(m$status == "swap"), n_strand = sum(grepl("strand", m$status)),
                     n_palin_drop = sum(m$status == "palindromic_dropped"),
                     n_mismatch = sum(m$status == "allele_mismatch"))
  list(data = ok, summary = summ, all = m)
}

wald_ratio <- function(bx, by, sx, sy) {
  b <- by / bx
  se <- sqrt(sy^2 / bx^2 + (by^2 * sx^2) / bx^4)   # second-order delta method
  c(beta = b, se = se)
}

ivw_fixed <- function(bx, by, sx, sy, label = "") {
  w <- 1 / sy^2
  b <- sum(w * bx * by) / sum(w * bx^2)
  se <- sqrt(1 / sum(w * bx^2))
  n <- length(bx)
  Q <- sum(w * (by - b * bx)^2); df <- n - 1
  Qp <- if (df > 0) pchisq(Q, df, lower.tail = FALSE) else NA_real_
  phi <- if (df > 0) max(1, Q / df) else 1
  data.table(method = "IVW", label = label, nsnp = n, beta = b, se = se, se_re = se * sqrt(phi),
             pval = 2 * pnorm(-abs(b / se)),
             pval_re = 2 * pnorm(-abs(b / (se * sqrt(phi)))),
             Q = Q, Q_df = df, Q_pval = Qp, phi = phi)
}

mr_egger <- function(bx, by, sx, sy, label = "") {
  n <- length(bx)
  if (n < 3) return(data.table(method = "MR-Egger", label = label, nsnp = n))
  s <- sign(bx); s[s == 0] <- 1
  bxs <- bx * s; bys <- by * s; sys <- sy
  w <- 1 / sys^2
  fit <- lm(bys ~ bxs, weights = w)
  cf <- summary(fit)$coefficients
  b1 <- cf[2, 1]; se1 <- cf[2, 2]; b0 <- cf[1, 1]; se0 <- cf[1, 2]
  sig <- summary(fit)$sigma
  if (is.finite(sig) && sig > 1) { b1 <- b1 / sig; se1 <- se1 / sig; b0 <- b0 / sig; se0 <- se0 / sig }
  ## MR-Egger is a regression, so inference must use t with n-2 df (TwoSampleMR convention).
  ## Using the normal approximation badly overstates significance when n is small.
  df <- n - 2
  data.table(method = "MR-Egger", label = label, nsnp = n, beta = b1, se = se1,
             pval = 2 * pt(-abs(b1 / se1), df),
             intercept = b0, intercept_se = se0,
             intercept_pval = 2 * pt(-abs(b0 / se0), df))
}

wm_median_w <- function(ratio, w) {
  ord <- order(ratio); ratio <- ratio[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  if (any(cw == 0.5)) return(ratio[which(cw == 0.5)[1]])
  if (!length(which(cw < 0.5))) return(ratio[1])
  if (!length(which(cw > 0.5))) return(ratio[length(ratio)])
  lo <- max(which(cw < 0.5)); hi <- min(which(cw > 0.5))
  ratio[lo] + (ratio[hi] - ratio[lo]) * (0.5 - cw[lo]) / (cw[hi] - cw[lo])
}

mr_wmedian <- function(bx, by, sx, sy, label = "", B = 2000, seed = 20251029) {
  n <- length(bx)
  if (n < 3) return(data.table(method = "Weighted median", label = label, nsnp = n))
  ratio <- by / bx
  w <- bx^2 / sy^2
  b <- wm_median_w(ratio, w)
  set.seed(seed)
  bb <- numeric(B)
  for (i in seq_len(B)) {
    bxs <- rnorm(n, bx, sx); bys <- rnorm(n, by, sy)
    rr <- bys / bxs
    bb[i] <- wm_median_w(rr, bxs^2 / sy^2)
  }
  se <- sd(bb)
  data.table(method = "Weighted median", label = label, nsnp = n, beta = b, se = se,
             pval = 2 * pnorm(-abs(b / se)))
}

mr_loo <- function(bx, by, sx, sy, label = "") {
  n <- length(bx)
  if (n < 3) return(NULL)
  res <- data.table()
  for (i in seq_len(n)) {
    r <- ivw_fixed(bx[-i], by[-i], sx[-i], sy[-i], label = label)
    res <- rbind(res, data.table(label = label, dropped = i, beta = r$beta, se = r$se,
                                 pval = r$pval, nsnp = n - 1))
  }
  full <- ivw_fixed(bx, by, sx, sy, label = label)
  rbind(res, data.table(label = label, dropped = 0, beta = full$beta, se = full$se,
                        pval = full$pval, nsnp = n))
}

f_stat <- function(z) mean(z^2)

steiger <- function(p_exp, n_exp, lor, f_alt, ncase, ncontrol, prevalence = 0.10) {
  ## r for continuous exposure
  tval <- qt(p_exp / 2, df = n_exp - 2, lower.tail = FALSE)
  r_exp <- tval / sqrt(tval^2 + n_exp - 2)
  ## r on the liability scale for the binary outcome
  K <- prevalence
  thr <- qnorm(1 - K)
  b_liab <- lor * K * (1 - K) / dnorm(thr)
  r_out <- sqrt(pmin(1, 2 * f_alt * (1 - f_alt) * b_liab^2))
  z1 <- atanh(r_exp); z2 <- atanh(r_out)
  se_diff <- sqrt(1 / (n_exp - 3) + 1 / (ncase + ncontrol - 3))
  zstat <- (z1 - z2) / se_diff
  c(r_exp = r_exp, r_out = r_out, z = zstat, pval = 2 * pnorm(-abs(zstat)))
}
