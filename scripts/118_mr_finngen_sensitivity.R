## 118_mr_finngen_sensitivity.R -- is the FinnGen AQP3 signal robust?
##  1. leave-one-out for the FinnGen AQP3 instrument set
##  2. multiplicity across the gene-level tests, within each outcome
##  3. side-by-side HERMES vs FinnGen table with direction concordance
##  4. diagnosis of the Egger instability (bimodal instrument strength)
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
source("E:/another data file/GEObulk analysis project/scripts/mr_core.R")
use("data.table")

MDIR <- file.path(DIR$results, "MR")
hm <- fread(file.path(MDIR, "01_tables", "T08_mr_main_results.csv"))
ff <- fread(file.path(MDIR, "01_tables", "T15_finngen_mr_results.csv"))
hd <- fread(file.path(MDIR, "01_data", "T10_harmonised_snps.csv"))
fd <- fread(file.path(MDIR, "01_tables", "T16_finngen_harmonised.csv"))

## ---------- 1. leave-one-out, FinnGen AQP3 ----------
d <- fd[exposure == "eQTLGen_blood" & gene == "AQP3"]
setorder(d, rsid)
loo <- mr_loo(d$bx, d$by, d$sx, d$sy, label = "AQP3_FinnGen")
loo[, reason := "none"]
loo[dropped > 0, reason := paste0("drop ", d$rsid[dropped])]
fwrite(loo, file.path(MDIR, "01_tables", "T18_finngen_loo_AQP3.csv"))
cat("== FinnGen AQP3 leave-one-out ==\n"); print(loo[, .(reason, beta = round(beta,4), se = round(se,4), pval = signif(pval,3))])

## ---------- 2. Egger instability: instrument-strength geometry ----------
geo <- d[, .(rsid, bx_oriented = abs(bx), by_oriented = by * sign(bx), ratio = by/bx, fstat)]
setorder(geo, bx_oriented)
cat("\n== FinnGen AQP3 instrument geometry (bx oriented positive) ==\n"); print(geo)
cat("correlation between |bx| and the Wald ratio: ", round(cor(abs(d$bx), d$by/d$bx), 3), "\n")
cat("number of SNPs with |bx| in the lower half: ", sum(abs(d$bx) <= median(abs(d$bx))), "\n")

## drop-the-outliers sensitivity (post hoc, labelled as such)
d2 <- d[sign(bx) < 0]   # the three SNPs whose expression effect is negative
r2 <- ivw_fixed(d2$bx, d2$by, d2$sx, d2$sy, label = "AQP3_FinnGen_negSNPsOnly")
cat("post-hoc IVW using the three negative-bx SNPs only: OR ",
    round(exp(r2$beta), 3), " (", round(exp(r2$beta - 1.96*r2$se), 3), "-",
    round(exp(r2$beta + 1.96*r2$se), 3), "), p = ", signif(r2$pval, 3), "\n", sep = "")

## ---------- 3. multiplicity within each outcome ----------
pick <- function(dt) dt[method %in% c("IVW", "Wald ratio"), .(exposure, gene, nsnp, beta, se, pval)]
prim <- rbind(pick(hm)[, outcome := "HERMES"], pick(ff)[, outcome := "FinnGen"])
prim[, `:=`(OR = exp(beta), lo = exp(beta - 1.96*se), hi = exp(beta + 1.96*se))]
prim[, fdr := p.adjust(pval, method = "BH"), by = outcome]
prim[, bonf := pmin(1, pval * .N), by = outcome]
setorder(prim, outcome, pval)
fwrite(prim, file.path(MDIR, "01_tables", "T20_two_outcome_multiplicity.csv"))
cat("\n== gene-level tests, multiplicity within each outcome ==\n")
print(prim[, .(outcome, exposure, gene, nsnp, OR = round(OR,3),
               lo = round(lo,3), hi = round(hi,3), pval = signif(pval,3),
               fdr = signif(fdr,3), bonf = signif(bonf,3))])

## ---------- 4. side-by-side comparison ----------
cmp <- merge(
  prim[outcome == "HERMES", .(exposure, gene, or_hm = OR, lo_hm = lo, hi_hm = hi, p_hm = pval, n_hm = nsnp)],
  prim[outcome == "FinnGen", .(exposure, gene, or_ff = OR, lo_ff = lo, hi_ff = hi, p_ff = pval, n_ff = nsnp)],
  by = c("exposure", "gene"), all = TRUE)
cmp[, `:=`(same_direction = sign(log(or_hm)) == sign(log(or_ff)),
           fold_change_between_cohorts = or_ff / or_hm)]
setorder(cmp, exposure, gene)
fwrite(cmp, file.path(MDIR, "01_tables", "T19_hermes_vs_finngen.csv"))
cat("\n== HERMES vs FinnGen ==\n")
print(cmp[, .(exposure, gene, or_hm = round(or_hm,3), or_ff = round(or_ff,3),
              same_direction, fold = round(fold_change_between_cohorts,2),
              p_hm = signif(p_hm,3), p_ff = signif(p_ff,3))])
cat("\nconcordant directions: ", sum(cmp$same_direction, na.rm=TRUE), "/", nrow(cmp), "\n")

## ---------- 5. cross-cohort meta-analysis of the two MR estimates ----------
## NOTE: the two outcome GWAS are not fully independent (shared cohorts), so this is
## reported as a descriptive synthesis, not as an independent replication test.
mm <- cmp[!is.na(or_hm) & !is.na(or_ff)]
meta <- mm[, {
  y <- c(log(or_hm), log(or_ff))
  v <- c(((log(hi_hm) - log(lo_hm)) / (2 * 1.96))^2,
         ((log(hi_ff) - log(lo_ff)) / (2 * 1.96))^2)
  w <- 1 / v
  b_fe <- sum(w * y) / sum(w); se_fe <- sqrt(1 / sum(w))
  Q <- sum(w * (y - b_fe)^2); df <- length(y) - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  w_re <- 1 / (v + tau2)
  b_re <- sum(w_re * y) / sum(w_re); se_re <- sqrt(1 / sum(w_re))
  .(k = length(y), beta_fe = b_fe, se_fe = se_fe, p_fe = 2 * pnorm(-abs(b_fe / se_fe)),
    beta_re = b_re, se_re = se_re, p_re = 2 * pnorm(-abs(b_re / se_re)),
    Q = Q, Q_pval = pchisq(Q, df, lower.tail = FALSE),
    I2 = max(0, (Q - df) / Q) * 100)
}, by = .(exposure, gene)]
meta[, `:=`(OR_fe = exp(beta_fe), OR_fe_lo = exp(beta_fe - 1.96*se_fe), OR_fe_hi = exp(beta_fe + 1.96*se_fe),
            OR_re = exp(beta_re), OR_re_lo = exp(beta_re - 1.96*se_re), OR_re_hi = exp(beta_re + 1.96*se_re))]
fwrite(meta, file.path(MDIR, "01_tables", "T21_cross_cohort_meta.csv"))
cat("\n== cross-cohort synthesis (descriptive; cohorts overlap) ==\n")
print(meta[, .(exposure, gene, OR_fe = round(OR_fe,3), p_fe = signif(p_fe,3),
               OR_re = round(OR_re,3), p_re = signif(p_re,3),
               Q_pval = signif(Q_pval,3), I2 = round(I2,1))])


write_session(file.path(MDIR, "03_logs"), "118_finngen_sensitivity")
