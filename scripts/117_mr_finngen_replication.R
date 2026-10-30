## 117_mr_finngen_replication.R -- replicate the MR analysis with FinnGen R11 heart failure as the outcome
##
## STATUS: parsing and statistics are smoke-tested on a synthetic file, but the script has NOT
## been run on the real 780 MB FinnGen release (the download could not be completed on the
## analysis machine). Run it once the file exists and check the log for the SNP match rate.
##
## usage:
##   Rscript 117_mr_finngen_replication.R
##   Rscript 117_mr_finngen_replication.R <finngen.gz> <out_dir>
##
## FinnGen R11 I9_HEARTFAIL: 33,250 cases / 420,483 controls
## columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
source("E:/another data file/GEObulk analysis project/scripts/mr_core.R")
use("data.table")

args <- commandArgs(trailingOnly = TRUE)
FF   <- if (length(args) >= 1) args[1] else file.path(DIR$raw, "mr", "finngen_R11_I9_HEARTFAIL.gz")
ODIR <- if (length(args) >= 2) args[2] else new_dir(DIR$results, "MR", "01_tables")
stopifnot(file.exists(FF))
NCASE <- 33250; NCONTROL <- 420483

say <- function(...) cat(..., "\n", sep = "")

## ---- instruments (same tables as the HERMES analysis) ----------------------
## eQTLGen only reports Z; per-allele beta/SE are reconstructed exactly as in 113_mr_run.R
eq <- fread(file.path(DIR$results, "MR", "01_data", "T02_eqtlgen_instruments_pruned.csv"))
gt <- fread(file.path(DIR$results, "MR", "01_data", "T05_gtex_instruments_pruned.csv"))
mf <- fread(file.path(DIR$results, "MR", "01_data", "T07_snp_maf.csv"))
mf[, `:=`(maf = as.numeric(maf), n_alleles = as.integer(n_alleles))]
eq <- merge(eq, mf[, .(rsid, maf, n_alleles)], by = "rsid", all.x = TRUE)
## same exclusion as 113_mr_run.R: no MAF, or a multi-allelic record
eq <- eq[!is.na(maf) & (is.na(n_alleles) | n_alleles <= 2)]
eq[, denom := sqrt(2 * maf * (1 - maf) * (n + z^2))]
eq[, `:=`(beta = z / denom, se = 1 / denom)]

inst <- rbind(
  eq[, .(exposure = "eQTLGen_blood", gene, rsid, effect_allele, other_allele, beta, se, z, pval, n, fstat)],
  gt[, .(exposure = "GTEx_heartLV", gene, rsid, effect_allele, other_allele, beta, se,
         z = NA_real_, pval, n = 386, fstat)]
)
inst[, palin := paste0(toupper(effect_allele), toupper(other_allele)) %in% c("AT","TA","CG","GC")]
inst <- inst[palin == FALSE]
say("instruments carried over: ", nrow(inst), " SNPs for ", uniqueN(inst$gene), " genes")

## ---- stream the FinnGen file, keep only rows whose rsids field mentions our SNPs
pat <- paste0("(^|,)(", paste(inst$rsid, collapse = "|"), ")($|,)")
say("scanning ", basename(FF), " (", round(file.size(FF)/1e6, 1), " MB)")
ff <- fread(FF, sep = "\t", header = TRUE, showProgress = FALSE,
            select = c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval","af_alt"))
setnames(ff, c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval","af_alt"),
             c("chr","pos38","ref","alt","rsids","beta","se","pval","eaf"))
say("FinnGen rows: ", nrow(ff))
hit <- ff[grepl(pat, rsids)]
say("rows whose rsids field contains an instrument: ", nrow(hit))

## one instrument per rsid: explode comma separated ids and keep the first match
long <- hit[, .(rsid = unlist(strsplit(rsids, ","))), by = .(chr, pos38, ref, alt, beta, se, pval, eaf)]
long <- long[rsid %in% inst$rsid]
long <- unique(long, by = c("rsid"))
setnames(long, c("ref","alt"), c("other_allele","effect_allele"))
long <- long[, .(rsid, chr, pos38, effect_allele, other_allele, beta, se, pval, eaf)]
say("unique instrument SNPs found in FinnGen: ", nrow(long), " / ", uniqueN(inst$rsid))

## ---- harmonise and run -----------------------------------------------------
res <- list(); snapl <- list(); snp <- list()
for (ex in unique(inst$exposure)) {
  for (g in unique(inst$gene)) {
    e <- inst[exposure == ex & gene == g]
    if (!nrow(e)) next
    o <- long[rsid %in% e$rsid]
    if (!nrow(o)) { say(ex, "/", g, ": no outcome data"); next }
    hm <- harmonise_mr(e, o, label = paste(ex, g, sep = "|"))
    d <- hm$data
    say(ex, "/", g, " -> matched ", hm$summary$n_matched, ", used ", hm$summary$n_used)
    if (!nrow(d)) next
    snapl[[paste(ex,g)]] <- d[, .(exposure=ex, gene=g, rsid, effect_allele, other_allele, bx, sx, by, sy, fstat)]
    w <- d[, {r <- wald_ratio(bx, by, sx, sy); .(beta = r["beta"], se = r["se"])}, by = rsid]
    w[, `:=`(exposure = ex, gene = g, pval = 2 * pnorm(-abs(beta/se)))]
    snp[[paste(ex,g)]] <- w
    if (nrow(d) == 1) {
      r <- wald_ratio(d$bx, d$by, d$sx, d$sy)
      res[[paste(ex,g)]] <- data.table(method="Wald ratio", exposure=ex, gene=g, nsnp=1,
                                       beta=r["beta"], se=r["se"], pval=2*pnorm(-abs(r["beta"]/r["se"])))
    } else {
      iv <- ivw_fixed(d$bx, d$by, d$sx, d$sy, label=g); setnames(iv, "label", "gene"); iv[, exposure := ex]
      out <- list(iv)
      if (nrow(d) >= 3) {
        eg <- mr_egger(d$bx, d$by, d$sx, d$sy, label=g); setnames(eg, "label", "gene"); eg[, exposure := ex]
        wm <- mr_wmedian(d$bx, d$by, d$sx, d$sy, label=g, seed=20251029); setnames(wm, "label", "gene"); wm[, exposure := ex]
        out <- c(out, list(eg, wm))
      }
      res[[paste(ex,g)]] <- rbindlist(out, fill = TRUE)
    }
  }
}
out <- rbindlist(res, fill = TRUE)
out[, `:=`(OR = exp(beta), OR_lower = exp(beta - 1.96*se), OR_upper = exp(beta + 1.96*se))]
fwrite(out, file.path(ODIR, "T15_finngen_mr_results.csv"))
if (length(snapl)) fwrite(rbindlist(snapl), file.path(ODIR, "T16_finngen_harmonised.csv"))
if (length(snp)) fwrite(rbindlist(snp), file.path(ODIR, "T17_finngen_per_snp.csv"))
say("---- FinnGen replication ----")
print(out[method %in% c("IVW","Wald ratio"), .(exposure, gene, nsnp, OR = round(OR,3),
      lower = round(OR_lower,3), upper = round(OR_upper,3), pval = signif(pval,3))])
say("done; compare with T08 (HERMES) for agreement in direction and magnitude")
