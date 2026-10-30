## 113_mr_run.R -- two-sample MR: genetically predicted expression of AQP3/GGT5/ISLR/LAD1 -> heart failure
## exposures: eQTLGen whole blood (primary) and GTEx v8 heart left ventricle (tissue-matched sensitivity)
## outcome  : HERMES heart failure GWAS (47,309 cases / 930,014 controls)
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({library(data.table)})
source("E:/another data file/GEObulk analysis project/scripts/mr_core.R")

ROOT <- "E:/another data file/GEObulk analysis project"
OUT  <- file.path(ROOT, "results/MR")
LOG <- file(file.path(OUT,"03_logs","113_mr_run.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || is.na(a)) b else a
GENES <- c("AQP3","GGT5","ISLR","LAD1")
NCASE <- 47309; NCONTROL <- 930014; PREV <- 0.02
say("== 113 MR run ==")

hf <- fread(file.path(OUT,"01_data","T06_hermes_hf_instrument_subset.csv"))
maf <- fread(file.path(OUT,"01_data","T07_snp_maf.csv"))
maf[, `:=`(maf=as.numeric(maf), n_alleles=as.integer(n_alleles))]
maf <- maf[!is.na(maf)]
say("HERMES subset rows: ", nrow(hf), " ; maf rows: ", nrow(maf))

## clean outcome: single-base alleles only, unique rsid (keep lowest p)
hf <- hf[nchar(effect_allele)==1 & nchar(other_allele)==1]
hf <- hf[!is.na(beta) & !is.na(se) & se > 0]
setorder(hf, rsid, pval); hf <- hf[!duplicated(rsid)]
hf <- merge(hf, maf[, .(rsid, maf, n_alleles)], by="rsid", all.x=TRUE)
say("outcome after cleaning: ", nrow(hf), " ; with maf: ", sum(!is.na(hf$maf)))

## ---------- build exposure tables ----------
eq <- fread(file.path(OUT,"01_data","T02_eqtlgen_instruments_pruned.csv"))
eq <- merge(eq, maf[, .(rsid, maf, n_alleles)], by="rsid", all.x=TRUE)
say("eQTLGen pruned: ", nrow(eq), " ; lacking maf: ", sum(is.na(eq$maf)), " ; multiallelic: ", sum(eq$n_alleles>2, na.rm=TRUE))
eq <- eq[!is.na(maf) & (is.na(n_alleles) | n_alleles <= 2) & maf > 0 & maf < 1]
## z -> per-allele beta (two-sampleMR convention for z-only eQTL summary data)
eq[, denom := sqrt(2*maf*(1-maf)*(n + z^2))]
eq[, `:=`(beta = z/denom, se = 1/denom)]
say("eQTLGen usable instruments: ", nrow(eq), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(eq$gene==g))), collapse=" "), ")")

gt <- fread(file.path(OUT,"01_data","T05_gtex_instruments_pruned.csv"))
gt <- merge(gt, maf[, .(rsid, maf_all=maf, n_alleles)], by="rsid", all.x=TRUE)
## GTEx betas are measured directly, so multiallelic records are kept (MAF only used for Steiger)
say("GTEx usable instruments: ", nrow(gt), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(gt$gene==g))), collapse=" "), ")")

## allele-convention sanity check between the two exposure resources (GTEx NES is per ALT allele)
cm <- merge(eq[, .(gene, rsid, ea=effect_allele, oa=other_allele, z)],
            gt[, .(gene, rsid, ea2=effect_allele, oa2=other_allele, nes=beta)], by=c("gene","rsid"))
if (nrow(cm)) {
  cm[, conc := fifelse(ea==ea2 & oa==oa2, sign(z)==sign(nes),
              fifelse(ea==oa2 & oa==ea2, sign(z)==-sign(nes), NA))]
  say("allele-convention check on shared SNPs: ", paste(sapply(GENES, function(g) {s <- cm[gene==g]; paste0(g,"=",sum(s$conc,na.rm=TRUE),"/",nrow(s))}), collapse=" "))
  fwrite(cm, file.path(OUT,"01_tables","T13_allele_concordance.csv"))
}

exposures <- list(eQTLGen_blood = eq[, .(gene, rsid, effect_allele, other_allele, beta, se, pval, fstat, n)],
                  GTEx_heartLV = gt[, .(gene, rsid, effect_allele, other_allele, beta, se, pval, fstat, n=386)])

res_main <- list(); res_snp <- list(); res_loo <- list(); res_harm <- list(); res_steig <- list()

for (exname in names(exposures)) {
  ex <- exposures[[exname]]
  for (g in GENES) {
    e <- ex[gene == g]
    if (!nrow(e)) { say(exname, " / ", g, ": no instrument"); next }
    hm <- harmonise_mr(e, hf, label = paste(exname, g, sep="|"))
    say(exname, " / ", g, " harmonise: ", paste(names(hm$summary)[-1], unlist(hm$summary[1,-1]), sep="=", collapse=" "))
    d <- hm$data
    if (!nrow(d)) next
    res_harm[[paste(exname,g)]] <- d[, .(exposure=exname, gene=g, rsid, effect_allele, other_allele,
                                         bx, sx, by, sy, fstat)]
    ## per-SNP Wald ratios
    wr <- d[, { w <- wald_ratio(bx, by, sx, sy); .(beta=w["beta"], se=w["se"]) }, by=.(rsid)]
    wr[, `:=`(exposure=exname, gene=g, pval=2*pnorm(-abs(beta/se)))]
    res_snp[[paste(exname,g)]] <- wr
    ## gene level
    rows <- list()
    if (nrow(d) == 1) {
      w <- wald_ratio(d$bx, d$by, d$sx, d$sy)
      rows[[1]] <- data.table(method="Wald ratio", exposure=exname, gene=g, nsnp=1,
                              beta=w["beta"], se=w["se"], pval=2*pnorm(-abs(w["beta"]/w["se"])))
    } else {
      iv <- ivw_fixed(d$bx, d$by, d$sx, d$sy, label=g)
      setnames(iv, "label", "gene"); iv[, exposure := exname]
      rows[[1]] <- iv
      if (nrow(d) >= 3) {
        eg <- mr_egger(d$bx, d$by, d$sx, d$sy, label=g); setnames(eg,"label","gene"); eg[, exposure := exname]; rows[[2]] <- eg
        wm <- mr_wmedian(d$bx, d$by, d$sx, d$sy, label=g); setnames(wm,"label","gene"); wm[, exposure := exname]; rows[[3]] <- wm
        lo <- mr_loo(d$bx, d$by, d$sx, d$sy, label=g); lo[, `:=`(exposure=exname, gene=g)]
        lo[, rsid := c(d$rsid, "All SNPs")]; res_loo[[paste(exname,g)]] <- lo
      }
    }
    r <- rbindlist(rows, fill=TRUE)
    r[, `:=`(F_mean = mean(d$fstat), F_min = min(d$fstat), n_used = nrow(d))]
    ord <- r[method=="Wald ratio"]
    if (!nrow(ord)) ord <- r[method=="IVW"][1]
    res_main[[paste(exname,g)]] <- r
    ## Steiger on the lead SNP
    lead <- d[which.max(fstat)]
    pm <- maf[rsid==lead$rsid, maf]
    st <- steiger(e[rsid==lead$rsid, pval][1], e[rsid==lead$rsid, n][1] %||% 31684,
                  lead$by, ifelse(length(pm), pm, NA), NCASE, NCONTROL, PREV)
    res_steig[[paste(exname,g)]] <- data.table(exposure=exname, gene=g, lead_snp=lead$rsid,
                                               r_exp=st["r_exp"], r_out=st["r_out"],
                                               steiger_z=st["z"], steiger_p=st["pval"],
                                               direction=ifelse(st["r_exp"] > st["r_out"], "exposure>outcome", "outcome>exposure"))
  }
}

main <- rbindlist(res_main, fill=TRUE)
main[, `:=`(OR = exp(beta), OR_lower = exp(beta-1.96*se), OR_upper = exp(beta+1.96*se),
            pval_re = fifelse(is.na(pval_re), pval, pval_re))]
fwrite(main, file.path(OUT,"01_tables","T08_mr_main_results.csv"))
fwrite(rbindlist(res_snp), file.path(OUT,"01_tables","T09_mr_per_snp_wald.csv"))
fwrite(rbindlist(res_harm), file.path(OUT,"01_data","T10_harmonised_snps.csv"))
if (length(res_loo)) fwrite(rbindlist(res_loo, fill=TRUE), file.path(OUT,"01_tables","T11_mr_leave_one_out.csv"))
if (length(res_steig)) fwrite(rbindlist(res_steig, fill=TRUE), file.path(OUT,"01_tables","T12_steiger.csv"))

say("---- headline (IVW / Wald) ----")
print(main[method %in% c("IVW","Wald ratio"), .(exposure, gene, nsnp, beta=round(beta,3), se=round(se,3),
                                                pval=signif(pval,3), OR=round(OR,2), F_mean=round(F_mean,1))])
if (length(warnings())) say("WARNINGS: ", paste(unique(names(warnings())), collapse=" | "))
writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_113.txt"))
close(LOG)
cat("113 done\n")
