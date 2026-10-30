## 110_mr_instruments.R -- build cis-eQTL instruments for AQP3/GGT5/ISLR/LAD1
## exposure 1: eQTLGen whole blood cis-eQTL (n ~ 31,684)
## exposure 2: GTEx v8 Heart_Left_Ventricle cis-eQTL (n ~ 386)
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))
set.seed(20251029)

ROOT <- "E:/another data file/GEObulk analysis project"
RAW  <- file.path(ROOT, "raw_data/mr")
OUT  <- file.path(ROOT, "results/MR")
for (d in c("01_data","01_tables","02_figures","03_logs")) dir.create(file.path(OUT,d), recursive=TRUE, showWarnings=FALSE)
LOG <- file(file.path(OUT,"03_logs","110_mr_instruments.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

GENES   <- c("AQP3","GGT5","ISLR","LAD1")
PTHR    <- 5e-8
PRUNE_BP <- 100000   # distance based pruning window (no LD reference panel available offline)

say("== 110 build instruments ==")
say("R ", getRversion(), " ; genes: ", paste(GENES, collapse=","))

## ---------------- 1. eQTLGen ----------------
f_eq <- file.path(RAW, "eqtlgen_cis_eqtl_FDR005.txt.gz")
say("reading ", f_eq, " (", round(file.size(f_eq)/1e6,1), " MB gz)")
eq <- fread(f_eq, sep="\t", header=TRUE, showProgress=FALSE)
say("eQTLGen rows read: ", nrow(eq))
say("eQTLGen columns: ", paste(names(eq), collapse=","))

need <- c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","Zscore","GeneSymbol","NrSamples")
miss <- setdiff(need, names(eq))
if (length(miss)) stop("missing columns: ", paste(miss, collapse=","))

eq <- eq[GeneSymbol %chin% GENES]
say("rows for 4 genes (FDR<0.05 file): ", nrow(eq))
say("per gene: ", paste(sapply(GENES, function(g) paste0(g,"=",sum(eq$GeneSymbol==g))), collapse=" "))
say("significance summary (p<5e-8): ", paste(sapply(GENES, function(g) paste0(g,"=",sum(eq$GeneSymbol==g & eq$Pvalue< PTHR))), collapse=" "))

eqi <- eq[Pvalue < PTHR, .(gene=GeneSymbol, rsid=SNP, chr=SNPChr, pos=SNPPos,
                           effect_allele=AssessedAllele, other_allele=OtherAllele,
                           z=Zscore, pval=Pvalue, n=NrSamples)]
eqi[, fstat := z^2]
eqi[, palin := paste0(toupper(effect_allele), toupper(other_allele)) %in% c("AT","TA","CG","GC")]
say("palindromic SNPs excluded from pruning: ", sum(eqi$palin))
say("instruments before pruning: ", nrow(eqi), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(eqi$gene==g))), collapse=" "), ")")

## distance based pruning: keep lowest-p SNP, drop any within PRUNE_BP
prune_dist <- function(dt, win) {
  dt <- copy(dt); setorder(dt, gene, pval)
  keep <- logical(nrow(dt))
  for (g in unique(dt$gene)) {
    idx <- which(dt$gene == g); kept <- numeric(0)
    for (i in idx) {
      if (!length(kept) || min(abs(dt$pos[i] - kept)) > win) { keep[i] <- TRUE; kept <- c(kept, dt$pos[i]) }
    }
  }
  dt[keep]
}
eqp <- prune_dist(eqi[palin==FALSE], PRUNE_BP)
say("instruments after ", PRUNE_BP/1000, "kb distance pruning: ", nrow(eqp),
    " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(eqp$gene==g))), collapse=" "), ")")
say("F statistics (after pruning): ", paste(sapply(GENES, function(g) { s <- eqp[gene==g]; if (!nrow(s)) paste0(g,"=NA") else paste0(g,"=",round(mean(s$fstat),1)) }), collapse=" "))

fwrite(eqi, file.path(OUT,"01_data","T01_eqtlgen_instruments_raw.csv"))
fwrite(eqp, file.path(OUT,"01_data","T02_eqtlgen_instruments_pruned.csv"))

## ---------------- 2. GTEx v8 heart LV ----------------
f_g <- file.path(RAW, "gtex_hlv_eqtl_4genes.csv")
if (file.exists(f_g)) {
  gx <- fread(f_g)
  say("GTEx rows: ", nrow(gx), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(gx$gene==g))), collapse=" "), ")")
  ## nes is the slope of inverse-normalised expression on the ALT allele of variantId
  gx[, z := qnorm(pval/2, lower.tail=FALSE) * sign(nes)]
  gx[, se := abs(nes)/abs(z)]
  gx[, fstat := z^2]
  gi <- gx[pval < PTHR, .(gene, rsid=snp, chr=as.integer(chrom), pos=as.integer(pos),
                          effect_allele=alt, other_allele=ref, beta=nes, se=se, z=z, pval=pval, fstat=fstat)]
  gi[, palin := paste0(toupper(effect_allele), toupper(other_allele)) %in% c("AT","TA","CG","GC")]
  say("GTEx genome-wide significant cis-eQTL: ", nrow(gi), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(gi$gene==g))), collapse=" "), ")")
  gip <- prune_dist(gi[palin==FALSE], PRUNE_BP)
  say("GTEx after pruning: ", nrow(gip), " (", paste(sapply(GENES, function(g) paste0(g,"=",sum(gip$gene==g))), collapse=" "), ")")
  fwrite(gx, file.path(OUT,"01_data","T03_gtex_hlv_all.csv"))
  fwrite(gi, file.path(OUT,"01_data","T04_gtex_instruments_raw.csv"))
  fwrite(gip, file.path(OUT,"01_data","T05_gtex_instruments_pruned.csv"))
} else say("GTEx file not found, skipped")

writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_110.txt"))
close(LOG)
cat("110 done\n")
