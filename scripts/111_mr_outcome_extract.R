## 111_mr_outcome_extract.R -- subset HERMES heart failure GWAS to the instrument SNPs
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

ROOT <- "E:/another data file/GEObulk analysis project"
RAW  <- file.path(ROOT, "raw_data/mr")
OUT  <- file.path(ROOT, "results/MR")
LOG <- file(file.path(OUT,"03_logs","111_mr_outcome.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

say("== 111 HERMES HF outcome extraction ==")
inst <- unique(rbind(
  fread(file.path(OUT,"01_data","T02_eqtlgen_instruments_pruned.csv"))[, .(rsid, chr, pos)],
  if (file.exists(file.path(OUT,"01_data","T05_gtex_instruments_pruned.csv")))
    fread(file.path(OUT,"01_data","T05_gtex_instruments_pruned.csv"))[, .(rsid, chr, pos)] else data.table(rsid=character(), chr=integer(), pos=integer())
))
say("instrument rsIDs: ", nrow(inst), " unique ", uniqueN(inst$rsid))

f <- file.path(RAW, "hermes_HF_GCST009541_h.tsv.gz")
say("reading ", f, " (", round(file.size(f)/1e6,1), " MB gz)")
sel <- c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele",
         "hm_beta","standard_error","p_value","hm_effect_allele_frequency","hm_code")
hf <- fread(f, sep="\t", select=sel, showProgress=FALSE)
say("HERMES rows: ", nrow(hf))
say("columns used: ", paste(sel, collapse=","))

sub <- hf[hm_rsid %chin% inst$rsid]
say("matched instrument SNPs in HERMES: ", nrow(sub), " / ", uniqueN(inst$rsid), " (",
    round(100*nrow(sub)/uniqueN(inst$rsid),1), "%)")
say("hm_code distribution: ", paste(names(table(sub$hm_code)), table(sub$hm_code), sep=":", collapse=" "))
say("eAF available: ", sum(!is.na(sub$hm_effect_allele_frequency)), " / ", nrow(sub))

setnames(sub, c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele",
                "hm_beta","standard_error","p_value","hm_effect_allele_frequency","hm_code"),
              c("rsid","chr38","pos38","other_allele","effect_allele","beta","se","pval","eaf","hm_code"))
fwrite(sub, file.path(OUT,"01_data","T06_hermes_hf_instrument_subset.csv"))

## genome-wide summary of the source GWAS for the Methods section
say("HERMES total SNPs: ", nrow(hf), "; median SE: ", round(median(hf$standard_error, na.rm=TRUE),4))

writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_111.txt"))
close(LOG)
cat("111 done\n")
