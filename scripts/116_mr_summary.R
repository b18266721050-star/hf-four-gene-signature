## 116_mr_summary.R -- concordance with the bulk direction, power, and a short numeric summary
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

ROOT <- "E:/another data file/GEObulk analysis project"
OUT  <- file.path(ROOT, "results/MR")
LOG <- file(file.path(OUT,"03_logs","116_mr_summary.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

GENES <- c("AQP3","GGT5","ISLR","LAD1")
main <- fread(file.path(OUT,"01_tables","T08_mr_main_results.csv"))
snp  <- fread(file.path(OUT,"01_tables","T09_mr_per_snp_wald.csv"))

deg <- fread(file.path(ROOT,"results/DEG/GSE57338_DEGall.csv"))
deg <- deg[GeneSymbol %chin% GENES, .(gene=GeneSymbol, log2FC=log2FoldChange, padj)]

prim <- main[exposure=="eQTLGen_blood" & method %in% c("IVW","Wald ratio"),
             .(gene, exposure, method, nsnp, beta, se, pval, F_mean)]
prim <- merge(GENES_dt <- data.table(gene=GENES), prim, by="gene", all.x=TRUE)
prim <- merge(prim, deg, by="gene", all.x=TRUE)
prim[, `:=`(OR = exp(beta), OR_lower = exp(beta-1.96*se), OR_upper = exp(beta+1.96*se),
            min_detectable_OR_80pwr = exp(2.802*se))]
prim[, concordant_with_bulk := fifelse(is.na(beta), NA,
      fifelse(sign(beta)==sign(log2FC), "yes", "no"))]

## per-SNP direction consistency for the multi-SNP gene
for (g in unique(snp$gene)) {
  s <- snp[exposure=="eQTLGen_blood" & gene==g]
  if (nrow(s) >= 3) {
    k <- sum(sign(s$beta) == sign(prim[gene==g, beta]))
    say(g, ": ", k, "/", nrow(s), " per-SNP ratios share the IVW sign; binomial p = ",
        signif(binom.test(k, nrow(s))$p.value, 3))
  }
}
say("primary table:")
print(prim[, .(gene, log2FC, nsnp, beta=round(beta,4), pval=signif(pval,3),
               OR=round(OR,3), OR_lower=round(OR_lower,3), OR_upper=round(OR_upper,3),
               min_detectable_OR_80pwr=round(min_detectable_OR_80pwr,3), concordant_with_bulk)])

fwrite(prim, file.path(OUT,"01_tables","T14_mr_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_116.txt"))
close(LOG)
cat("116 done\n")
