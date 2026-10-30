## 119_mr_locus_pleiotropy.R -- are the AQP3 instruments also cis-eQTLs for neighbouring genes?
## Uses the local eQTLGen significant-pairs table (FDR < 0.05) to list every gene each
## instrument SNP regulates in cis, i.e. a direct read-out of locus-level pleiotropy.
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
use("data.table")

MDIR <- file.path(DIR$results, "MR")
eq <- fread(file.path(DIR$raw, "mr", "eqtlgen_cis_eqtl_FDR005.txt.gz"),
            sep = "\t", header = TRUE, showProgress = FALSE,
            select = c("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                       "OtherAllele", "Zscore", "GeneSymbol", "NrSamples"))
setnames(eq, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","Zscore","GeneSymbol","NrSamples"),
             c("pval","rsid","chr","pos","effect_allele","other_allele","z","gene","n"))

## instrument sets actually used in the MR
inst <- rbind(
  fread(file.path(MDIR, "01_data", "T02_eqtlgen_instruments_pruned.csv"))[, .(exposure="eQTLGen_blood", gene_assayed=gene, rsid, effect_allele, other_allele, z_instrument=z)],
  fread(file.path(MDIR, "01_data", "T05_gtex_instruments_pruned.csv"))[, .(exposure="GTEx_heartLV", gene_assayed=gene, rsid, effect_allele, other_allele, z_instrument=NA_real_)]
)
inst[, palin := paste0(toupper(effect_allele), toupper(other_allele)) %in% c("AT","TA","CG","GC")]
inst <- inst[palin == FALSE]

rs <- unique(inst$rsid)
cat("instrument SNPs:", length(rs), "\n")

## every significant cis target of those SNPs
tg <- eq[rsid %chin% rs & pval < 5e-8]
setorder(tg, rsid, pval)
cat("significant cis targets (p < 5e-8):", nrow(tg), "\n")

## align the target-gene Z to the instrument's effect allele
m <- merge(inst[, .(rsid, gene_assayed, exposure, inst_effect=effect_allele, inst_other=other_allele, z_instrument)],
           tg[, .(rsid, target_gene=gene, t_effect=effect_allele, t_other=other_allele, z_target=z, p_target=pval, n_target=n)],
           by = "rsid", allow.cartesian = TRUE)
m[, same_allele := toupper(inst_effect) == toupper(t_effect)]
m[, z_target_aligned := ifelse(same_allele, z_target, -z_target)]
m[, is_assayed_gene := target_gene == gene_assayed]

cat("\n== per instrument: how many genes does it regulate in cis? ==\n")
print(m[, .(n_targets = uniqueN(target_gene),
            targets = paste(sort(unique(target_gene)), collapse = ",")), by = .(exposure, rsid)])

cat("\n== targets of the AQP3 instruments that are NOT AQP3 ==\n")
print(m[gene_assayed == "AQP3" & is_assayed_gene == FALSE,
        .(exposure, rsid, target_gene, z_target_aligned = round(z_target_aligned, 2),
          p_target = signif(p_target, 3))][order(rsid, p_target)])

cat("\n== AQP7 specifically ==\n")
aqp7 <- m[target_gene == "AQP7"]
if (nrow(aqp7)) {
  print(aqp7[, .(exposure, rsid, gene_assayed, z_target_aligned = round(z_target_aligned, 2),
                 p_target = signif(p_target, 3))])
  cat("instruments of AQP3 that are also AQP7 eQTLs: ",
      uniqueN(aqp7[gene_assayed == "AQP3"]$rsid), " / ", uniqueN(inst[gene_assayed == "AQP3"]$rsid), "\n")
} else cat("no AQP7 signal among the instruments\n")

fwrite(m, file.path(MDIR, "01_tables", "T22_locus_cis_targets.csv"))
write_session(file.path(MDIR, "03_logs"), "119_locus_pleiotropy")
cat("\nwrote T22_locus_cis_targets.csv\n")
