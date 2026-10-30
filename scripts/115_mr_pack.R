## 115_mr_pack.R -- pack MR tables into one xlsx + write English figure legends
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({library(data.table); library(openxlsx)})

ROOT <- "E:/another data file/GEObulk analysis project"
OUT  <- file.path(ROOT, "results/MR")
LOG <- file(file.path(OUT,"03_logs","115_mr_pack.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

rd <- function(p) if (file.exists(p)) fread(p) else NULL
t01 <- rd(file.path(OUT,"01_data","T01_eqtlgen_instruments_raw.csv"))
t02 <- rd(file.path(OUT,"01_data","T02_eqtlgen_instruments_pruned.csv"))
t04 <- rd(file.path(OUT,"01_data","T04_gtex_instruments_raw.csv"))
t05 <- rd(file.path(OUT,"01_data","T05_gtex_instruments_pruned.csv"))
t08 <- rd(file.path(OUT,"01_tables","T08_mr_main_results.csv"))
t09 <- rd(file.path(OUT,"01_tables","T09_mr_per_snp_wald.csv"))
t10 <- rd(file.path(OUT,"01_data","T10_harmonised_snps.csv"))
t11 <- rd(file.path(OUT,"01_tables","T11_mr_leave_one_out.csv"))
t12 <- rd(file.path(OUT,"01_tables","T12_steiger.csv"))
t13 <- rd(file.path(OUT,"01_tables","T13_allele_concordance.csv"))
t14 <- rd(file.path(OUT,"01_tables","T14_mr_summary.csv"))
t18 <- rd(file.path(OUT,"01_tables","T18_finngen_loo_AQP3.csv"))
t19 <- rd(file.path(OUT,"01_tables","T19_hermes_vs_finngen.csv"))
t20 <- rd(file.path(OUT,"01_tables","T20_two_outcome_multiplicity.csv"))
t21 <- rd(file.path(OUT,"01_tables","T21_cross_cohort_meta.csv"))
t15 <- rd(file.path(OUT,"01_tables","T15_finngen_mr_results.csv"))
t16 <- rd(file.path(OUT,"01_tables","T16_finngen_harmonised.csv"))
t17 <- rd(file.path(OUT,"01_tables","T17_finngen_per_snp.csv"))
t22 <- rd(file.path(OUT,"01_tables","T22_locus_cis_targets.csv"))

wb <- createWorkbook()
addWorksheet(wb, "README")
rmk <- data.frame(x = c(
  "Mendelian randomisation: genetically predicted expression of AQP3 / GGT5 / ISLR / LAD1 and heart failure",
  "",
  "Exposure 1 (primary): eQTLGen whole-blood cis-eQTL, n = 31,684 (Vosa et al. 2021).",
  "Exposure 2 (tissue-matched sensitivity): GTEx v8 Heart_Left_Ventricle cis-eQTL, n = 386.",
  "Outcome: HERMES heart failure GWAS meta-analysis (Shah et al. 2020), 47,309 cases / 930,014 controls,",
  "         summary statistics obtained from GWAS Catalog GCST009541 (harmonised, GRCh38).",
  "",
  "Instrument selection: cis (within 1 Mb of the gene), p < 5e-8, distance-based pruning at 100 kb",
  "(no LD reference panel was available offline; residual LD between retained SNPs is a stated limitation).",
  "eQTLGen reports Z scores only; per-allele beta and SE were reconstructed as",
  "  se = 1 / sqrt(2*MAF*(1-MAF)*(N + Z^2)),  beta = Z * se,  with MAF from Ensembl (1000G/gnomAD).",
  "Palindromic SNPs (A/T, C/G) were dropped because HERMES harmonised files do not carry allele frequency.",
  "",
  "Estimators: Wald ratio (1 SNP), inverse-variance weighted (fixed and multiplicative random effects),",
  "MR-Egger (with residual-SE correction and intercept test), weighted median (2,000-parameter bootstrap SE),",
  "Cochran Q, leave-one-out, mean F statistic, Steiger directionality.",
  "All estimators were implemented in-house (scripts/mr_core.R) because TwoSampleMR could not be",
  "installed in this offline R environment (no Rtools, no OpenGWAS JWT token).",
  "",
  "Sheets: M01 main results | M02 per-SNP Wald ratios | M03 harmonised SNPs | M04 leave-one-out |",
  "        M05 Steiger | M06 eQTLGen instruments (raw) | M07 eQTLGen instruments (pruned) |",
  "        M08 GTEx instruments (raw) | M09 GTEx instruments (pruned)"), stringsAsFactors = FALSE)
writeData(wb, "README", rmk, colNames = FALSE)
setColWidths(wb, "README", 1, 120)

add_sheet <- function(name, dt) {
  if (is.null(dt) || !nrow(dt)) return(invisible(NULL))
  addWorksheet(wb, name); writeData(wb, name, dt)
  setColWidths(wb, name, 1:ncol(dt), "auto")
}
add_sheet("M01_summary", t14); add_sheet("M01b_main", t08); add_sheet("M02_per_snp", t09); add_sheet("M03_harmonised", t10)
add_sheet("M04_loo", t11); add_sheet("M05_steiger", t12)
add_sheet("M06_eqtlgen_raw", t01); add_sheet("M07_eqtlgen_pruned", t02)
add_sheet("M08_gtex_raw", t04); add_sheet("M09_gtex_pruned", t05)
add_sheet("M10_allele_check", t13)
add_sheet("M11_finngen_main", t15); add_sheet("M12_finngen_snps", t16)
add_sheet("M13_finngen_persnp", t17); add_sheet("M14_finngen_loo", t18)
add_sheet("M15_two_outcome", t20); add_sheet("M16_cohort_compare", t19)
add_sheet("M17_cross_cohort_meta", t21)
add_sheet("M18_locus_targets", t22)
saveWorkbook(wb, file.path(OUT,"MR_supplementary_tables.xlsx"), overwrite=TRUE)
say("wrote MR_supplementary_tables.xlsx with ", length(names(wb)), " sheets")

leg <- c(
"Figure 8. Mendelian randomisation of genetically predicted expression of the four candidate genes and heart failure.",
"(A) Gene-level causal estimates (log odds ratio of heart failure per 1-SD higher genetically predicted expression)",
"with 95% confidence intervals, shown separately for the primary blood-based exposure (eQTLGen, n = 31,684) and the",
"tissue-matched sensitivity exposure (GTEx v8 left ventricle, n = 386), against the HERMES heart failure GWAS.",
"Numbers next to the intervals are the number of cis-eQTL instruments retained after pruning.",
"(B, C) Scatter plots of per-SNP effects on expression (x axis) versus per-SNP effects on heart failure (y axis) for",
"AQP3 (eQTLGen instruments) and GGT5 (GTEx instruments). Point size is proportional to the inverse variance of the",
"SNP-heart failure association. Solid line, inverse-variance weighted slope; dashed line, MR-Egger slope; dotted line,",
"weighted median slope; each line spans the observed range of instrument effects.",
"(D) Leave-one-out analysis for AQP3: each point is the inverse-variance weighted estimate obtained after omitting one",
"SNP; the last row (All SNPs) is the estimate using all instruments. Horizontal bars are 95% confidence intervals.",
"(E) The same gene-level estimates repeated in an independent outcome GWAS (FinnGen release 11, I9_HEARTFAIL,",
"33,250 cases and 420,483 controls; blue, HERMES; red, FinnGen). AQP3 was associated with heart failure in FinnGen",
"(odds ratio 0.79 per 1-SD, 95% confidence interval 0.69-0.90), but the effect was smaller and non-significant in HERMES",
"and the between-cohort heterogeneity was substantial.",
"MR, Mendelian randomisation; SD, standard deviation.")
writeLines(leg, file.path(OUT,"Figure8_legend.txt"))
say("wrote Figure8_legend.txt")
writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_115.txt"))
close(LOG)
cat("115 done\n")
