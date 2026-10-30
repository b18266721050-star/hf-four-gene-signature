## 120_table1.R -- Table 1 (cohort characteristics) for the manuscript
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
use("data.table")

MDIR <- new_dir(DIR$results, "Manuscript")

## ---- discovery cohort -----------------------------------------------------
md <- fread(file.path(DIR$results, "Batch", "01_tables", "T01_sample_batch_metadata.csv"))
md[, group := factor(group, levels = c("NF", "DCM", "ICM"))]
cat("discovery n =", nrow(md), "\n"); print(table(md$group))

agg <- function(d) {
  data.table(n = nrow(d),
             age_mean = round(mean(d$age_num, na.rm = TRUE), 1),
             age_sd = round(sd(d$age_num, na.rm = TRUE), 1),
             age_median = round(median(d$age_num, na.rm = TRUE), 1),
             age_q1 = round(quantile(d$age_num, .25, na.rm = TRUE), 1),
             age_q3 = round(quantile(d$age_num, .75, na.rm = TRUE), 1),
             n_male = sum(d$sex == "Male", na.rm = TRUE),
             n_female = sum(d$sex == "Female", na.rm = TRUE),
             n_age_missing = sum(is.na(d$age_num)))
}
disc <- rbind(agg(md)[, group := "All"], agg(md[group == "NF"])[, group := "NF"],
              agg(md[group == "DCM"])[, group := "DCM"], agg(md[group == "ICM"])[, group := "ICM"])
disc[, male_pct := round(100 * n_male / n, 1)]
setcolorder(disc, c("group", "n", "age_mean", "age_sd", "age_median", "age_q1", "age_q3",
                    "n_male", "n_female", "male_pct", "n_age_missing"))
cat("\n== discovery ==\n"); print(disc)

## ---- external cohorts -----------------------------------------------------
ext0 <- fread(file.path(DIR$results, "ExternalValidation", "01_tables", "T10_cohort_baseline.csv"))
auc  <- fread(file.path(DIR$results, "ExternalValidation", "01_tables", "T02_external_auc.csv"))
auc2 <- auc[subset == "all disease vs NF" & version == "V2", .(cohort, n_analysed = n,
                                                               n_case, n_ctrl, auc = round(auc, 3))]
ext <- merge(ext0, auc2, by = "cohort", all.x = TRUE)
ext[, age_text := fifelse(is.na(age_median), "n/a",
                          paste0(age_median, " (", age_range, ")"))]
ext[, male_text := fifelse(is.na(n_male), "n/a", paste0(n_male, " (", round(100*n_male/(n_male+n_female)), "%)"))]
cat("\n== external ==\n"); print(ext[, .(cohort, platform, n, NF, DCM, ICM, n_analysed, auc, age_text, male_text)])

## ---- single-nucleus cohort ------------------------------------------------
sg <- fread(file.path(DIR$results, "Singlecel", "fourgene", "tables", "64_sample_group.csv"))
comp <- fread(file.path(DIR$results, "Singlecel", "fourgene", "tables", "68_sc_celltype_composition_donor.csv"))
cat("\n== snRNA-seq ==", "\n"); print(table(sg$group))
n_nuclei <- sum(comp$n_nuclei)
n_donor <- uniqueN(comp$sample)
cat("donors =", n_donor, " nuclei =", n_nuclei, " cell types =", uniqueN(comp$celltype), "\n")

## ---- manuscript-ready table ----------------------------------------------
age_txt <- function(g) sprintf("%.1f (%.1f)", disc[group == g, age_mean], disc[group == g, age_sd])
t1 <- data.table(
  Cohort = c("GSE57338 (discovery)", "  NF", "  DCM", "  ICM",
             "GSE5406", "GSE116250", "GSE55296", "GSE292067 (snRNA-seq)",
             "HERMES (GWAS)", "FinnGen R11 (GWAS)"),
  Platform = c("Affymetrix HuGene-1.0-ST (GPL11532)", "", "", "",
               "Affymetrix HG-U133A (GPL96)", "Illumina HiSeq 2500 (RNA-seq)",
               "RNA-seq (normalised counts)", "10x snRNA-seq (3' GEX)",
               "GWAS meta-analysis (29 studies)", "GWAS (Finnish founder population)"),
  n = c(nrow(md), disc[group=="NF", n], disc[group=="DCM", n], disc[group=="ICM", n],
        ext[cohort=="GSE5406", n_analysed], ext[cohort=="GSE116250", n],
        ext[cohort=="GSE55296", n_analysed], n_donor, 977323, 453733),
  Cases = c(177, NA, 82, 95, ext[cohort=="GSE5406", n_case], ext[cohort=="GSE116250", n_case],
            ext[cohort=="GSE55296", n_case], 17, 47309, 33250),
  Controls = c(136, NA, NA, NA, ext[cohort=="GSE5406", n_ctrl], ext[cohort=="GSE116250", n_ctrl],
               ext[cohort=="GSE55296", n_ctrl], 8, 930014, 420483),
  Age = c(age_txt("All"), age_txt("NF"), age_txt("DCM"), age_txt("ICM"),
          "n/a", ext[cohort=="GSE116250", age_text], "n/a", "n/a", "n/a", "n/a"),
  Male = c(paste0(disc[group=="All", n_male], " (", disc[group=="All", male_pct], "%)"),
           paste0(disc[group=="NF", n_male], " (", disc[group=="NF", male_pct], "%)"),
           paste0(disc[group=="DCM", n_male], " (", disc[group=="DCM", male_pct], "%)"),
           paste0(disc[group=="ICM", n_male], " (", disc[group=="ICM", male_pct], "%)"),
           "n/a", ext[cohort=="GSE116250", male_text], "n/a", "n/a", "n/a", "n/a"),
  Note = c("left ventricular myocardium", "non-failing donor hearts", "dilated cardiomyopathy",
           "ischaemic cardiomyopathy", "no age/sex in GEO record", "2 of 14 controls lacked usable expression and were excluded from scoring (62 scored)",
           "3 samples without disease annotation excluded", "25 donors: 8 NF / 6 NICM / 11 DoxCM; 150,750 nuclei, 12 cell types",
           "Shah et al. 2020; GCST009541", "strict I9_HEARTFAIL endpoint")
)
fwrite(t1, file.path(MDIR, "T01_baseline.csv"))

t1d <- copy(t1)
for (cl in c("n", "Cases", "Controls")) {
  t1d[[cl]] <- ifelse(is.na(t1[[cl]]), "NA", format(t1[[cl]], big.mark = ",", trim = TRUE))
}
md_lines <- c("# Table 1. Characteristics of the cohorts used in this study", "",
              paste0("| ", paste(names(t1d), collapse = " | "), " |"),
              paste0("| ", paste(rep("---", ncol(t1d)), collapse = " | "), " |"))
for (i in seq_len(nrow(t1d))) {
  md_lines <- c(md_lines, paste0("| ", paste(unlist(t1d[i]), collapse = " | "), " |"))
}
writeLines(md_lines, file.path(MDIR, "Table1.md"))
cat("\nwrote", file.path(MDIR, "Table1.md"), "\n")
print(t1[, .(Cohort, n, Cases, Controls, Age, Male)])
write_session(file.path(MDIR, "logs"), "120_table1")
