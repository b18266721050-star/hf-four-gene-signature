###############################################################################
## 95_candidate70_composition_link.R
## Why do all 70 candidate genes move together? Quantify how tightly the
## WGCNA-yellow x DEG intersection tracks the myeloid compartment.
##
## Output: results/Composition/01_tables/T17_candidate70_vs_composition.csv
##
## Run:  Rscript scripts/95_candidate70_composition_link.R
###############################################################################

.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages(library(data.table))

proj_dir <- "E:/another data file/GEObulk analysis project"
out_dir <- file.path(proj_dir, "results", "Composition", "01_tables")

comp <- fread(file.path(proj_dir, "results", "Composition", "01_tables",
                        "T04_composition_merged.csv"))
e <- fread(file.path(proj_dir, "clean_data", "GSE57338_clean_expression_matrix.csv"),
           data.table = FALSE)
m <- as.matrix(e[, -1]); rownames(m) <- e[[1]]
m <- m[, comp$sample]

cand <- fread(file.path(proj_dir, "results", "HubGene", "GSE57338_candidate_genes.csv"))
g70 <- cand[in_WGCNA_set == TRUE & in_DEG_set == TRUE, GeneSymbol]

feats <- c(Macrophages = "xcell_Macrophages", Fibroblasts = "xcell_Fibroblasts",
           Endothelial = "xcell_Endothelial_cells", Myocytes = "xcell_Myocytes",
           CD8T = "xcell_CD8_T_cells", ImmuneScore = "xcell_ImmuneScore")
out <- rbindlist(lapply(g70, function(g) {
  rbindlist(lapply(names(feats), function(nm) {
    data.table(GeneSymbol = g, compartment = nm,
               rho = round(cor(m[g, ], comp[[feats[[nm]]]], method = "spearman"), 3))
  }))
}))
write.csv(out, file.path(out_dir, "T17_candidate70_vs_composition.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

summ <- out[, .(median_rho = round(median(rho), 2),
                q25 = round(quantile(rho, .25), 2),
                q75 = round(quantile(rho, .75), 2),
                n_pos_above_0.3 = sum(rho > 0.3),
                n_neg_below_minus0.3 = sum(rho < -0.3),
                n = .N), by = compartment]
print(summ)
