## 134_locus_ld_relations.R -- robust LD relations at 9p13.3 (no fine-mapping assumptions)
## For each gene's lead cis-eQTL SNP: r2 with the six MR instruments, and r2 with
## the lead SNP of every other gene. This quantifies "the instruments do not tag
## the AQP3 signal" without relying on credible sets, which are resolution-limited
## with a 503-sample reference panel.
## ASCII only.

lib <- "E:/paper test/.workbuddy/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))

ROOT <- "E:/another data file/GEObulk analysis project"
LD   <- file.path(ROOT, "results", "Coloc", "05_ld")
DAT  <- file.path(ROOT, "results", "Coloc", "01_data")
TAB  <- file.path(ROOT, "results", "Coloc", "01_tables")
LOG  <- file.path(ROOT, "results", "Coloc", "03_logs")
sink(file.path(LOG, "134_locus_ld_relations.log"), split = TRUE)

GENES <- c("AQP3", "RP11-176F3.7", "GALT", "DCTN3", "NOL6", "RP11-384P7.7", "BAG1")
INSTR <- c("rs10971815", "rs116972937", "rs12551364", "rs3814511", "rs573538", "rs626069")

align <- read.csv(file.path(LD, "align.csv"), stringsAsFactors = FALSE)
flipmap <- setNames(align$flip, align$snp)

## one big LD matrix is not stored, so rebuild pairwise r from the per-gene matrices
mats <- list()
for (g in GENES) {
  fn <- file.path(LD, sprintf("ld_%s.csv", gsub("[.]", "_", g)))
  if (!file.exists(fn)) next
  s <- read.csv(file.path(LD, sprintf("ld_%s_snps.csv", gsub("[.]", "_", g))),
                stringsAsFactors = FALSE)$snp
  R <- as.matrix(read.csv(fn, header = FALSE)); dimnames(R) <- list(s, s)
  mats[[g]] <- R
}

## lead SNP per gene (max |z| among its own SNPs)
leads <- data.frame()
for (g in names(mats)) {
  d <- read.csv(file.path(DAT, sprintf("coloc_%s_FinnGen.csv", gsub("[.]", "_", g))),
                stringsAsFactors = FALSE)
  s <- rownames(mats[[g]])
  m <- match(s, d$snp)
  z <- rep(NA_real_, length(s)); z[!is.na(m)] <- d$z1[m[!is.na(m)]]
  i <- which.max(abs(z))
  leads <- rbind(leads, data.frame(gene = g, lead = s[i], z = z[i], stringsAsFactors = FALSE))
}
print(leads)

## r2 between gene leads: need both SNPs in one matrix; use the matrix of the gene
## that contains both, otherwise fall back to computing from the union matrix
getR <- function(a, b) {
  for (g in names(mats)) {
    R <- mats[[g]]
    if (a %in% rownames(R) && b %in% rownames(R)) return(R[a, b])
  }
  NA_real_
}

pairs <- data.frame()
for (i in 1:(nrow(leads) - 1)) for (j in (i + 1):nrow(leads)) {
  pairs <- rbind(pairs, data.frame(
    gene_a = leads$gene[i], snp_a = leads$lead[i],
    gene_b = leads$gene[j], snp_b = leads$lead[j],
    r2 = getR(leads$lead[i], leads$lead[j])^2, stringsAsFactors = FALSE))
}
pairs <- pairs[order(-pairs$r2), ]
write.csv(pairs, file.path(TAB, "T10_gene_lead_r2.csv"), row.names = FALSE)
cat("\n== r2 between gene lead SNPs ==\n"); print(pairs, row.names = FALSE, digits = 3)

## instruments vs each gene's lead
out <- data.frame()
for (g in names(mats)) {
  R <- mats[[g]]; s <- rownames(R)
  lp <- leads$lead[leads$gene == g]
  for (r in INSTR) {
    if (!(r %in% s)) next
    out <- rbind(out, data.frame(gene = g, lead = lp, instrument = r,
                                 r2 = R[r, lp]^2, stringsAsFactors = FALSE))
  }
}
write.csv(out, file.path(TAB, "T11_instrument_vs_lead_r2.csv"), row.names = FALSE)
cat("\n== instruments vs gene lead SNPs ==\n")
print(reshape(out[, c("instrument", "gene", "r2")], idvar = "instrument",
              timevar = "gene", direction = "wide"), row.names = FALSE, digits = 3)

## r2 among the instruments themselves (are they independent?)
R <- NULL
for (g in names(mats)) {
  s <- rownames(mats[[g]])
  if (sum(INSTR %in% s) >= 3) { R <- mats[[g]]; break }
}
if (!is.null(R)) {
  ins <- INSTR[INSTR %in% rownames(R)]
  M <- R[ins, ins]^2
  write.csv(M, file.path(TAB, "T12_instrument_r2_matrix.csv"))
  cat("\n== r2 among instruments ==\n"); print(round(M, 3))
}

## is the AQP3 lead separated from the top SNPs by real LD or by panel noise?
if ("AQP3" %in% names(mats)) {
  R <- mats[["AQP3"]]
  top <- names(sort(-abs(R[leads$lead[leads$gene == "AQP3"], ]))[1:8])
  cat("\n== top 8 SNPs by |r| with the AQP3 lead ==\n")
  print(data.frame(snp = top, r = round(R[leads$lead[leads$gene == "AQP3"], top], 3),
                   r2 = round(R[leads$lead[leads$gene == "AQP3"], top]^2, 3)), row.names = FALSE)
}
sink()
