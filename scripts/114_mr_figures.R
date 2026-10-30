## 114_mr_figures.R -- figures for the MR module (Figure 8 composite + single panels)
.libPaths(c("D:/R/library/4.6", .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

ROOT <- "E:/another data file/GEObulk analysis project"
OUT  <- file.path(ROOT, "results/MR")
FIG  <- file.path(OUT,"02_figures")
LOG <- file(file.path(OUT,"03_logs","114_mr_figures.log"), open="wt")
say <- function(...) { m <- paste0(...); cat(m,"\n"); writeLines(m, LOG) }

main <- fread(file.path(OUT,"01_tables","T08_mr_main_results.csv"))
snp  <- fread(file.path(OUT,"01_tables","T09_mr_per_snp_wald.csv"))
harm <- fread(file.path(OUT,"01_data","T10_harmonised_snps.csv"))
loo  <- if (file.exists(file.path(OUT,"01_tables","T11_mr_leave_one_out.csv")))
          fread(file.path(OUT,"01_tables","T11_mr_leave_one_out.csv")) else data.table()

GENES <- c("AQP3","GGT5","ISLR","LAD1")
gcol <- c(AQP3="#2E5C8A", GGT5="#C0504D", ISLR="#4F9D69", LAD1="#8064A2")
exlab <- c(eQTLGen_blood="eQTLGen blood (n=31,684)", GTEx_heartLV="GTEx heart LV (n=386)")
main[, exposure_lab := factor(exlab[exposure], levels=unname(exlab))]
say("main rows: ", nrow(main), " ; per-snp rows: ", nrow(snp), " ; harmonised: ", nrow(harm))

save2 <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name,".pdf")), p, width=w, height=h, device="pdf")
  ggsave(file.path(FIG, paste0(name,".png")), p, width=w, height=h, dpi=600)
  say("wrote ", name)
}
th <- theme_bw(base_size=8) +
  theme(panel.grid.minor=element_blank(), strip.background=element_blank(),
        plot.title=element_text(size=8, face="bold"), legend.key.size=unit(8,"pt"),
        legend.position="bottom")

## ---------- A: gene level forest ----------
a_dt <- main[method %in% c("IVW","Wald ratio")]
setorder(a_dt, exposure, gene)
a_dt[, lab := factor(gene, levels=rev(GENES))]
pA <- ggplot(a_dt, aes(x=beta, y=lab, colour=gene)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), height=0.18) +
  geom_point(size=1.9) +
  geom_text(aes(x=beta+1.96*se, label=paste0("n=", nsnp)), hjust=-0.35, size=2.1,
            colour="grey35", show.legend=FALSE) +
  facet_wrap(~exposure_lab) +
  scale_colour_manual(values=gcol, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.12, 0.22))) +
  labs(x="log OR of heart failure per 1-SD genetically predicted expression", y=NULL,
       title="A  Gene-level causal estimates") + th

## ---------- scatter panels ----------
scatter_panel <- function(g, ex, letter) {
  d <- harm[exposure==ex & gene==g]
  if (!nrow(d)) return(ggplot() + theme_void())
  fits <- main[exposure==ex & gene==g]
  gb <- function(mm) { v <- fits[method==mm, beta]; if (length(v)==1) v else NA_real_ }
  b_ivw <- gb("IVW"); if (is.na(b_ivw)) b_ivw <- gb("Wald ratio")
  b_eg <- gb("MR-Egger"); b_wm <- gb("Weighted median")
  lns <- data.table(slope=numeric(), intercept=numeric(), method=character())
  if (is.finite(b_ivw)) lns <- rbind(lns, data.table(slope=b_ivw, intercept=0, method="IVW"))
  if (is.finite(b_eg))  lns <- rbind(lns, data.table(slope=b_eg,
                        intercept=fits[method=="MR-Egger", intercept], method="MR-Egger"))
  if (is.finite(b_wm))  lns <- rbind(lns, data.table(slope=b_wm, intercept=0, method="Weighted median"))
  if (nrow(lns)) {
    lns[, method := factor(method, levels=c("IVW","MR-Egger","Weighted median"))]
    xr <- range(d$bx)
    lns[, `:=`(x = xr[1], xend = xr[2],
               y = intercept + slope * xr[1], yend = intercept + slope * xr[2])]
  }
  p <- ggplot(d, aes(x=bx, y=by)) +
    geom_hline(yintercept=0, colour="grey85") + geom_vline(xintercept=0, colour="grey85") +
    geom_point(aes(size=1/sy^2), colour=gcol[g], alpha=.75) +
    scale_size_continuous(guide="none")
  if (nrow(lns)) p <- p + geom_segment(data=lns, aes(x=x, xend=xend, y=y, yend=yend, linetype=method))
  p <- p + scale_x_continuous(expand=expansion(mult=0.22)) +
    scale_y_continuous(expand=expansion(mult=0.22))
  p + labs(x=paste0("SNP effect on ", g, " expression (SD per allele)"),
           y="SNP effect on heart failure (log OR)", linetype=NULL,
           title=paste0(letter, "  ", g, ", ", ifelse(ex=="eQTLGen_blood","eQTLGen","GTEx"),
                        " (n=", nrow(d), ")")) +
    th + theme(legend.position="bottom")
}
pB <- scatter_panel("AQP3", "eQTLGen_blood", "B")
pC <- scatter_panel("GGT5", "GTEx_heartLV", "C")

## ---------- D: leave-one-out ----------
if (nrow(loo)) {
  dl <- loo[exposure=="eQTLGen_blood"]
  dl[, term := fifelse(dropped==0, "All SNPs", paste0("drop ", rsid))]
  lv <- c("All SNPs", dl[dropped!=0][order(beta), term])
  dl[, term := factor(term, levels=lv)]
  pD <- ggplot(dl, aes(x=beta, y=term)) +
    geom_vline(xintercept=0, linetype=2, colour="grey50") +
    geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), height=.15, colour="grey25") +
    geom_point(aes(colour = dropped==0), size=1.4) +
    scale_colour_manual(values=c(`TRUE`="#C0504D", `FALSE`="grey25"), guide="none") +
    facet_wrap(~gene, nrow=1) +
    labs(x="IVW log OR after dropping the SNP", y=NULL,
         title="D  Leave-one-out (eQTLGen)") + th
} else { pD <- ggplot() + theme_void() }

## ---------- E: funnel (AQP3) ----------
fun_dt <- snp[exposure=="eQTLGen_blood" & gene=="AQP3"]
pE <- ggplot(fun_dt, aes(x=beta, y=1/se)) +
  geom_point(colour=gcol["AQP3"], size=1.2) +
  geom_vline(xintercept=main[exposure=="eQTLGen_blood" & gene=="AQP3" & method=="IVW", beta],
             linetype=2, colour="grey30") +
  labs(x="Wald ratio (log OR per SD)", y="1 / SE of ratio", title="Funnel plot, AQP3 (eQTLGen)") + th

## ---------- F: per-SNP forest, AQP3 ----------
pf <- snp[exposure=="eQTLGen_blood" & gene=="AQP3"]
if (nrow(pf)) {
  pf[, lab := factor(rsid, levels=unique(rsid[order(beta)]))]
  pF <- ggplot(pf, aes(x=beta, y=lab)) +
    geom_vline(xintercept=0, linetype=2, colour="grey50") +
    geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), height=.2, colour=gcol["AQP3"]) +
    geom_point(size=1.2, colour=gcol["AQP3"]) +
    labs(x="Wald ratio (log OR per SD)", y=NULL, title="Per-SNP estimates, AQP3 (eQTLGen)") + th
} else pF <- ggplot() + theme_void()

## ---------- G: two-outcome comparison (HERMES vs FinnGen) ----------
f20 <- file.path(OUT,"01_tables","T20_two_outcome_multiplicity.csv")
if (file.exists(f20)) {
  t20 <- fread(f20)
  t20[, lab := paste0(gene, "  ", ifelse(exposure=="eQTLGen_blood","blood","heart"))]
  t20[, outcome := factor(outcome, levels=c("HERMES","FinnGen"))]
  lv <- t20[outcome=="FinnGen"][order(gene), unique(lab)]
  t20[, lab := factor(lab, levels=rev(lv))]
  pG <- ggplot(t20, aes(x=OR, y=lab, colour=outcome)) +
    geom_vline(xintercept=1, linetype=2, colour="grey50") +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=.16,
                   position=position_dodge(width=.55)) +
    geom_point(size=1.7, position=position_dodge(width=.55)) +
    scale_x_log10() +
    scale_colour_manual(values=c(HERMES="#4A6FA5", FinnGen="#C0504D"), name=NULL) +
    labs(x="OR of heart failure per 1-SD expression (log scale)", y=NULL,
         title="E  Two outcome GWAS") + th
} else pG <- ggplot() + theme_void()

## ---------- save ----------
save2(pA, "F01_gene_level_forest", 7.2, 2.6)
save2(pB, "F02_scatter_AQP3", 3.4, 3.6)
save2(pC, "F03_scatter_GGT5_GTEx", 3.4, 3.6)
save2(pD, "F04_leave_one_out", 7.2, 2.8)
save2(pE, "F05_funnel_AQP3", 3.4, 3.0)
save2(pF, "F06_per_snp_AQP3", 3.4, 3.4)
save2(pG, "F08_two_outcome_forest", 4.0, 3.2)

composite <- (pA) / (pB | pC) / (pD | pG) + plot_layout(heights=c(1.0, 1.45, 1.35))
save2(composite, "Figure8_mr_composite", 7.2, 8.2)
save2((pB | pF) / (pE | pC), "F07_panels_BCEF", 6.8, 6.6)

writeLines(capture.output(sessionInfo()), file.path(OUT,"03_logs","sessionInfo_114.txt"))
close(LOG)
cat("114 done\n")
