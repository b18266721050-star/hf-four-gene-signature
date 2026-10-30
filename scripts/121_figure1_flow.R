## 121_figure1_flow.R -- Figure 1: study design flow chart (ASCII only: R runs in a C locale here)
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
use("data.table", "ggplot2")

MDIR <- new_dir(DIR$results, "Manuscript")
FIG  <- new_dir(DIR$results, "Manuscript", "figures")

boxes <- data.table(
  y = c(95, 81.5, 65.5, 50.5, 36.5, 22.5, 8.5),
  h = c(9.5, 9.5, 12.5, 9.5, 10.5, 9.5, 10.5),
  w = c(78, 78, 88, 62, 78, 70, 84),
  fill = c("#E9EEF4", "#E9EEF4", "#FBF0E2", "#FBF0E2", "#E8F1EA", "#E8F1EA", "#F3EAF1"),
  tag = c("1", "2", "3", "4", "5", "6", "7"),
  label = c(
    "Discovery cohort\nGSE57338 left ventricular myocardium; 313 samples\n(NF 136 / DCM 82 / ICM 95)",
    "Differential expression (limma)\n269 genes: 135 up-regulated / 134 down-regulated\n(|log2FC| >= 0.585, FDR < 0.05)",
    "Two candidate gene sets\nWGCNA HF-associated module intersect DEG = 70 genes (used for the composition-adjusted re-analysis)\nSubtype-consistent DEGs (DCM intersect ICM) = 212 genes (input to feature selection)",
    "Feature selection: LASSO / SVM-RFE / random forest\nintersection = 4 genes: AQP3, GGT5, ISLR, LAD1",
    "Logistic models M1 (clinical) / M2 (4 genes) / M3 (combined)\n7:3 stratified split; 1,000-bootstrap optimism; calibration; decision-curve analysis",
    "External validation of the frozen model\nGSE5406 AUC 0.92; GSE116250 0.93; GSE55296 0.84",
    "Mechanism and causality\nsnRNA-seq (25 donors, 150,750 nuclei); cell-cell communication\ntwo-sample MR against HERMES and FinnGen")
)
boxes[, x := 50]

arrows <- data.table(
  x = 50, xend = 50,
  y = boxes$y[-nrow(boxes)] - boxes$h[-nrow(boxes)] / 2,
  yend = boxes$y[-1] + boxes$h[-1] / 2 + 0.6
)

p <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey30", linewidth = 0.32) +
  scale_fill_identity() +
  geom_text(data = boxes, aes(x = x - w/2 + 4.2, y = y, label = label),
            size = 2.45, lineheight = 1.05, colour = "grey8", hjust = 0) +
  geom_text(data = boxes, aes(x = x - w/2 + 2.0, y = y + h/2 - 2.9, label = tag),
            size = 2.45, fontface = "bold", colour = "grey45") +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = grid::unit(1.5, "mm"), type = "closed"),
               colour = "grey35", linewidth = 0.35) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(3, 4, 3, 4))

save_pub(p, "Figure1_study_design", dir = FIG, width = 7.1, height = 6.4)
cat("wrote Figure1_study_design\n")

leg <- c(
"Figure 1. Study design.",
"Bulk transcriptomes from left-ventricular myocardium were used to discover differentially expressed genes (GSE57338;",
"313 samples: 136 non-failing donor hearts, 82 with dilated cardiomyopathy and 95 with ischaemic cardiomyopathy) and to",
"define two candidate gene sets: the intersection of the heart-failure-associated WGCNA module with the differentially",
"expressed genes (70 genes), and the set of genes differentially expressed in both cardiomyopathy subtypes (212 genes).",
"Three feature-selection algorithms were applied to the 212 genes and their intersection (AQP3, GGT5, ISLR and LAD1)",
"was carried forward into logistic models M1 (clinical), M2 (four genes) and M3 (combined). The frozen M2 model was",
"validated in three independent cohorts (GSE5406, GSE116250 and GSE55296). Mechanisms were addressed with single-nucleus",
"RNA-sequencing (GSE292067; 25 donors, 150,750 nuclei), a ligand-receptor communication analysis, and two-sample",
"Mendelian randomisation against two heart-failure GWAS (HERMES and FinnGen release 11).",
"AUC, area under the receiver-operating characteristic curve; DEG, differentially expressed gene; MR, Mendelian",
"randomisation; NF, non-failing; WGCNA, weighted gene co-expression network analysis.")
writeLines(leg, file.path(FIG, "Figure1_legend.txt"))
cat("wrote Figure1_legend.txt\n")
write_session(file.path(MDIR, "logs"), "121_figure1")
