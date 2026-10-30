# =============================================================================
# 53_sc_annotation_review.R
# Stage 3 - Step 3: build the annotation review table.
#
# Inputs : results/Singlecel/annotation/tables/cluster_qc.csv
#          results/Singlecel/annotation/tables/panel_pct_expr.csv
#          results/Singlecel/annotation/tables/doublet_like_index.csv
#          results/annotation/tables/{cluster_auto_module,module_mean_pct}.csv
# Outputs: results/annotation/annotation review.csv
#          results/annotation/tables/doublet_like_check.csv
#
# NOTE 1: GSE292067 deposited metadata carries sample / donor / condition only.
#         There are NO author-provided per-cluster cell-type labels, so the
#         `reference_annotation` column cites published human heart atlases and
#         says so explicitly instead of pretending to match author labels.
# NOTE 2: `codet_*` co-detection is evaluated against the product of the two
#         marginal detection rates. Co-detection that matches the product is
#         what independent (ambient) co-occurrence produces; only a clear
#         excess over the product argues for true doublets.
# =============================================================================

proj  <- "E:/another data file/GEObulk analysis project"
sdir  <- file.path(proj, "results", "Singlecel", "annotation", "tables")
adr   <- file.path(proj, "results", "annotation")
tdir  <- file.path(adr, "tables")

qc   <- read.csv(file.path(sdir, "cluster_qc.csv"), stringsAsFactors = FALSE)
pct  <- read.csv(file.path(sdir, "panel_pct_expr.csv"), stringsAsFactors = FALSE)
dl   <- read.csv(file.path(sdir, "doublet_like_index.csv"), stringsAsFactors = FALSE)
auto <- read.csv(file.path(tdir, "cluster_auto_module.csv"), stringsAsFactors = FALSE)
modp <- read.csv(file.path(tdir, "module_mean_pct.csv"), check.names = FALSE)

# per-cluster co-detection summary (68 tests vs the independence expectation)
codet <- NULL
f_codet <- file.path(tdir, "codet_recheck.csv")
if (file.exists(f_codet)) {
  codet <- read.csv(f_codet, stringsAsFactors = FALSE)
  cat("codet_recheck loaded:", nrow(codet), "tests; ratio range",
      round(min(codet$ratio, na.rm = TRUE), 2), "-",
      round(max(codet$ratio, na.rm = TRUE), 2), "\n")
}

rownames(pct) <- pct$gene
N <- sum(qc$n_cells)

# ------------------------------------------- doublet check: obs vs expected --
pairs <- list(
  COL1A1_CD3D   = c("COL1A1", "CD3D"),
  PECAM1_MYH11  = c("PECAM1", "MYH11"),
  DCN_CD3D      = c("DCN",    "CD3D"),
  TTN_PTPRC     = c("TTN",    "PTPRC"),
  LYZ_TNNT2     = c("LYZ",    "TNNT2"))

cl_cols <- paste0("C", dl$cluster)
chk <- do.call(rbind, lapply(names(pairs), function(p) {
  col <- paste0("codet_", p)
  if (!col %in% colnames(dl)) return(NULL)
  g <- pairs[[p]]
  a <- as.numeric(unlist(pct[g[1], cl_cols]))
  b <- as.numeric(unlist(pct[g[2], cl_cols]))
  obs <- as.numeric(unlist(dl[[col]]))
  data.frame(pair = p, lineage = paste(g, collapse = " x "),
             cluster = dl$cluster,
             observed_pct = obs,
             expected_pct = a * b * 100,
             stringsAsFactors = FALSE)
}))
chk <- chk[!is.na(chk$observed_pct) & !is.na(chk$expected_pct), ]
chk$ratio <- chk$observed_pct / chk$expected_pct
chk$verdict <- ifelse(chk$observed_pct >= 5 & chk$ratio >= 2,
                      "enriched beyond independence - doublet candidate",
                      "consistent with independent co-detection (ambient)")
chk <- chk[order(-chk$observed_pct), ]
write.csv(chk, file.path(tdir, "doublet_like_check.csv"), row.names = FALSE)
big <- chk[chk$observed_pct >= 1, ]
cat("doublet check (3 populated pairs, observed >= 1%): n =", nrow(big),
    " ratio range", round(min(big$ratio), 2), "-", round(max(big$ratio), 2),
    " | max observed", round(max(chk$observed_pct), 2), "%\n")

# ------------------------------------------------ curated review content ---
REF_ATLAS <- paste("human heart atlas: Litvinukova 2020 Nature;",
                   "Koenig 2022 Nat Cardiovasc Res; Kanemaru 2023 Nature.",
                   "GSE292067 metadata has NO author per-cluster labels (sample/donor/condition only).")

r <- list()

r[["0"]] <- c(
  label = "Cardiomyocyte",
  sub   = "atrial-like (MYH6/CHRM2/ALPK2 high); broad label kept",
  supp  = "TTN 1.00; MYH7 1.00; TNNT2 1.00; TNNI3 0.97; RYR2 1.00; MYBPC3 0.99; RBM20 0.99; MYL2 0.96; MYL3 0.98; CHRM2 0.94; ALPK2 0.91; MYH6 0.72",
  evid  = "CM module 0.89 = highest in dataset; 17 of the 20 CM panel genes exceed 0.91 (TTN/MYH7/RYR2 at 1.00), the only exceptions being NPPA and NPPB (0.27 each, secreted natriuretic peptides depleted in nuclei) and MYH6 (0.72); sarcomere + calcium-handling programme complete",
  contra= "MYOCD 0.88, ABCC9 0.81, PROX1 0.65 are co-detected but are also genuine cardiomyocyte genes, not foreign lineages; FB module 0.08 (DCN 0.31) and MYE 0.02 = ambient",
  conf  = "High",
  alt   = "none at major-type level; optional atrial vs ventricular split later",
  dbl   = "No - TTN x PTPRC 2.79% equals 2.79% expected (ratio 1.00) and all four cross-lineage pairs sit within 0.87-1.05 of independence; 46.4% of nuclei carry >=2 CM anchor genes, which is normal for a sarcomere-rich lineage",
  stat  = "accepted",
  act   = "none for major typing; revisit only if atrial/ventricular CM is a study aim",
  note  = "mean_nCount 12,740 = deepest cluster, so ambient/background detection is highest here. IMPORTANT: TTN is detected in 75-93% of nuclei in EVERY cluster (mast C14 93.1%, T/NK C10 84.6%, FB C1 81.4%), i.e. sarcomere transcripts dominate the ambient pool - TTN ranks CM correctly within the dataset but is NOT CM-specific and must not be used as a sole CM anchor"
)

r[["1"]] <- c(
  label = "Fibroblast",
  sub   = "quiescent / matrix-homeostatic (PDGFRA+ABCA9+, MFAP5 low)",
  supp  = "DCN 0.90; COL6A3 0.91; FBN1 0.87; ABCA9 0.85; COL1A2 0.75; COL3A1 0.74; PDGFRA 0.69; LUM 0.44; CFD 0.44; COL1A1 0.46",
  evid  = "FB module 0.59 with 2nd module CM only 0.23; PDGFRA + DCN + ABCA9 co-detected in the majority of nuclei",
  contra= "MFAP5 0.30 and SERPINE2 0.32 are well below the activated-FB level seen in C4 (0.66 / 0.64); LUM 0.44 lower than expected for fibroblast",
  conf  = "High",
  alt   = "Fibroblast (activated) - assigned to C4 instead",
  dbl   = "No - pct_ge2_anchor 77.9 is expected because COL1A1/COL1A2/DCN/FN1 are co-expressed in one lineage; cross-lineage pair DCN x CD3D only 0.05%",
  stat  = "accepted",
  act   = "none",
  note  = "the pct_ge2_anchor metric counts >=2 anchor genes within a lineage and must not be read as a doublet rate"
)

r[["2"]] <- c(
  label = "Endothelial cell",
  sub   = "capillary/arterial (BTNL9+CA4+, ACKR1-)",
  supp  = "VWF 0.92; FLT1 0.91; PECAM1 0.80; EMCN 0.77; ADGRL4 0.77; BTNL9 0.75; F8 0.68; CDH5 0.67; ERG 0.66; KDR 0.40; CA4 0.20",
  evid  = "EC module 0.57 = highest of the three EC clusters; BTNL9/CA4/ADGRL4 capillary-arterial programme; ERG (EC TF) 0.66",
  contra= "ACKR1 0.004 excludes a venous identity; WT1 0.19 and MESO module 0.05 are low but nonzero - mesothelial contamination not supported",
  conf  = "High",
  alt   = "Endothelial cell (venous) = C7; Endothelial cell (endocardial-like) = C8",
  dbl   = "No - PECAM1 x MYH11 1.35% vs 1.28% expected under independence (ratio 1.05)",
  stat  = "accepted",
  act   = "none",
  note  = "C2/C7/C8 share PECAM1/CDH5/VWF and are therefore listed as each other's alternatives, in line with the no-subtyping rule"
)

r[["3"]] <- c(
  label = "Pericyte",
  sub   = "mural / pericyte (RGS5+HIGD1B+EGFLAM+, MYH11-)",
  supp  = "RGS5 0.94; ABCC9 0.92; GUCY1A2 0.92; EGFLAM 0.88; PDGFRB 0.83; CARMN 0.82; NOTCH3 0.51; KCNJ8 0.44; HIGD1B 0.30",
  evid  = "PC module 0.67 = highest in dataset; HIGD1B and EGFLAM are pericyte-restricted and clearly above C11 (0.08 / 0.50); capillary-adjacent mural position on UMAP",
  contra= "MUC3A 0.66 was the only discordant gene and is resolved as a MUC3A-specific signal rather than a lineage conflict: Ensembl confirms ENSG00000169894 is canonical MUC3A (7:100,949,414-100,968,348, protein_coding); it is present in 25/25 samples (49.6-88.0% of C3 nuclei per sample), 45.5% of positive nuclei carry >=3 counts, no mucin paralogue co-detection excess (MUC12 1.37x, MUC6 0.94x, MUC16 0.50x, MUC17/MUC5B absent) and no epithelial programme (KRT8 0.007, EPCAM 0.000, CDH1 0.000). SMC programme still excluded (MYH11 0.02, MYOCD 0.13)",
  conf  = "High",
  alt   = "Smooth muscle cell (rejected: MYH11 0.02 / MYOCD 0.13); mural cell, unspecified",
  dbl   = "No - all cross-lineage pairs < 0.1%",
  stat  = "accepted",
  act   = "keep Pericyte; do not include MUC3A in the ambient/contamination score for C3, and do not use it as a marker downstream",
  note  = "MUC3A provenance checked in scripts 55/56: detection is not depth-proportional beyond a mild effect (rho 0.36 vs RGS5 0.46), not sample driven (25/25 samples), and spatially structured within C3 (UMAP binomial dispersion 47.3). MUC3A+ nuclei also show uniformly higher pericyte-marker detection (KCNJ8 +11.4, PDGFRB +10.5, NOTCH3 +10.2 points) with nCount +27% and unchanged DCN/KRT8, i.e. better-captured pericyte nuclei rather than a foreign lineage. LIMITATION: read-level verification (BAM / multi-mapping) is not possible from the deposited count matrix, so an unannotated transcript at the same locus cannot be formally excluded"
)

r[["4"]] <- c(
  label = "Fibroblast",
  sub   = "activated / matrix-remodelling (MFAP5+SERPINE2+FBLN1 high)",
  supp  = "DCN 0.97; FBN1 0.96; COL6A3 0.92; COL1A2 0.89; COL3A1 0.85; PDGFRA 0.68; MFAP5 0.66; SERPINE2 0.64; FBLN1 0.54; POSTN 0.28; HAPLN1 0.21",
  evid  = "FB module 0.71 = highest FB value; MFAP5 and SERPINE2 roughly double C1, which is the canonical activated vs quiescent contrast",
  contra= "POSTN only 0.28, so a full myofibroblast identity is not supported; no immune or endothelial co-expression",
  conf  = "High",
  alt   = "Fibroblast (quiescent) = C1",
  dbl   = "No - PECAM1 x MYH11 0.25%; DCN x CD3D 0.08%",
  stat  = "accepted",
  act   = "none; keep broad Fibroblast rather than splitting quiescent/activated in the final labels",
  note  = "C1 and C4 are both Fibroblast by design - the quiescent/activated distinction is recorded here as evidence only"
)

r[["5"]] <- c(
  label = "Macrophage",
  sub   = "resident / MRC1+CD163+LYVE1+ (M2-like)",
  supp  = "F13A1 0.91; MRC1 0.88; STAB1 0.82; CD163 0.81; PTPRC 0.79; MS4A4A 0.75; MS4A7 0.72; CSF1R 0.70; LILRB5 0.64; LYVE1 0.59",
  evid  = "MYE module 0.52 = highest; the resident-macrophage gene set (MRC1/CD163/F13A1/LYVE1/STAB1/MS4A4A/MS4A7/LILRB5) is internally coherent and PTPRC confirms haematopoietic origin",
  contra= "LYZ 0.13, CD68 0.25, CD14 0.28, ITGAX 0.12, FCN1 0.02 are low; this is the expected nuclear-transcript profile of resident macrophages (cytoplasmic myeloid transcripts are depleted in nuclei) rather than a lineage conflict",
  conf  = "High",
  alt   = "Monocyte/DC = C9; Macrophage (proliferating) = C15",
  dbl   = "No - PECAM1 x MYH11 0.83%; no CM/TNK co-expression",
  stat  = "accepted",
  act   = "none",
  note  = "the LYZ-low / MRC1-high pattern is a nuclei-vs-cells effect; do not use whole-cell myeloid gates on this object. TNNT2 is ambient at 34.1% here, so LYZ x TNNT2 co-detection is 4.18% vs 4.36% expected and carries no doublet signal"
)

r[["6"]] <- c(
  label = "Cardiomyocyte",
  sub   = "ventricular (MYBPC3/RBM20/MYL2 high, MYH6 lower than C0)",
  supp  = "TTN 1.00; MYH7 0.98; RYR2 0.98; MYBPC3 0.94; RBM20 0.90; MYL2 0.89; TNNI3 0.88; DES 0.87; MYL3 0.85; CHRM2 0.70; ALPK2 0.60; MYH6 0.43",
  evid  = "CM module 0.77; complete sarcomere set plus the ventricular MLC programme (MYL2/MYL3) and the ventricular splicing factor RBM20",
  contra= "MYH6 0.43 and CHRM2 0.70 still present, so the atrial programme is only partially absent; NPR3 0.20 and POSTN 0.08 are ambient. mean_nCount 5,356 vs 12,740 in C0 - the depth gap alone lowers every detection rate, which must not be read as a different lineage",
  conf  = "High",
  alt   = "Cardiomyocyte (atrial-like) = C0",
  dbl   = "No - PECAM1 x MYH11 0.11%",
  stat  = "accepted",
  act   = "none; C0 and C6 both stay Cardiomyocyte",
  note  = "C0 vs C6 differ in MYH6/CHRM2/ALPK2 vs MYBPC3/RBM20/MYL2 and in sequencing depth; a later atrial/ventricular split would need depth-matched comparison"
)

r[["7"]] <- c(
  label = "Endothelial cell",
  sub   = "venous (ACKR1+NR2F2+, BTNL9/CA4 low)",
  supp  = "VWF 0.93; PECAM1 0.89; FLT1 0.85; EMCN 0.82; ADGRL4 0.76; ERG 0.71; CDH5 0.50; F8 0.37; NR2F2 0.28; ACKR1 0.23",
  evid  = "EC module 0.50; ACKR1 0.23 is the highest ACKR1 value in the dataset (all other clusters <= 0.02), which is the defining venous-EC transcript",
  contra= "CDH5 0.50 below C2 (0.67); BTNL9 0.10 and CA4 0.04 rule out capillary/arterial identity",
  conf  = "High",
  alt   = "Endothelial cell (capillary/arterial) = C2",
  dbl   = "No - PECAM1 x MYH11 1.58% vs 1.48% expected (ratio 1.07)",
  stat  = "accepted",
  act   = "none",
  note  = "ACKR1 is the only clean discriminator between C7 and C2, and it is a single gene - reason enough to keep both under the broad Endothelial cell label"
)

r[["8"]] <- c(
  label = "Endothelial cell",
  sub   = "endocardial (NRG1/SMOC1/SALL1/GATA5/NPR3) - subtype resolved in script 62",
  supp  = "PECAM1 0.90; EMCN 0.93; VWF 0.88; CDH5 0.58; ERG 0.78; POSTN 0.87; SMOC1 0.78; COL3A1 0.68; NPR3 0.62; NRG1 0.46; NFATC1 0.39; SALL1 0.17; GATA5 0.13",
  evid  = "98.9% of C8 nuclei are positive for >=2 of PECAM1/CDH5/VWF/EMCN/ERG (C2 96.2%, C7 97.3%), i.e. C8 is as endothelial as the other EC clusters. Against C2+C7 the top up genes are NRG1 +8.90 (46.0% vs 1.9%), SMOC1 +8.65 (77.6% vs 0.5%), SALL1 +8.23, SALL3 +8.09, GATA5 +7.78, NPR3 +5.92, POSTN +5.56 - a coherent endocardial programme in which SALL1/GATA5/NFATC1 are endocardial transcription factors, not a single-gene artefact",
  contra= "blood-EC markers are absent (BTNL9 0.010, CA4 0.003, ACKR1 0.017, FLT1 0.78 vs 0.91 in C2) and lymphatic markers are absent (FLT4 0.006, PROX1 0.061, CCL21 0.029), so neither blood EC nor LEC; HAPLN1 0.024 is low but is not required for an endocardial call; TTN 0.75 / MYH7 0.47 reflect cardiomyocyte ambient (endocardium borders myocardium)",
  conf  = "High",
  alt   = "Endothelial cell (capillary/arterial) = C2; Endothelial cell (venous) = C7; Doublet-like (EC x CM) rejected (see doublet_assessment)",
  dbl   = "No - PECAM1 x MYH11 21.7% looks alarming but equals 0.905 x 0.219 = 19.8% expected under independence (ratio 1.10); it is ambient MYH11, not doublet enrichment. The CM module itself is 0.19, the lowest of all clusters, which is incompatible with a genuine CM-containing doublet population",
  stat  = "accepted",
  act   = "none; keep the broad Endothelial cell label and record the subtype in celltype_note",
  note  = "originally withheld as Uncertain at the user request; resolved with C8 vs C2+C7 differential expression plus cell-level co-expression (scripts 62/63). This cluster also exemplifies why co-detection must be compared with the product of marginals - a raw 21.7% would have produced a false Doublet-like call"
)

r[["9"]] <- c(
  label = "Monocyte/DC",
  sub   = "non-resident myeloid (LYZ/FCN1/ITGAX+), monocyte vs DC unresolved",
  supp  = "PTPRC 0.79; KYNU 0.62; CSF1R 0.51; ITGAX 0.50; CD163 0.44; MRC1 0.37; LYZ 0.36; FCN1 0.28; CD68 0.16; CD14 0.16",
  evid  = "MYE module 0.35; LYZ + FCN1 + ITGAX is the classical non-resident myeloid combination and is essentially absent from C5 (0.13 / 0.02 / 0.12), which separates the two myeloid clusters",
  contra= "no single coherent identity: MRC1 0.37 and CD163 0.44 sit between the resident macrophage level (0.88 / 0.81) and background; the panel lacks CD1C/CLEC9A/LILRA4, so monocyte vs cDC vs pDC cannot be decided here",
  conf  = "Moderate",
  alt   = "Macrophage = C5; Monocyte; Dendritic cell",
  dbl   = "No - PECAM1 x MYH11 0.81%; DCN x CD3D 0.18%",
  stat  = "accepted as broad myeloid class",
  act   = "keep Monocyte/DC; add CD1C/CLEC9A/LILRA4/FCGR3A to the panel before any monocyte-DC split",
  note  = "cluster markers from FindMarkers are KYNU/TBXAS1/DOCK2, all shared with C5 or C15 - consistent with a mixed myeloid identity. TNNT2 ambient is 31.2% here, giving LYZ x TNNT2 11.6% vs 11.4% expected (ratio 1.02), i.e. ambient not doublets"
)

r[["10"]] <- c(
  label = "T/NK cell",
  sub   = "T and NK lymphocytes combined, not separated",
  supp  = "PTPRC 0.93; SKAP1 0.86; CD247 0.66; ITK 0.64; BCL11B 0.60; THEMIS 0.56; IL7R 0.40; CD2 0.36; KLRD1 0.28; GNLY 0.25; CD3E 0.24; NKG7 0.18; CD8A 0.13; CD3D 0.10",
  evid  = "TNK module 0.38 with a coordinated T-cell TF set (THEMIS, BCL11B, ITK, SKAP1) plus NK effectors (NKG7, GNLY, KLRD1) and near-absent myeloid/stromal modules",
  contra= "CD3D 0.10 and CD3E 0.24 are low in absolute terms; CD3 transcripts are partly cytoplasmic and under-detected in nuclei, and the coordinated TF panel is what carries the call. BPL module 0.005 excludes B cells",
  conf  = "High",
  alt   = "T cell; NK cell (both subsets of this cluster)",
  dbl   = "No - DCN x CD3D 2.67% is the highest value of that pair in the dataset but only 1.35x the independence expectation (1.98%)",
  stat  = "accepted",
  act   = "none; do not split T vs NK at this stage",
  note  = "pct.1 in the FindMarkers table is computed within the cluster and was 0.095 for CD3D - the panel detection rate, not pct.1, is the correct readout here"
)

r[["11"]] <- c(
  label = "Smooth muscle cell",
  sub   = "vascular SMC (MYH11/MYOCD/LMOD1+), possible pericyte admixture",
  supp  = "CARMN 0.91; SLIT3 0.85; MYH11 0.73; ACTA2 0.69; NTRK3 0.65; LMOD1 0.60; MYOCD 0.55; TAGLN 0.44; CNN1 0.11",
  evid  = "SMC module 0.59 with the contractile TF MYOCD and the SMC-restricted markers MYH11/LMOD1/NTRK3; clearly separated from C3 by MYH11 (0.73 vs 0.02) and MYOCD (0.55 vs 0.13)",
  contra= "PC module 0.39 - RGS5 0.62, PDGFRB 0.76, EGFLAM 0.50, GUCY1A2 0.59, NOTCH3 0.56 are co-detected. These are shared mural-cell genes rather than foreign lineage genes, and the pericyte-restricted markers stay low (HIGD1B 0.08, KCNJ8 0.08), so this is not a pericyte cluster",
  conf  = "High",
  alt   = "Pericyte = C3; mural cell, unspecified",
  dbl   = "No - PECAM1 x MYH11 2.84% vs 3.26% expected (ratio 0.87), i.e. below independence",
  stat  = "accepted",
  act   = "none; if mural subtypes matter later, sub-cluster C3+C11 jointly rather than relabelling C11",
  note  = "the C3 vs C11 distinction rests on MYH11/MYOCD/LMOD1 vs HIGD1B/KCNJ8/EGFLAM; RGS5, PDGFRB, GUCY1A2 and NOTCH3 are shared and must not be used alone"
)

r[["12"]] <- c(
  label = "Neural/glial cell",
  sub   = "peripheral glia / Schwann-like (SOX10/PLP1/CDH19+)",
  supp  = "NRXN1 0.99; XKR4 0.95; GPM6B 0.78; PTPRZ1 0.74; SOX10 0.41; PLP1 0.38; S100B 0.38; EGFLAM 0.30; MPZ 0.22",
  evid  = "NEUR module 0.65; SOX10 plus XKR4/PTPRZ1/GPM6B/NRXN1 is a coherent peripheral-glia programme and no competing lineage module exceeds 0.26",
  contra= "SOX10 0.41 and PLP1 0.38 are moderate rather than high - myelin and cytoplasmic glial transcripts are partly lost in nuclei; MPZ 0.22 is low. CDH19 was dropped from the evidence after the final check: it is 0.96 in C12 but also 0.82 in C1 and 0.90 in C4, so it is NOT glial-specific in this dataset and must not be cited as a Schwann marker",
  conf  = "High",
  alt   = "Schwann cell (more specific, same cluster); Glial cell",
  dbl   = "No - all cross-lineage pairs <= 0.3%",
  stat  = "accepted",
  act   = "none; keep the broad Neural/glial label, Schwann cell is recorded only as a note",
  note  = "cardiac nerves are a normal constituent of ventricular tissue; no evidence of a neuronal (as opposed to glial) population in this panel"
)

r[["13"]] <- c(
  label = "Lymphatic endothelial cell",
  sub   = "LEC (PROX1/CCL21/FLT4/MMRN1+)",
  supp  = "PKHD1L1 0.82; FLT4 0.81; MMRN1 0.80; PROX1 0.70; CCL21 0.69; PECAM1 0.67; CDH5 0.55; STAB2 0.51; KDR 0.49; LYVE1 0.35; TFF3 0.30",
  evid  = "LEC module 0.62 and the lymphatic-restricted set PROX1/CCL21/FLT4/MMRN1/STAB2/TFF3 is co-detected, with a 2nd module (EC) at only 0.28",
  contra= "none at lineage level; ACKR1 0.002 excludes venous EC. PKHD1L1 is also 0.97 in C8, so this gene alone does not define lymphatics",
  conf  = "High",
  alt   = "Endothelial cell (blood) = C2/C7",
  dbl   = "No - PECAM1 x MYH11 0.77%",
  stat  = "accepted",
  act   = "none",
  note  = "PROX1 0.37 in C6 and 0.65 in C0 shows PROX1 is not LEC-exclusive in heart nuclei; the four-gene combination is what supports the call"
)

r[["14"]] <- c(
  label = "Mast cell",
  sub   = "mast cell (KIT/CPA3/TPSB2/MS4A2/HDC+)",
  supp  = "KIT 0.96; CPA3 0.84; HDC 0.77; TPSAB1 0.59; MS4A2 0.56; TPSB2 0.54; PTPRC 0.56",
  evid  = "MAST module 0.62 with KIT plus three granule genes (CPA3, TPSAB1/TPSB2, MS4A2) and the histamine synthesis enzyme HDC - a complete mast-cell programme",
  contra= "CMA1 only 0.05 (MCTC vs MCT phenotype unresolved); SKAP1 0.70 and ERG 0.43 are haematopoietic/endothelial genes that are also expressed in mast cells and do not indicate another lineage",
  conf  = "High",
  alt   = "none",
  dbl   = "No cross-lineage co-detection, but mean_dbl_score 0.068 is the second highest of all clusters (after C16) - worth a local scDblFinder re-check if C14 is analysed on its own",
  stat  = "accepted",
  act   = "none",
  note  = "911 nuclei across 25 samples, so the label is well replicated despite the small size"
)

r[["15"]] <- c(
  label = "Proliferating cell",
  sub   = "cycling, myeloid/macrophage-biased; not a lineage",
  supp  = "MKI67 0.76; BUB1B 0.65; RRM2 0.64; CENPF 0.64; TOP2A 0.62; PTPRC 0.80; CD163 0.80; F13A1 0.80; MRC1 0.77; STAB1 0.74; MS4A4A 0.72; LYVE1 0.70; KYNU 0.63; CD14 0.39; CD68 0.34",
  evid  = "PROL module 0.66 (five independent cell-cycle genes) with MYE 0.51 as the only substantial lineage module; the myeloid genes are the same resident-macrophage set seen in C5",
  contra= "no lineage identity beyond proliferation: FB 0.13, CM 0.22, TNK 0.02, EC 0.09; pct_disease 23.1% is by far the lowest of all clusters, which is the opposite of what a disease-driven expansion would show and suggests a dissociation/repair-associated state",
  conf  = "Moderate",
  alt   = "Proliferating macrophage; Cycling cell; merge into Macrophage (C5) only after sub-clustering",
  dbl   = "No - PECAM1 x MYH11 1.89%; LYZ x TNNT2 12.28% vs 11.28% expected (ratio 1.09), explained by 34.4% ambient TNNT2 rather than by doublets",
  stat  = "accepted as a broad, non-lineage class",
  act   = "keep Proliferating cell as its own label; do not merge into C5 without a cell-cycle-regressed comparison",
  note  = "847 nuclei; because a cycling class is defined by state rather than lineage, it should be reported separately and not counted as a cell type in composition analyses without a footnote"
)

r[["16"]] <- c(
  label = "Adipocyte",
  sub   = "adipocyte, likely epicardial adipose; small and sample-restricted",
  supp  = "PLIN1 0.94; GPAM 0.94; GPD1 0.86; ADIPOQ 0.86; CIDEA 0.61; LEP 0.27",
  evid  = "ADIPO module 0.75 with PLIN1/ADIPOQ/GPD1/GPAM/CIDEA - a complete adipocyte lipid programme, and 2nd module (FB) only 0.26",
  contra= "CARMN 0.84 and SLIT3 0.81 are unexpectedly high for adipocytes; COL3A1 0.64 / FBN1 0.65 indicate ECM ambient. Only 257 nuclei from 18 of 25 samples, so the cluster may reflect epicardial fat captured in a subset of biopsies rather than ventricular adipocytes",
  conf  = "Moderate",
  alt   = "Adipocyte (epicardial adipose); Fibroblast/adipocyte doublet (rejected: codet pairs 0.39%)",
  dbl   = "No - PECAM1 x MYH11 0.39%",
  stat  = "accepted with caveat",
  act   = "keep the label but report n and sample coverage; do not build composition statistics on 257 nuclei without stating the 18/25 sample restriction",
  note  = "mean_nCount 8,349 and mean_pct_mt 0.11 - high depth, low mitochondrial content, consistent with genuine adipocyte nuclei rather than debris"
)

# ------------------------------------------------------------- assemble -----
stopifnot(setequal(names(r), as.character(qc$cluster)))

rev <- do.call(rbind, lapply(as.character(qc$cluster), function(k) {
  x <- r[[k]]
  data.frame(
    cluster            = k,
    resolution_column  = "harmony_snn_res.0.6",
    n_cells            = qc$n_cells[qc$cluster == k],
    pct_of_all_nuclei  = round(qc$n_cells[qc$cluster == k] / N * 100, 2),
    pct_disease        = qc$pct_disease[qc$cluster == k],
    n_samples          = qc$n_samples[qc$cluster == k],
    mean_nCount        = qc$mean_nCount[qc$cluster == k],
    candidate_label    = x[["label"]],
    subtype_note       = x[["sub"]],
    supporting_markers = x[["supp"]],
    supporting_evidence= x[["evid"]],
    contradictory_markers = x[["contra"]],
    auto_annotation    = auto$auto_module_rank[auto$cluster == k],
    reference_annotation = REF_ATLAS,
    confidence         = x[["conf"]],
    alternative_labels = x[["alt"]],
    doublet_assessment = x[["dbl"]],
    review_status      = x[["stat"]],
    recommended_action = x[["act"]],
    notes              = x[["note"]],
    stringsAsFactors   = FALSE, check.names = FALSE)
}))

rev <- rev[order(as.integer(rev$cluster)), ]

# ------------------------------------- attach co-detection summary column ----
if (!is.null(codet)) {
  agg <- do.call(rbind, lapply(unique(codet$cluster), function(k) {
    x  <- codet[codet$cluster == k, ]
    x1 <- x[x$observed_pct >= 1, ]          # ratios below 1% rest on a handful of nuclei
    if (!nrow(x1)) x1 <- x
    data.frame(cluster = k,
               max_codet_pct = max(x$observed_pct, na.rm = TRUE),
               max_codet_ratio = max(x1$ratio, na.rm = TRUE),
               codet_enrichment = ifelse(any(x1$observed_pct >= 5 & x1$ratio >= 2),
                                         "yes", "none"), stringsAsFactors = FALSE)
  }))
  rev <- merge(rev, agg, by = "cluster", all.x = TRUE, sort = FALSE)
  rev <- rev[order(as.integer(rev$cluster)), ]
  rev$codet_note <- sprintf(
    "max cross-lineage co-detection %.2f%% (obs/exp %.2f) across 4 pairs; %s",
    rev$max_codet_pct, rev$max_codet_ratio,
    ifelse(rev$codet_enrichment == "yes", "ENRICHED - doublet candidate",
           "no enrichment beyond independence"))
  rev$codet_enrichment <- NULL
}

# ------------------------------------------- user decisions (2026-10-29) ----
f_dec <- file.path(tdir, "annotation_decisions.csv")
if (file.exists(f_dec)) {
  dec <- read.csv(f_dec, stringsAsFactors = FALSE, check.names = FALSE)
  names(dec)[1] <- "cluster"          # guard against a UTF-8 BOM in the header
  dec <- dec[, c("cluster", "decision", "decided_on", "note")]
  names(dec)[names(dec) == "note"] <- "user_decision_note"
  rev <- merge(rev, dec, by = "cluster", all.x = TRUE, sort = FALSE)
  rev <- rev[order(as.integer(rev$cluster)), ]
  rev$in_final_label_set <- ifelse(rev$decision == "exclude", FALSE, TRUE)
  rev$review_status[rev$decision == "exclude"] <-
    "EXCLUDED by user decision - drop cluster from downstream analysis"
  rev$review_status[rev$decision == "pending"] <-
    "label pending user decision - must be raised again at conclusion stage"
  cat("decisions applied:", paste(names(table(rev$decision)), table(rev$decision),
                                 sep = "=", collapse = " "), "\n")
  cat("nuclei in final label set:", sum(rev$n_cells[rev$in_final_label_set]),
      "of", sum(rev$n_cells), "\n")
}

# move the co-detection note next to doublet_assessment for readability
ord_cols <- c("cluster","resolution_column","n_cells","pct_of_all_nuclei","pct_disease",
              "n_samples","mean_nCount","candidate_label","subtype_note",
              "supporting_markers","supporting_evidence","contradictory_markers",
              "auto_annotation","reference_annotation","confidence","alternative_labels",
              "doublet_assessment","codet_note","review_status","recommended_action","notes",
              "decision","in_final_label_set","user_decision_note","decided_on")
if (all(ord_cols %in% colnames(rev))) rev <- rev[, ord_cols]

out <- file.path(adr, "annotation review.csv")
write.csv(rev, out, row.names = FALSE, fileEncoding = "UTF-8")

cat("written:", out, "\n")

# ------------------------------------- final label set (after user decisions) --
if ("in_final_label_set" %in% colnames(rev)) {
  fin <- data.frame(cluster = rev$cluster,
                    candidate_label = rev$candidate_label,
                    include_in_analysis = rev$in_final_label_set,
                    n_cells = rev$n_cells,
                    confidence = rev$confidence,
                    status = rev$review_status,
                    stringsAsFactors = FALSE)
  write.csv(fin, file.path(tdir, "cluster_labels_final.csv"), row.names = FALSE)
  cat("final label set written; included clusters:",
      paste(fin$cluster[fin$include_in_analysis], collapse = ","), "\n")
}
cat("rows:", nrow(rev), " cols:", ncol(rev), "\n")
print(table(rev$candidate_label))
print(table(rev$confidence))
print(table(rev$review_status))
cat("\nsum pct =", sum(rev$pct_of_all_nuclei), " sum n =", sum(rev$n_cells), "\n")
writeLines(capture.output(sessionInfo()), file.path(adr, "sessionInfo_53.txt"))
