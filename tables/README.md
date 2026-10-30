# Supplementary tables

Numeric source of every claim in the manuscript. Files are named `TableS<n>_<k>.csv` when a single
table was split across several sheets; `.xlsx` files keep the original multi-sheet layout.

| Table | Content | Produced by |
| --- | --- | --- |
| S1 | Sample metadata and group assignment for the 313 discovery samples | `70_batch_diagnosis.R` |
| S2 | Full differential-expression table (18,726 genes) | `04_deg.R` |
| S3 | GO / KEGG over-representation and Hallmark GSEA results | `07_enrichment.R`, `08_gsea.R` |
| S4 | WGCNA module–trait relationships and the 70-gene module | `09_wgcna.R` |
| S5 | Genes retained by LASSO / SVM-RFE / random forest, with coefficients and importances | `09_ml_feature_selection.R` |
| S6 | Regression coefficients, VIF, events per variable, bootstrap optimism | `10_nomogram.R` |
| S7 | Immune infiltration: 28 cell types, pairwise contrasts, gene–cell correlations | `13_immune_stats.R` |
| S8 | Batch diagnostics and ComBat before/after comparison | `70`–`73` |
| S9 | External validation: AUC, calibration, DCA, direction of effect in three cohorts | `81_ext_validation.R` |
| S10 | Sample-overlap audit between GSE57338 and GSE5406 | `83_sample_overlap_audit.R` |
| S11 | xCell/EPIC deconvolution and the composition-adjusted differential expression | `90`–`92` |
| S12 | CellChat pathways, information flow and edge tables | `99`–`101` |
| S13 | Full Mendelian randomisation output (HERMES and FinnGen) | `113`–`118` |
| S14 | STRING node degrees for the 68 mapped genes | `17_ppi_network.R` |
| S15 | Feature selection restricted to the 219 training samples (sensitivity analysis) | `09b_ml_trainonly.R` |
| S16 | 9p13.3 colocalisation, SuSiE fine-mapping and instrument specificity (12 sheets) | `124`–`135` |

`TableS15_*` is an independent re-implementation used only to judge robustness: re-running the
original three algorithms on the full cohort with that implementation returns the same four-gene
intersection, but the per-algorithm counts differ (23/21/29 versus 20/43/29) because the
cross-validation folds and the "optimal number of features" rule differ.
