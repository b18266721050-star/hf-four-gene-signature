# A four-gene expression signature indexes myocardial remodelling rather than driving heart failure

Analysis code, package versions and source data tables for the manuscript:

> **A four-gene expression signature indexes myocardial remodelling rather than driving heart
> failure: a multi-cohort transcriptomic, single-nucleus and Mendelian randomisation study**

This repository contains everything needed to regenerate the tables and figures of the paper from
public data. It does **not** contain the raw expression matrices or the single-nucleus object
(≈37 GB); those are downloaded from GEO and the GWAS consortia by the scripts in `scripts/`.

---

## What is in here

| Path | Contents |
| --- | --- |
| `scripts/` | 110 R / Python scripts, numbered by analysis stage (see below) |
| `docs/package_versions.tsv` | Exact version of every R package used (R 4.6.1, 61 packages recorded) |
| `docs/CODE_AVAILABILITY.md` | Ready-to-paste data and code availability statement |
| `docs/zenodo_metadata.json` | Metadata template used for the archived release |
| `tables/` | Supplementary Tables S1–S16 (the numeric source of every claim in the paper) |
| `figures/` | Main Figures 1–8 as vector PDF |

## Data sources (all public, none redistributed here)

| Dataset | Source | Use |
| --- | --- | --- |
| GSE57338 | GEO, Affymetrix HuGene-1.0-ST (GPL11532) | discovery cohort, 313 left-ventricular samples |
| GSE5406 | GEO, Affymetrix HG-U133A (GPL96) | external validation, n = 210 |
| GSE116250 | GEO, Illumina RNA-seq | external validation, 64 deposited / 62 scored |
| GSE55296 | GEO, RNA-seq | external validation, n = 33 |
| GSE292067 | GEO, 10x Chromium 5' v2 | single-nucleus RNA-seq, 25 donors, 150,750 nuclei |
| HERMES (GCST009541) | GWAS Catalog | primary heart-failure outcome |
| FinnGen release 11 (`I9_HEARTFAIL`) | FinnGen | replication outcome |
| eQTLGen (n = 31,684) | eQTLGen consortium | primary cis-eQTL exposure |
| GTEx v8 left ventricle (n = 386) | GTEx | tissue-matched exposure |
| 1000 Genomes phase 3 EUR (n = 503) | 1000G | LD reference for fine-mapping and pairwise r² |

`scripts/01_download.R` retrieves the GEO series; the GWAS summary statistics, eQTL resources and
the 1000G VCF are downloaded manually (links and file sizes are listed in the manuscript's
Table S1 and in `docs/CODE_AVAILABILITY.md`).

## Requirements

- R 4.6.1 (developed on `x86_64-w64-mingw32`); package versions in `docs/package_versions.tsv`
- Python 3.11+ for the five `.py` helpers (standard library only, plus `pandas`/`numpy` where noted)
- No API keys are required. No network access is needed after the data are in place, except for
  `01_download.R` and the Ensembl minor-allele-frequency lookups in `112_…`/`125_…`.

## Running it

Two environment variables replace all hard-coded paths:

```r
Sys.setenv(GEOBULK_ROOT = "/path/to/hf_project")   # working directory
Sys.setenv(GEOBULK_LIB  = "/path/to/R/library")    # optional, only if your library is non-standard
```

Every script begins with `source("scripts/00_setup.R")`, which resolves `raw_data/`, `clean_data/`,
`results/` and `_tmp/` beneath `GEOBULK_ROOT` and sets the plotting theme and seeds. The expected
directory layout is created by `00_setup.R` on first run.

Run in numeric order; each stage writes into `results/<module>/01_tables` and `results/<module>/02_figures`
and logs its `sessionInfo`:

| Block | Stage |
| --- | --- |
| `00`–`08` | setup, download, cleaning, differential expression, GO/KEGG, Hallmark GSEA |
| `09`–`19` | WGCNA, three feature-selection algorithms, logistic models (M1–M3), immune infiltration (ssGSEA), STRING network |
| `20`–`69` | single-nucleus pipeline: QC, doublet detection, Harmony integration, annotation, four-gene localisation |
| `70`–`73` | batch and confounding diagnostics, ComBat sensitivity analysis |
| `80`–`83` | external validation, calibration, decision curves, sample-overlap audit |
| `90`–`96e` | xCell/EPIC deconvolution, composition adjustment, 70-gene axis replication |
| `97`–`101` | CellChat mass-action communication, donor-level statistics |
| `110`–`119` | Mendelian randomisation: instruments, harmonisation, IVW/Egger/weighted median, FinnGen replication, locus pleiotropy |
| `120`–`123` | Table 1, Figure 1 flow diagram, composite figure assembly |
| `124`–`135` | colocalisation, SuSiE fine-mapping with the 1000G panel, instrument-specificity and LD relations |

### Two components that are re-implemented rather than called from a package

Neither `TwoSampleMR` nor the CellChat C++ toolchain could be installed in the analysis
environment, so both were written out explicitly and are included, so that the estimators can be
checked line by line:

- `scripts/mr_core.R` — Wald ratio, inverse-variance weighting (fixed and multiplicative random
  effects), MR-Egger with intercept test and t-based inference, weighted median, Cochran's Q,
  leave-one-out, Steiger, following Bowden et al. (2017).
- `scripts/98_sc_cellchat_network.R` — the mass-action/Hill model (K = 0.5) applied to the official
  CellChatDB human resource, with the donor, not the nucleus, as the unit of replication.

## Reproducibility notes

- Random seeds are fixed inside each script: `20251029`, `20261027`, `123`.
- Figures are written as vector PDF (600 dpi PNG rasters are produced alongside in the submission
  package, not stored here).
- Two analyses are reported as sensitivity analyses and *not* as corrections: the ComBat run
  (`71_combat_reanalysis.R`) and the train-only feature selection (`09b_ml_trainonly.R`).
- The internal discrimination estimate is optimistically biased because feature selection preceded
  the training/test split; the manuscript says so and points to `Table S15`.

## Licence

Code: MIT (`LICENSE`). Tables and figures: CC-BY-4.0.

## Contact

Corresponding author details are withheld until the manuscript is accepted.
