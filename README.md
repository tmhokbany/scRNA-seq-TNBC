# Single-Cell Resolution Mapping of Cellular States Driving Antigen Presentation Suppression in Triple-Negative Breast Cancer

[![R 4.5.2](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.0-green.svg)](https://satijalab.org/seurat/)
[![scRNA-seq](https://img.shields.io/badge/scRNA--seq-10x_Genomics-orange.svg)](https://www.10xgenomics.com/)
[![GEO](https://img.shields.io/badge/Dataset-GSE161529-red.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529)

A single-cell transcriptomic atlas of the TNBC tumor microenvironment identifying cellular populations driving antigen presentation suppression, characterizing their functional states, and discovering candidate markers for therapeutic intervention.

---

## Problem & System Solution

Triple-negative breast cancer (TNBC) is an aggressive subtype with poor prognosis and limited therapeutic options. Immune evasion via suppression of antigen presentation (AP) represents a key mechanism by which TNBC tumors escape immune surveillance. However, the specific cell populations driving this suppression and their underlying molecular characteristics remain poorly understood.

The system transforms single-cell transcriptomic data into actionable biological insights:

1. **Single-Cell Transcriptomic Atlas**: Comprehensive UMAP and t-SNE visualization of 54,820 cells from TNBC tissue samples, revealing 10 transcriptionally distinct clusters encompassing malignant, stromal, and immune populations.

2. **Antigen Presentation Scoring**: 11-gene signature (HLA-A, HLA-B, HLA-C, HLA-DRA, HLA-DRB1, HLA-DPB1, TAP1, TAP2, PSMB8, PSMB9, B2M) enabling quantitative measurement of AP capacity at single-cell resolution.

3. **Cell State & Lineage Classification**: Systematic characterization of tumor cell states (proliferative, epithelial-like, basal-like, mesenchymal-like, stressed, cycling), immune cell types, and stromal populations.

4. **Dominant Suppressor Cluster Identification**: Unbiased ranking of clusters by AP suppression percentage, identifying Tumor_Epithelial cells (Cluster 0) as the primary pathological suppressor population (38.7% AP-suppressed).

5. **Candidate Marker Discovery**: Overlap analysis identifying FXYD3, PERP, KRT8, and KRT19 as high-confidence markers converging on ER stress and metabolic reprogramming pathways.

---

## System Architecture
```

                                     TNBC scRNA-seq Dataset (GSE161529)
                                     54,820 cells, 10x Genomics Platform
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │ Seurat QC & Data Normalization │ (Quality Control, Log-Normalization, Scaling) │
                └─────────────────────────────────────┬──────────────────────────────────────────┘
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │ Dimensionality Reduction & Clustering │ (PCA, UMAP, t-SNE, Louvain Clustering) │
                └─────────────────────────────────────┬──────────────────────────────────────────┘
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │ Gene Signature Scoring & AP Suppression │ (Cell State, Lineage & AP Scoring)   │
                │ - Cell State Score Calculation                                                 │
                │ - AP Score Calculation                                                         │
                │ - AP-Suppressed Cell Classification (25th %ile)                                │
                └─────────────────────────────────────┬──────────────────────────────────────────┘
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │ Cluster Characterization & Annotation │ (Lineage Assignment, Functional State  │
                │ - Lineage Assignment (Tumor/Immune/Stromal) │ Classification, Cluster Renaming)│
                │ - Functional State Classification                                              │
                └─────────────────────────────────────┬──────────────────────────────────────────┘
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │ Candidate Marker Discovery │ (Differential Expression, Correlation,            │
                │ - Differential Expression (Wilcoxon) │ Overlap Filtering)                      │
                │ - Spearman Correlation                                                         │ 
                │ - Overlap Filtering                                                            │
                └─────────────────────────────────────┬──────────────────────────────────────────┘
                                                      │
                                                      ▼
                ┌────────────────────────────────────────────────────────────────────────────────┐
                │     Visualization & Interpretation │ UMAPs, t-SNE, Feature Plots, Heatmaps     │
                └────────────────────────────────────────────────────────────────────────────────┘
```
---

## Key Findings

• **AP Suppression is Cell-Type Specific** : Pathological AP suppression is confined to tumor epithelial cells (38.7%), not immune cells (< 9%) or CAFs (lineage-intrinsic) 

• **Immune Exclusion is Driven by AP Suppression** : Strong negative correlation (r = -0.906) indicates AP suppression precedes and prevents T cell infiltration 

• **Distinct Tumor States Influence Suppression** : Epithelial-differentiated tumors show highest suppression; proliferative tumors show lowest 

• **ER Stress & Metabolic Reprogramming** : Candidate markers converge on UPR/ER stress pathways, linking metabolic adaptation to immune evasion 

• **Inter-Patient Heterogeneity** : AP scores and immune infiltration vary substantially across patients, informing stratification strategies 

---

## Start and Execution Guide

### 1. Prerequisites

- **R version 4.5.2** or higher
- RStudio (recommended) or R command-line interface

### 2. Installation

```bash
# Clone the repository
git clone https://github.com/tmhokbany/scRNA-seq-TNBC-analysis.git
cd scRNA-seq-TNBC-analysis
```

Install required R packages:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")
install.packages(c("ggplot2", "dplyr", "pheatmap", "viridis"))

# Or, for a reproducible install pinned to the exact versions this was
# developed with, use renv instead (see "Reproducible environment" below):
# install.packages("renv")
# renv::restore()
```

#### Reproducible environment (recommended)

This repo includes an [`renv`](https://rstudio.github.io/renv/) lockfile
(`renv.lock`) pinning the exact R and package versions the analysis was
developed with. Open `scRNA-seq-TNBC-analysis.Rproj` in RStudio (renv
activates automatically), then run:

```r
renv::restore()
```

This installs matching versions of Seurat, ggplot2, dplyr, pheatmap, and
viridis (plus their dependencies) into a project-local library, so results
don't silently drift as packages update over time.

### 3. Download Dataset

Download the pre-processed Seurat object from Figshare:

```bash
# Option 1: Using wget
wget -O SeuratObject_TNBC.rds https://figshare.com/ndownloader/files/36286674

# Option 2: Manual download from DOI: 10.6084/m9.figshare.17058077
```

### 4. Run the analysis

Open `TNBC_analysis.R` in RStudio (or run it with `Rscript TNBC_analysis.R`),
set your working directory to this folder as instructed at the top of the
script, and run it. Gene marker lists, cluster labels, and other constants
live in `config.R`, which is sourced automatically. Output plots and tables
are written to `analysis_output_tnbc/`.

---

### Academic References & Research Citations
Single-Cell Atlas of Human Breast: Pal B, Chen Y, Vaillant F, et al. (2021). A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. EMBO Journal, 40(11). DOI: 10.15252/embj.2020107333

Seurat Single-Cell Analysis: Hao Y, Hao S, Andersen-Nissen E, et al. (2021). Integrated analysis of multimodal single-cell data. Cell, 184(13):3573-3587.e29. DOI: 10.1016/j.cell.2021.04.048

MHC Class I Downregulation in Cancer: Cornel AM, Mimpen IL, Nierkens S. (2020). MHC class I downregulation in cancer: Underlying mechanisms and potential targets for cancer immunotherapy. Cancers, 12(7):1-33. DOI: 10.3390/cancers12071760

Cancer Immune Evasion through MHC Loss: Dhatchinamoorthy K, Colbert JD, Rock KL. (2021). Cancer immune evasion through loss of MHC class I antigen presentation. Frontiers in Immunology, 12:636568. DOI: 10.3389/fimmu.2021.636568

Tumor-Infiltrating Lymphocytes in Breast Cancer: Denkert C, von Minckwitz G, Darb-Esfahani S, et al. (2018). Tumour-infiltrating lymphocytes and prognosis in different subtypes of breast cancer: a pooled analysis of 3771 patients treated with neoadjuvant therapy. Lancet Oncology, 19(1):40-50. DOI: 10.1016/S1470-2045(17)30904-X

Antigen Processing Pathways: Blum JS, Wearsch PA, Cresswell P. (2013). Pathways of antigen processing. Annual Review of Immunology, 31:443-473. DOI: 10.1146/annurev-immunol-032712-095910

Cancer-Associated Fibroblasts: Kalluri R. (2016). The biology and function of fibroblasts in cancer. Nature Reviews Cancer, 16(9):582-598. DOI: 10.1038/nrc.2016.73

CAF Heterogeneity in Breast Cancer: Costa A, Kieffer Y, Scholer-Dahirel A, et al. (2018). Fibroblast heterogeneity and immunosuppressive environment in human breast cancer. Cancer Cell, 33(3):463-479.e10. DOI: 10.1016/j.ccell.2018.01.011

Cytokeratin 19 and ER Stress: Bambang IF, Lu D, Li H, et al. (2009). Cytokeratin 19 regulates endoplasmic reticulum stress and inhibits ERp29 expression via p38 MAPK/XBP-1 signaling in breast cancer cells. Experimental Cell Research, 315(11):1964-1974. DOI: 10.1016/j.yexcr.2009.02.017

PERP as p53 Target: Attardi LD, Reczek EE, Cosmas C, et al. (2000). PERP, an apoptosis-associated target of p53, is a novel member of the PMP-22/gas3 family. Genes & Development, 14(6):704-718. DOI: 10.1101/gad.14.6.704

# License
This project is licensed under the [MIT License](LICENSE). If you use this
code, please cite it — see [`CITATION.cff`](CITATION.cff).

Author: Tala Alhokbany

Supervisor: Dr. Hannah Al-Mubarak

Department: Biochemistry, College of Science, King Saud University

Date: February 2026
