# Metabolomic Discovery of Multi-fluid Biomarkers of Food Intake in a Randomized Acute Feeding Trial


## 📌 Overview

This repository contains the simulated sample data and main analysis code used for investigation of candidate biomarkers of 8 test foods and estimation of their kinetic properties. The analysis utilizes diet, metabolomics data from food, plasma, and urine samples.


## 📂 Repository Structure
```text
├── data/
│   ├── plasma_sample.RData
│   ├── urine_sample.RData
│
├── code/
│   ├── Figure 2.R
│   ├── Figure 3.R
│   ├── Figure 4.R
│   ├── Figure 5.R
│   └── Figure 6.R
│
└── README.md

Note: Individual-level data are not included in this repository because of data-use and participant confidentiality restrictions.
```

## 🔬 Sample Data
- plasma_sample.Rdata: A simulated plasma dataset with the same sample size and metabolite variables as the original dataset.

- urine_sample.RData: A simulated urine dataset with the same sample size and metabolite variables as the original dataset.


## 🧬 Main Analysis Code
- Figure 2.R: Time-course analysis based on sequential postprandial metabolomics data.

- Figure 3.R: Identification of plasma metabolites in response to food intake using linear mixed model.

- Figure 4.R: Identification of urinary metabolites in response to food intake using linear mixed model.

- Figure 5.R: Classification performance of multi-metabolite signature for differentiating test foods. LASSO model was used to select metabolites. 

- Figure 6.R: Estimation of kinetic parameters of annotated metabolites uisng non-compartmental model.


## 💻 Environment & Dependencies
- R version 4.3.3 (2024-02-29)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Rocky Linux 9.8 (Blue Onyx)

- Matrix products: default
- BLAS:   /app/R-4.3.3@i86-rhel9.0/lib64/R/lib/libRblas.so 
- LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

The following core packages are required (full list available in the /functions folder):
- **Data processing**: data.table_1.18.2.1, dplyr_1.2.0, tidyr_1.3.2 
- **Data analysis**:   Mfuzz_2.62.0, lme4_2.0-1, emmeans_2.0.4, PKNCA_0.12.1, glmnet_4.1-10, cvTools_0.3.3, PredictABEL_1.2-4
- **Visualization**:   ggplot2_4.0.2, ComplexUpset_1.3.3, ComplexHeatmap_2.25.3, circlize_0.4.15, gridExtra_2.3, cowplot_1.2.0


## 📧 Contact
For questions regarding the analysis or code, please contact:

Huan Yun, Harvard T.H. Chan School of Public Health
huanyun@hsph.harvard.edu
