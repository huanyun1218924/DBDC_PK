# Metabolomic Discovery of Multi-fluid Biomarkers of Food Intake in a Randomized Acute Feeding Trial


## 📌 Overview

This repository contains the simulated sample data and main analysis code for a integrative study investigating the associations between diet, circulating metabolome in relation to incident type 2 diabetes. The analysis utilizes diet, metabolomics data from food, plasma, and urine samples.


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

## 🔬 Sample data
- plasma_sample.Rdata: A simulated plasma dataset with the same sample size and metabolite variables as the original dataset

- urine_sample.RData: A simulated urine dataset with the same sample size and metabolite variables as the original dataset


## 🧬 Main analysis code
- Figure 2.R: Time-course metabolites change patterns after food intake;

- Figure 3.R: Identification of plasma metabolites in response to food intake;

- Figure 4.R: Identification of urinary metabolites in response to food intake

- Figure 5.R: Multi-metabolite signature differentiating different foods;

- Figure 6.R: Kinetic parameters of annotated metabolites.


## 💻 Environment & Dependencies
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 9.7 (Blue Onyx)

Matrix products: default
BLAS:   /app/R-4.3.3@i86-rhel9.0/lib64/R/lib/libRblas.so 
LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

The following core packages are required (full list available in the /functions folder):
- Data processing: tidyr_1.3.2    plyr_1.8.9  dplyr_1.2.0  data.table_1.18.2.1
- Model training:  cvTools_0.3.3  glmnet_4.1-10
- Data analysis:   survival_3.8-3 metafor_4.8-0  coxme_2.2-22  lmerTest_3.2-1  Maaslin2_1.16.0
- Visualization:   ComplexHeatmap_2.25.3  RColorBrewer_1.1-3  circlize_0.4.15  gridExtra_2.3  ComplexHeatmap_2.25.3   circlize_0.4.15  venn_1.12  cowplot_1.2.0    


## 📧 Contact
For questions regarding the analysis or code, please contact:

Huan Yun, Harvard T.H. Chan School of Public Health
huanyun@hsph.harvard.edu
