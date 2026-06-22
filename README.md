# A Disease-Agnostic Nasal Microbiome Wellness Index for Standardized Assessment of Upper-Airway Respiratory Health. 
Barrera-Suarez, M.A., _et al_, manuscript under review

## Overview

**NMWI (Nasal Microbiome Wellness Index)** is a machine learning model that produces a single, continuous, disease-agnostic score of nasal microbiome health from 16S rRNA gene taxonomic abundance data. Rather than counting taxa or summarizing diversity, it learns which taxa—and in what balance—characterize a healthy nose, and returns the predicted log-odds that a profile resembles a healthy state. The model is built on LASSO-penalized (L1-regularized) logistic regression, which selects a parsimonious signature of 24 taxa whose weighted relative abundances determine the score. Positive scores indicate resemblance to healthy configurations and negative scores a shift toward non-healthy states, with the magnitude reflecting confidence; binary classification is an optional downstream step obtained by thresholding, not the model's primary output.
    
This repository implements the full pipeline behind the index: standardized 16S rRNA gene sequence processing, taxonomic annotation, model training and hyperparameter tuning, multi-framework validation (leave-one-study-out, leave-one-disease-out, and external cohorts), and visualization. The index was trained on 1,654 nasal microbiome samples (589 healthy, 1,065 non-healthy) pooled from 27 publicly available studies spanning seven chronic disease conditions across four continents, and is computed directly from existing 16S rRNA gene data without any new sequencing.

## Features

- **Data Integration**: Merges genus-level abundance tables from multiple NCBI BioProjects
- **Machine Learning Model**: Logistic Regression classifier with optimized hyperparameters
- **Cross-Validation**: 10x10 repeated stratified k-fold cross-validation for robust performance evaluation
- **Visualization**: Jupyter notebooks for generating publication-ready figures
- **Preprocessing**: Automated filtering and normalization of microbiome data

## Project Structure

```
NMWI/
├── Codes/
│   ├── Auxiliary_files/
│   │   ├── Generate_data_index.py      # Training data preparation
│   │   ├── Generate_data_val.py        # Validation data preparation
│   │   ├── Generate_data_rarefied.py   # Training-rarefied data preparation
│   │   ├── LOSO.py                     # Leave-one-study-out cross validation
│   │   ├── NMWI.py                     # Full training of NMWI L1-logistic regression
│   │   ├── PCOA_analysis.R             # PCoA analysis of the training data
│   │   └── Diff_abundance_testing.R    # Differential abundance training for genera and family levels
│   └── Notebooks/                      # Jupyter notebooks for figures
│       ├── Figure1.ipynb
│       ├── Figure2.ipynb
│       ├── Figure3.ipynb
│       ├── Figure4.ipynb
│       └── Figure5.ipynb
├── Data/
│   └── Supplementary_table1.xlsx       # Metadata and sample information
├── Features_tables/                    # Processed genus abundance tables
│   └── Final__*_final_table_depth.tsv
└── README.md
```

## Requirements

- Python 3.7+
- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn
- openpyxl (for Excel file reading)
- joblib

## Data Description

The project uses microbiome sequencing data from the NCBI Sequence Read Archive (SRA):

- **Training Dataset**: Multiple BioProjects with case/control samples
- **Validation Dataset**: Independent BioProjects for model validation
- **Features**: Genus-level taxonomic abundance (relative abundance)
- **Metadata**: Sample information including BioProject IDs, categories, and disease status

## Usage

Run the notebooks sequentially to generate the figures and the corresponding data visualization.

## Model Details

- **Algorithm**: Logistic Regression with L1 penalty
- **Hyperparameters**:
  - C = 0.9375 (regularization strength)
  - class_weight = "balanced"
  - solver = "liblinear"
  - max_iter = 2000
- **Evaluation**: Balanced accuracy score
- **Cross-validation**: 10x10 Repeated Stratified K-Fold

## Output Files

- `NMWI_coefficients.csv`: Feature importance (genus coefficients)
- Various figure files from notebooks
