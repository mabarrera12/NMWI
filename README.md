# A Disease-Agnostic Nasal Microbiome Wellness Index for Standardized Assessment of Upper-Airway Respiratory Health

## Overview

NMWI (Nasal Microbiome Wellness Index) is a machine learning model developed for classifying microbiome samples into non-healthy and healthy groups using genus-level taxonomic abundance data. The model utilizes Logistic Regression with L1 regularization to identify microbial signatures associated with disease states.

This project implements a comprehensive pipeline for microbiome data analysis, including data preprocessing, model training, validation, and visualization of results.

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
  - C = 0.973 (regularization strength)
  - class_weight = "balanced"
  - solver = "liblinear"
  - max_iter = 2000
- **Evaluation**: Balanced accuracy score
- **Cross-validation**: 10x10 Repeated Stratified K-Fold

## Output Files

- `NMWI_coefficients.csv`: Feature importance (genus coefficients)
- Various figure files from notebooks
