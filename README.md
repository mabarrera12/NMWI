## A Health-Referenced Nasal Microbiome Wellness Index Derived from Multi-Study 16S rRNA Data

**Barrera-Suarez, M.A., et al. — manuscript under review**

The **Nasal Microbiome Wellness Index (NMWI)** is a continuous, health-referenced score that summarizes how strongly a nasal microbiome profile resembles microbiome configurations observed in healthy versus non-healthy individuals.

The NMWI was developed from multi-study 16S rRNA gene amplicon data using **L1-regularized (LASSO) logistic regression**. The model uses taxonomic relative abundances and returns the fitted **log-odds of the healthy state**:

- **NMWI > 0:** more healthy-like microbiome configuration
- **NMWI < 0:** more non-healthy-like microbiome configuration
- **|NMWI|:** strength of the model prediction

A score of zero is the model-derived decision boundary, not a universal biological or clinical cutoff.

> **Repository scope:** this repository contains the processed taxonomic feature tables, curated metadata, derived analysis data, model-analysis code, and Manuscript-figures used for the NMWI study. 

---

### Study overview

The NMWI development cohort contains:

- **1,654 nasal microbiome samples**
- **589 healthy/control samples**
- **1,065 non-healthy/case samples**
- **27 studies**
- **7 chronic disease categories**
- **184 retained taxonomic features** after filtering

The independent validation cohort contains **350 samples from 11 studies** and is processed in `Codes/Manuscript-Figures/Figure5.ipynb`.

The final NMWI model retains **24 non-zero taxonomic coefficients** from the 184-feature input space.

---

### Repository structure

The scientifically relevant contents of the repository are organized as follows:

```text
NMWI/
├── Codes/
│   ├── Auxiliary_files/
│   │   ├── Cohort_summary.ipynb
│   │   ├── Generate_data_index.py
│   │   ├── Lasso_TableS6.py
│   │   ├── NMWI.py
│   │   └── PCOA_analysis.R
│   │
│   └── Manuscript-Figures/
│       ├── Figure1.ipynb
│       ├── Figure2.ipynb
│       ├── Figure3.ipynb
│       ├── Figure4.ipynb
│       └── Figure5.ipynb
│
├── Data/
│   ├── Additional file 1.xlsx
│   ├── AlphaDiversity_metrics.csv
│   ├── Genus_level_raw_count_featureTable.tsv
│   ├── loso.tsv
│   └── nested_lasso_NMWI_OOF_predictions.csv
│
├── Feature_tables/
│   └── Final__<BioProject>_final_table_depth.tsv
│       # 47 processed BioProject-level feature tables
│
├── LICENSE
└── README.md
```
---

### NMWI model

The final NMWI is implemented with scikit-learn logistic regression using:

```python
LogisticRegression(
    penalty="l1",
    solver="liblinear",
    class_weight="balanced",
    C=0.9375,
    max_iter=2000,
    random_state=42,
)
```

The model is fit to the **184 relative-abundance features without additional transformation or standardization**.

The NMWI score for a sample is obtained using the model's `decision_function`, i.e. the fitted log-odds for the healthy/control class.

The final model contains:

```text
24 non-zero coefficients
14 positive coefficients
10 negative coefficients
160 coefficients shrunk to zero
```
---

## Environment and software requirements
The Python analyses in this repository were developed and tested using the versions specified in environment.yml. Using the provided environment is recommended to ensure reproducibility of the manuscript analyses and figures.

The environment includes:

```text
Python: 3.12.2
pandas: 2.3.2
numpy: 2.2.6
scipy: 1.16.2
scikit-learn: 1.7.2
matplotlib: 3.10.6
seaborn: 0.13.2
openpyxl: 3.1.5
jupyter: 1.1.1
```

**Create the Conda environment**

From the repository root, run and activate the conda environment:

```text
conda env create -f environment.yml
conda activate nmwi
jupyter notebook
```

For the R analysis:

```r
remotes::install_version("vegan", version = "2.7.1")
remotes::install_version("ade4", version = "1.7.23")
remotes::install_version("ggplot2", version = "3.5.2")
```

---

### Notes:
The Manuscript-figures contain setup cells that change the working directory to the repository root. These cells assume the initial working directory corresponds to the Manuscript-figure's own directory; if your Jupyter/IDE configuration behaves differently, set the working directory to the repository root before running data-loading cells.

---

## License

This repository is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.

---
