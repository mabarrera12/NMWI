# A Health-Referenced Nasal Microbiome Wellness Index Derived from Multi-Study 16S rRNA Data

**Barrera-Suarez, M.A., et al. — manuscript under review**

The **Nasal Microbiome Wellness Index (NMWI)** is a continuous, health-referenced score that summarizes how strongly a nasal microbiome profile resembles microbiome configurations observed in healthy versus non-healthy individuals.

The NMWI was developed from multi-study 16S rRNA gene amplicon data using **L1-regularized (LASSO) logistic regression**. The model uses taxonomic relative abundances and returns the fitted **log-odds of the healthy/control state**:

- **NMWI > 0:** more healthy-like microbiome configuration
- **NMWI < 0:** more non-healthy-like microbiome configuration
- **|NMWI|:** strength of the model prediction

A score of zero is the model-derived decision boundary and should not be interpreted as a universal biological or clinical cutoff.

> **Repository scope:** this repository contains curated metadata, processed taxonomic feature tables, derived analysis outputs, model-analysis code, and notebooks used to reproduce the NMWI manuscript analyses and figures.

---

## Study overview

The NMWI development cohort contains:

- **1,654 nasal microbiome samples**
- **589 healthy/control samples**
- **1,065 non-healthy/case samples**
- **27 studies**
- **7 chronic disease categories**
- **184 retained taxonomic features**

The independent validation cohort contains **350 samples from 11 studies**.

The final NMWI model retains **24 non-zero taxonomic coefficients** from the 184-feature input space:

- **14 positive coefficients**
- **10 negative coefficients**
- **160 coefficients shrunk to zero**

---

## Repository structure

```text
NMWI/
├── Codes/
│   ├── Auxiliary_files/
│   │   ├── 10x10cv.py
│   │   ├── Generate_data_index.py
│   │   ├── nested_LOSO.py
│   │   ├── NMWI.py
│   │   └── PCoA_analysis_.R
│   │
│   └── Manuscript-Figures/
│       ├── Table1.ipynb
│       ├── Figure1.ipynb
│       ├── Figure2.ipynb
│       ├── Figure3.ipynb
│       ├── Figure4.ipynb
│       └── Figure5.ipynb
│
├── Data/
│   ├── Input/
│   │   ├── Additional file 1.xlsx
│   │   └── AlphaDiversity_metrics.csv
│   │
│   └── Output/
│       ├── 10x10_cutoff_fold_results.csv
│       ├── 10x10_cutoff_summary.csv
│       ├── 10x10_oof_predictions.csv
│       ├── NMWI_coefficients.csv
│       ├── nested_loso_coefficients.tsv
│       ├── nested_loso_oof_predictions.tsv
│       └── nested_loso_summary.tsv
│
├── Feature_tables/
│   └── Final__<BioProject>_final_table_depth.tsv
│
├── .gitignore
├── CITATION.cff
├── environment.yml
├── LICENSE
└── README.md
```

`Data/Input/` contains source files used directly by the analysis code. `Data/Output/` contains derived analysis results used by downstream notebooks. `Feature_tables/` contains the processed BioProject-level taxonomic feature tables.

---

## NMWI model

The final NMWI is implemented with scikit-learn logistic regression using:

```python
LogisticRegression(
    penalty="l1",
    solver="liblinear",
    class_weight="balanced",
    C=0.9375,
    max_iter=2000,
    random_state=42
)
```

The model is fit to the **184 relative-abundance features without additional transformation or standardization**.The NMWI score for a sample is obtained using the model's `decision_function`, corresponding to the fitted log-odds of the healthy/control class.

---

## Environment and software requirements

The Python analyses were developed and tested using the versions specified in [`environment.yml`](environment.yml):

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

The PCoA/PERMANOVA analysis for Figure 1g uses R with:

```text
vegan: 2.7-1
ade4: 1.7-23
ggplot2: 3.5.2
```

Specific R package versions can be installed with:

```r
install.packages("remotes")

remotes::install_version("vegan", version = "2.7-1")
remotes::install_version("ade4", version = "1.7-23")
remotes::install_version("ggplot2", version = "3.5.2")
```

---

## How to run

### 1. Clone the repository

```bash
git clone https://github.com/mabarrera12/NMWI.git
cd NMWI
```

All commands below assume the repository root (`NMWI/`) is the working directory.

### 2. Create and activate the Python environment

```bash
conda env create -f environment.yml
conda activate nmwi
```

Start Jupyter with:

```bash
jupyter notebook
```

The manuscript notebooks contain setup code that identifies the `NMWI` repository root and sets it as the working directory.

## How to run

Run the manuscript notebooks in:

```text
Codes/Manuscript-Figures/
```

Each notebook reproduces the corresponding manuscript figure or table and calls the required auxiliary code when needed. Computationally intensive analyses, including nested LOSO and repeated 10×10 cross-validation, have already been run. Their outputs are stored in:

```text
Data/Output/
```

To reduce runtime, the notebooks read these precomputed CSV/TSV files rather than rerunning the analyses. The notebooks indicate which auxiliary scripts generated the corresponding outputs. Users who wish to reproduce the intensive analyses from scratch can run the scripts in:

```text
Codes/Auxiliary_files/
```
---

## Data and outputs

### `Data/Input/`

- `Additional file 1.xlsx` — study- and sample-level metadata used to define the development and validation cohorts.
- `AlphaDiversity_metrics.csv` — alpha-diversity metrics computed using QIIME 2 and exported to CSV for downstream analyses.

### `Data/Output/`

Contains derived results used by manuscript notebooks, including:

- Nested LOSO held-out predictions and model coefficients
- Repeated 10×10 CV predictions and cutoff analyses
- Final NMWI model coefficients

---

## Reproducibility notes

- Run notebooks **top to bottom in a fresh kernel**.
- Use the package versions specified in `environment.yml`.
- Run scripts from the repository root unless the script explicitly identifies the project root itself.
- Precomputed cross-validation outputs are included to avoid requiring users to rerun computationally intensive analyses before viewing manuscript figures.
- The final NMWI model uses relative abundances directly, without additional transformation or feature standardization.

---

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). GitHub's **Cite this repository** option can also be used.

The preferred citation is the associated manuscript:

> Barrera-Suarez, M.A., et al. *A Health-Referenced Nasal Microbiome Wellness Index Derived from Multi-Study 16S rRNA Data.* Manuscript under review.

The citation metadata will be updated with the final journal and DOI after publication.

---

## License

This repository is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.
