# Bias Analysis in Skeletal Stature Estimation
### A Comparative Study of Osteometric Methods Using the Goldman Dataset

---

## Overview

This repository contains the R script used for the statistical analysis of stature estimation bias across multiple osteometric methods applied to skeletal remains. The study evaluates the agreement, accuracy, and systematic differences between eleven widely used forensic anthropology formulas, using femur, tibia, humerus, and radius measurements from the Goldman Skeletal Collection.

---

## Repository Structure

```
/
├── analysis_stature_estimation.R   # Main analysis script
├── goldman_dataset.xlsx            # Input dataset (not included — see Data Access)
├── outputs/
│   ├── Estadistica_Descriptiva_Definitiva.csv
│   ├── Matriz_BlandAltman_Total.csv
│   ├── Grafico_5_BlandAltman_Extremos.png
│   ├── Grafico_6a_LMM_Metodos.png
│   ├── Grafico_6b_LMM_Paises.png
│   ├── Grafico_7a_Bayes_Metodos.png
│   ├── Grafico_7b_Bayes_Paises.png
│   └── Grafico_7c_Bayes_Huesos.png
└── README.md
```

---

## Workflow Description

The script follows a sequential, modular workflow organized into eight analytical steps.

### Step 0 — Data Loading and Outlier Cleaning

The script loads the Goldman dataset from an `.xlsx` file and cleans outliers by replacing any bone measurement exceeding 600 mm with `NA`. This threshold is used to remove digitization errors and missing-data codes (e.g., `999`).

### Step 1 — Definition of Stature Estimation Functions

Eleven stature estimation methods are defined as individual R functions. Each function accepts a bone identifier, a biological sex code (`0` = male, `1` = female), and a bone length in millimeters. All functions return stature estimates in centimeters.

The methods implemented are:

| Author(s) | Year | Bones Covered |
|---|---|---|
| Pearson | 1899 | Femur, Tibia, Humerus, Radius |
| Telkka | 1950 | Femur, Tibia, Humerus, Radius |
| Trotter & Gleser (European) | 1952 | Femur |
| Trotter & Gleser (African-American) | 1958 | Femur |
| Genovés | 1967 | Femur, Tibia |
| Del Angel & Cisneros | 2004 | Femur, Tibia, Humerus, Radius |
| Béguelin | 2011 | Femur, Tibia, Humerus, Radius |
| Belmonte | 2011 | Tibia (females only) |
| Ross | 2011 | Femur, Tibia, Humerus |
| Abarca | 2013 | Femur |
| Menéndez et al. | 2018 | Femur, Tibia, Humerus |

### Step 2 — Application of Formulas to the Dataset

All eleven methods are applied to each individual in the dataset using `dplyr::mutate()`. The result is a wide-format dataframe where each column represents the stature estimate produced by a specific method-bone combination (e.g., `Est_Pearson_Femur`, `Est_Ross_Tibia`).

### Step 3 — Descriptive Statistics

The wide dataframe is pivoted to long format. A descriptive summary table is calculated per method-bone combination, including sample size (N), mean, standard deviation, minimum, and maximum. Results are exported to `Estadistica_Descriptiva_Definitiva.csv`.

### Step 4 — Parametric Assumption Testing

Two assumption tests are conducted before inferential analysis:

- **Kolmogorov-Smirnov test**: evaluates whether the pooled stature estimates follow a normal distribution.
- **Levene's test** (`car` package): evaluates the homogeneity of variances across method-bone groups.

### Step 5 — Inferential Analysis

Because Levene's test detects heterogeneous variances across groups, the script applies:

- **Welch's one-way ANOVA** (`oneway.test(..., var.equal = FALSE)`): a robust alternative to classical ANOVA that does not assume equal variances.
- **Games-Howell post-hoc test** (`rstatix` package): identifies which specific pairs of methods differ significantly (adjusted p < 0.05).

### Step 6 — Bland-Altman Agreement Analysis

The script performs an exhaustive pairwise Bland-Altman analysis across all possible combinations of method-bone estimators. For each pair with more than 30 shared observations, it computes:

- Mean bias (systematic difference between methods)
- Standard deviation of differences
- 95% limits of agreement (mean ± 1.96 SD)
- Total amplitude of agreement interval

All results are exported to `Matriz_BlandAltman_Total.csv`. A forest plot (`Grafico_5_BlandAltman_Extremos.png`) visualizes the top 10 most concordant and the top 10 most discordant pairs.

### Step 7 — Linear Mixed Model (LMM)

A frequentist linear mixed model is fitted using the `lme4` package. The model estimates the effect of the estimation method, skeletal element (bone), and country of origin on stature estimates, while controlling for repeated measurements per individual (random intercept by `ID_Individuo`). Coefficients are visualized using `sjPlot::plot_model()` and saved as `Grafico_6a_LMM_Metodos.png` and `Grafico_6b_LMM_Paises.png`.

### Step 8 — Bayesian Generalized Linear Mixed Model (BGLMM)

A Bayesian version of the mixed model is fitted using the `brms` package with four MCMC chains (2,000 iterations each, 1,000 warmup). The model structure mirrors Step 7. Posterior means and 95% credibility intervals for each predictor type (method, bone, country) are extracted and plotted in three separate forest plots:

- `Grafico_7a_Bayes_Metodos.png`
- `Grafico_7b_Bayes_Paises.png`
- `Grafico_7c_Bayes_Huesos.png`

A practical equivalence band of ±2 cm is displayed on each plot to contextually interpret whether posterior differences are forensically meaningful.

---

## Requirements

### R Packages

```r
install.packages(c(
  "conflicted", "dplyr", "tidyr", "readr", "tidyverse",
  "readxl", "car", "rstatix", "ggplot2", "stringr",
  "lme4", "sjPlot", "brms", "tidybayes", "scales"
))
```

> **Note:** The `brms` package requires a working installation of Stan. See [mc-stan.org](https://mc-stan.org/users/interfaces/rstan) for setup instructions.

### Input Data

The script expects a file named `goldman_dataset.xlsx` containing at minimum the following columns:

| Column | Description |
|---|---|
| `sexo_biológico` | Biological sex (`0` = male, `1` = female) |
| `p_LMF` | Maximum length of the femur (mm) |
| `p_LMT` | Maximum length of the tibia (mm) |
| `p_LMH` | Maximum length of the humerus (mm) |
| `p_LMR` | Maximum length of the radius (mm) |
| `país` | Country of origin of the individual |

---

## Citation

If this script or analysis is used in academic work, please cite the original thesis and the corresponding osteometric references listed in Step 1.

---

## License

This code is shared for academic and reproducibility purposes. See `LICENSE` for details.
