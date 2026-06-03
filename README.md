# TCGA Survival Analysis

**Reproducible pipeline for clinical and transcriptomic survival analysis on TCGA data.**

> Template: swap `CANCER_TYPE` in `config.R` to reuse for any TCGA cancer cohort.

---

## What this does

1. Downloads TCGA-BRCA clinical + RNA-seq data via the GDC API
2. Cleans survival endpoints and filters low-count genes
3. Runs Kaplan-Meier curves and multivariable Cox regression on clinical variables
4. Screens 5,000 genes with univariate Cox models → volcano plot + top-gene KM curves
5. Produces interactive Plotly figures (HTML) and a self-contained Quarto report

---

## Example outputs

| Figure | Description |
|--------|-------------|
| `km_stage.png` | Kaplan-Meier by AJCC tumor stage |
| `km_age.png` | Kaplan-Meier by age group (<60 vs ≥60) |
| `cox_forest.png` | Multivariable Cox forest plot |
| `volcano_cox.png` | Gene-level survival volcano plot |
| `km_top_genes.png` | KM curves for top 3 survival genes |
| `volcano_interactive.html` | Hover-enabled volcano (Plotly) |
| `km_interactive.html` | Interactive KM with CI bands (Plotly + lifelines) |

---

## Requirements

**R (≥ 4.2):**
```
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2"))
install.packages(c("survival", "survminer", "ggplot2", "dplyr", "tibble", "cowplot"))
```

**Python (≥ 3.9):**
```
pip install pyreadr plotly pandas numpy lifelines kaleido
```

**Quarto:** [quarto.org](https://quarto.org)

---

## Quickstart

```r
# 1. Clone and set working directory to repo root
# 2. (Optional) Edit config.R to change cancer type

source("R/01_download.R")           # ~20 min — downloads ~2 GB
source("R/02_preprocess.R")         # ~2 min
source("R/03_survival_analysis.R")  # ~1 min
source("R/04_gene_survival.R")      # ~10 min (5000 Cox models)
```

```bash
python python/01_interactive_viz.py
quarto render report/report.qmd
```

---

## Adapting to other cancer types

Edit `config.R`:
```r
CANCER_TYPE <- "LUAD"   # e.g., lung adenocarcinoma
```
Then re-run the pipeline. All TCGA project codes work
(LUAD, COAD, GBM, OV, KIRC, …).

---

## Methods

- **Normalization**: DESeq2 variance-stabilizing transform (VST)
- **Survival**: Kaplan-Meier + log-rank test; multivariable Cox PH regression
- **Gene screening**: univariate Cox per gene, Benjamini-Hochberg FDR correction
- **Visualization**: ggplot2 / survminer (static); Plotly + lifelines (interactive)

---

## License

MIT — free to use, adapt, and share.
