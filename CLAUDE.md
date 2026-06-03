# TCGA Survival Analysis — Project Context for Claude

## Project Goal
Reproducible survival analysis pipeline on TCGA breast cancer (BRCA) data.
Published as a GitHub template others can reuse for any TCGA cancer type.

## Stack
- **R**: data download (TCGAbiolinks), preprocessing, KM curves, Cox regression
- **Python**: interactive visualizations (plotly), optional ML (scikit-survival)
- **Report**: Quarto (.qmd) rendering to HTML

## Workflow Rules (Phase 2 — Plan → Execute → Review)
- Always plan before writing code; confirm plan with user before executing
- One script = one atomic task; no multi-task scripts
- After each script is written: run it, check output, review before moving on
- Never silently skip errors — surface them immediately

## File Conventions
- R scripts: `R/01_download.R`, `R/02_preprocess.R`, etc. (numbered, snake_case)
- Python scripts: `python/01_visualize.py`, etc.
- Raw data: `data/raw/` — never modified after download
- Processed data: `data/processed/` — outputs of R preprocessing
- Figures: `results/figures/`
- Report: `report/report.qmd`

## Dataset
- TCGA-BRCA (Breast Invasive Carcinoma)
- Clinical data + RNA-seq (HTSeq counts) via TCGAbiolinks or cBioPortal fallback
- Target: Overall Survival (OS_STATUS, OS_MONTHS)

## Key R Packages
- TCGAbiolinks, SummarizedExperiment
- survival, survminer
- DESeq2, limma (if differential expression added later)
- ggplot2, dplyr, tidyr

## Key Python Packages
- pandas, numpy
- plotly, matplotlib
- lifelines (survival), scikit-survival (ML extension)

## GitHub Repo Goal
- Clean README with badges, methods summary, example output figures
- `.gitignore` excludes `data/raw/` (large files) — users re-run download script
- MIT License
- Reusable: cancer type swappable via one config variable
