# Project-wide configuration — change CANCER_TYPE to reuse for other TCGA projects
CANCER_TYPE  <- "BRCA"          # TCGA project code
DATA_DIR     <- "data/raw"
PROCESSED_DIR <- "data/processed"
RESULTS_DIR  <- "results"
FIGURES_DIR  <- "results/figures"

# Survival endpoint
SURV_TIME   <- "OS.time"        # column name after preprocessing
SURV_STATUS <- "OS"             # 1 = event (death), 0 = censored

# RNA-seq filtering thresholds
MIN_COUNT   <- 10               # minimum count across samples
MIN_SAMPLES <- 0.2              # fraction of samples that must meet MIN_COUNT
