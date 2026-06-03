# 01_download.R
# Downloads TCGA-BRCA clinical data and RNA-seq counts via TCGAbiolinks.
# Run once — outputs saved to data/raw/. Re-running is safe (uses cache).

source("config.R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("TCGAbiolinks", "SummarizedExperiment")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Clinical data ---
message("Downloading clinical data for TCGA-", CANCER_TYPE, "...")
clinical_query <- GDCquery(
  project      = paste0("TCGA-", CANCER_TYPE),
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)
GDCdownload(clinical_query, directory = DATA_DIR)
clinical_data <- GDCprepare(clinical_query, directory = DATA_DIR)

# BCR Biotab returns a list; the patient-level table is what we need
clinical_df <- clinical_data[["clinical_patient_brca"]]
saveRDS(clinical_df, file.path(PROCESSED_DIR, "clinical_raw.rds"))
message("Clinical data saved: ", nrow(clinical_df), " patients")

# --- RNA-seq (HTSeq counts, STAR) ---
message("Querying RNA-seq data (this may take a few minutes)...")
rna_query <- GDCquery(
  project           = paste0("TCGA-", CANCER_TYPE),
  data.category     = "Transcriptome Profiling",
  data.type         = "Gene Expression Quantification",
  workflow.type     = "STAR - Counts"
)
GDCdownload(rna_query, directory = DATA_DIR)
rna_se <- GDCprepare(rna_query, directory = DATA_DIR)

saveRDS(rna_se, file.path(PROCESSED_DIR, "rna_se_raw.rds"))
message("RNA-seq SummarizedExperiment saved: ",
        ncol(rna_se), " samples, ", nrow(rna_se), " genes")

message("\nDownload complete. Next: run R/02_preprocess.R")
