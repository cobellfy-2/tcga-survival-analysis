# 02_preprocess.R
# Cleans clinical data, extracts survival endpoint, filters RNA-seq counts.
# Input:  data/processed/clinical_raw.rds, rna_se_raw.rds
# Output: data/processed/survival_df.rds, counts_filtered.rds

source("config.R")
library(SummarizedExperiment)
library(dplyr)

dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Clinical / Survival ---
clinical_raw <- readRDS(file.path(PROCESSED_DIR, "clinical_raw.rds"))

survival_df <- clinical_raw |>
  as.data.frame() |>
  # BCR Biotab uses character columns; convert key fields
  mutate(
    bcr_patient_barcode = toupper(bcr_patient_barcode),
    OS.time = suppressWarnings(as.numeric(days_to_death)),
    # For censored patients, use days_to_last_followup
    OS.time = ifelse(is.na(OS.time),
                     suppressWarnings(as.numeric(days_to_last_followup)),
                     OS.time),
    OS = ifelse(!is.na(suppressWarnings(as.numeric(days_to_death))), 1L, 0L)
  ) |>
  filter(!is.na(OS.time), OS.time > 0) |>
  select(bcr_patient_barcode, OS.time, OS,
         age_at_initial_pathologic_diagnosis,
         ajcc_pathologic_tumor_stage,
         histological_type)

message("Survival data: ", nrow(survival_df), " patients after filtering")
message("Events (deaths): ", sum(survival_df$OS), " (",
        round(mean(survival_df$OS) * 100, 1), "%)")

saveRDS(survival_df, file.path(PROCESSED_DIR, "survival_df.rds"))

# --- RNA-seq filtering ---
rna_se <- readRDS(file.path(PROCESSED_DIR, "rna_se_raw.rds"))

counts <- assay(rna_se, "unstranded")

# Keep genes with >= MIN_COUNT in at least MIN_SAMPLES fraction of samples
keep <- rowSums(counts >= MIN_COUNT) >= (MIN_SAMPLES * ncol(counts))
counts_filtered <- counts[keep, ]

message("RNA-seq: ", nrow(counts_filtered), " genes retained after filtering (",
        nrow(counts) - nrow(counts_filtered), " removed)")

# Trim sample barcodes to patient level (first 12 chars) for merging with clinical
colnames(counts_filtered) <- substr(colnames(counts_filtered), 1, 12)

# Remove duplicate patient entries (keep first occurrence)
counts_filtered <- counts_filtered[, !duplicated(colnames(counts_filtered))]

saveRDS(counts_filtered, file.path(PROCESSED_DIR, "counts_filtered.rds"))

# Overlap between clinical and RNA-seq
overlap <- intersect(survival_df$bcr_patient_barcode, colnames(counts_filtered))
message("Patients with both clinical and RNA-seq data: ", length(overlap))

message("\nPreprocessing complete. Next: run R/03_survival_analysis.R")
