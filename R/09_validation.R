# 09_validation.R
# Model validation: C-index, Brier Score, train/test comparison.
#
# BACKGROUND:
# Building a model on the full dataset and reporting its p-values is
# CIRCULAR — the model was optimized on the same data it's evaluated on.
# Proper validation requires splitting data into:
#   Training set (70%)  → fit the model
#   Test set    (30%)   → evaluate on unseen data
#
# Metrics:
#   C-index (Concordance index): probability that a randomly chosen patient
#     who died earlier has a higher predicted risk than one who survived longer.
#     0.5 = random, 1.0 = perfect. >0.65 is considered good for survival models.
#
#   Brier Score: mean squared error between predicted survival probability
#     and actual outcome. Lower = better. Null model (no predictors) ≈ 0.25.
#     Integrated Brier Score (IBS) summarizes across all time points.
#
# We compare three models:
#   1. Clinical only  (age + stage)
#   2. Lasso signature (risk score from script 06)
#   3. Combined       (clinical + risk score)
#
# Input:  data/processed/survival_df.rds, lasso_results.rds, counts_filtered.rds
# Output: results/validation_results.rds
#         results/figures/validation_cindex.png
#         results/figures/validation_brier.png
#         results/figures/calibration.png

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "dplyr",
              "tibble", "pec", "prodlim")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(survival)
library(ggplot2)
library(dplyr)
library(tibble)
library(pec)
library(prodlim)

set.seed(42)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw  <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))
lasso_res   <- readRDS(file.path(RESULTS_DIR,   "lasso_results.rds"))

# --- Prepare data ---
common <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
surv_df <- survival_df[match(common, survival_df$bcr_patient_barcode), ]
counts  <- counts_raw[, common]

# log2(CPM+1)
lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm   <- log2(cpm + 1)

# Lasso risk score
sig_genes  <- lasso_res$sig_df$ensembl
sig_coefs  <- lasso_res$sig_df$coef
valid_genes <- intersect(sig_genes, rownames(log_cpm))
sig_coefs_v <- sig_coefs[sig_genes %in% valid_genes]

risk_score <- as.vector(t(log_cpm[valid_genes, ]) %*% sig_coefs_v)

df <- surv_df |>
  mutate(
    age_num     = as.numeric(age_at_diagnosis),
    stage_clean = case_when(
      grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
      grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
      grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
      grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
      TRUE ~ NA_character_
    ),
    risk_score  = risk_score,
    OS.time_m   = OS.time / 30.44
  ) |>
  filter(!is.na(stage_clean), !is.na(age_num), OS.time_m > 0)

# --- 70/30 Train/Test split ---
n       <- nrow(df)
train_i <- sample(seq_len(n), size = floor(0.7 * n))
train   <- df[train_i,  ]
test    <- df[-train_i, ]
message("Train: ", nrow(train), " | Test: ", nrow(test))

# --- Fit three Cox models on TRAINING set ---
cox_clinical <- coxph(Surv(OS.time_m, OS) ~ age_num + stage_clean,
                      data = train, x = TRUE)
cox_lasso    <- coxph(Surv(OS.time_m, OS) ~ risk_score,
                      data = train, x = TRUE)
cox_combined <- coxph(Surv(OS.time_m, OS) ~ age_num + stage_clean + risk_score,
                      data = train, x = TRUE)

# --- C-index on TEST set ---
c_clinical <- concordance(cox_clinical, newdata = test)$concordance
c_lasso    <- concordance(cox_lasso,    newdata = test)$concordance
c_combined <- concordance(cox_combined, newdata = test)$concordance
c_null     <- 0.5

message("\nC-index on test set:")
message("  Null model:       0.500")
message("  Clinical only:    ", round(c_clinical, 3))
message("  Lasso signature:  ", round(c_lasso,    3))
message("  Combined:         ", round(c_combined, 3))

# --- Brier Score via pec ---
# Use fixed time points up to 80th percentile of follow-up
max_time <- quantile(df$OS.time_m, 0.80)
times    <- seq(12, max_time, by = 6)

brier_res <- tryCatch(
  pec(list(
    "Clinical"  = cox_clinical,
    "Lasso"     = cox_lasso,
    "Combined"  = cox_combined),
    formula = Surv(OS.time_m, OS) ~ 1,
    data    = test,
    times   = times,
    exact   = FALSE),
  error = function(e) { message("Brier score error: ", e$message); NULL })

# --- Visualizations ---
# C-index bar chart
cindex_df <- tibble(
  Model   = c("Null", "Clinical\n(age+stage)", "Lasso\nsignature", "Combined"),
  Cindex  = c(c_null, c_clinical, c_lasso, c_combined),
  fill    = c("grey70", "#4CAF50", "#2196F3", "#F44336")
)

p_cindex <- ggplot(cindex_df, aes(x = reorder(Model, Cindex), y = Cindex, fill = fill)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(Cindex, 3)), hjust = -0.2, size = 4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_fill_identity() +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(cindex_df$Cindex, na.rm=TRUE) * 1.12)) +
  labs(title    = "Model Comparison: C-index on Test Set (30%)",
       subtitle = "C-index > 0.5 = better than random | > 0.65 = good",
       x = NULL, y = "C-index (concordance)") +
  theme_bw(base_size = 13)

ggsave(file.path(FIGURES_DIR, "validation_cindex.png"),
       plot = p_cindex, width = 8, height = 4, dpi = 150)
message("Saved: validation_cindex.png")

# Brier score plot
if (!is.null(brier_res)) {
  brier_df <- as.data.frame(brier_res$AppErr) |>
    mutate(time = brier_res$time) |>
    tidyr::pivot_longer(-time, names_to = "model", values_to = "brier")

  p_brier <- ggplot(brier_df, aes(x = time, y = brier, color = model)) +
    geom_line(linewidth = 1) +
    labs(title    = "Integrated Brier Score over Time (Test Set)",
         subtitle = "Lower = better prediction | Reference = Kaplan-Meier (no covariates)",
         x = "Time (months)", y = "Brier Score", color = "Model") +
    scale_color_manual(values = c(
      "Reference" = "grey60", "Clinical" = "#4CAF50",
      "Lasso"     = "#2196F3", "Combined" = "#F44336")) +
    theme_bw(base_size = 13)

  ggsave(file.path(FIGURES_DIR, "validation_brier.png"),
         plot = p_brier, width = 9, height = 5, dpi = 150)
  message("Saved: validation_brier.png")

  ibs <- crps(brier_res, times = times)
  message("\nIntegrated Brier Score (lower = better):")
  print(round(ibs, 4))
}

saveRDS(list(
  cindex = list(clinical = c_clinical, lasso = c_lasso, combined = c_combined),
  brier  = brier_res,
  models = list(clinical = cox_clinical, lasso = cox_lasso, combined = cox_combined),
  train_n = nrow(train), test_n = nrow(test)
), file.path(RESULTS_DIR, "validation_results.rds"))

message("\nValidation complete. Next: source('R/10_nomogram.R')")
