# 10_nomogram.R
# Clinical nomogram for individualized survival prediction.
#
# BACKGROUND:
# A nomogram is a graphical tool that translates a statistical model into
# a practical clinical instrument. For each predictor, a point score is
# assigned based on its value. The total score maps to a predicted
# survival probability at specific time points (e.g. 3-year, 5-year OS).
#
# Why this matters: a Cox model gives you a p-value and a HR.
# A nomogram gives a doctor something actionable:
#   "This 65-year-old patient with Stage III disease and high CCDC9B
#    expression has an estimated 5-year survival probability of 58%."
#
# We use the `rms` package (Regression Modeling Strategies by Frank Harrell)
# which is the gold standard for clinical prediction models in R.
#
# Input:  data/processed/survival_df.rds, lasso_results.rds, counts_filtered.rds
# Output: results/figures/nomogram.png
#         results/figures/calibration_3yr.png

source("config.R")

for (pkg in c("rms", "survival", "ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(rms)
library(survival)
library(ggplot2)
library(dplyr)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw  <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))
lasso_res   <- readRDS(file.path(RESULTS_DIR,   "lasso_results.rds"))

# --- Prepare data ---
common <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
surv_df <- survival_df[match(common, survival_df$bcr_patient_barcode), ]
counts  <- counts_raw[, common]

lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm   <- log2(cpm + 1)

sig_genes   <- lasso_res$sig_df$ensembl
sig_coefs   <- lasso_res$sig_df$coef
valid_genes <- intersect(sig_genes, rownames(log_cpm))
sig_coefs_v <- sig_coefs[sig_genes %in% valid_genes]
risk_score  <- as.vector(t(log_cpm[valid_genes, ]) %*% sig_coefs_v)

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
  filter(!is.na(stage_clean), !is.na(age_num), OS.time_m > 0) |>
  mutate(stage_clean = factor(stage_clean, levels = c("I","II","III","IV")))

# rms requires datadist for automatic knot placement
dd <- datadist(df)
options(datadist = "dd")

# --- Fit survival model with rms::cph ---
fit_nomo <- cph(
  Surv(OS.time_m, OS) ~ age_num + stage_clean + risk_score,
  data     = df,
  x        = TRUE,
  y        = TRUE,
  surv     = TRUE,
  time.inc = 12,
  control  = coxph.control(iter.max = 500)  # prevent non-convergence warnings
)

message("rms model summary:")
print(fit_nomo)

# --- Build nomogram ---
# Survival() returns a function f(times, lp); wrap for nomogram's fun= argument
surv_fn <- Survival(fit_nomo)
fun_3yr <- function(lp) surv_fn(36, lp)
fun_5yr <- function(lp) surv_fn(60, lp)

nomo <- nomogram(
  fit_nomo,
  fun      = list(fun_3yr, fun_5yr),
  funlabel = c("3-year OS", "5-year OS"),
  lp       = FALSE
)

# Save as PNG
png(file.path(FIGURES_DIR, "nomogram.png"),
    width = 1800, height = 900, res = 150)
par(mar = c(2, 2, 4, 2), bg = "white")
plot(nomo,
     xfrac       = 0.32,
     label.every = 1,
     col.grid    = gray(c(0.88, 0.96)),
     lplabel     = "Linear Predictor",
     fun.side    = c(1, 3),
     cex.axis    = 0.85,
     cex.var     = 1.0,
     col.conf    = c("#2196F3", "#F44336"),
     main        = "Nomogram — Predicted 3-year & 5-year Overall Survival (TCGA-BRCA)",
     cex.main    = 1.2)
dev.off()
message("Saved: nomogram.png")

# --- Calibration plot (3-year) ---
# Compares predicted vs. observed survival
png(file.path(FIGURES_DIR, "calibration_3yr.png"),
    width = 1200, height = 900, res = 150)
par(mar = c(5, 5, 5, 3), bg = "white", cex.lab = 1.2, cex.axis = 1.1)
cal <- suppressWarnings(
  calibrate(fit_nomo, cmethod = "KM", method = "boot",
            u = 36, B = 100, cuts = 5)
)
plot(cal,
     xlab  = "Predicted 3-year Survival Probability",
     ylab  = "Observed 3-year Survival Probability (Kaplan-Meier)",
     main  = "Calibration Plot — 3-year Overall Survival",
     col   = "#2196F3",
     lwd   = 2,
     subtitles = FALSE)
abline(0, 1, col = "#F44336", lty = 2, lwd = 2)
legend("topleft",
       legend = c("Calibration curve", "Ideal (perfect calibration)"),
       col    = c("#2196F3", "#F44336"),
       lty    = c(1, 2), lwd = 2, bty = "n", cex = 1.1)
dev.off()
message("Saved: calibration_3yr.png")

saveRDS(list(fit = fit_nomo, nomo = nomo),
        file.path(RESULTS_DIR, "nomogram_results.rds"))

message("\nNomogram complete. Next: source('R/11_mutations.R')")
