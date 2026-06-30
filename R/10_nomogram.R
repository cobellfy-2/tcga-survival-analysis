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
     xfrac    = 0.35,
     cex.axis = 0.8,
     cex.var  = 1.0,
     lwd      = 2,
     col.grid = gray(c(0.9, 0.97)))
title(main = "Nomogram — Predicted 3-year & 5-year Overall Survival (TCGA-BRCA)",
      cex.main = 1.2)
dev.off()
message("Saved: nomogram.png")

# --- Calibration plot (3-year) ---
# Custom calibration: group patients by predicted 3-year survival into quartiles,
# then compare mean predicted vs. observed Kaplan-Meier survival per group.
# This is far more robust than rms::calibrate for low-event cohorts like BRCA,
# where predicted probabilities cluster near 1.0 and bootstrap binning collapses.
df$lp_cal   <- predict(fit_nomo, type = "lp")
df$pred_3yr <- surv_fn(36, df$lp_cal)

# Quartile groups (drop duplicate breaks if predictions are tightly clustered)
brks <- unique(quantile(df$pred_3yr, probs = seq(0, 1, 0.25), na.rm = TRUE))
df$cal_grp <- cut(df$pred_3yr, breaks = brks, include.lowest = TRUE, labels = FALSE)

cal_list <- lapply(sort(unique(df$cal_grp)), function(g) {
  sub <- df[df$cal_grp == g & !is.na(df$cal_grp), ]
  km  <- survfit(Surv(OS.time_m, OS) ~ 1, data = sub)
  s36 <- summary(km, times = 36, extend = TRUE)
  data.frame(
    predicted = mean(sub$pred_3yr, na.rm = TRUE),
    observed  = s36$surv,
    lower     = s36$lower,
    upper     = s36$upper,
    n         = nrow(sub)
  )
})
cal_df <- do.call(rbind, cal_list)

rng <- range(c(cal_df$predicted, cal_df$lower, cal_df$upper), na.rm = TRUE)
pad <- 0.03
lim <- c(max(0, rng[1] - pad), min(1, rng[2] + pad))

p_cal <- ggplot(cal_df, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#ff4fa3", linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = (lim[2] - lim[1]) * 0.02, color = "#b983ff", linewidth = 0.9) +
  geom_line(color = "#ff85c0", linewidth = 1) +
  geom_point(size = 4, color = "#ff4fa3") +
  coord_cartesian(xlim = lim, ylim = lim) +
  labs(
    title    = paste0("Calibration Plot — 3-year Overall Survival (~",
                      round(mean(cal_df$n)), " patients per quartile)"),
    subtitle = "Pink dashed line = perfect calibration · error bars = 95% CI",
    x = "Predicted 3-year Survival Probability",
    y = "Observed 3-year Survival (Kaplan-Meier)"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", color = "#d6357f"))

ggsave(file.path(FIGURES_DIR, "calibration_3yr.png"),
       plot = p_cal, width = 8.5, height = 6, dpi = 150)
message("Saved: calibration_3yr.png")

saveRDS(list(fit = fit_nomo, nomo = nomo),
        file.path(RESULTS_DIR, "nomogram_results.rds"))

message("\nNomogram complete. Next: source('R/11_mutations.R')")
