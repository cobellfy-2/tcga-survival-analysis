# 03_survival_analysis.R
# Kaplan-Meier curves and multivariable Cox regression on TCGA-BRCA.
# Input:  data/processed/survival_df.rds
# Output: results/figures/km_stage.png, km_age.png, cox_forest.png
#         results/cox_results.rds

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))

# --- Prepare covariates ---
df <- survival_df |>
  mutate(
    age_group = ifelse(
      suppressWarnings(as.numeric(age_at_diagnosis)) >= 60,
      ">=60", "<60"),
    stage_clean = case_when(
      grepl("Stage I[^V]|Stage I$", ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "I/II",
      grepl("Stage II",             ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "I/II",
      grepl("Stage III",            ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "III/IV",
      grepl("Stage IV",             ajcc_pathologic_tumor_stage, ignore.case = TRUE) ~ "III/IV",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(stage_clean), !is.na(age_group))

surv_obj <- Surv(df$OS.time / 30.44, df$OS)  # convert days -> months

# --- KM: Tumor stage ---
km_stage <- survfit(surv_obj ~ stage_clean, data = df)

p_stage <- ggsurvplot(
  km_stage, data = df,
  palette      = c("#2196F3", "#F44336"),
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = TRUE,
  xlab         = "Time (months)",
  ylab         = "Overall Survival probability",
  title        = "TCGA-BRCA: Overall Survival by Tumor Stage",
  legend.title = "Stage",
  ggtheme      = theme_bw()
)
ggsave(file.path(FIGURES_DIR, "km_stage.png"),
       plot = p_stage$plot, width = 8, height = 5, dpi = 150)
message("Saved: km_stage.png")

# --- KM: Age group ---
km_age <- survfit(surv_obj ~ age_group, data = df)

p_age <- ggsurvplot(
  km_age, data = df,
  palette      = c("#4CAF50", "#FF9800"),
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = TRUE,
  xlab         = "Time (months)",
  ylab         = "Overall Survival probability",
  title        = "TCGA-BRCA: Overall Survival by Age Group",
  legend.title = "Age",
  ggtheme      = theme_bw()
)
ggsave(file.path(FIGURES_DIR, "km_age.png"),
       plot = p_age$plot, width = 8, height = 5, dpi = 150)
message("Saved: km_age.png")

# --- Multivariable Cox regression ---
cox_model <- coxph(
  Surv(OS.time / 30.44, OS) ~ age_group + stage_clean,
  data = df
)

cox_summary <- summary(cox_model)
message("\nCox model summary:")
print(cox_summary)

saveRDS(list(model = cox_model, summary = cox_summary),
        file.path(RESULTS_DIR, "cox_results.rds"))

# Forest plot
p_forest <- ggforest(cox_model, data = df,
                     main = "TCGA-BRCA Multivariable Cox Regression")
ggsave(file.path(FIGURES_DIR, "cox_forest.png"),
       plot = p_forest, width = 8, height = 4, dpi = 150)
message("Saved: cox_forest.png")

message("\nSurvival analysis complete. Next: run R/04_gene_survival.R")
