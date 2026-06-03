# 05_subtypes.R
# Defines molecular subtypes from receptor status (ER/PR/HER2)
# and runs survival analysis per subtype.
#
# BACKGROUND:
# Breast cancer is not one disease — it has 4 main molecular subtypes
# with very different prognoses and treatment responses:
#   Luminal A   : ER+/PR+, HER2-  → best prognosis, hormone therapy
#   Luminal B   : ER+/PR+, HER2+  → intermediate, hormone + targeted therapy
#   HER2-enrich : ER-, PR-, HER2+ → aggressive, but HER2-targeted therapy helps
#   Triple-neg  : ER-, PR-, HER2- → worst prognosis, only chemotherapy
#
# Input:  data/processed/survival_df.rds
# Output: data/processed/survival_subtypes.rds
#         results/figures/km_subtypes.png
#         results/figures/subtype_distribution.png

source("config.R")
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))

# --- Define subtypes ---
# Normalize receptor status values to Positive / Negative / Unknown
norm_status <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("positive")                          ~ "Positive",
    x %in% c("negative")                          ~ "Negative",
    x %in% c("indeterminate", "[not evaluated]",
              "[not available]", "equivocal", "")  ~ NA_character_,
    TRUE                                           ~ NA_character_
  )
}

df <- survival_df |>
  mutate(
    er  = norm_status(er_status_by_ihc),
    pr  = norm_status(pr_status_by_ihc),
    # HER2: FISH is gold standard; fall back to IHC if FISH missing
    her2_fish = norm_status(her2_fish_status),
    her2_ihc  = norm_status(her2_status_by_ihc),
    her2 = case_when(
      !is.na(her2_fish) ~ her2_fish,
      !is.na(her2_ihc)  ~ her2_ihc,
      TRUE              ~ NA_character_
    ),
    subtype = case_when(
      (er == "Positive" | pr == "Positive") & her2 == "Negative" ~ "Luminal A",
      (er == "Positive" | pr == "Positive") & her2 == "Positive" ~ "Luminal B",
      er == "Negative" & pr == "Negative"  & her2 == "Positive"  ~ "HER2-enriched",
      er == "Negative" & pr == "Negative"  & her2 == "Negative"  ~ "Triple-Negative",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(subtype))

message("Subtype distribution:")
print(table(df$subtype))
message("Total patients with subtype: ", nrow(df))

saveRDS(df, file.path(PROCESSED_DIR, "survival_subtypes.rds"))

# --- Bar chart: subtype distribution ---
subtype_colors <- c(
  "Luminal A"      = "#2196F3",
  "Luminal B"      = "#4CAF50",
  "HER2-enriched"  = "#FF9800",
  "Triple-Negative"= "#F44336"
)

p_bar <- df |>
  count(subtype) |>
  mutate(pct = round(n / sum(n) * 100, 1),
         label = paste0(n, "\n(", pct, "%)")) |>
  ggplot(aes(x = reorder(subtype, -n), y = n, fill = subtype)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = subtype_colors) +
  labs(title = "TCGA-BRCA: Molecular Subtype Distribution",
       x = NULL, y = "Number of patients") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(FIGURES_DIR, "subtype_distribution.png"),
       plot = p_bar, width = 7, height = 5, dpi = 150)
message("Saved: subtype_distribution.png")

# --- Kaplan-Meier by subtype ---
df$subtype <- factor(df$subtype,
  levels = c("Luminal A", "Luminal B", "HER2-enriched", "Triple-Negative"))

km_fit <- survfit(Surv(OS.time / 30.44, OS) ~ subtype, data = df)

p_km <- ggsurvplot(
  km_fit, data = df,
  palette    = unname(subtype_colors),
  risk.table = TRUE,
  pval       = TRUE,
  conf.int   = TRUE,
  xlab       = "Time (months)",
  ylab       = "Overall Survival probability",
  title      = "TCGA-BRCA: Overall Survival by Molecular Subtype",
  legend.title = "Subtype",
  ggtheme    = theme_bw()
)

ggsave(file.path(FIGURES_DIR, "km_subtypes.png"),
       plot = p_km$plot, width = 10, height = 6, dpi = 150)
message("Saved: km_subtypes.png")

# --- Cox regression with subtype ---
# Remove subtypes with < 30 patients (too small for stable Cox estimation)
subtype_counts <- table(df$subtype)
valid_subtypes  <- names(subtype_counts[subtype_counts >= 30])
message("Subtypes included in Cox (n >= 30): ", paste(valid_subtypes, collapse = ", "))

df_cox <- df |>
  filter(subtype %in% valid_subtypes) |>
  mutate(
    subtype          = factor(subtype, levels = valid_subtypes),
    age_at_diagnosis = as.numeric(age_at_diagnosis)
  )

cox_subtype <- coxph(
  Surv(OS.time / 30.44, OS) ~ subtype + age_at_diagnosis,
  data = df_cox
)

message("\nCox regression with subtype:")
print(summary(cox_subtype)$conf.int)

p_forest <- ggforest(cox_subtype, data = df_cox,
                     main = "Cox Regression: Molecular Subtypes")
ggsave(file.path(FIGURES_DIR, "cox_subtypes_forest.png"),
       plot = p_forest, width = 9, height = 5, dpi = 150)
message("Saved: cox_subtypes_forest.png")

message("\nSubtype analysis complete. Next: source('R/06_lasso_cox.R')")
