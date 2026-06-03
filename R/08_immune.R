# 08_immune.R
# Immune landscape estimation via ESTIMATE + survival correlation.
#
# BACKGROUND:
# Tumor tissue is not just cancer cells — it contains immune cells,
# fibroblasts, endothelial cells, etc. The composition of this
# "tumor microenvironment" strongly influences patient outcomes.
#
# ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumors
# using Expression data) uses gene expression signatures to compute:
#
#   ImmuneScore  — level of immune cell infiltration
#   StromalScore — level of stromal (connective tissue) cell infiltration
#   ESTIMATEScore — combined tumor purity estimate
#   TumorPurity  — fraction of the sample that is actual tumor cells
#
# Key insight: high ImmuneScore = more immune cells in tumor → in many
# cancers this correlates with better survival (immune surveillance).
# In BRCA the relationship is subtype-dependent.
#
# Input:  data/processed/counts_filtered.rds, survival_df.rds
# Output: results/immune_scores.rds
#         results/figures/immune_boxplot.png
#         results/figures/immune_survival_scatter.png
#         results/figures/km_immune.png

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "dplyr", "tibble", "ggpubr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("estimate", quietly = TRUE)) {
  # ESTIMATE is not on CRAN — install from GitHub via a tar.gz
  install.packages(
    "https://bioinformatics.mdanderson.org/estimate/rpackage/estimate_1.0.13.tar.gz",
    repos = NULL, type = "source"
  )
}

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw  <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))

# ESTIMATE needs a tab-separated gene expression file
# Use log2(CPM+1) and map to gene symbols first
message("Preparing expression matrix for ESTIMATE...")

common  <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
counts  <- counts_raw[, common]

lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm   <- log2(cpm + 1)

# Map Ensembl IDs to gene symbols (ESTIMATE needs symbols)
library(org.Hs.eg.db)
ensembl_ids <- sub("\\..*", "", rownames(log_cpm))
symbols <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = ensembl_ids,
         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
)

# Keep only genes with valid symbols, remove duplicates
valid   <- !is.na(symbols)
log_cpm <- log_cpm[valid, ]
rownames(log_cpm) <- symbols[valid]
log_cpm <- log_cpm[!duplicated(rownames(log_cpm)), ]

message("Expression matrix: ", nrow(log_cpm), " genes × ", ncol(log_cpm), " samples")

# Write temp file for ESTIMATE
tmp_expr <- file.path(PROCESSED_DIR, "expr_for_estimate.gct")
tmp_out  <- file.path(PROCESSED_DIR, "estimate_scores.gct")

# GCT format: 2 header lines, then NAME | Description | samples...
gct_mat <- cbind(NAME = rownames(log_cpm),
                 Description = rownames(log_cpm),
                 as.data.frame(log_cpm))

con <- file(tmp_expr, "w")
writeLines("#1.2", con)
writeLines(paste(nrow(log_cpm), ncol(log_cpm), sep = "\t"), con)
write.table(gct_mat, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

# Run ESTIMATE
library(estimate)
filterCommonGenes(input.f  = tmp_expr,
                  output.f = file.path(PROCESSED_DIR, "expr_filtered.gct"),
                  id       = "GeneSymbol")
estimateScore(file.path(PROCESSED_DIR, "expr_filtered.gct"),
              tmp_out, platform = "illumina")

# Parse output
scores_raw <- read.table(tmp_out, skip = 2, header = TRUE,
                          sep = "\t", check.names = FALSE)
scores_t <- as.data.frame(t(scores_raw[, -(1:2)]))
colnames(scores_t) <- scores_raw$NAME
scores_t$patient <- substr(rownames(scores_t), 1, 12)

immune_df <- scores_t |>
  dplyr::select(patient, StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity) |>
  mutate(across(c(StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity),
                as.numeric)) |>
  inner_join(survival_df, by = c("patient" = "bcr_patient_barcode"))

message("Immune scores computed for ", nrow(immune_df), " patients")
saveRDS(immune_df, file.path(RESULTS_DIR, "immune_scores.rds"))

# --- Boxplot: immune score by stage ---
stage_df <- immune_df |>
  mutate(stage_clean = case_when(
    grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
    grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
    grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
    grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
    TRUE ~ NA_character_
  )) |> filter(!is.na(stage_clean))

p_box <- ggboxplot(stage_df, x = "stage_clean", y = "ImmuneScore",
                   fill = "stage_clean",
                   palette = c("#2196F3","#4CAF50","#FF9800","#F44336"),
                   add = "jitter", add.params = list(alpha = 0.3, size = 0.8)) +
  stat_compare_means(label = "p.signif", ref.group = "I") +
  labs(title = "TCGA-BRCA: Immune Score by Tumor Stage",
       x = "Stage", y = "ESTIMATE Immune Score") +
  theme(legend.position = "none")

ggsave(file.path(FIGURES_DIR, "immune_boxplot.png"),
       plot = p_box, width = 7, height = 5, dpi = 150)
message("Saved: immune_boxplot.png")

# --- Scatter: ImmuneScore vs. TumorPurity ---
p_scatter <- ggplot(immune_df,
                    aes(x = ImmuneScore, y = TumorPurity)) +
  geom_point(alpha = 0.4, size = 1, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#F44336", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = 0.6) +
  labs(title   = "Immune Infiltration vs. Tumor Purity",
       subtitle = "Higher immune score → lower tumor purity (more immune cells)",
       x = "Immune Score", y = "Tumor Purity") +
  theme_bw()

ggsave(file.path(FIGURES_DIR, "immune_survival_scatter.png"),
       plot = p_scatter, width = 7, height = 5, dpi = 150)
message("Saved: immune_survival_scatter.png")

# --- KM: high vs. low immune score ---
immune_df <- immune_df |>
  mutate(immune_group = ifelse(ImmuneScore >= median(ImmuneScore),
                               "High immune", "Low immune"))

km_immune <- survfit(Surv(OS.time / 30.44, OS) ~ immune_group, data = immune_df)

p_km <- ggsurvplot(km_immune, data = immune_df,
                   palette      = c("#2196F3", "#F44336"),
                   risk.table   = TRUE,
                   pval         = TRUE,
                   conf.int     = TRUE,
                   xlab         = "Time (months)",
                   title        = "TCGA-BRCA: Survival by Immune Infiltration",
                   legend.title = "Immune Score",
                   ggtheme      = theme_bw())

ggsave(file.path(FIGURES_DIR, "km_immune.png"),
       plot = p_km$plot, width = 8, height = 5, dpi = 150)
message("Saved: km_immune.png")

message("\nImmune analysis complete. Next: Shiny app (app/)")
