# 11_mutations.R
# Somatic mutation landscape for TCGA-BRCA.
#
# BACKGROUND:
# RNA-seq tells us WHAT is being expressed, but not WHY. Somatic mutations
# are the drivers: specific DNA changes that initiate or reshape tumor biology.
# TCGA provides Mutation Annotation Format (MAF) files with every called
# somatic variant per patient.
#
# We compute:
#   1. Mutation frequency per gene (top 20 ranked bar chart)
#   2. Mutation type breakdown (missense, nonsense, frameshift, etc.)
#   3. Tumor mutational burden (TMB) by stage
#   4. Survival by mutation status for top driver genes
#   5. OncoPrint (co-mutation landscape) via maftools
#
# Input:  survival_df.rds  |  TCGA GDC API
# Output: results/mutation_results.rds
#         results/figures/mutation_freq.png
#         results/figures/mutation_tmb.png
#         results/figures/oncoprint.png
#         results/figures/mutation_cooccurrence.png

source("config.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in c("maftools", "survival", "survminer",
              "ggplot2", "dplyr", "tidyr", "ggpubr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "maftools") BiocManager::install(pkg, ask = FALSE)
    else install.packages(pkg)
  }
}

library(maftools)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))

# ── 1. Load mutations via maftools::tcgaLoad ─────────────────────────────────
# tcgaLoad() downloads pre-processed TCGA mutation data from GitHub
# (PoisonAlien/TCGAmutations) — no GDC tar.gz extraction needed, works on Windows.
maf_cache <- file.path(PROCESSED_DIR, "maf_obj.rds")

if (!file.exists(maf_cache)) {
  message("Downloading TCGA-BRCA mutations via maftools::tcgaLoad ...")
  # TCGAmutations is a GitHub-only data package; use BiocManager which handles GitHub installs
  if (!requireNamespace("TCGAmutations", quietly = TRUE)) {
    message("Installing TCGAmutations from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github("PoisonAlien/TCGAmutations", upgrade = "never")
  }
  maf_obj <- maftools::tcgaLoad(study = "BRCA")
  saveRDS(maf_obj, maf_cache)
  message("Cached MAF object.")
} else {
  message("Loading cached MAF object...")
  maf_obj <- readRDS(maf_cache)
}

# Extract the underlying data frame for custom analyses
maf_df <- as.data.frame(maf_obj@data)
names(maf_df) <- make.unique(names(maf_df))   # fix duplicate column names
maf_df$Tumor_Sample_Barcode <- substr(maf_df$Tumor_Sample_Barcode, 1, 12)
message("MAF: ", nrow(maf_df), " mutations")

n_patients <- length(unique(maf_df$Tumor_Sample_Barcode))
message("Patients in MAF: ", n_patients)

# ── 4. Mutation frequency table ──────────────────────────────────────────────
top_genes_mut <- getGeneSummary(maf_obj)
message("Top 10 mutated genes:")
print(head(top_genes_mut[, c("Hugo_Symbol", "MutatedSamples", "AlteredSamples")], 10))

# Mutation type breakdown for top 20
type_cols <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
               "Translation_Start_Site", "Nonstop_Mutation")
available_cols <- intersect(type_cols, colnames(top_genes_mut))
top20 <- head(top_genes_mut, 20)
top20$pct <- round(100 * top20$MutatedSamples / n_patients, 1)

type_long <- top20 |>
  as.data.frame() |>
  dplyr::select(Hugo_Symbol, pct, all_of(available_cols)) |>
  pivot_longer(cols = all_of(available_cols),
               names_to  = "mut_type",
               values_to = "n") |>
  filter(n > 0) |>
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(top20$Hugo_Symbol)))

# ── 5. TMB by stage ──────────────────────────────────────────────────────────
tmb_df <- maf_df |>
  filter(Variant_Classification %in%
           c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
             "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins")) |>
  group_by(Tumor_Sample_Barcode) |>
  summarise(tmb = n(), .groups = "drop") |>
  left_join(survival_df |>
              dplyr::select(bcr_patient_barcode, OS.time, OS,
                            ajcc_pathologic_tumor_stage),
            by = c("Tumor_Sample_Barcode" = "bcr_patient_barcode")) |>
  mutate(
    OS.time_m   = OS.time / 30.44,
    stage_clean = case_when(
      grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
      grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
      grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
      grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
      TRUE ~ NA_character_
    ),
    tmb_group = ifelse(tmb >= median(tmb, na.rm = TRUE), "High TMB", "Low TMB")
  )
message("TMB: median=", round(median(tmb_df$tmb, na.rm=TRUE), 1),
        " range=", min(tmb_df$tmb, na.rm=TRUE), "-", max(tmb_df$tmb, na.rm=TRUE))

# ── 6. Survival by mutation status (top driver genes) ────────────────────────
driver_genes <- head(top_genes_mut$Hugo_Symbol, 10)

mut_matrix_raw <- mutCountMatrix(maf_obj, removeNonMutated = FALSE)
mut_binary     <- (mut_matrix_raw > 0) * 1L

common_pat <- intersect(survival_df$bcr_patient_barcode, colnames(mut_binary))
surv_mut   <- survival_df |>
  filter(bcr_patient_barcode %in% common_pat) |>
  mutate(OS.time_m = OS.time / 30.44)

for (g in driver_genes) {
  if (g %in% rownames(mut_binary)) {
    surv_mut[[g]] <- as.integer(mut_binary[g, surv_mut$bcr_patient_barcode])
  } else {
    surv_mut[[g]] <- 0L
  }
}

# Fill patients not in MAF as WT (0)
for (g in driver_genes) {
  surv_mut[[g]][is.na(surv_mut[[g]])] <- 0L
}

# ── 7. Save figures ──────────────────────────────────────────────────────────
mut_colors <- c(
  "Missense_Mutation"      = "#2196F3",
  "Nonsense_Mutation"      = "#F44336",
  "Frame_Shift_Del"        = "#FF9800",
  "Frame_Shift_Ins"        = "#9C27B0",
  "Splice_Site"            = "#4CAF50",
  "In_Frame_Del"           = "#00BCD4",
  "In_Frame_Ins"           = "#8BC34A",
  "Translation_Start_Site" = "#FFC107",
  "Nonstop_Mutation"       = "#E91E63"
)

# Figure 1: mutation frequency
pct_labels <- top20 |>
  dplyr::select(Hugo_Symbol, pct) |>
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = rev(top20$Hugo_Symbol)),
         total_n = top20$MutatedSamples)

p_freq <- ggplot(type_long, aes(x = Hugo_Symbol, y = n, fill = mut_type)) +
  geom_col() +
  geom_text(data = pct_labels,
            aes(x = Hugo_Symbol, y = total_n + max(top20$MutatedSamples) * 0.02,
                label = paste0(pct, "%")),
            inherit.aes = FALSE, size = 3, hjust = 0, color = "#333") +
  coord_flip() +
  scale_fill_manual(values = mut_colors, na.value = "grey70") +
  labs(title    = "Top 20 Mutated Genes — TCGA-BRCA",
       subtitle = paste0("n = ", n_patients,
                         " patients | % = fraction with ≥1 somatic mutation"),
       x = NULL, y = "Number of mutations", fill = "Mutation type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        axis.text.y = element_text(face = "bold", size = 10))

png(file.path(FIGURES_DIR, "mutation_freq.png"),
    width = 1400, height = 900, res = 150)
print(p_freq)
dev.off()
message("Saved: mutation_freq.png")

# Figure 2: TMB by stage
p_tmb <- ggplot(tmb_df |> filter(!is.na(stage_clean)),
                aes(x = stage_clean, y = log10(tmb + 1), fill = stage_clean)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.25, size = 0.6) +
  scale_fill_manual(values = c("I"="#2196F3","II"="#4CAF50",
                               "III"="#FF9800","IV"="#F44336")) +
  ggpubr::stat_compare_means(label = "p.signif", ref.group = "I") +
  labs(title    = "Tumor Mutational Burden (TMB) by Tumor Stage",
       subtitle = "TMB = somatic non-synonymous mutations per patient | log10 scale",
       x = "Tumor Stage", y = "log10(TMB + 1)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

png(file.path(FIGURES_DIR, "mutation_tmb.png"),
    width = 900, height = 700, res = 150)
print(p_tmb)
dev.off()
message("Saved: mutation_tmb.png")

# Figure 3: OncoPrint (top 15 genes)
png(file.path(FIGURES_DIR, "oncoprint.png"),
    width = 1600, height = 800, res = 130)
tryCatch(
  oncoplot(maf = maf_obj, top = 15L,
           titleText = "TCGA-BRCA: OncoPrint — Top 15 Mutated Genes"),
  error = function(e) message("Oncoplot error: ", e$message)
)
dev.off()
message("Saved: oncoprint.png")

# Figure 4: Co-occurrence / mutual exclusivity
png(file.path(FIGURES_DIR, "mutation_cooccurrence.png"),
    width = 900, height = 800, res = 130)
tryCatch(
  somaticInteractions(maf = maf_obj, top = 15, fontSize = 0.7),
  error = function(e) message("somaticInteractions error: ", e$message)
)
dev.off()
message("Saved: mutation_cooccurrence.png")

# ── 8. Save RDS ──────────────────────────────────────────────────────────────
saveRDS(
  list(
    gene_freq    = as.data.frame(top_genes_mut),
    top20        = top20,
    type_long    = as.data.frame(type_long),
    tmb_df       = tmb_df,
    surv_mut     = surv_mut,
    driver_genes = driver_genes,
    n_patients   = n_patients,
    mut_colors   = mut_colors
  ),
  file.path(RESULTS_DIR, "mutation_results.rds")
)
message("Saved: mutation_results.rds")
message("\nMutation analysis complete. Restart Shiny app to see the new tab.")
