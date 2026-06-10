# 11_mutations.R
# Somatic mutation analysis: OncoPrint + survival impact of top mutated genes.
#
# BACKGROUND:
# RNA-seq tells us WHAT is being expressed — but not WHY.
# Somatic mutations are the drivers: specific DNA changes that cause
# cancer or alter its behavior. TCGA provides Mutation Annotation Format
# (MAF) files with every somatic mutation per patient.
#
# OncoPrint: compact visual showing which patients have mutations in
# which genes. A standard figure in cancer genomics papers.
#
# We then ask: do patients with mutations in top genes survive differently?
# This combines mutational and survival data — true multi-omics.
#
# Input:  TCGA-BRCA MAF file (downloaded via TCGAbiolinks)
#         data/processed/survival_df.rds
# Output: data/processed/mutation_matrix.rds
#         results/figures/oncoprint.png
#         results/figures/km_mutations.png

source("config.R")

for (pkg in c("TCGAbiolinks", "maftools", "survival",
              "survminer", "ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("TCGAbiolinks", "maftools")) BiocManager::install(pkg)
    else install.packages(pkg)
  }
}

library(TCGAbiolinks)
library(maftools)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))

# --- Download MAF (somatic mutations) ---
message("Downloading somatic mutation data...")
maf_query <- GDCquery(
  project           = paste0("TCGA-", CANCER_TYPE),
  data.category     = "Simple Nucleotide Variation",
  data.type         = "Masked Somatic Mutation",
  workflow.type     = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query, directory = DATA_DIR)
maf_df <- GDCprepare(maf_query, directory = DATA_DIR)

saveRDS(maf_df, file.path(PROCESSED_DIR, "maf_raw.rds"))
message("MAF downloaded: ", nrow(maf_df), " mutations")

# --- Load into maftools ---
# Need clinical data in maftools format
clin_maf <- survival_df |>
  mutate(
    Tumor_Sample_Barcode = bcr_patient_barcode,
    OS_days  = OS.time,
    OS_event = OS
  ) |>
  as.data.frame()

maf_obj <- read.maf(
  maf            = maf_df,
  clinicalData   = clin_maf,
  verbose        = FALSE
)

# --- Summary stats ---
message("\nTop 20 mutated genes:")
top_genes_mut <- getGeneSummary(maf_obj)
print(head(top_genes_mut[, c("Hugo_Symbol","MutatedSamples","AlteredSamples")], 20))

# --- OncoPrint (top 15 genes) ---
top15 <- head(top_genes_mut$Hugo_Symbol, 15)

png(file.path(FIGURES_DIR, "oncoprint.png"),
    width = 1400, height = 700, res = 120)
oncoplot(
  maf         = maf_obj,
  top         = 15,
  clinicalFeatures = c("ajcc_pathologic_tumor_stage"),
  sortByAnnotation = TRUE,
  title       = "TCGA-BRCA: OncoPrint — Top 15 Mutated Genes"
)
dev.off()
message("Saved: oncoprint.png")

# --- Mutation co-occurrence / exclusivity ---
png(file.path(FIGURES_DIR, "mutation_cooccurrence.png"),
    width = 800, height = 700, res = 120)
somaticInteractions(
  maf   = maf_obj,
  top   = 15,
  pvalue = c(0.05, 0.01),
  fontSize = 0.7
)
dev.off()
message("Saved: mutation_cooccurrence.png")

# --- KM: survival by mutation status of top 3 genes ---
mut_matrix <- mutCountMatrix(maf_obj, removeNonMutated = FALSE)
mut_binary <- (mut_matrix > 0) * 1

# Align with survival data
common_pat <- intersect(survival_df$bcr_patient_barcode, colnames(mut_binary))
surv_mut   <- survival_df[survival_df$bcr_patient_barcode %in% common_pat, ]
mut_sub    <- mut_binary[top15[1:3], common_pat, drop = FALSE]

saveRDS(list(mut_matrix = mut_binary, top_genes = top_genes_mut),
        file.path(PROCESSED_DIR, "mutation_matrix.rds"))

km_plots <- lapply(top15[1:3], function(gene) {
  if (!gene %in% rownames(mut_sub)) return(NULL)
  df_g <- data.frame(
    OS.time = surv_mut$OS.time / 30.44,
    OS      = surv_mut$OS,
    mutated = factor(ifelse(mut_sub[gene, surv_mut$bcr_patient_barcode] == 1,
                            "Mutated", "Wild-type"))
  )
  n_mut <- sum(df_g$mutated == "Mutated")
  if (n_mut < 5) return(NULL)

  km <- survfit(Surv(OS.time, OS) ~ mutated, data = df_g)
  ggsurvplot(km, data = df_g,
             palette      = c("#F44336", "#2196F3"),
             pval         = TRUE, conf.int = TRUE,
             risk.table   = FALSE,
             title        = paste0(gene, " (n mut=", n_mut, ")"),
             xlab         = "Time (months)",
             legend.title = "Status",
             ggtheme      = theme_bw(base_size = 11))$plot
})

km_plots <- Filter(Negate(is.null), km_plots)
if (length(km_plots) > 0) {
  library(cowplot)
  combined <- plot_grid(plotlist = km_plots, ncol = length(km_plots))
  ggsave(file.path(FIGURES_DIR, "km_mutations.png"),
         plot = combined,
         width = length(km_plots) * 6, height = 5, dpi = 150)
  message("Saved: km_mutations.png")
}

message("\nMutation analysis complete. Next: update Shiny app tabs.")
