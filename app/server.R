# server.R — loads all data once at startup, passes to modules
library(shiny)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(plotly)
library(enrichplot)
library(clusterProfiler)
library(tidyr)

message("Loading data...")
survival_df  <- readRDS("../data/processed/survival_df.rds")
survival_sub <- readRDS("../data/processed/survival_subtypes.rds")
cox_genes    <- readRDS("../results/gene_cox_results.rds")
lasso_res    <- readRDS("../results/lasso_results.rds")
enrich_res   <- readRDS("../results/enrichment_results.rds")
immune_raw   <- readRDS("../results/immune_scores.rds")
val_res      <- readRDS("../results/validation_results.rds")
counts_filt  <- readRDS("../data/processed/counts_filtered.rds")
mut_res      <- tryCatch(
  readRDS("../results/mutation_results.rds"),
  error = function(e) { message("mutation_results.rds not yet available"); NULL }
)

# Add subtype to immune data
subtype_map  <- setNames(survival_sub$subtype, survival_sub$bcr_patient_barcode)
immune_df    <- immune_raw
immune_df$subtype <- subtype_map[immune_df$patient]

FIGURES_DIR  <- "../results/figures"
message("Data loaded.")

server <- function(input, output, session) {
  mod_overview_server("overview",    survival_df, survival_sub)
  mod_clinical_server("clinical",    survival_df, survival_sub)
  mod_gene_server("gene",            survival_df, counts_filt, cox_genes)
  mod_lasso_server("lasso",          lasso_res)
  mod_validation_server("validation", val_res, FIGURES_DIR)
  mod_pathway_server("pathway",      enrich_res)
  mod_immune_server("immune",        immune_df)
  mod_mutations_server("mutations",  mut_res, FIGURES_DIR)
}
