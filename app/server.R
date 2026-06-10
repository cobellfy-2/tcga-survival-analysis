# server.R — loads all data once, passes to modules

library(shiny)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(plotly)
library(enrichplot)
library(clusterProfiler)

# ── Load pre-computed results (once at startup) ──────────────────────────────
message("Loading data...")
survival_df   <- readRDS("../data/processed/survival_df.rds")
survival_sub  <- readRDS("../data/processed/survival_subtypes.rds")

immune_df_raw <- readRDS("../results/immune_scores.rds")
cox_genes     <- readRDS("../results/gene_cox_results.rds")
lasso_res     <- readRDS("../results/lasso_results.rds")
enrich_res    <- readRDS("../results/enrichment_results.rds")
subtype_lookup <- setNames(survival_sub$subtype, survival_sub$bcr_patient_barcode)
immune_df <- immune_df_raw
immune_df$subtype <- subtype_lookup[immune_df$patient]
counts_filt   <- readRDS("../data/processed/counts_filtered.rds")
message("Data loaded.")

server <- function(input, output, session) {
  mod_overview_server("overview", survival_df, survival_sub)
  mod_clinical_server("clinical", survival_df, survival_sub)
  mod_gene_server("gene",         survival_df, counts_filt, cox_genes)
  mod_lasso_server("lasso",       lasso_res)
  mod_pathway_server("pathway",   enrich_res)
  mod_immune_server("immune",     immune_df)
}
