# server.R — loads all data once, passes to modules

library(shiny)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(plotly)

# ── Load pre-computed results (once at startup) ──────────────────────────────
message("Loading data...")
survival_df   <- readRDS("../data/processed/survival_df.rds")
survival_sub  <- readRDS("../data/processed/survival_subtypes.rds")
cox_genes     <- readRDS("../results/gene_cox_results.rds")
lasso_res     <- readRDS("../results/lasso_results.rds")
enrich_res    <- readRDS("../results/enrichment_results.rds")
immune_df     <- readRDS("../results/immune_scores.rds")
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
