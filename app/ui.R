# ui.R — overall layout and tab structure

library(shiny)
library(bslib)
library(plotly)

# Load all modules
for (f in list.files("modules", full.names = TRUE)) source(f)

ui <- page_navbar(
  title          = "TCGA-BRCA Survival Analysis",
  theme          = bs_theme(bootswatch = "flatly", primary = "#2196F3"),
  navbar_options = navbar_options(bg = "#1a1a2e"),

  # ── Tab 1: Overview ──────────────────────────────────────────────────────
  nav_panel("📊 Overview",      mod_overview_ui("overview")),

  # ── Tab 2: Clinical Survival ─────────────────────────────────────────────
  nav_panel("🏥 Clinical",      mod_clinical_ui("clinical")),

  # ── Tab 3: Gene Explorer ─────────────────────────────────────────────────
  nav_panel("🧬 Gene Explorer", mod_gene_ui("gene")),

  # ── Tab 4: Lasso Signature ───────────────────────────────────────────────
  nav_panel("🎯 Lasso",         mod_lasso_ui("lasso")),

  # ── Tab 5: Pathway Enrichment ────────────────────────────────────────────
  nav_panel("🔬 Pathways",      mod_pathway_ui("pathway")),

  # ── Tab 6: Immune Landscape ──────────────────────────────────────────────
  nav_panel("🛡️ Immune",        mod_immune_ui("immune")),

  nav_spacer(),
  nav_item(tags$small(style = "color:#aaa; padding:10px;",
                      "Data: TCGA-BRCA | Tools: R + Shiny"))
)
