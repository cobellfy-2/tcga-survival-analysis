# ui.R — overall layout and tab structure
library(shiny)
library(bslib)
library(plotly)

for (f in list.files("modules", full.names = TRUE)) source(f)

ui <- page_navbar(
  title          = "TCGA-BRCA Survival Analysis",
  theme          = bs_theme(bootswatch = "flatly", primary = "#2196F3"),
  navbar_options = navbar_options(bg = "#1a1a2e"),

  nav_panel("📊 Overview",      mod_overview_ui("overview")),
  nav_panel("🏥 Clinical",      mod_clinical_ui("clinical")),
  nav_panel("🧬 Gene Explorer", mod_gene_ui("gene")),
  nav_panel("🎯 Lasso",         mod_lasso_ui("lasso")),
  nav_panel("✅ Validation",    mod_validation_ui("validation")),
  nav_panel("🔬 Pathways",      mod_pathway_ui("pathway")),
  nav_panel("🛡️ Immune",        mod_immune_ui("immune")),
  nav_panel("🧪 Mutations",      mod_mutations_ui("mutations")),

  nav_spacer(),
  nav_item(tags$small(style = "color:#aaa; padding:10px;",
                      "TCGA-BRCA | R + Python + Shiny"))
)
