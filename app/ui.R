# ui.R — overall layout and tab structure
library(shiny)
library(bslib)
library(plotly)

for (f in list.files("modules", full.names = TRUE)) source(f)

# ── Girly-pop theme ──────────────────────────────────────────────────────────
girly_theme <- bs_theme(
  version      = 5,
  bg           = "#fff7fb",
  fg           = "#3d2c3a",
  primary      = "#ff4fa3",
  secondary    = "#b983ff",
  success      = "#ff85c0",
  info         = "#ffa6d2",
  base_font    = font_google("Quicksand"),
  heading_font = font_google("Poppins"),
  "border-radius"      = "1rem",
  "navbar-bg"          = "#ff4fa3",
  "navbar-dark-color"  = "#ffffff"
)

# Custom pastel polish for cards, alert boxes and headings
girly_css <- tags$head(tags$style(HTML("
  body { background: #fff7fb; }
  h1,h2,h3,h4,h5,h6 { color: #d6357f; font-weight: 600; }
  .nav-link { font-weight: 600 !important; }
  .card { border-radius: 18px; border: none;
          box-shadow: 0 4px 16px rgba(255,79,163,.14); }
  .well, .wellPanel { border-radius: 16px; background: #fdeef7;
                      border: 1px solid #ffd2e8; }
  .alert { border-radius: 16px; border-width: 1px; }
  .alert-info      { background:#fff0f7; border-color:#ffc2dd; color:#7a2b52; }
  .alert-success   { background:#fdeafd; border-color:#f3b6f0; color:#6a2a66; }
  .alert-secondary { background:#f5f0ff; border-color:#d9c7ff; color:#4a3b6b; }
  .alert-warning   { background:#fff6e9; border-color:#ffd9a0; color:#7a5320; }
  .alert-light     { background:#fffafd; border-color:#ffe0ef; }
  .nav-tabs .nav-link.active { color:#d6357f; border-bottom:3px solid #ff4fa3; }
  .navbar .nav-link, .navbar-brand { color:#fff !important; }
  .btn-primary { background:#ff4fa3; border-color:#ff4fa3; }
  table.table-sm td { padding:.3rem .5rem; }
  .fa, .fas, .far { margin-right:.35rem; }
")))

# Helper: pictogram + label for nav titles
nav_title <- function(icon_name, label) {
  tagList(icon(icon_name), label)
}

ui <- page_navbar(
  title  = tagList(icon("ribbon"), "BRCA Survival Studio"),
  theme  = girly_theme,
  header = girly_css,

  nav_panel(nav_title("chart-pie",            "Overview"),      mod_overview_ui("overview")),
  nav_panel(nav_title("heart-pulse",          "Clinical"),      mod_clinical_ui("clinical")),
  nav_panel(nav_title("dna",                  "Gene Explorer"), mod_gene_ui("gene")),
  nav_panel(nav_title("wand-magic-sparkles",  "Lasso"),         mod_lasso_ui("lasso")),
  nav_panel(nav_title("circle-check",         "Validation"),    mod_validation_ui("validation")),
  nav_panel(nav_title("diagram-project",      "Pathways"),      mod_pathway_ui("pathway")),
  nav_panel(nav_title("shield-heart",         "Immune"),        mod_immune_ui("immune")),
  nav_panel(nav_title("vial",                 "Mutations"),     mod_mutations_ui("mutations")),

  nav_spacer(),
  nav_item(tags$small(style = "color:#ffe3f1; padding:10px;",
                      "TCGA-BRCA · R + Python + Shiny"))
)
