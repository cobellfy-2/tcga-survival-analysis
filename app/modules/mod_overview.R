# mod_overview.R — landing page with project summary and data stats

mod_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(12,
        h2("TCGA-BRCA Survival Analysis Dashboard"),
        p(style = "font-size:1.1em; color:#555;",
          "This dashboard presents a reproducible survival analysis of ",
          strong("1,035 breast cancer patients"), " from The Cancer Genome Atlas (TCGA-BRCA). ",
          "It combines clinical variables, RNA-seq gene expression, ",
          "pathway biology, and immune infiltration to identify predictors of overall survival.")
      )
    ),
    hr(),
    # Key stats row
    fluidRow(
      column(3, div(class = "card text-center p-3 shadow-sm",
        h1(textOutput(ns("n_patients")), style = "color:#2196F3; font-weight:bold;"),
        p("Patients"))),
      column(3, div(class = "card text-center p-3 shadow-sm",
        h1(textOutput(ns("n_events")), style = "color:#F44336; font-weight:bold;"),
        p("Deaths (events)"))),
      column(3, div(class = "card text-center p-3 shadow-sm",
        h1(textOutput(ns("n_genes")), style = "color:#4CAF50; font-weight:bold;"),
        p("Genes analysed"))),
      column(3, div(class = "card text-center p-3 shadow-sm",
        h1(textOutput(ns("median_fu")), style = "color:#FF9800; font-weight:bold;"),
        p("Median follow-up (months)")))
    ),
    br(),
    # Pipeline overview
    fluidRow(
      column(6,
        h4("Analysis Pipeline"),
        tags$ol(
          tags$li(strong("Clinical Survival:"), " Kaplan-Meier curves, log-rank test, Cox regression by stage, age, subtype"),
          tags$li(strong("Gene Explorer:"), " Search any gene — instant KM curve + hazard ratio"),
          tags$li(strong("Lasso Signature:"), " Penalized Cox regression selects a sparse prognostic gene signature"),
          tags$li(strong("Pathway Enrichment:"), " ORA + GSEA identifies biological processes enriched in survival genes"),
          tags$li(strong("Immune Landscape:"), " ESTIMATE scores quantify immune infiltration and correlate with survival")
        )
      ),
      column(6,
        h4("Dataset"),
        tags$ul(
          tags$li(strong("Source:"), " TCGA-BRCA via GDC API (TCGAbiolinks)"),
          tags$li(strong("Molecular:"), " RNA-seq HTSeq counts, STAR workflow"),
          tags$li(strong("Normalization:"), " VST (DESeq2) for gene analysis; log2(CPM+1) for Lasso + ESTIMATE"),
          tags$li(strong("Endpoint:"), " Overall Survival (OS) — time from diagnosis to death or last contact"),
          tags$li(strong("Subtypes:"), " ER/PR/HER2 receptor status → Luminal A/B, HER2-enriched, Triple-Negative")
        )
      )
    ),
    hr(),
    fluidRow(
      column(12,
        h4("How to use this dashboard"),
        p("Each tab focuses on one analysis method. Every tab contains an ",
          strong("'About this method'"), " box that explains ",
          "what the method does, why it is used, and how to interpret the results. ",
          "Use the controls on each tab to explore the data interactively.")
      )
    )
  )
}

mod_overview_server <- function(id, survival_df, survival_sub) {
  moduleServer(id, function(input, output, session) {
    output$n_patients <- renderText(nrow(survival_df))
    output$n_events   <- renderText(sum(survival_df$OS))
    output$n_genes    <- renderText("23,735")
    output$median_fu  <- renderText(
      round(median(survival_df$OS.time / 30.44, na.rm = TRUE), 0))
  })
}
