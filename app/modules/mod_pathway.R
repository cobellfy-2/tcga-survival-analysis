# mod_pathway.R — ORA and GSEA enrichment results

mod_pathway_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method"),
      p(strong("ORA (Over-Representation Analysis):"),
        " Takes the significant survival genes (FDR < 5%) and tests whether any GO Biological Process terms ",
        "are more common in this list than expected by chance (hypergeometric test). ",
        "Simple, but requires a hard p-value cutoff."),
      p(strong("GSEA (Gene Set Enrichment Analysis):"),
        " Ranks ALL genes by their survival association (sign × -log10p), then checks whether ",
        "genes in a pathway cluster at the top (risk genes) or bottom (protective genes) of the ranking. ",
        "No threshold needed — uses the full signal. The ", strong("Normalized Enrichment Score (NES)"),
        " > 0 means the pathway is enriched in risk genes; NES < 0 means enriched in protective genes.")
    ),
    tabsetPanel(
      tabPanel("ORA — Significant genes",
        br(),
        plotOutput(ns("ora_plot"), height = "550px")
      ),
      tabPanel("GSEA — All genes ranked",
        br(),
        fluidRow(
          column(8, plotOutput(ns("gsea_plot"), height = "550px")),
          column(4,
            br(),
            h6("Top enriched pathway"),
            plotOutput(ns("gsea_enrich"), height = "350px")
          )
        )
      )
    )
  )
}

mod_pathway_server <- function(id, enrich_res) {
  moduleServer(id, function(input, output, session) {
    requireNamespace("enrichplot", quietly = TRUE)
    dotplot <- enrichplot::dotplot

    output$ora_plot <- renderPlot({
      ora <- enrich_res$ora_go
      if (is.null(ora) || nrow(ora) == 0) {
        ggplot() + annotate("text", x=0.5, y=0.5,
          label="No significant ORA terms found.\nThis is a valid result — the gene list may be\ntoo small or the signal too diffuse for ORA.\nCheck the GSEA tab for full-ranked results.",
          size=5, hjust=0.5) + theme_void()
      } else {
        enrichplot::dotplot(ora, showCategory = 20,
                title = "ORA: GO Biological Process (FDR < 5% survival genes)") +
          theme(axis.text.y = element_text(size = 9))
      }
    })

    output$gsea_plot <- renderPlot({
      gsea <- enrich_res$gsea_go
      if (is.null(gsea) || nrow(gsea) == 0) {
        ggplot() + annotate("text", x=0.5, y=0.5,
          label="No significant GSEA pathways found.", size=6) + theme_void()
      } else {
        enrichplot::dotplot(gsea, showCategory = 20, split = ".sign",
                title = "GSEA: GO Biological Process") +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 8))
      }
    })

    output$gsea_enrich <- renderPlot({
      gsea <- enrich_res$gsea_go
      if (is.null(gsea) || nrow(gsea) == 0) return(NULL)
      enrichplot::gseaplot2(gsea,
                            geneSetID = gsea@result$ID[1],
                            title     = gsea@result$Description[1])
    })
  })
}
