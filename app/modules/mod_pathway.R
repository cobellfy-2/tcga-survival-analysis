# mod_pathway.R — ORA and GSEA enrichment results

mod_pathway_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About Pathway Enrichment Analysis"),
      p("After identifying survival-associated genes, the next question is: ",
        strong("do these genes cluster in specific biological processes?"),
        " Two complementary methods answer this:"),
      tags$ul(
        tags$li(strong("ORA (Over-Representation Analysis):"),
          " Takes only the ", em("significant"), " survival genes (FDR < 5%) and asks whether any GO term
          appears more often than expected by chance (hypergeometric test). Simple and interpretable,
          but requires a ", strong("large enough gene list"), " — with too few genes, statistical power collapses."),
        tags$li(strong("GSEA (Gene Set Enrichment Analysis):"),
          " Ranks ", strong("all ~5,000 genes"), " by their survival association (log2HR × -log10 adj.p),
          then tests whether genes in a pathway cluster at the top (enriched in risk genes)
          or bottom (enriched in protective genes) of this ranked list. No hard threshold needed — uses the full signal.",
          " The ", strong("Normalized Enrichment Score (NES)"), ": NES > 0 = pathway genes tend toward
          risk (worse survival when highly expressed); NES < 0 = pathway genes tend toward protection.")
      )
    ),
    div(class = "alert alert-success",
      h5("Key findings — GSEA (687 pathways tested)"),
      tags$ul(
        tags$li(strong("Homophilic cell-cell adhesion (NES = 2.07, padj = 0.004):"),
          " The strongest significant hit. E-cadherin and related adhesion molecules are
          risk-enriched — tumors with high adhesion gene expression have worse survival.
          This reflects epithelial-mesenchymal transition (EMT) dynamics in aggressive BRCA."),
        tags$li(strong("Lipid metabolism cluster:"),
          " Response to fatty acids (NES = 2.10), triglyceride metabolism (NES = 1.97),
          and neutral lipid metabolism (NES = 1.95) are all enriched in risk genes.
          Reprogrammed lipid metabolism is a hallmark of aggressive breast cancer — tumor cells
          upregulate lipid synthesis to fuel rapid proliferation."),
        tags$li(strong("Steroid hormone metabolism (NES = 1.97):"),
          " C21-steroid hormone processing (glucocorticoids, progestins) enriched in risk genes,
          consistent with hormone-receptor biology driving aggressive BRCA subtypes."),
        tags$li(strong("ORA: 0 significant terms"),
          " — expected with only 52 Lasso-selected genes. GSEA is the more informative method here.")
      )
    ),
    tabsetPanel(
      tabPanel("ORA — Significant genes",
        br(),
        div(class = "alert alert-warning",
          h5("Why ORA found no significant pathways"),
          p(strong("ORA requires a large gene list."), " The hypergeometric test asks:
            if I randomly draw N genes from the genome, what is the probability that
            k or more fall in pathway P? With only ", strong("52 Lasso-selected genes"),
            " as input, statistical power is very low — even a pathway with 10 of its 50 genes
            present in our list would barely reach significance after multiple-testing correction."),
          p(strong("This is not a failure — it is a known limitation of ORA."), " The Lasso selects
            a minimal set of independently predictive genes, deliberately avoiding redundancy.
            The trade-off is a small final list. GSEA, which uses all 5,000 genes ranked by
            survival signal, is far more powerful in this scenario."),
          p("See the ", strong("GSEA tab"), " for 687 pathways analyzed with the full gene ranking —
            including 5 pathways significant at FDR < 5%.")
        ),
        plotOutput(ns("ora_plot"), height = "300px")
      ),
      tabPanel("GSEA — All genes ranked",
        br(),
        div(class = "alert alert-secondary",
          h6("How to read this GSEA dot plot"),
          tags$ul(style = "margin-bottom:0",
            tags$li(strong("Each dot = one GO Biological Process pathway.")),
            tags$li(strong("X-axis (NES):"),
              " Normalized Enrichment Score. Positive NES = pathway genes cluster among ",
              em("risk"), " genes (high expression → worse survival). Negative NES = protective."),
            tags$li(strong("Dot size = gene ratio:"),
              " fraction of the pathway's genes that appear in our ranked gene list."),
            tags$li(strong("Color = adjusted p-value:"),
              " darker/more saturated = more significant. Only dots with padj < 0.2 shown."),
            tags$li(strong("Bottom panel (enrichment plot):"),
              " shows the running enrichment score for the top-ranked pathway.
              The peak indicates where pathway genes are concentrated in the full ranked list.")
          )
        ),
        fluidRow(
          column(8, plotOutput(ns("gsea_plot"), height = "520px")),
          column(4,
            br(),
            h6("Top pathway — enrichment profile"),
            p(style="color:#555; font-size:0.88em;",
              "Running sum of enrichment score across all ranked genes.
               The peak = where pathway members cluster. A high peak far to the left
               means pathway genes dominate the risk (high-hazard) end of the ranking."),
            plotOutput(ns("gsea_enrich"), height = "360px")
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
        ggplot() +
          annotate("text", x = 0.5, y = 0.55,
            label = "0 significant ORA terms",
            size = 7, fontface = "bold", color = "#555") +
          annotate("text", x = 0.5, y = 0.42,
            label = "Input: 52 Lasso-selected genes  |  Threshold: FDR < 5%\nSee explanation above and the GSEA tab for full results.",
            size = 4.5, color = "#777", hjust = 0.5) +
          xlim(0, 1) + ylim(0, 1) + theme_void()
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
        tryCatch({
          g_df <- gsea@result
          has_sign <- ".sign" %in% colnames(g_df) && length(unique(g_df$.sign)) > 1
          if (has_sign) {
            enrichplot::dotplot(gsea, showCategory = 15, split = ".sign",
                    title = "GSEA: GO Biological Process (top 15 per direction)") +
              facet_grid(. ~ .sign) +
              theme(axis.text.y = element_text(size = 8))
          } else {
            enrichplot::dotplot(gsea, showCategory = 20,
                    title = "GSEA: GO Biological Process (ranked by NES)") +
              theme(axis.text.y = element_text(size = 9))
          }
        }, error = function(e) {
          ggplot() + annotate("text", x=0.5, y=0.5,
            label=paste("Plot error:", e$message), size=4, hjust=0.5) + theme_void()
        })
      }
    })

    output$gsea_enrich <- renderPlot({
      gsea <- enrich_res$gsea_go
      if (is.null(gsea) || nrow(gsea) == 0) return(NULL)
      tryCatch(
        enrichplot::gseaplot2(gsea,
          geneSetID = gsea@result$ID[1],
          title     = gsea@result$Description[1]),
        error = function(e) {
          ggplot() + annotate("text", x=0.5, y=0.5,
            label=paste("Enrichment plot error:", e$message), size=4) + theme_void()
        }
      )
    })
  })
}
