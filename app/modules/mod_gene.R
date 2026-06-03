# mod_gene.R — search any gene, get instant KM + HR

mod_gene_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method"),
      p("For each gene, we ran a ", strong("univariate Cox regression"),
        " where the predictor is the normalized expression level. ",
        "The ", strong("Hazard Ratio (HR)"), " tells you: if expression increases by 1 unit, ",
        "the death hazard changes by HR-fold. ",
        "Patients are split into high/low expression at the ", strong("median"),
        " to visualize the KM curve. ",
        "P-values are adjusted for multiple testing using ",
        strong("Benjamini-Hochberg (FDR)"), " correction across 5,000 genes.")
    ),
    fluidRow(
      column(4,
        wellPanel(
          h5("Search Gene"),
          selectizeInput(ns("gene_search"), "Gene symbol:",
                         choices  = NULL,
                         options  = list(placeholder  = "e.g. SUSD3, LEF1, TP53",
                                         create       = FALSE)),
          hr(),
          uiOutput(ns("gene_stats_box"))
        )
      ),
      column(8,
        plotOutput(ns("km_gene"), height = "420px"),
        br(),
        div(class = "alert alert-secondary",
          h6("Volcano plot — all 5,000 genes"),
          p("Red = high expression worsens survival | Blue = protective | Hover for gene name"),
          plotlyOutput(ns("volcano"), height = "320px")
        )
      )
    )
  )
}

mod_gene_server <- function(id, survival_df, counts_filt, cox_genes) {
  moduleServer(id, function(input, output, session) {

    # Populate gene dropdown with symbols that have valid data
    valid_genes <- cox_genes |>
      filter(!is.na(symbol), symbol != "") |>
      arrange(padj)

    updateSelectizeInput(session, "gene_search",
                         choices = setNames(valid_genes$gene, valid_genes$symbol),
                         server  = TRUE)

    # KM plot for selected gene
    output$km_gene <- renderPlot({
      req(input$gene_search)
      gene_id <- input$gene_search
      sym     <- cox_genes$symbol[cox_genes$gene == gene_id][1]
      if (is.na(sym) || sym == "") sym <- gene_id

      # Get expression
      common <- intersect(survival_df$bcr_patient_barcode, colnames(counts_filt))
      surv   <- survival_df[match(common, survival_df$bcr_patient_barcode), ]
      expr   <- log2(counts_filt[gene_id, common] + 1)

      df_g <- data.frame(
        OS.time = surv$OS.time / 30.44,
        OS      = surv$OS,
        expr    = as.numeric(expr),
        group   = ifelse(as.numeric(expr) >= median(as.numeric(expr)),
                         "High expression", "Low expression")
      )

      km  <- survfit(Surv(OS.time, OS) ~ group, data = df_g)
      row <- cox_genes[cox_genes$gene == gene_id, ]
      subtitle <- if (nrow(row) > 0)
        paste0("HR = ", round(row$hr[1], 3),
               " | adj.p = ", formatC(row$padj[1], digits = 3, format = "e"))
      else ""

      p <- ggsurvplot(km, data = df_g,
                      palette      = c("#F44336", "#2196F3"),
                      conf.int     = TRUE, pval = TRUE,
                      risk.table   = TRUE,
                      xlab         = "Time (months)",
                      title        = paste0(sym, " expression and Overall Survival"),
                      subtitle     = subtitle,
                      legend.title = "Expression",
                      ggtheme      = theme_bw(base_size = 13))
      gridExtra::grid.arrange(p$plot, p$table, ncol = 1, heights = c(3, 1))
    })

    # Stats box
    output$gene_stats_box <- renderUI({
      req(input$gene_search)
      row <- cox_genes[cox_genes$gene == input$gene_search, ]
      if (nrow(row) == 0) return(NULL)
      div(
        h6(row$symbol[1]),
        tags$table(class = "table table-sm",
          tags$tr(tags$td("Hazard Ratio"),  tags$td(strong(round(row$hr[1], 3)))),
          tags$tr(tags$td("95% CI"),
                  tags$td(paste0(round(row$lower[1],3), " – ", round(row$upper[1],3)))),
          tags$tr(tags$td("p-value"),        tags$td(formatC(row$pval[1], digits=3, format="e"))),
          tags$tr(tags$td("adj. p-value"),   tags$td(formatC(row$padj[1], digits=3, format="e"))),
          tags$tr(tags$td("Direction"),
                  tags$td(if (row$hr[1] > 1)
                    span(style="color:#F44336", "High = worse survival")
                  else
                    span(style="color:#2196F3", "High = better survival")))
        )
      )
    })

    # Volcano
    output$volcano <- renderPlotly({
      df_v <- cox_genes |>
        filter(!is.na(log2HR), !is.na(neg_log10p)) |>
        mutate(
          label = ifelse(is.na(symbol) | symbol == "", gene, symbol),
          sig   = case_when(
            padj < 0.05 & log2HR > 0 ~ "High = worse OS",
            padj < 0.05 & log2HR < 0 ~ "High = better OS",
            TRUE                      ~ "Not significant"
          ),
          color = case_when(
            sig == "High = worse OS"   ~ "#F44336",
            sig == "High = better OS"  ~ "#2196F3",
            TRUE                       ~ "#BDBDBD"
          )
        )

      plot_ly(df_v,
              x    = ~log2HR, y = ~neg_log10p,
              type = "scatter", mode = "markers",
              text = ~paste0("<b>", label, "</b><br>HR: ", round(hr,3),
                             "<br>adj.p: ", formatC(padj, digits=2, format="e")),
              hoverinfo = "text",
              marker = list(color = ~color, size = 5, opacity = 0.6)) |>
        layout(
          xaxis = list(title = "log2 Hazard Ratio",
                       zerolinecolor = "#888", zerolinewidth = 1),
          yaxis = list(title = "-log10(adj. p-value)"),
          shapes = list(
            list(type="line", x0=0, x1=0,
                 y0=0, y1=max(df_v$neg_log10p, na.rm=TRUE),
                 line=list(dash="dot", color="grey")),
            list(type="line", x0=min(df_v$log2HR, na.rm=TRUE),
                 x1=max(df_v$log2HR, na.rm=TRUE),
                 y0=-log10(0.05), y1=-log10(0.05),
                 line=list(dash="dot", color="grey"))
          ),
          paper_bgcolor = "white", plot_bgcolor = "white",
          showlegend = FALSE
        )
    })
  })
}
