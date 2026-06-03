# mod_lasso.R — Lasso Cox signature and risk score

mod_lasso_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method"),
      p(strong("Lasso-penalized Cox regression"), " solves the problem of having far more genes than patients. ",
        "It adds an L1 penalty (λ × Σ|β|) to the Cox likelihood, which shrinks weak gene coefficients to exactly zero. ",
        "Only genes with real prognostic signal survive. ",
        strong("Cross-validation (10-fold)"), " tests ~100 λ values and picks the one that best predicts survival on held-out data. ",
        "We use ", strong("lambda.1se"), " — the most parsimonious model within 1 standard error of the optimum. ",
        "The resulting ", strong("risk score"), " (= weighted sum of selected gene expressions) ",
        "stratifies patients into high/low risk groups.")
    ),
    fluidRow(
      column(5,
        h5("Selected gene signature"),
        p(textOutput(ns("n_genes_selected")), style = "color:#555;"),
        plotOutput(ns("sig_plot"), height = "450px")
      ),
      column(7,
        h5("Risk score distribution"),
        plotOutput(ns("risk_hist"), height = "200px"),
        br(),
        h5("Kaplan-Meier: High vs. Low risk"),
        plotOutput(ns("km_lasso"), height = "300px")
      )
    )
  )
}

mod_lasso_server <- function(id, lasso_res) {
  moduleServer(id, function(input, output, session) {

    sig_df  <- lasso_res$sig_df
    risk_df <- lasso_res$risk_df

    output$n_genes_selected <- renderText({
      paste0(nrow(sig_df), " genes selected by Lasso (lambda.1se = ",
             round(lasso_res$lambda, 4), ")")
    })

    output$sig_plot <- renderPlot({
      df_plot <- sig_df |>
        slice_max(abs(coef), n = min(30, nrow(sig_df))) |>
        mutate(symbol = factor(symbol, levels = symbol[order(coef)]))

      ggplot(df_plot, aes(x = coef, y = symbol, fill = direction)) +
        geom_col() +
        scale_fill_manual(values = c(
          "Risk (high = worse)"        = "#F44336",
          "Protective (high = better)" = "#2196F3")) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#333") +
        labs(x = "Lasso coefficient", y = NULL, fill = NULL,
             title = "Prognostic gene signature") +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")
    })

    output$risk_hist <- renderPlot({
      ggplot(risk_df, aes(x = risk_score, fill = risk_group)) +
        geom_histogram(bins = 50, alpha = 0.8, color = "white") +
        scale_fill_manual(values = c("High risk" = "#F44336", "Low risk" = "#2196F3")) +
        geom_vline(xintercept = median(risk_df$risk_score),
                   linetype = "dashed", color = "#333") +
        labs(title = "Risk score distribution (split at median)",
             x = "Risk score", y = "Count", fill = NULL) +
        theme_bw(base_size = 12)
    })

    output$km_lasso <- renderPlot({
      km <- survfit(Surv(OS.time, OS) ~ risk_group, data = risk_df)
      p  <- ggsurvplot(km, data = risk_df,
                       palette      = c("#F44336", "#2196F3"),
                       conf.int     = TRUE, pval = TRUE,
                       risk.table   = FALSE,
                       xlab         = "Time (months)",
                       title        = "Lasso risk score: Overall Survival",
                       legend.title = "Risk group",
                       ggtheme      = theme_bw(base_size = 12))
      print(p$plot)
    })
  })
}
