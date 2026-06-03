# mod_immune.R — immune landscape and survival

mod_immune_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method"),
      p(strong("ESTIMATE"), " (Estimation of STromal and Immune cells in MAlignant Tumors using Expression data) ",
        "uses known gene signatures of immune and stromal cells to deconvolute bulk RNA-seq. ",
        "Because immune cells always express certain marker genes (e.g. CD3, CD8 for T cells), ",
        "their relative abundance can be inferred from the overall expression profile. ",
        strong("ImmuneScore"), " reflects immune cell infiltration; ",
        strong("TumorPurity"), " is the estimated fraction of actual cancer cells. ",
        "A 'hot' tumor (high ImmuneScore) often responds better to immunotherapy. ",
        "In BRCA, the relationship with survival is subtype-dependent.")
    ),
    fluidRow(
      column(4,
        wellPanel(
          h5("Controls"),
          selectInput(ns("score_type"), "Score to display:",
                      choices = c("ImmuneScore", "StromalScore",
                                  "ESTIMATEScore", "TumorPurity"),
                      selected = "ImmuneScore"),
          hr(),
          radioButtons(ns("group_by"), "Stratify by:",
                       choices = c("Stage" = "stage", "Subtype" = "subtype"),
                       selected = "stage")
        )
      ),
      column(8,
        plotOutput(ns("box_plot"), height = "320px")
      )
    ),
    fluidRow(
      column(6,
        h5("Immune Score vs. Tumor Purity"),
        plotOutput(ns("scatter_plot"), height = "320px")
      ),
      column(6,
        h5("Survival: High vs. Low Immune Score"),
        plotOutput(ns("km_immune"), height = "320px")
      )
    )
  )
}

mod_immune_server <- function(id, immune_df) {
  moduleServer(id, function(input, output, session) {

    plot_df <- reactive({
      immune_df |>
        mutate(
          stage_clean = case_when(
            grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
            grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
            grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
            grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
            TRUE ~ NA_character_
          )
        )
    })

    output$box_plot <- renderPlot({
      df  <- plot_df()
      col <- input$score_type
      gv  <- if (input$group_by == "stage") "stage_clean" else "subtype"

      df2 <- df |> filter(!is.na(.data[[gv]]))
      pal <- if (input$group_by == "stage")
        c("#2196F3","#4CAF50","#FF9800","#F44336")
      else
        c("#2196F3","#4CAF50","#FF9800","#F44336")

      ggplot(df2, aes(x = .data[[gv]], y = .data[[col]],
                      fill = .data[[gv]])) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        geom_jitter(width = 0.2, alpha = 0.2, size = 0.6) +
        scale_fill_manual(values = setNames(pal, sort(unique(df2[[gv]])))) +
        ggpubr::stat_compare_means(label = "p.signif",
                                   ref.group = sort(unique(df2[[gv]]))[1]) +
        labs(title = paste(col, "by", if(input$group_by=="stage") "Tumor Stage" else "Subtype"),
             x = NULL, y = col) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none")
    })

    output$scatter_plot <- renderPlot({
      df <- plot_df()
      ggplot(df, aes(x = ImmuneScore, y = TumorPurity)) +
        geom_point(alpha = 0.35, size = 1, color = "#2196F3") +
        geom_smooth(method = "lm", color = "#F44336", se = TRUE) +
        ggpubr::stat_cor(method = "pearson") +
        labs(title   = "Immune Score vs. Tumor Purity",
             subtitle = "More immune infiltration → lower tumor purity",
             x = "Immune Score", y = "Tumor Purity") +
        theme_bw(base_size = 12)
    })

    output$km_immune <- renderPlot({
      df <- plot_df() |>
        mutate(immune_group = ifelse(
          ImmuneScore >= median(ImmuneScore, na.rm = TRUE),
          "High immune", "Low immune"))

      km <- survfit(Surv(OS.time / 30.44, OS) ~ immune_group, data = df)
      p  <- ggsurvplot(km, data = df,
                       palette      = c("#2196F3","#F44336"),
                       conf.int     = TRUE, pval = TRUE,
                       risk.table   = FALSE,
                       xlab         = "Time (months)",
                       title        = "Survival by Immune Infiltration",
                       legend.title = "Immune Score",
                       ggtheme      = theme_bw(base_size = 12))
      print(p$plot)
    })
  })
}
