# mod_clinical.R — interactive KM curves by clinical variables

mod_clinical_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    # Method explanation box
    div(class = "alert alert-info",
      h5("About this method"),
      p(strong("Kaplan-Meier (KM) estimator"), " estimates the survival probability over time for groups of patients. ",
        "The ", strong("log-rank test"), " (p-value on the plot) tests whether survival differs between groups. ",
        strong("Cox proportional hazards regression"), " quantifies how much each variable increases or decreases the death hazard, ",
        "expressed as a ", strong("Hazard Ratio (HR):"),
        " HR > 1 = higher risk, HR < 1 = protective. ",
        "The forest plot shows HR with 95% confidence intervals.")
    ),
    fluidRow(
      # Controls
      column(3,
        wellPanel(
          h5("Controls"),
          radioButtons(ns("group_by"), "Stratify by:",
            choices = c("Tumor Stage" = "stage",
                        "Age Group"   = "age",
                        "Subtype"     = "subtype"),
            selected = "stage"),
          hr(),
          sliderInput(ns("time_max"), "Max follow-up (months):",
                      min = 24, max = 240, value = 180, step = 12),
          checkboxInput(ns("show_ci"), "Show confidence intervals", TRUE),
          checkboxInput(ns("show_risk"), "Show risk table", TRUE)
        )
      ),
      # KM plot
      column(9,
        plotOutput(ns("km_plot"), height = "500px"),
        br(),
        plotOutput(ns("forest_plot"), height = "280px")
      )
    )
  )
}

mod_clinical_server <- function(id, survival_df, survival_sub) {
  moduleServer(id, function(input, output, session) {

    # Prepare data reactive
    plot_data <- reactive({
      df <- survival_df |>
        mutate(
          age_group = ifelse(as.numeric(age_at_diagnosis) >= 60, "≥60", "<60"),
          stage_clean = case_when(
            grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
            grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
            grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
            grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
            TRUE ~ NA_character_
          )
        ) |> filter(!is.na(stage_clean))

      if (input$group_by == "subtype") {
        df <- survival_sub |>
          filter(!is.na(subtype)) |>
          mutate(stage_clean = case_when(
            grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
            grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
            grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
            grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
            TRUE ~ NA_character_
          ),
          age_group = ifelse(as.numeric(age_at_diagnosis) >= 60, "≥60", "<60"))
      }
      df |> mutate(OS.time_m = pmin(OS.time / 30.44, input$time_max))
    })

    group_var <- reactive({
      switch(input$group_by,
        stage   = "stage_clean",
        age     = "age_group",
        subtype = "subtype")
    })

    palette_map <- list(
      stage   = c("#2196F3","#4CAF50","#FF9800","#F44336"),
      age     = c("#4CAF50","#FF9800"),
      subtype = c("#2196F3","#4CAF50","#FF9800","#F44336")
    )

    output$km_plot <- renderPlot({
      df  <- plot_data()
      gv  <- group_var()
      pal <- palette_map[[input$group_by]]

      form <- as.formula(paste0("Surv(OS.time_m, OS) ~ ", gv))
      km   <- survfit(form, data = df)

      p <- ggsurvplot(km, data = df,
                      palette    = pal,
                      conf.int   = input$show_ci,
                      risk.table = input$show_risk,
                      pval       = TRUE,
                      xlab       = "Time (months)",
                      ylab       = "Survival probability",
                      title      = paste("Overall Survival by",
                                         switch(input$group_by,
                                           stage="Tumor Stage", age="Age Group", subtype="Molecular Subtype")),
                      ggtheme    = theme_bw(base_size = 13))
      if (input$show_risk) {
        gridExtra::grid.arrange(p$plot, p$table, ncol = 1, heights = c(3, 1))
      } else {
        print(p$plot)
      }
    })

    output$forest_plot <- renderPlot({
      df  <- plot_data()
      gv  <- group_var()
      df2 <- df |>
        mutate(age_num = as.numeric(age_at_diagnosis)) |>
        filter(!is.na(age_num))

      form_cox <- as.formula(paste0("Surv(OS.time_m, OS) ~ ", gv, " + age_num"))
      tryCatch({
        cox <- coxph(form_cox, data = df2)
        ggforest(cox, data = df2, main = "Cox Regression — Hazard Ratios")
      }, error = function(e) {
        ggplot() + annotate("text", x=0.5, y=0.5,
          label="Cox model could not converge for this grouping.",
          size=5) + theme_void()
      })
    })
  })
}
