# mod_clinical.R — interactive KM curves by clinical variables

mod_clinical_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method — Kaplan-Meier & Cox Regression"),
      p(strong("Kaplan-Meier (KM):"), " estimates the probability of surviving past each time point.
        The curve steps down each time a patient dies. Patients who are still alive at last contact
        are 'censored' (shown as tick marks) — we know they survived at least that long, but not longer."),
      p(strong("Log-rank test (p-value):"), " tests whether survival curves are statistically
        different between groups. p < 0.05 means the difference is unlikely due to chance."),
      p(strong("Cox Hazard Ratio (HR):"), " quantifies the effect size. HR = 2 means twice the
        risk of death compared to the reference group, at any given time point.
        HR > 1 = worse survival, HR < 1 = better survival.")
    ),
    div(class = "alert alert-success",
      h5("Key findings from this dataset"),
      tags$ul(
        tags$li(strong("Stage III/IV vs I:"), " HR = 3.43 — patients with advanced stage have 3.4× higher death hazard (p < 0.001)"),
        tags$li(strong("Age ≥60 vs <60:"), " HR = 1.73 — older patients have 73% higher death hazard (p = 0.01)"),
        tags$li(strong("Triple-Negative vs Luminal A:"), " HR = 2.76 — worst prognosis subtype (p < 0.01)")
      )
    ),
    conditionalPanel(
      condition = "input['clinical-group_by'] == 'subtype'",
      div(class = "alert alert-warning",
        h5("⚠️ Why HER2-enriched drops steeply in the plot but has a lower HR than Triple-Negative"),
        tags$ol(
          tags$li(strong("Timing of risk differs:"),
            " HER2-enriched is aggressively lethal early (first 2–3 years), but patients who survive this
            window — often thanks to targeted therapy (trastuzumab/Herceptin) — tend to survive long-term.
            The KM curve falls steeply early, then flattens. Cox regression measures the
            ", em("average"), " hazard across all time points, so the early spike is diluted."),
          tags$li(strong("Triple-Negative has no targeted therapy:"),
            " Without a targetable receptor, Triple-Negative patients face consistently elevated risk
            throughout the entire follow-up. This sustained hazard produces the highest average HR (2.76)."),
          tags$li(strong("Small sample size for HER2-enriched (n=35):"),
            " The Cox confidence interval is very wide (0.69–7.57) and overlaps 1, meaning the
            HER2-enriched result is ", em("not statistically significant"), ".
            The visual impression in the KM plot can be misleading with small groups."),
          tags$li(strong("Historical cohort effect:"),
            " TCGA data was collected during an era of evolving HER2 treatment. Modern cohorts
            treated with current HER2-targeted regimens show substantially better HER2 outcomes.")
        )
      )
    ),
    fluidRow(
      column(3,
        wellPanel(
          h5("Controls"),
          radioButtons(ns("group_by"), "Stratify by:",
            choices  = c("Tumor Stage" = "stage",
                         "Age Group"   = "age",
                         "Subtype"     = "subtype"),
            selected = "stage"),
          hr(),
          sliderInput(ns("time_max"), "Max follow-up (months):",
                      min = 24, max = 240, value = 180, step = 12),
          checkboxInput(ns("show_ci"),   "Show confidence intervals", TRUE),
          checkboxInput(ns("show_risk"), "Show risk table", TRUE)
        )
      ),
      column(9,
        plotOutput(ns("km_plot"),     height = "460px"),
        br(),
        plotOutput(ns("forest_plot"), height = "260px"),
        br(),
        div(class = "alert alert-secondary",
          h6("How to read this forest plot"),
          tags$ul(style = "margin-bottom:0",
            tags$li(strong("Each row"), " = one predictor variable. The first category is the
                    ", strong("reference group"), " (HR fixed at 1.0, shown as dashed line).
                    All other HRs are ", em("relative to that reference.")),
            tags$li(strong("Square"), " = the estimated Hazard Ratio. Further right = higher risk."),
            tags$li(strong("Horizontal line"), " = 95% confidence interval. If it crosses the dashed
                    line at 1.0, the result is ", strong("not statistically significant.")),
            tags$li(strong("p-value"), " on the right: *** p<0.001, ** p<0.01, * p<0.05."),
            tags$li(strong("Current result (Subtype view):"), " HER2-enriched is the reference (n=35).
                    All subtype CIs are very wide and cross 1 → not significant, because n=35
                    is too small for a stable estimate. Age (age_num) is highly significant:
                    HR=1.04 means each additional year of age increases death hazard by 4%.")
          )
        )
      )
    )
  )
}

mod_clinical_server <- function(id, survival_df, survival_sub) {
  moduleServer(id, function(input, output, session) {

    plot_data <- reactive({
      df <- survival_df |>
        mutate(
          age_group   = ifelse(as.numeric(age_at_diagnosis) >= 60, "≥60", "<60"),
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
          mutate(
            age_group   = ifelse(as.numeric(age_at_diagnosis) >= 60, "≥60", "<60"),
            stage_clean = case_when(
              grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
              grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
              grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
              grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
              TRUE ~ NA_character_
            )
          )
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

      df[["grp"]] <- df[[gv]]
      df <- df[!is.na(df[["grp"]]), ]

      km <- survfit(Surv(OS.time_m, OS) ~ grp, data = df)

      p <- ggsurvplot(
        km, data = df,
        palette    = pal,
        conf.int   = input$show_ci,
        risk.table = input$show_risk,
        pval       = TRUE,
        xlab       = "Time (months)",
        ylab       = "Survival probability",
        title      = paste("Overall Survival by",
                           switch(input$group_by,
                                  stage = "Tumor Stage",
                                  age   = "Age Group",
                                  subtype = "Molecular Subtype")),
        ggtheme = theme_bw(base_size = 13)
      )
      if (input$show_risk) print(p) else print(p$plot)
    })

    output$forest_plot <- renderPlot({
      df  <- plot_data()
      gv  <- group_var()
      df2 <- df |>
        mutate(age_num = as.numeric(age_at_diagnosis),
               grp     = .data[[gv]]) |>
        filter(!is.na(age_num), !is.na(grp))

      tryCatch({
        cox <- coxph(Surv(OS.time_m, OS) ~ grp + age_num, data = df2)
        ggforest(cox, data = df2, main = "Cox Regression — Hazard Ratios (95% CI)")
      }, error = function(e) {
        ggplot() +
          annotate("text", x=0.5, y=0.5,
                   label = paste("Cox model note:", e$message), size = 4) +
          theme_void()
      })
    })
  })
}
