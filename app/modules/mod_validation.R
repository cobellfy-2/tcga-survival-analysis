# mod_validation.R — model validation: C-index, Brier score, nomogram, calibration

mod_validation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About Model Validation"),
      p(strong("Why validate?"),
        " Evaluating a model on the same data it was trained on is circular — the model has
        already seen those patients and will look artificially good. Instead we use a ",
        strong("70/30 train-test split"), ": train on 707 patients, evaluate on 303 patients
        the model has never seen. Any result on the test set is a genuine out-of-sample estimate."),
      p("We measure model quality two ways: ",
        strong("discrimination"), " (does the model rank patients correctly — C-index) and ",
        strong("calibration"), " (are the predicted probabilities accurate in absolute terms — Brier score + calibration plot).")
    ),
    div(class = "alert alert-success",
      h5("Key findings at a glance"),
      tags$ul(
        tags$li(strong("Lasso gene signature: C-index = 0.881"),
          " — in 88% of randomly chosen patient pairs, the model correctly identifies
          who has higher risk. Excellent (clinical benchmark is typically 0.65–0.75)."),
        tags$li(strong("Combined model: C-index = 0.893"),
          " — adding age and tumor stage on top of the 52-gene signature gives a marginal
          +1.2 percentage points. The gene signature already captures most clinical signal."),
        tags$li(strong("Clinical only: C-index = 0.752"),
          " — standard variables (age + stage) alone are good but the molecular signature
          adds ", strong("+12.9 points"), " of discrimination."),
        tags$li(strong("Brier score at 36 months"),
          ": Reference = 0.065 | Clinical = 0.060 | Lasso = 0.035 | Combined = 0.028.
          The gene model reduces prediction error by ", strong("46%"), " relative to the
          no-covariate baseline.")
      )
    ),
    tabsetPanel(
      # --- C-index ---
      tabPanel("C-index",
        br(),
        plotOutput(ns("cindex_plot"), height = "360px"),
        br(),
        div(class = "alert alert-secondary",
          h6("How to read this chart"),
          tags$ul(style = "margin-bottom:0",
            tags$li(strong("C-index (Concordance index)"),
              " is the probability that, for two randomly chosen patients, the model
              assigns a higher risk score to the one who dies earlier. It answers:
              'does the model rank patients in the right order?'"),
            tags$li(strong("0.5 = random (coin flip)"), " | 0.7 = good | 0.8 = excellent | 0.9 = near-perfect.
              Published survival models in oncology typically achieve 0.65–0.75."),
            tags$li(strong("Grey bar (Null = 0.5):"),
              " The floor. Any useful model must beat this."),
            tags$li(strong("Green (Clinical = 0.752):"),
              " Age and tumor stage alone already provide good discrimination —
              these are well-established prognostic factors."),
            tags$li(strong("Blue (Lasso = 0.881):"),
              " The 52-gene expression signature adds a large jump.
              The genes capture molecular variation that stage and age cannot see
              — e.g., a Stage II patient with an aggressive molecular subtype
              looks the same clinically as a Stage II patient with a dormant one."),
            tags$li(strong("Red (Combined = 0.893):"),
              " Combining everything gives the best discrimination,
              but the marginal gain over Lasso alone is small (+0.012),
              suggesting the gene signature is the primary driver of predictive power.")
          )
        )
      ),
      # --- Brier Score ---
      tabPanel("Brier Score",
        br(),
        plotOutput(ns("brier_plot"), height = "380px"),
        br(),
        div(class = "alert alert-secondary",
          h6("How to read this chart"),
          tags$ul(style = "margin-bottom:0",
            tags$li(strong("Brier score"),
              " = mean squared error between predicted survival probability and actual outcome
              (1 = died, 0 = alive). Lower is better. Unlike C-index, Brier score tests
              ", strong("absolute accuracy"), ", not just ranking. A model that says
              'everyone has 60% survival' can have a good C-index but a bad Brier score."),
            tags$li(strong("Reference (grey dashed):"),
              " The Kaplan-Meier curve with no covariates — predicts the same probability
              for every patient. Everything below this line is a genuine improvement."),
            tags$li(strong("At 36 months"),
              ": Reference = 0.065 | Clinical = 0.060 | Lasso = 0.035 | Combined = 0.028.",
              " The Lasso model cuts error by ", strong("46%"), " vs Reference;
              Combined cuts it by ", strong("57%"), "."),
            tags$li(strong("Why are all values so small (< 0.08)?"),
              " BRCA has a low event rate (~10% died in this cohort).
              Even a naive model that predicts 'everyone survives' will have a low Brier
              score because it's right 90% of the time. What matters is the ",
              em("relative"), " improvement over the reference line — and that is large."),
            tags$li(strong("Early time points (0–18 months):"),
              " Lasso and Combined are slightly above the reference here.
              Very few deaths occur in the first year of BRCA follow-up, so the model's
              gene-driven risk scores are not yet 'activated' — this is normal for a
              disease with a slow natural history.")
          )
        )
      ),
      # --- Nomogram ---
      tabPanel("Nomogram",
        br(),
        div(class = "alert alert-info",
          h5("What is a nomogram?"),
          p("A nomogram translates a statistical model into a ", strong("clinical decision tool"),
            " — a graphical calculator a doctor can use at the bedside without running code.
            For each predictor, you assign a point score. The total score maps to a predicted
            survival probability. The nomogram here uses: ",
            strong("Age, Tumor Stage"), ", and the ", strong("52-gene Lasso risk score.")),
          h6("How to use it — step by step"),
          tags$ol(
            tags$li("Find your patient's ", strong("Age"), " on the Age row.
              Draw a vertical line straight up to the ", strong("Points"), " scale. Note the value."),
            tags$li("Repeat for ", strong("Tumor Stage"), " (I / II / III / IV)."),
            tags$li("Repeat for the ", strong("Risk Score"),
              " (the numeric output of the 52-gene Lasso model)."),
            tags$li(strong("Sum all three point values"), " → locate the total on the ",
              strong("Total Points"), " row."),
            tags$li("Draw a vertical line down to read predicted ",
              strong("3-year survival"), " and ", strong("5-year survival"), " probabilities.")
          ),
          p(style = "margin-bottom:0; color:#555;",
            strong("Example:"), " a 70-year-old with Stage III disease and high risk score
            accumulates ~180 total points → predicted 5-year survival ≈ 45%.")
        ),
        div(style = "text-align:center; background:#f8f9fa; padding:12px; border-radius:6px;",
          imageOutput(ns("nomogram_img"), height = "520px", width = "100%")
        )
      ),
      # --- Calibration ---
      tabPanel("Calibration",
        br(),
        div(class = "alert alert-info",
          h5("What is a calibration plot?"),
          p("C-index tells you whether the model ", em("ranks"), " patients correctly.
            Calibration tells you whether the model's ", strong("absolute probabilities are accurate."),
            " If the model says 70% 3-year survival, do 70% of those patients actually survive 3 years?"),
          tags$ul(
            tags$li(strong("X-axis — Predicted probability:"),
              " Model's 3-year survival forecast. Patients are grouped into 5 bins of similar predicted risk."),
            tags$li(strong("Y-axis — Observed probability (Kaplan-Meier):"),
              " Actual fraction alive at 3 years within each bin."),
            tags$li(strong("Red dashed diagonal = perfect calibration"),
              ": predicted = observed. Points above the diagonal = model underestimates survival
              (patients do better than predicted). Points below = overestimates."),
            tags$li(strong("Blue line = your model's calibration curve."),
              " Points with ", strong("error bars"), " = 95% CI from 100 bootstrap samples.
              Wide bars = few patients in that bin, high uncertainty."),
            tags$li(strong("Your result:"),
              " Points track close to the diagonal, especially in the mid-range (0.6–0.85),
              confirming the model produces trustworthy absolute survival predictions.
              Minor deviation at extremes is expected with small bin sizes (n ≈ 60 per bin).")
          )
        ),
        div(style = "text-align:center; background:#f8f9fa; padding:12px; border-radius:6px;",
          imageOutput(ns("calibration_img"), height = "550px", width = "100%")
        )
      )
    )
  )
}

mod_validation_server <- function(id, validation_res, figures_dir) {
  moduleServer(id, function(input, output, session) {

    output$cindex_plot <- renderPlot({
      req(validation_res)
      ci <- validation_res$cindex
      df <- data.frame(
        Model  = c("Null\n(random)", "Clinical\n(age + stage)",
                   "Lasso\nsignature", "Combined"),
        Cindex = c(0.5, ci$clinical, ci$lasso, ci$combined),
        fill   = c("grey70", "#4CAF50", "#2196F3", "#F44336")
      )
      ggplot(df, aes(x = reorder(Model, Cindex), y = Cindex, fill = fill)) +
        geom_col(width = 0.55) +
        geom_text(aes(label = sprintf("%.3f", Cindex)), hjust = -0.15,
                  size = 5.5, fontface = "bold") +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
        annotate("text", x = 0.6, y = 0.515, label = "Random chance (0.5)",
                 color = "grey40", size = 3.5) +
        scale_fill_identity() +
        coord_flip() +
        scale_y_continuous(limits = c(0, max(df$Cindex) * 1.18)) +
        labs(title    = "Model Discrimination: C-index on Held-out Test Set (n = 303)",
             subtitle = "Higher = better | Trained on 70% (n = 707), evaluated on unseen 30%",
             x = NULL, y = "C-index (Concordance)") +
        theme_bw(base_size = 14) +
        theme(panel.grid.major.y = element_blank())
    })

    output$brier_plot <- renderPlot({
      req(validation_res)
      br <- validation_res$brier
      if (is.null(br)) {
        return(ggplot() +
          annotate("text", x=0.5, y=0.5, label="Brier score not available.", size=6) +
          theme_void())
      }
      brier_df <- as.data.frame(br$AppErr) |>
        mutate(time = br$time) |>
        tidyr::pivot_longer(-time, names_to = "model", values_to = "brier") |>
        filter(!is.na(brier))

      ggplot(brier_df, aes(x = time, y = brier, color = model, linetype = model)) +
        geom_line(linewidth = 1.3) +
        geom_point(size = 2.5) +
        scale_color_manual(values = c(
          "Reference" = "grey50", "Clinical" = "#4CAF50",
          "Lasso"     = "#2196F3", "Combined" = "#F44336")) +
        scale_linetype_manual(values = c(
          "Reference" = "dashed", "Clinical" = "solid",
          "Lasso"     = "solid",  "Combined" = "solid")) +
        annotate("text", x = 40, y = 0.071, label = "Reference\n(KM only)",
                 color = "grey50", size = 3.2, hjust = 1) +
        labs(title    = "Brier Score over Follow-up Time (Test Set, n = 303)",
             subtitle = "Lower = better | Below the grey reference = genuine predictive improvement",
             x = "Follow-up time (months)", y = "Brier Score",
             color = "Model", linetype = "Model") +
        theme_bw(base_size = 13) +
        theme(legend.position = "right")
    })

    output$nomogram_img <- renderImage({
      path <- file.path(figures_dir, "nomogram.png")
      list(src = path, contentType = "image/png", width = "100%", alt = "Nomogram")
    }, deleteFile = FALSE)

    output$calibration_img <- renderImage({
      path <- file.path(figures_dir, "calibration_3yr.png")
      list(src = path, contentType = "image/png", width = "100%", alt = "Calibration plot")
    }, deleteFile = FALSE)
  })
}
