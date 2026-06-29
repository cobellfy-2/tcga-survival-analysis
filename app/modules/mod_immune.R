# mod_immune.R — immune landscape and survival

mod_immune_ui <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    div(class = "alert alert-info",
      h5("About this method — GSVA-based Immune Deconvolution"),
      p(strong("What are we measuring?"), " Bulk RNA-seq measures the combined expression of all cells in a tumor biopsy —
        cancer cells, immune cells, fibroblasts, endothelium. We cannot directly count cell types,
        but we can infer their relative abundance using ", strong("gene expression signatures."),
        " Cell types always express certain marker genes (e.g. CD3E, CD8A for T cells;
        FAP, COL1A1 for stromal fibroblasts). If a tumor has many T cells, those marker genes
        will be highly expressed in the bulk profile."),
      p(strong("Method:"), " We use ", strong("GSVA (Gene Set Variation Analysis)"),
        " with the immune and stromal gene sets from ",
        strong("Yoshihara et al. 2013 (Nature Communications)"),
        " — the same gene sets underlying the ESTIMATE algorithm.
        GSVA computes per-sample enrichment scores ranging roughly ", strong("-1 to +1"),
        " (not absolute cell counts — relative enrichment compared to the cohort average).
        A score near ", strong("0"), " = cohort-average level; positive = above average; negative = below average."),
      tags$ul(
        tags$li(strong("ImmuneScore:"), " Relative immune cell infiltration (T cells, B cells, NK cells)."),
        tags$li(strong("StromalScore:"), " Relative stromal content (fibroblasts, endothelial cells, ECM)."),
        tags$li(strong("ESTIMATEScore:"), " Combined immune + stromal (inverse of tumor purity)."),
        tags$li(strong("TumorPurity:"), " Estimated fraction of malignant cells (1 − ESTIMATE signal).")
      )
    ),
    div(class = "alert alert-success",
      h5("Key findings — TCGA-BRCA immune landscape (n = 1,033)"),
      tags$ul(
        tags$li(strong("Tumor purity median = 82.2%:"),
          " The typical BRCA tumor biopsy is ~82% cancer cells, ~18% microenvironment.
          This is moderately high purity — consistent with BRCA being less infiltrated than
          'hot' tumors like melanoma or NSCLC."),
        tags$li(strong("ImmuneScore is near zero for most patients"),
          " (median ≈ 0, IQR: −0.1 to +0.1 on the GSVA scale),
          indicating most tumors show close-to-average immune infiltration.
          A minority of outliers have strongly positive or negative scores."),
        tags$li(strong("High vs. Low immune infiltration does NOT significantly predict survival"),
          " (log-rank p = 0.90). This is a biologically meaningful null result: unlike melanoma
          or NSCLC, where immune infiltration strongly predicts immunotherapy response,
          BRCA is generally considered an 'immune-cold' cancer. The survival impact of
          immune cells is subtype-dependent — it matters most in Triple-Negative BRCA,
          where immune checkpoint inhibitors are increasingly used."),
        tags$li(strong("Stromal score median = 0.5"),
          " (on a 0–1 scale after transformation): moderate stromal content,
          consistent with the desmoplastic (fibrotic) microenvironment characteristic of invasive BRCA.")
      )
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
                       selected = "stage"),
          hr(),
          div(style = "font-size:0.85em; color:#555;",
            p(strong("Score scale:"), " GSVA enrichment scores are relative to the cohort mean.
              A score of 0 = average; ±0.5 = strong deviation. TumorPurity is converted
              to a 0–1 fraction (not on the GSVA scale).")
          )
        )
      ),
      column(8,
        plotOutput(ns("box_plot"), height = "340px"),
        br(),
        div(class = "alert alert-secondary",
          h6("How to read the box plot"),
          tags$ul(style = "margin-bottom:0",
            tags$li("Each box = distribution of the selected score across patients in that group.
              The horizontal line = median; box edges = 25th and 75th percentiles; whiskers = 1.5× IQR."),
            tags$li("Significance labels (ns / * / ** / ***) compare each group to the leftmost group
              (reference) using the Wilcoxon rank-sum test — a non-parametric test appropriate for
              non-normally distributed scores."),
            tags$li(strong("For Stage:"),
              " If ImmuneScore increases with stage, it suggests more advanced tumors recruit more
              immune cells — potentially as a failed anti-tumor response."),
            tags$li(strong("For Subtype:"),
              " Triple-Negative BRCA typically shows the highest immune infiltration
              (immune-hot subtype), while Luminal A/B tend to be immune-cold.")
          )
        )
      )
    ),
    fluidRow(
      column(6,
        h5("Immune Score vs. Tumor Purity"),
        plotOutput(ns("scatter_plot"), height = "320px"),
        div(class = "alert alert-secondary",
          style = "margin-top:8px; font-size:0.88em;",
          p(strong("Expected anti-correlation:"),
            " A tumor with more immune cells (high ImmuneScore) should have fewer cancer cells
            (lower TumorPurity). A strong negative correlation here validates that our
            GSVA-based scores behave as expected. The Pearson r and p-value are shown on the plot.")
        )
      ),
      column(6,
        h5("Survival: High vs. Low Immune Score"),
        plotOutput(ns("km_immune"), height = "320px"),
        div(class = "alert alert-secondary",
          style = "margin-top:8px; font-size:0.88em;",
          p(strong("Log-rank p = 0.90 — not significant."),
            " Splitting patients at the median immune score does not separate their survival curves.
            This does not mean immune biology is irrelevant in BRCA — it means the signal is
            subtype-specific (strongest in Triple-Negative) and may require longer follow-up
            or continuous-score modeling rather than a median split.")
        )
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

      if (!gv %in% colnames(df)) {
        return(
          ggplot() + annotate("text", x=0.5, y=0.5,
            label=paste("Column", gv, "not available in this dataset.\n",
                        "Run R/05_subtypes.R to add subtype annotations."),
            size=4.5, hjust=0.5) + theme_void()
        )
      }
      df2 <- df |> filter(!is.na(.data[[gv]]))
      n_groups <- length(unique(df2[[gv]]))
      pal <- c("#2196F3","#4CAF50","#FF9800","#F44336")[seq_len(n_groups)]

      ggplot(df2, aes(x = .data[[gv]], y = .data[[col]],
                      fill = .data[[gv]])) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        geom_jitter(width = 0.2, alpha = 0.2, size = 0.6) +
        scale_fill_manual(values = setNames(pal, sort(unique(df2[[gv]])))) +
        ggpubr::stat_compare_means(label = "p.signif",
                                   ref.group = sort(unique(df2[[gv]]))[1]) +
        labs(title = paste(col, "by", if(input$group_by=="stage") "Tumor Stage" else "Molecular Subtype"),
             subtitle = "Significance vs. leftmost group (Wilcoxon). GSVA scores: 0 = cohort mean.",
             x = NULL, y = col) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none")
    })

    output$scatter_plot <- renderPlot({
      df <- plot_df()
      ggplot(df, aes(x = ImmuneScore, y = TumorPurity)) +
        geom_point(alpha = 0.3, size = 0.9, color = "#2196F3") +
        geom_smooth(method = "lm", color = "#F44336", se = TRUE) +
        ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
        labs(title   = "Immune Score vs. Tumor Purity",
             subtitle = "Higher immune infiltration → lower tumor cell fraction (anti-correlation validates method)",
             x = "Immune Score (GSVA)", y = "Tumor Purity (estimated)") +
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
                       pval.method  = TRUE,
                       risk.table   = FALSE,
                       xlab         = "Time (months)",
                       title        = "Survival by Immune Infiltration (median split)",
                       subtitle     = "Log-rank p = 0.90 — not significant in overall BRCA cohort",
                       legend.title = "Immune Score",
                       ggtheme      = theme_bw(base_size = 12))
      print(p$plot)
    })
  })
}
