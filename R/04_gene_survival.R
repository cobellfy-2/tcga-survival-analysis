# 04_gene_survival.R
# Univariate Cox regression for each gene — identifies survival-associated genes.
# Input:  data/processed/survival_df.rds, counts_filtered.rds
# Output: results/gene_cox_results.rds, results/figures/volcano_cox.png
#         results/figures/km_top_genes.png

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "dplyr", "tibble", "DESeq2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "DESeq2") BiocManager::install("DESeq2") else install.packages(pkg)
  }
}

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tibble)
library(DESeq2)

survival_df  <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw   <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))

# Intersect patients
common <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
surv   <- survival_df[survival_df$bcr_patient_barcode %in% common, ]
counts <- counts_raw[, common]
surv   <- surv[match(colnames(counts), surv$bcr_patient_barcode), ]

# VST normalization via DESeq2 (variance-stabilizing transform)
message("Running VST normalization on ", ncol(counts), " samples...")
col_data <- data.frame(row.names = colnames(counts), condition = rep("tumor", ncol(counts)))
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = col_data,
                               design    = ~ 1)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

# --- Univariate Cox: one model per gene ---
# Limit to top 5000 most variable genes for speed
gene_var   <- apply(vst_mat, 1, var)
top_genes  <- names(sort(gene_var, decreasing = TRUE))[1:5000]
vst_top    <- vst_mat[top_genes, ]

message("Running univariate Cox for ", nrow(vst_top), " genes...")

run_cox <- function(expr_vec) {
  df_tmp <- data.frame(
    OS.time = surv$OS.time / 30.44,
    OS      = surv$OS,
    expr    = as.numeric(expr_vec)
  )
  tryCatch({
    fit <- coxph(Surv(OS.time, OS) ~ expr, data = df_tmp)
    s   <- summary(fit)
    c(hr   = s$conf.int[1, "exp(coef)"],
      lower = s$conf.int[1, "lower .95"],
      upper = s$conf.int[1, "upper .95"],
      pval  = s$coefficients[1, "Pr(>|z|)"])
  }, error = function(e) c(hr = NA, lower = NA, upper = NA, pval = NA))
}

cox_results <- t(apply(vst_top, 1, run_cox)) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  mutate(
    padj      = p.adjust(pval, method = "BH"),
    log2HR    = log2(hr),
    neg_log10p = -log10(padj)
  ) |>
  arrange(padj)

saveRDS(cox_results, file.path(RESULTS_DIR, "gene_cox_results.rds"))
message("Top 10 survival-associated genes:")
print(head(cox_results[, c("gene", "hr", "pval", "padj")], 10))

# --- Volcano plot ---
cox_plot <- cox_results |>
  mutate(
    sig = case_when(
      padj < 0.05 & log2HR > 0  ~ "High expr = worse OS",
      padj < 0.05 & log2HR < 0  ~ "High expr = better OS",
      TRUE                       ~ "Not significant"
    )
  )

top_labels <- bind_rows(
  filter(cox_plot, sig == "High expr = worse OS")  |> slice_min(padj, n = 5),
  filter(cox_plot, sig == "High expr = better OS") |> slice_min(padj, n = 5)
)

p_volcano <- ggplot(cox_plot, aes(x = log2HR, y = neg_log10p, color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text(data = top_labels, aes(label = gene),
            size = 3, vjust = -0.5, hjust = 0.5, color = "black") +
  scale_color_manual(values = c(
    "High expr = worse OS"  = "#F44336",
    "High expr = better OS" = "#2196F3",
    "Not significant"       = "grey70")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title  = "TCGA-BRCA: Gene-level Survival Association (Univariate Cox)",
       x      = "log2 Hazard Ratio",
       y      = "-log10 adjusted p-value",
       color  = NULL) +
  theme_bw()

ggsave(file.path(FIGURES_DIR, "volcano_cox.png"),
       plot = p_volcano, width = 9, height = 6, dpi = 150)
message("Saved: volcano_cox.png")

# --- KM plot for top 3 genes ---
top3 <- head(cox_results$gene, 3)

plots <- lapply(top3, function(g) {
  med   <- median(vst_top[g, ])
  group <- ifelse(vst_top[g, ] >= med, "High", "Low")
  df_g  <- data.frame(OS.time = surv$OS.time / 30.44,
                       OS      = surv$OS,
                       group   = group)
  km <- survfit(Surv(OS.time, OS) ~ group, data = df_g)
  ggsurvplot(km, data = df_g, pval = TRUE, conf.int = TRUE,
             palette = c("#F44336", "#2196F3"),
             title = paste0(g, " expression"),
             xlab = "Time (months)", legend.title = "Expression",
             ggtheme = theme_bw())$plot
})

combined <- cowplot::plot_grid(plotlist = plots, ncol = 3)
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
library(cowplot)
combined <- plot_grid(plotlist = plots, ncol = 3)
ggsave(file.path(FIGURES_DIR, "km_top_genes.png"),
       plot = combined, width = 18, height = 5, dpi = 150)
message("Saved: km_top_genes.png")

message("\nGene survival analysis complete. Next: python/01_interactive_viz.py")
