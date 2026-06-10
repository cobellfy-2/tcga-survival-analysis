# 06_lasso_cox.R
# Lasso-penalized Cox regression to identify a prognostic gene signature.
#
# BACKGROUND:
# Standard Cox regression fails when predictors (genes) >> patients.
# Lasso (Least Absolute Shrinkage and Selection Operator) adds an L1 penalty
# that shrinks weak coefficients to exactly zero → automatic gene selection.
# Cross-validation picks the optimal penalty strength (lambda).
# Result: a sparse gene signature + a continuous risk score per patient.
#
# Input:  data/processed/survival_df.rds, counts_filtered.rds
# Output: results/lasso_results.rds
#         results/figures/lasso_cv.png
#         results/figures/lasso_signature.png
#         results/figures/km_lasso_risk.png

source("config.R")

for (pkg in c("glmnet", "survival", "survminer", "ggplot2", "dplyr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(glmnet)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tibble)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw  <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))

# Load gene symbols from script 04 results
cox_gene_results <- readRDS(file.path(RESULTS_DIR, "gene_cox_results.rds"))

# --- Prepare matrix ---
common  <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
surv_df <- survival_df[match(common, survival_df$bcr_patient_barcode), ]
counts  <- counts_raw[, common]

# VST normalization (reuse from script 04 logic, quick version via log2 CPM)
# For Lasso we use log2(CPM+1) — faster than full DESeq2 VST
lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm   <- log2(cpm + 1)

# Use top 3000 most variable genes (balance speed vs. signal)
gene_var  <- apply(log_cpm, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:3000]
X         <- t(log_cpm[top_genes, ])   # samples × genes matrix
Y         <- Surv(surv_df$OS.time / 30.44, surv_df$OS)

message("Lasso input: ", nrow(X), " samples × ", ncol(X), " genes")

# --- Cross-validated Lasso Cox ---
set.seed(42)
message("Running 10-fold cross-validation (may take ~2 min)...")
cv_fit <- cv.glmnet(
  X, Y,
  family   = "cox",
  alpha    = 1,
  nfolds   = 10,
  parallel = FALSE,
  cox.ties = "efron",                    # efron is more accurate; silences v5.1 warning
  control  = list(maxit = 250000)       # more iterations for convergence at small lambdas
)

# Plot CV curve
png(file.path(FIGURES_DIR, "lasso_cv.png"), width = 800, height = 500)
plot(cv_fit, main = "Lasso Cox: Cross-validation — choosing lambda")
abline(v = log(cv_fit$lambda.min),  col = "blue",  lty = 2)
abline(v = log(cv_fit$lambda.1se),  col = "red",   lty = 2)
legend("topright",
       legend = c("lambda.min (best CV)", "lambda.1se (more sparse)"),
       col    = c("blue", "red"), lty = 2)
dev.off()
message("Saved: lasso_cv.png")

# --- Extract signature genes ---
# Try lambda.1se first (parsimonious); fall back to lambda.min if 0 genes selected
coef_fit  <- coef(cv_fit, s = "lambda.1se")
sig_genes <- rownames(coef_fit)[coef_fit[, 1] != 0]

if (length(sig_genes) == 0) {
  message("lambda.1se selected 0 genes — falling back to lambda.min")
  coef_fit  <- coef(cv_fit, s = "lambda.min")
  sig_genes <- rownames(coef_fit)[coef_fit[, 1] != 0]
}

sig_coefs <- coef_fit[sig_genes, 1]
message("\nLasso signature: ", length(sig_genes), " genes selected")

# Map to gene symbols
ensembl_clean <- sub("\\..*", "", sig_genes)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(org.Hs.eg.db)
  symbols <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = ensembl_clean,
           column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
} else {
  symbols <- setNames(sig_genes, sig_genes)
}

sig_df <- tibble(
  ensembl  = sig_genes,
  symbol   = ifelse(is.na(symbols), ensembl_clean, symbols),
  coef     = sig_coefs,
  direction = ifelse(sig_coefs > 0, "Risk (high = worse)", "Protective (high = better)")
) |> arrange(desc(abs(coef)))

message("\nTop signature genes:")
print(sig_df, n = 20)

# --- Signature bar plot ---
p_sig <- sig_df |>
  slice_max(abs(coef), n = min(30, nrow(sig_df))) |>
  mutate(symbol = factor(symbol, levels = symbol[order(coef)])) |>
  ggplot(aes(x = coef, y = symbol, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Risk (high = worse)"        = "#F44336",
                               "Protective (high = better)" = "#2196F3")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title    = "TCGA-BRCA: Lasso Cox Prognostic Gene Signature",
       subtitle = paste0(length(sig_genes), " genes selected by Lasso (lambda.1se)"),
       x = "Lasso coefficient", y = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "lasso_signature.png"),
       plot = p_sig, width = 8, height = max(5, length(sig_genes) * 0.3 + 2),
       dpi = 150)
message("Saved: lasso_signature.png")

# --- Compute risk score and KM: high vs. low risk ---
risk_score <- as.vector(X[, sig_genes, drop = FALSE] %*% sig_coefs)
risk_df    <- data.frame(
  patient    = rownames(X),
  risk_score = risk_score,
  OS.time    = surv_df$OS.time / 30.44,
  OS         = surv_df$OS
) |>
  mutate(risk_group = ifelse(risk_score >= median(risk_score), "High risk", "Low risk"))

km_risk <- survfit(Surv(OS.time, OS) ~ risk_group, data = risk_df)
p_km_risk <- ggsurvplot(
  km_risk, data = risk_df,
  palette      = c("#F44336", "#2196F3"),
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = TRUE,
  xlab         = "Time (months)",
  ylab         = "Overall Survival probability",
  title        = "TCGA-BRCA: Lasso Risk Score — High vs. Low Risk",
  legend.title = "Risk group",
  ggtheme      = theme_bw()
)

ggsave(file.path(FIGURES_DIR, "km_lasso_risk.png"),
       plot = p_km_risk$plot, width = 8, height = 5, dpi = 150)
message("Saved: km_lasso_risk.png")

# Save all results
used_lambda <- if (length(sig_genes) > 0) cv_fit$lambda.min else cv_fit$lambda.1se
saveRDS(list(cv_fit    = cv_fit,
             sig_df    = sig_df,
             risk_df   = risk_df,
             lambda    = used_lambda),
        file.path(RESULTS_DIR, "lasso_results.rds"))

message("\nLasso analysis complete. Next: source('R/07_enrichment.R')")
