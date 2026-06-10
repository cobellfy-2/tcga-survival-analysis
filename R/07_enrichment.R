# 07_enrichment.R
# Pathway enrichment analysis on survival-associated genes.
#
# BACKGROUND:
# From script 04 we have hundreds of genes associated with survival.
# But a list of gene names is hard to interpret biologically.
# Enrichment analysis asks: which biological PROCESSES or PATHWAYS
# are overrepresented among our survival genes?
#
# Two complementary approaches:
#
# 1. ORA (Over-Representation Analysis)
#    - Take the top significant genes (hard threshold, e.g. padj < 0.05)
#    - Ask: are any GO terms / KEGG pathways more common in this list
#      than expected by chance? (Hypergeometric test)
#    - Simple and fast, but loses rank information
#
# 2. GSEA (Gene Set Enrichment Analysis)
#    - Rank ALL genes by their association with survival (log2HR * -log10p)
#    - Ask: do genes in a pathway cluster at the top or bottom of this ranking?
#    - No hard threshold needed, uses the full signal
#    - More powerful and less sensitive to cutoff choice
#
# Input:  results/gene_cox_results.rds (with symbol column)
# Output: results/enrichment_results.rds
#         results/figures/ora_dotplot.png
#         results/figures/gsea_dotplot.png
#         results/figures/gsea_top_pathway.png

source("config.R")

for (pkg in c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "ggplot2", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "enrichplot", "org.Hs.eg.db")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

cox_results <- readRDS(file.path(RESULTS_DIR, "gene_cox_results.rds"))

# Clean Ensembl IDs (remove version suffix e.g. ".11")
cox_results <- cox_results |>
  mutate(ensembl_clean = sub("\\..*", "", gene))

# Map Ensembl → Entrez ID (needed for KEGG)
entrez_map <- suppressMessages(
  mapIds(org.Hs.eg.db,
         keys    = cox_results$ensembl_clean,
         column  = "ENTREZID",
         keytype = "ENSEMBL",
         multiVals = "first")
)
cox_results$entrez <- entrez_map

# ── 1. ORA ────────────────────────────────────────────────────────────────────
# Significant survival genes (FDR < 5%)
sig_genes <- cox_results |>
  filter(padj < 0.05, !is.na(entrez)) |>
  pull(entrez)

background <- cox_results |>
  filter(!is.na(entrez)) |>
  pull(entrez)

message("ORA input: ", length(sig_genes), " significant genes / ",
        length(background), " background")

ora_go <- enrichGO(
  gene          = sig_genes,
  universe      = background,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",          # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (nrow(ora_go) > 0) {
  p_ora <- dotplot(ora_go, showCategory = 20, title = "ORA: GO Biological Process") +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(FIGURES_DIR, "ora_dotplot.png"),
         plot = p_ora, width = 10, height = 8, dpi = 150)
  message("Saved: ora_dotplot.png — ", nrow(ora_go), " enriched GO terms")
} else {
  message("No significant ORA terms found")
}

# ── 2. GSEA ───────────────────────────────────────────────────────────────────
# Rank statistic: signed -log10(p) — captures both direction and significance
# Positive = high expression → worse survival (risk gene)
# Negative = high expression → better survival (protective gene)
gsea_input <- cox_results |>
  filter(!is.na(entrez), !is.na(pval)) |>
  mutate(rank_stat = sign(log2HR) * (-log10(pval))) |>
  arrange(desc(rank_stat))

gene_list <- setNames(gsea_input$rank_stat, gsea_input$entrez)
gene_list <- gene_list[!duplicated(names(gene_list))]

message("GSEA input: ", length(gene_list), " ranked genes")

set.seed(42)
gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  nPermSimple   = 10000,   # more permutations → fewer NA p-values for unbalanced pathways
  verbose       = FALSE
)

if (nrow(gsea_go) > 0) {
  # Dotplot: top enriched + top depleted pathways
  p_gsea <- dotplot(gsea_go, showCategory = 20,
                    split = ".sign", title = "GSEA: GO Biological Process") +
    facet_grid(. ~ .sign) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(FIGURES_DIR, "gsea_dotplot.png"),
         plot = p_gsea, width = 14, height = 8, dpi = 150)
  message("Saved: gsea_dotplot.png — ", nrow(gsea_go), " enriched pathways")

  # Enrichment plot for top pathway
  top_pathway <- gsea_go@result$ID[1]
  p_enrich <- gseaplot2(gsea_go, geneSetID = top_pathway,
                        title = gsea_go@result$Description[1])
  ggsave(file.path(FIGURES_DIR, "gsea_top_pathway.png"),
         plot = p_enrich, width = 9, height = 6, dpi = 150)
  message("Saved: gsea_top_pathway.png")
  message("Top enriched pathway: ", gsea_go@result$Description[1])
} else {
  message("No significant GSEA pathways found — try pvalueCutoff = 0.1")
}

# Save results
saveRDS(list(ora_go = ora_go, gsea_go = gsea_go),
        file.path(RESULTS_DIR, "enrichment_results.rds"))

message("\nEnrichment analysis complete. Next: source('R/08_immune.R')")
