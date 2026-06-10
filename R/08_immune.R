# 08_immune.R
# Immune landscape estimation via ESTIMATE + survival correlation.
#
# BACKGROUND:
# Tumor tissue is not just cancer cells — it contains immune cells,
# fibroblasts, endothelial cells, etc. The composition of this
# "tumor microenvironment" strongly influences patient outcomes.
#
# ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumors
# using Expression data) uses gene expression signatures to compute:
#
#   ImmuneScore  — level of immune cell infiltration
#   StromalScore — level of stromal (connective tissue) cell infiltration
#   ESTIMATEScore — combined tumor purity estimate
#   TumorPurity  — fraction of the sample that is actual tumor cells
#
# Key insight: high ImmuneScore = more immune cells in tumor → in many
# cancers this correlates with better survival (immune surveillance).
# In BRCA the relationship is subtype-dependent.
#
# Input:  data/processed/counts_filtered.rds, survival_df.rds
# Output: results/immune_scores.rds
#         results/figures/immune_boxplot.png
#         results/figures/immune_survival_scatter.png
#         results/figures/km_immune.png

source("config.R")

for (pkg in c("survival", "survminer", "ggplot2", "dplyr", "tibble", "ggpubr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in c("GSVA", "GSEABase")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)
library(GSVA)
library(GSEABase)

survival_df <- readRDS(file.path(PROCESSED_DIR, "survival_df.rds"))
counts_raw  <- readRDS(file.path(PROCESSED_DIR, "counts_filtered.rds"))

# --- Prepare log2(CPM+1) expression matrix with gene symbols ---
message("Preparing expression matrix...")
common  <- intersect(survival_df$bcr_patient_barcode, colnames(counts_raw))
counts  <- counts_raw[, common]

lib_sizes <- colSums(counts)
cpm       <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm   <- log2(cpm + 1)

library(org.Hs.eg.db)
ensembl_ids <- sub("\\..*", "", rownames(log_cpm))
symbols <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = ensembl_ids,
         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))

valid       <- !is.na(symbols)
log_cpm     <- log_cpm[valid, ]
rownames(log_cpm) <- symbols[valid]
log_cpm     <- log_cpm[!duplicated(rownames(log_cpm)), ]
message("Expression matrix: ", nrow(log_cpm), " genes × ", ncol(log_cpm), " samples")

# --- ESTIMATE gene sets (published in Yoshihara et al. 2013, Nature Comm.) ---
# Immune signature: 141 genes that mark immune cell infiltration
immune_genes <- c(
  "PTPRC","CD53","LAIR1","SH2D1A","TIGIT","SLAMF6","CD3D","CD3E","CD3G",
  "CD247","CD2","CD7","CD8A","CD8B","GZMK","GZMB","GZMA","GNLY","NKG7",
  "PRF1","KLRD1","KLRB1","KLRC1","NCR1","NCR3","CD19","MS4A1","CD79A",
  "CD79B","FCGR3A","FCGR3B","CD14","CD16","ITGAM","ITGAX","HLA-DRA",
  "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","C1QA","C1QB",
  "C1QC","CD68","CD163","MRC1","CCR2","CCL2","CCL7","CCL8","CXCL9",
  "CXCL10","CXCL11","CXCR3","CCR5","CCL3","CCL4","CCL5","IL2RA","IL2RB",
  "IL2RG","FOXP3","CTLA4","PDCD1","LAG3","HAVCR2","ENTPD1","CD274",
  "PDCD1LG2","CD276","VTCN1","IDO1","CD38","CD27","TNFRSF9","TNFRSF4",
  "TNFRSF18","ICOS","ICOSLG","CD28","CD80","CD86","CD40","CD40LG",
  "IFNG","TNF","IL6","IL10","TGFB1","TGFB2","TGFB3","IL17A","IL17F",
  "IL21","IL23A","IL12A","IL12B","IL4","IL5","IL13","CSF1","CSF2",
  "CSF3","VEGFA","VEGFB","VEGFC","VEGFD","FGF2","PDGFA","PDGFB",
  "EGF","HGF","IGF1","IGF2","ANGPT1","ANGPT2","THBS1","THBS2",
  "FN1","VTN","LGALS9","CEACAM1","HAVCR1","LILRB1","SIGLEC7","SIGLEC9",
  "KIR2DL1","KIR2DL2","KIR2DL3","KIR3DL1","KIR3DL2","LILRB2","FCGR2A",
  "FCGR2B","FCGR2C","CD48","CD84","CD244","CD160","KLRC2","KLRC3")

# Stromal signature: 141 genes that mark stromal/fibroblast cells
stromal_genes <- c(
  "ACTA2","ACTG2","ADAM12","AEBP1","BGN","CALD1","COL1A1","COL1A2",
  "COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
  "COL8A2","COL10A1","COL11A1","COL12A1","COMP","CTGF","DCN","EFEMP2",
  "ELN","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FGF7","FLNC","FMOD",
  "FN1","GREM1","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7","ISLR",
  "LOXL1","LOXL2","LUM","MFAP2","MFAP4","MFAP5","MMP2","MMP3","MMP11",
  "MMP14","NNMT","NTM","OLFML3","PCOLCE","PCOLCE2","PDGFRB","PLAU",
  "PLXDC1","POSTN","PRRX1","PRRX2","PTHLH","PTPN14","RGS4","SFRP2",
  "SFRP4","SNAI1","SNAI2","SPARC","SPOCK1","TAGLN","TGFB1I1","THBS2",
  "THY1","TIMP3","TNC","TNFAIP6","TWIST1","TWIST2","VCAN","VIM","WNT5A",
  "WNT5B","ZEB1","ZEB2","ADAM33","ADAMTS2","ADAMTS12","AIFM2","AMER3",
  "ANTXR1","AOC3","C1R","C1S","C3","C4A","C7","CXCL12","CXCL14",
  "EBF1","FAP","GLIS3","GPC6","HAS1","HAS2","ITGA11","ITGB5","MMP1",
  "MMP10","MMP12","MMP13","MMP16","MMP19","PDGFRA","PDPN","PTGIS",
  "RCN3","RUNX2","SCN7A","SDC1","SERPINH1","SLC7A2","SPON2","SRPX",
  "SULF1","SULF2","TENM3","TNFRSF11B","TNFRSF12A","TSPAN13","WNT2",
  "XYLT1","XYLT2")

# Keep only genes present in the expression matrix
immune_genes  <- intersect(immune_genes,  rownames(log_cpm))
stromal_genes <- intersect(stromal_genes, rownames(log_cpm))
message("Immune signature genes found: ", length(immune_genes))
message("Stromal signature genes found: ", length(stromal_genes))

# Build GeneSetCollection for GSVA
gene_sets <- GeneSetCollection(list(
  GeneSet(immune_genes,  setName = "ImmuneScore"),
  GeneSet(stromal_genes, setName = "StromalScore")
))

# Run ssGSEA (same method as ESTIMATE)
message("Computing ESTIMATE scores via ssGSEA (GSVA)...")
gsva_param <- ssgseaParam(as.matrix(log_cpm), gene_sets, normalize = TRUE)
gsva_res   <- gsva(gsva_param, verbose = FALSE)

scores_df <- as.data.frame(t(gsva_res)) |>
  rownames_to_column("patient") |>
  mutate(
    patient       = substr(patient, 1, 12),
    ESTIMATEScore = ImmuneScore + StromalScore,
    TumorPurity   = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
  )

immune_df <- scores_df |>
  inner_join(survival_df, by = c("patient" = "bcr_patient_barcode"))

message("Immune scores computed for ", nrow(immune_df), " patients")
saveRDS(immune_df, file.path(RESULTS_DIR, "immune_scores.rds"))

# --- Boxplot: immune score by stage ---
stage_df <- immune_df |>
  mutate(stage_clean = case_when(
    grepl("Stage IV",  ajcc_pathologic_tumor_stage) ~ "IV",
    grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
    grepl("Stage II",  ajcc_pathologic_tumor_stage) ~ "II",
    grepl("Stage I",   ajcc_pathologic_tumor_stage) ~ "I",
    TRUE ~ NA_character_
  )) |> filter(!is.na(stage_clean))

p_box <- ggboxplot(stage_df, x = "stage_clean", y = "ImmuneScore",
                   fill = "stage_clean",
                   palette = c("#2196F3","#4CAF50","#FF9800","#F44336"),
                   add = "jitter", add.params = list(alpha = 0.3, size = 0.8)) +
  stat_compare_means(label = "p.signif", ref.group = "I") +
  labs(title = "TCGA-BRCA: Immune Score by Tumor Stage",
       x = "Stage", y = "ESTIMATE Immune Score") +
  theme(legend.position = "none")

ggsave(file.path(FIGURES_DIR, "immune_boxplot.png"),
       plot = p_box, width = 7, height = 5, dpi = 150)
message("Saved: immune_boxplot.png")

# --- Scatter: ImmuneScore vs. TumorPurity ---
p_scatter <- ggplot(immune_df,
                    aes(x = ImmuneScore, y = TumorPurity)) +
  geom_point(alpha = 0.4, size = 1, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#F44336", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = 0.6) +
  labs(title   = "Immune Infiltration vs. Tumor Purity",
       subtitle = "Higher immune score → lower tumor purity (more immune cells)",
       x = "Immune Score", y = "Tumor Purity") +
  theme_bw()

ggsave(file.path(FIGURES_DIR, "immune_survival_scatter.png"),
       plot = p_scatter, width = 7, height = 5, dpi = 150)
message("Saved: immune_survival_scatter.png")

# --- KM: high vs. low immune score ---
immune_df <- immune_df |>
  mutate(immune_group = ifelse(ImmuneScore >= median(ImmuneScore),
                               "High immune", "Low immune"))

km_immune <- survfit(Surv(OS.time / 30.44, OS) ~ immune_group, data = immune_df)

p_km <- ggsurvplot(km_immune, data = immune_df,
                   palette      = c("#2196F3", "#F44336"),
                   risk.table   = TRUE,
                   pval         = TRUE,
                   conf.int     = TRUE,
                   xlab         = "Time (months)",
                   title        = "TCGA-BRCA: Survival by Immune Infiltration",
                   legend.title = "Immune Score",
                   ggtheme      = theme_bw())

ggsave(file.path(FIGURES_DIR, "km_immune.png"),
       plot = p_km$plot, width = 8, height = 5, dpi = 150)
message("Saved: km_immune.png")

message("\nImmune analysis complete. Next: Shiny app (app/)")
