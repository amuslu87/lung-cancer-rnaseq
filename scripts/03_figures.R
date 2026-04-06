# =============================================================================
# Script 03: Visualization
# Project: Lung Cancer RNA-seq - Osimertinib Treatment Response
# Dataset: GSE140941 (NCI-H1975, safe genotype, osimertinib vs DMSO)
# Author: Matt Muslu
# Date: 2026-04-06
# =============================================================================

# NOTE: Run once in console if needed:
BiocManager::install("EnhancedVolcano")
install.packages("pheatmap")

# --- Load libraries -----------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# --- Set paths ----------------------------------------------------------------
project_dir <- "C:/Users/amusl/Desktop/DrugDiscov/lung-cancer-rnaseq"
proc_dir    <- file.path(project_dir, "data", "processed")
tables_dir  <- file.path(project_dir, "results", "tables")
figures_dir <- file.path(project_dir, "results", "figures")

# --- Load data ----------------------------------------------------------------
dds    <- readRDS(file.path(proc_dir, "dds_object.rds"))
res_df <- read_csv(file.path(tables_dir, "deseq2_all_results.csv"))
meta   <- read_csv(file.path(proc_dir, "metadata_clean.csv")) %>%
  mutate(condition = factor(condition, levels = c("control", "treated")))

# =============================================================================
# FIGURE 1: PCA Plot
# Purpose: Show that treatment groups separate cleanly — validates the experiment
# =============================================================================

# Variance-stabilizing transformation for PCA (better than raw counts)
vsd <- vst(dds, blind = FALSE)

# Compute PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(vjust = -0.8, size = 3.2, show.legend = FALSE) +
  scale_color_manual(
    values = c("control" = "#4575b4", "treated" = "#d73027"),
    labels = c("control" = "DMSO (Control)", "treated" = "Osimertinib")
  ) +
  labs(
    title    = "PCA: Osimertinib vs DMSO (NCI-H1975)",
    subtitle = "Variance-stabilized counts",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)"),
    color    = "Condition"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "01_pca_plot.pdf"), pca_plot, width = 7, height = 5)
ggsave(file.path(figures_dir, "01_pca_plot.png"), pca_plot, width = 7, height = 5, dpi = 150)
cat("Saved: 01_pca_plot\n")

# =============================================================================
# FIGURE 2: Volcano Plot
# Purpose: Show the full landscape of differential expression
# =============================================================================

volcano_plot <- EnhancedVolcano(
  res_df,
  lab            = res_df$gene,
  x              = "log2FoldChange",
  y              = "padj",
  pCutoff        = 0.05,
  FCcutoff       = 1,
  title          = "Osimertinib vs DMSO — NCI-H1975",
  subtitle       = "DESeq2 | padj < 0.05 | |LFC| > 1",
  xlab           = "Log2 Fold Change",
  ylab           = "-Log10 Adjusted P-value",
  pointSize      = 1.5,
  labSize        = 3.0,
  col            = c("grey70", "grey70", "#4575b4", "#d73027"),
  colAlpha       = 0.6,
  legendPosition = "right",
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  max.overlaps   = 20
)

ggsave(file.path(figures_dir, "02_volcano_plot.pdf"), volcano_plot, width = 9, height = 7)
ggsave(file.path(figures_dir, "02_volcano_plot.png"), volcano_plot, width = 9, height = 7, dpi = 150)
cat("Saved: 02_volcano_plot\n")

# =============================================================================
# FIGURE 3: Heatmap of Top DE Genes
# Purpose: Show expression patterns across samples for the most significant genes
# =============================================================================

# Select top 50 DE genes by adjusted p-value (|LFC| > 1)
top_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene)

# Extract normalized counts for top genes
norm_counts <- counts(dds, normalized = TRUE)
heatmap_mat <- norm_counts[top_genes, ]

# Z-score scale each gene (row) so patterns are visible across samples
heatmap_scaled <- t(scale(t(heatmap_mat)))

# Annotation for columns (samples)
col_annotation <- data.frame(
  Condition = meta$condition,
  row.names = meta$sample_id
)

# Color scheme
ann_colors <- list(
  Condition = c("control" = "#4575b4", "treated" = "#d73027")
)

# Save heatmap as PDF
pdf(file.path(figures_dir, "03_heatmap_top50_genes.pdf"), width = 8, height = 12)
pheatmap(
  heatmap_scaled,
  annotation_col  = col_annotation,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 7,
  fontsize_col    = 9,
  main            = "Top 50 DE Genes: Osimertinib vs DMSO (Z-score)",
  border_color    = NA
)
dev.off()

# Save as PNG
png(file.path(figures_dir, "03_heatmap_top50_genes.png"), width = 800, height = 1200, res = 120)
pheatmap(
  heatmap_scaled,
  annotation_col  = col_annotation,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 7,
  fontsize_col    = 9,
  main            = "Top 50 DE Genes: Osimertinib vs DMSO (Z-score)",
  border_color    = NA
)
dev.off()
cat("Saved: 03_heatmap_top50_genes\n")

# =============================================================================
# FIGURE 4: Sample Distance Heatmap
# Purpose: Show overall similarity between samples — confirms replicates cluster
# =============================================================================

sample_dists <- dist(t(assay(vsd)))
dist_matrix  <- as.matrix(sample_dists)

# Color: white = similar, dark blue = dissimilar
dist_colors <- colorRampPalette(brewer.pal(9, "Blues"))(100)

pdf(file.path(figures_dir, "04_sample_distance_heatmap.pdf"), width = 6, height = 5)
pheatmap(
  dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color       = dist_colors,
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  main        = "Sample-to-Sample Distances (VST counts)",
  border_color = NA,
  fontsize     = 9
)
dev.off()

png(file.path(figures_dir, "04_sample_distance_heatmap.png"), width = 700, height = 600, res = 120)
pheatmap(
  dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color       = dist_colors,
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  main        = "Sample-to-Sample Distances (VST counts)",
  border_color = NA,
  fontsize     = 9
)
dev.off()
cat("Saved: 04_sample_distance_heatmap\n")

cat("\nAll figures saved to:", figures_dir, "\n")
cat("\nScript 03 complete.\n")
