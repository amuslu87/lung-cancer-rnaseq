# =============================================================================
# Script 02: Differential Expression Analysis with DESeq2
# Project: Lung Cancer RNA-seq - Osimertinib Treatment Response
# Dataset: GSE140941 (NCI-H1975, safe genotype, osimertinib vs DMSO)
# Author: Matt Muslu
# Date: 2026-04-06
# =============================================================================

# NOTE: Run once in console if needed:
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db"))

# --- Load libraries -----------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- Set paths ----------------------------------------------------------------
project_dir <- "C:/Users/amusl/Desktop/DrugDiscov/lung-cancer-rnaseq"
proc_dir    <- file.path(project_dir, "data", "processed")
tables_dir  <- file.path(project_dir, "results", "tables")

# --- Load data ----------------------------------------------------------------
count_matrix <- read_csv(file.path(proc_dir, "count_matrix_raw.csv")) %>%
  column_to_rownames("gene")

meta <- read_csv(file.path(proc_dir, "metadata_clean.csv")) %>%
  mutate(condition = factor(condition, levels = c("control", "treated")))

# Confirm sample order matches between matrix and metadata
stopifnot(all(colnames(count_matrix) == meta$sample_id))
cat("Sample order check passed.\n")

# --- Filter low-count genes ---------------------------------------------------
# Keep genes with at least 10 total counts across all samples
# This removes noise from very lowly expressed genes before DESeq2
keep <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep, ]
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")

# --- Create DESeq2 object -----------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData   = meta,
  design    = ~ condition
)

# --- Run DESeq2 ---------------------------------------------------------------
dds <- DESeq(dds)
cat("\nDESeq2 size factors (library normalization):\n")
print(sizeFactors(dds))

# --- Extract results ----------------------------------------------------------
# contrast: treated vs control (positive logFC = higher in osimertinib)
res <- results(dds,
               contrast  = c("condition", "treated", "control"),
               alpha     = 0.05)

cat("\n=== DESeq2 Results Summary ===\n")
summary(res)

# --- Apply log fold change shrinkage (apeglm) ---------------------------------
# Shrinkage reduces noise for low-count genes — makes volcano plots cleaner
res_shrunk <- lfcShrink(dds,
                        coef = "condition_treated_vs_control",
                        type = "apeglm")

cat("\n=== Results after LFC shrinkage ===\n")
summary(res_shrunk)

# --- Annotate with gene descriptions ------------------------------------------
# Map gene symbols to Entrez IDs and gene descriptions using org.Hs.eg.db
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Add Entrez ID
res_df$entrez_id <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$gene,
  column  = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Add gene description
res_df$description <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$gene,
  column  = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)

# --- Summarize DE genes -------------------------------------------------------
sig_up   <- res_df %>% filter(padj < 0.05, log2FoldChange > 1)
sig_down <- res_df %>% filter(padj < 0.05, log2FoldChange < -1)

cat("\n=== Differentially Expressed Genes (padj < 0.05, |LFC| > 1) ===\n")
cat("Upregulated in osimertinib:  ", nrow(sig_up), "\n")
cat("Downregulated in osimertinib:", nrow(sig_down), "\n")

cat("\nTop 10 upregulated genes:\n")
print(sig_up %>% dplyr::select(gene, log2FoldChange, padj, description) %>% head(10))

cat("\nTop 10 downregulated genes:\n")
print(sig_down %>% dplyr::select(gene, log2FoldChange, padj, description) %>% head(10))

# --- Save outputs -------------------------------------------------------------
# Full results table (all genes)
write_csv(res_df, file.path(tables_dir, "deseq2_all_results.csv"))

# Significant genes only
write_csv(
  res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1),
  file.path(tables_dir, "deseq2_significant_genes.csv")
)

# Save the DESeq2 object for use in downstream scripts
saveRDS(dds, file.path(proc_dir, "dds_object.rds"))

cat("\nSaved:\n")
cat(" -", file.path(tables_dir, "deseq2_all_results.csv"), "\n")
cat(" -", file.path(tables_dir, "deseq2_significant_genes.csv"), "\n")
cat(" -", file.path(proc_dir, "dds_object.rds"), "\n")
cat("\nScript 02 complete.\n")
