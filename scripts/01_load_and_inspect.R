# =============================================================================
# Script 01: Data Loading and Inspection
# Project: Lung Cancer RNA-seq - Osimertinib Treatment Response
# Dataset: GSE140941 (NCI-H1975, safe genotype, osimertinib vs DMSO)
# Author: Matt Muslu
# Date: 2026-04-05
# =============================================================================

# NOTE: Run this once in console if needed:

# install.packages("rlang", type = "binary")

install.packages(c("tidyverse", "readr", "janitor"))
# --- Load libraries -----------------------------------------------------------
library(tidyverse)
library(readr)
library(janitor)
library(rlang)
# --- Set paths ----------------------------------------------------------------
project_dir <- "C:/Users/amusl/Desktop/DrugDiscov/lung-cancer-rnaseq"
raw_dir     <- file.path(project_dir, "data", "raw")
meta_dir    <- file.path(project_dir, "data", "metadata")
proc_dir    <- file.path(project_dir, "data", "processed")

# --- Load metadata ------------------------------------------------------------
meta <- read_csv(file.path(meta_dir, "metadata.csv")) %>%
  clean_names() %>%
  mutate(condition = factor(condition, levels = c("control", "treated")))

cat("=== Metadata ===\n")
print(meta)
cat("\nSample counts per condition:\n")
print(table(meta$condition))

# --- Load count files and merge into one matrix -------------------------------
# Each file has 7 columns: Geneid, Chr, Start, End, Strand, Length, counts
# We only need columns 1 (Geneid = gene name) and 7 (counts)

load_counts <- function(filepath, sample_id) {
  read_tsv(filepath, comment = "#", col_types = cols(.default = "c")) %>%
    select(gene = 1, count = last_col()) %>%
    mutate(count = as.integer(count),
           sample = sample_id)
}

cat("\nLoading count files...\n")
counts_long <- map2_dfr(
  file.path(raw_dir, meta$filename),
  meta$sample_id,
  load_counts
)

# Pivot to wide matrix: rows = genes, columns = samples
count_matrix <- counts_long %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  column_to_rownames("gene")

cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# --- Basic QC -----------------------------------------------------------------

# 1. Library sizes (total counts per sample)
lib_sizes <- colSums(count_matrix)
cat("\n=== Library Sizes (total counts per sample) ===\n")
print(lib_sizes)

# 2. How many genes have zero counts in all samples?
all_zero <- sum(rowSums(count_matrix) == 0)
cat("\nGenes with zero counts in all samples:", all_zero, "\n")

# 3. Quick check: are sample names in matrix matching metadata order?
stopifnot(all(colnames(count_matrix) == meta$sample_id))
cat("Sample order check passed.\n")

# --- Save processed count matrix and metadata ---------------------------------
write_csv(
  count_matrix %>% rownames_to_column("gene"),
  file.path(proc_dir, "count_matrix_raw.csv")
)
write_csv(meta, file.path(proc_dir, "metadata_clean.csv"))

cat("\nSaved:\n")
cat(" -", file.path(proc_dir, "count_matrix_raw.csv"), "\n")
cat(" -", file.path(proc_dir, "metadata_clean.csv"), "\n")
cat("\nScript 01 complete.\n")

