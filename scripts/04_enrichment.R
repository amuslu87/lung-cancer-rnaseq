# =============================================================================
# Script 04: Pathway Enrichment Analysis (GSEA + ORA)
# Project: Lung Cancer RNA-seq - Osimertinib Treatment Response
# Dataset: GSE140941 (NCI-H1975, safe genotype, osimertinib vs DMSO)
# Author: Matt Muslu
# Date: 2026-04-06
# =============================================================================

# NOTE: Run once in console if needed:
# BiocManager::install(c("clusterProfiler", "fgsea", "enrichplot"))
# install.packages("msigdbr")

# --- Load libraries -----------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(org.Hs.eg.db)

# --- Set paths ----------------------------------------------------------------
project_dir  <- "C:/Users/amusl/Desktop/DrugDiscov/lung-cancer-rnaseq"
tables_dir   <- file.path(project_dir, "results", "tables")
figures_dir  <- file.path(project_dir, "results", "figures")
enrich_dir   <- file.path(project_dir, "results", "enrichment")

# --- Load DE results ----------------------------------------------------------
res_df <- read_csv(file.path(tables_dir, "deseq2_all_results.csv")) %>%
  mutate(entrez_id = as.character(as.integer(entrez_id)))  # clean: "7157" not "7157.0"

# =============================================================================
# PART 1: GSEA — Hallmark (fgsea + msigdbr)
# =============================================================================

# --- Build ranked gene list (Entrez IDs as names) -----------------------------
ranked_genes <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(entrez_id)) %>%
  mutate(rank_score = sign(log2FoldChange) * -log10(padj + 1e-300)) %>%
  arrange(desc(rank_score)) %>%
  dplyr::select(entrez_id, rank_score) %>%
  deframe()

cat("Ranked gene list length:", length(ranked_genes), "\n")

# --- Load Hallmark gene sets (msigdbr v10 API) --------------------------------
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, ncbi_gene) %>%
  mutate(ncbi_gene = as.character(ncbi_gene))

hallmark_list <- split(hallmark_sets$ncbi_gene, hallmark_sets$gs_name)
cat("Hallmark gene sets loaded:", length(hallmark_list), "\n")

# --- Run GSEA on Hallmark -----------------------------------------------------
set.seed(42)
gsea_hallmark <- fgsea(
  pathways    = hallmark_list,
  stats       = ranked_genes,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = 10000,
  nproc       = 1  # disable parallel processing (Windows compatibility)
)

gsea_hallmark_clean <- gsea_hallmark %>%
  arrange(padj) %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed"))

cat("\n=== Significant Hallmark Pathways (padj < 0.05) ===\n")
cat("Activated: ", sum(gsea_hallmark_clean$direction == "Activated"), "\n")
cat("Suppressed:", sum(gsea_hallmark_clean$direction == "Suppressed"), "\n")
print(gsea_hallmark_clean %>% dplyr::select(pathway, NES, padj, direction) %>% head(20))

# =============================================================================
# PART 2: GSEA — KEGG (gseKEGG from clusterProfiler)
# =============================================================================

BiocParallel::register(BiocParallel::SerialParam())  # Windows fix
set.seed(42)
gsea_kegg <- gseKEGG(
  geneList     = ranked_genes,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose      = FALSE
)

gsea_kegg_clean <- gsea_kegg@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  mutate(direction = ifelse(NES > 0, "Activated", "Suppressed"))

cat("\n=== Significant KEGG Pathways (padj < 0.05) ===\n")
cat("Activated: ", sum(gsea_kegg_clean$direction == "Activated"), "\n")
cat("Suppressed:", sum(gsea_kegg_clean$direction == "Suppressed"), "\n")
print(gsea_kegg_clean %>% dplyr::select(Description, NES, p.adjust, direction) %>% head(20))

# =============================================================================
# PART 3: ORA (Over-Representation Analysis on significant DE genes)
# =============================================================================

sig_entrez_up <- res_df %>%
  filter(padj < 0.05, log2FoldChange > 1, !is.na(entrez_id)) %>%
  pull(entrez_id)

sig_entrez_down <- res_df %>%
  filter(padj < 0.05, log2FoldChange < -1, !is.na(entrez_id)) %>%
  pull(entrez_id)

universe_entrez <- res_df %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id)

ora_up <- enrichKEGG(
  gene          = sig_entrez_up,
  universe      = universe_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

ora_down <- enrichKEGG(
  gene          = sig_entrez_down,
  universe      = universe_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

cat("\n=== ORA: KEGG pathways enriched in upregulated genes ===\n")
if (!is.null(ora_up) && nrow(ora_up@result) > 0) {
  print(ora_up@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% head(10))
} else {
  cat("No significant pathways found.\n")
}

cat("\n=== ORA: KEGG pathways enriched in downregulated genes ===\n")
if (!is.null(ora_down) && nrow(ora_down@result) > 0) {
  print(ora_down@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% head(10))
} else {
  cat("No significant pathways found.\n")
}

# =============================================================================
# FIGURES
# =============================================================================

# --- Figure 5: Hallmark GSEA dot plot -----------------------------------------
top_hallmark <- gsea_hallmark %>%
  filter(padj < 0.05) %>%
  arrange(NES) %>%
  mutate(
    pathway_clean = str_replace(pathway, "HALLMARK_", "") %>% str_replace_all("_", " "),
    direction     = ifelse(NES > 0, "Activated", "Suppressed")
  )

if (nrow(top_hallmark) > 0) {
  hallmark_plot <- ggplot(top_hallmark,
    aes(x = NES, y = reorder(pathway_clean, NES), fill = direction, size = -log10(padj))) +
    geom_point(shape = 21, color = "grey30") +
    scale_fill_manual(values = c("Activated" = "#d73027", "Suppressed" = "#4575b4")) +
    scale_size_continuous(range = c(3, 9), name = "-log10(padj)") +
    labs(
      title    = "Hallmark Pathway Enrichment (GSEA)",
      subtitle = "Osimertinib vs DMSO | NCI-H1975",
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      fill     = "Direction"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(figures_dir, "05_gsea_hallmark_dotplot.pdf"), hallmark_plot,
         width = 10, height = max(5, nrow(top_hallmark) * 0.35 + 2))
  ggsave(file.path(figures_dir, "05_gsea_hallmark_dotplot.png"), hallmark_plot,
         width = 10, height = max(5, nrow(top_hallmark) * 0.35 + 2), dpi = 150)
  cat("Saved: 05_gsea_hallmark_dotplot\n")
}

# --- Figure 6: KEGG GSEA bar plot ---------------------------------------------
if (nrow(gsea_kegg_clean) > 0) {
  top_kegg_plot <- gsea_kegg_clean %>%
    arrange(NES) %>%
    slice(c(1:min(10, sum(NES < 0)), (n() - min(10, sum(NES > 0)) + 1):n()))

  kegg_plot <- ggplot(top_kegg_plot,
    aes(x = NES, y = reorder(Description, NES), fill = direction)) +
    geom_col(color = "grey30", width = 0.7) +
    scale_fill_manual(values = c("Activated" = "#d73027", "Suppressed" = "#4575b4")) +
    labs(
      title    = "Top KEGG Pathway Enrichment (GSEA)",
      subtitle = "Osimertinib vs DMSO | NCI-H1975",
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      fill     = "Direction"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(figures_dir, "06_gsea_kegg_barplot.pdf"), kegg_plot,
         width = 10, height = max(5, nrow(top_kegg_plot) * 0.4 + 2))
  ggsave(file.path(figures_dir, "06_gsea_kegg_barplot.png"), kegg_plot,
         width = 10, height = max(5, nrow(top_kegg_plot) * 0.4 + 2), dpi = 150)
  cat("Saved: 06_gsea_kegg_barplot\n")
} else {
  cat("No significant KEGG pathways — skipping figure 6.\n")
}

# =============================================================================
# SAVE RESULTS TABLES
# =============================================================================

gsea_hallmark %>%
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) %>%
  arrange(padj) %>%
  write_csv(file.path(enrich_dir, "gsea_hallmark_results.csv"))

gsea_kegg@result %>%
  arrange(p.adjust) %>%
  write_csv(file.path(enrich_dir, "gsea_kegg_results.csv"))

if (!is.null(ora_up) && nrow(ora_up@result) > 0)
  write_csv(ora_up@result,   file.path(enrich_dir, "ora_kegg_upregulated.csv"))
if (!is.null(ora_down) && nrow(ora_down@result) > 0)
  write_csv(ora_down@result, file.path(enrich_dir, "ora_kegg_downregulated.csv"))

cat("\nSaved enrichment tables to:", enrich_dir, "\n")
cat("\nScript 04 complete.\n")
