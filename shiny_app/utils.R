# =============================================================================
# utils.R — Helper functions for Shiny app
# Project: Lung Cancer RNA-seq - Osimertinib Treatment Response
# =============================================================================

# Paths relative to shiny_app/ directory
data_dir   <- "../data/processed"
tables_dir <- "../results/tables"
enrich_dir <- "../results/enrichment"

# --- Load all data once at startup -------------------------------------------
load_all_data <- function() {
  list(
    res_df = read_csv(file.path(tables_dir, "deseq2_all_results.csv"),
                      show_col_types = FALSE) %>%
      mutate(
        neg_log10_padj = -log10(padj),
        significance   = case_when(
          padj < 0.05 & log2FoldChange >  1 ~ "Up",
          padj < 0.05 & log2FoldChange < -1 ~ "Down",
          TRUE                               ~ "NS"
        )
      ),

    meta = read_csv(file.path(data_dir, "metadata_clean.csv"),
                    show_col_types = FALSE) %>%
      mutate(condition = factor(condition, levels = c("control", "treated"))),

    hallmark = read_csv(file.path(enrich_dir, "gsea_hallmark_results.csv"),
                        show_col_types = FALSE) %>%
      mutate(
        pathway_clean = str_replace(pathway, "HALLMARK_", "") %>%
          str_replace_all("_", " ") %>%
          str_to_title(),
        direction = ifelse(NES > 0, "Activated", "Suppressed"),
        neg_log10_padj = -log10(padj)
      ) %>%
      arrange(padj),

    kegg = read_csv(file.path(enrich_dir, "gsea_kegg_results.csv"),
                    show_col_types = FALSE) %>%
      mutate(
        direction = ifelse(NES > 0, "Activated", "Suppressed"),
        neg_log10_padj = -log10(p.adjust)
      ) %>%
      arrange(p.adjust),

    norm_counts = read_csv(file.path(tables_dir, "deseq2_all_results.csv"),
                           show_col_types = FALSE) %>%
      dplyr::select(gene)  # placeholder — heatmap uses dds object
  )
}

# --- Volcano plot (plotly) ----------------------------------------------------
make_volcano <- function(res_df, padj_cut = 0.05, lfc_cut = 1, label_n = 15) {
  colors <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")

  top_genes <- res_df %>%
    filter(significance != "NS") %>%
    arrange(padj) %>%
    slice_head(n = label_n) %>%
    pull(gene)

  plot_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      label     = ifelse(gene %in% top_genes, gene, ""),
      hover_txt = paste0(
        "<b>", gene, "</b><br>",
        "Log2FC: ", round(log2FoldChange, 3), "<br>",
        "padj: ", signif(padj, 3), "<br>",
        "Status: ", significance
      )
    )

  plot_ly(plot_df,
    x    = ~log2FoldChange,
    y    = ~neg_log10_padj,
    type = "scatter",
    mode = "markers",
    color     = ~significance,
    colors    = colors,
    text      = ~hover_txt,
    hoverinfo = "text",
    marker    = list(size = 5, opacity = 0.6)
  ) %>%
    add_segments(x = -lfc_cut, xend = -lfc_cut,
                 y = 0, yend = max(plot_df$neg_log10_padj, na.rm = TRUE),
                 line = list(dash = "dash", color = "grey40", width = 1),
                 showlegend = FALSE, hoverinfo = "none") %>%
    add_segments(x = lfc_cut, xend = lfc_cut,
                 y = 0, yend = max(plot_df$neg_log10_padj, na.rm = TRUE),
                 line = list(dash = "dash", color = "grey40", width = 1),
                 showlegend = FALSE, hoverinfo = "none") %>%
    add_segments(x = min(plot_df$log2FoldChange, na.rm = TRUE),
                 xend = max(plot_df$log2FoldChange, na.rm = TRUE),
                 y = -log10(padj_cut), yend = -log10(padj_cut),
                 line = list(dash = "dash", color = "grey40", width = 1),
                 showlegend = FALSE, hoverinfo = "none") %>%
    layout(
      title  = list(text = "Volcano Plot: Osimertinib vs DMSO", font = list(size = 16)),
      xaxis  = list(title = "Log\u2082 Fold Change"),
      yaxis  = list(title = "-Log\u2081\u2080 Adjusted P-value"),
      legend = list(title = list(text = "Significance")),
      plot_bgcolor  = "white",
      paper_bgcolor = "white"
    )
}

# --- Hallmark dot plot (plotly) -----------------------------------------------
make_hallmark_plot <- function(hallmark_df, padj_cut = 0.05) {
  df <- hallmark_df %>%
    filter(padj < padj_cut) %>%
    arrange(NES)

  if (nrow(df) == 0) return(plotly_empty())

  plot_ly(df,
    x         = ~NES,
    y         = ~reorder(pathway_clean, NES),
    type      = "scatter",
    mode      = "markers",
    color     = ~direction,
    colors    = c("Activated" = "#d73027", "Suppressed" = "#4575b4"),
    size      = ~neg_log10_padj,
    sizes     = c(8, 30),
    text      = ~paste0("<b>", pathway_clean, "</b><br>NES: ", round(NES, 3),
                        "<br>padj: ", signif(padj, 3)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title  = list(text = "Hallmark Pathway Enrichment (GSEA)", font = list(size = 15)),
      xaxis  = list(title = "Normalized Enrichment Score (NES)",
                    zeroline = TRUE, zerolinecolor = "grey80"),
      yaxis  = list(title = ""),
      legend = list(title = list(text = "Direction")),
      plot_bgcolor  = "white",
      paper_bgcolor = "white",
      margin = list(l = 220)
    )
}

# --- KEGG bar plot (plotly) ---------------------------------------------------
make_kegg_plot <- function(kegg_df, padj_cut = 0.05, top_n = 20) {
  df <- kegg_df %>%
    filter(p.adjust < padj_cut) %>%
    arrange(NES) %>%
    dplyr::slice(c(1:min(top_n %/% 2, sum(NES < 0, na.rm = TRUE)),
                   (n() - min(top_n %/% 2, sum(NES > 0, na.rm = TRUE)) + 1):n()))

  if (nrow(df) == 0) return(plotly_empty())

  plot_ly(df,
    x         = ~NES,
    y         = ~reorder(Description, NES),
    type      = "bar",
    orientation = "h",
    color     = ~direction,
    colors    = c("Activated" = "#d73027", "Suppressed" = "#4575b4"),
    text      = ~paste0("<b>", Description, "</b><br>NES: ", round(NES, 3),
                        "<br>padj: ", signif(p.adjust, 3)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title  = list(text = "Top KEGG Pathway Enrichment (GSEA)", font = list(size = 15)),
      xaxis  = list(title = "Normalized Enrichment Score (NES)",
                    zeroline = TRUE, zerolinecolor = "grey80"),
      yaxis  = list(title = ""),
      legend = list(title = list(text = "Direction")),
      plot_bgcolor  = "white",
      paper_bgcolor = "white",
      margin = list(l = 250)
    )
}
