# =============================================================================
# Shiny Dashboard: Lung Cancer RNA-seq Osimertinib Treatment Response
# Dataset: GSE140941 | NCI-H1975 | Osimertinib vs DMSO
# Author: Matt Muslu
# =============================================================================

## Only run this a thte console below: shiny::runApp("C:/Users/amusl/Desktop/DrugDiscov/lung-cancer-rnaseq/shiny_app")

# NOTE: Run once in console if needed:
# install.packages(c("shiny", "shinydashboard", "DT", "plotly", "bslib"))

library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(plotly)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

source("utils.R")

# --- Load data once at startup ------------------------------------------------
dat <- load_all_data()

# Summary stats for Overview tab
n_total     <- nrow(dat$res_df)
n_up        <- sum(dat$res_df$significance == "Up",   na.rm = TRUE)
n_down      <- sum(dat$res_df$significance == "Down", na.rm = TRUE)
n_hallmark  <- sum(dat$hallmark$padj < 0.05, na.rm = TRUE)

# Normalized counts for heatmap (load dds)
dds <- readRDS("../data/processed/dds_object.rds")
norm_counts <- counts(dds, normalized = TRUE)

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span("Osimertinib Response", style = "font-size:16px; font-weight:bold;")
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview",          tabName = "overview",  icon = icon("house")),
      menuItem("DE Gene Explorer",  tabName = "volcano",   icon = icon("dna")),
      menuItem("Pathway Analysis",  tabName = "pathways",  icon = icon("diagram-project")),
      menuItem("Heatmap Viewer",    tabName = "heatmap",   icon = icon("table-cells")),
      menuItem("About",             tabName = "about",     icon = icon("circle-info"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f9f9f9; }
      .box { border-top: 3px solid #2c3e7a; }
      .value-box .inner h3 { font-size: 28px; }
      .skin-blue .main-header .logo { background-color: #1a2a5e; }
      .skin-blue .main-header .navbar { background-color: #1a2a5e; }
      .skin-blue .main-sidebar { background-color: #2c3e7a; }
    "))),

    tabItems(

      # -----------------------------------------------------------------------
      # TAB 1: OVERVIEW
      # -----------------------------------------------------------------------
      tabItem(tabName = "overview",
        fluidRow(
          valueBox(n_total,    "Genes Analyzed",      icon = icon("microscope"),  color = "navy"),
          valueBox(n_up,       "Upregulated Genes",   icon = icon("arrow-up"),    color = "red"),
          valueBox(n_down,     "Downregulated Genes", icon = icon("arrow-down"),  color = "blue")
        ),
        fluidRow(
          valueBox(n_hallmark, "Significant Hallmark Pathways", icon = icon("route"),      color = "purple"),
          valueBox("NCI-H1975","Cell Line",                     icon = icon("flask"),      color = "teal"),
          valueBox("GSE140941","GEO Accession",                 icon = icon("database"),   color = "olive")
        ),
        fluidRow(
          box(width = 12, title = "Project Summary", status = "primary", solidHeader = TRUE,
            h4("Transcriptomic Response to Osimertinib in EGFR-Mutant Lung Cancer"),
            p(HTML(
              "This dashboard presents results from a bulk RNA-seq reanalysis of NCI-H1975 lung
              adenocarcinoma cells (EGFR L858R + T790M) treated with <b>osimertinib</b>
              (1 µM, 12 hours) compared to vehicle control (DMSO). Data are from the publicly
              available dataset <b>GSE140941</b> (NCBI GEO)."
            )),
            br(),
            h4("Key Findings"),
            tags$ul(
              tags$li(HTML("<b>4,308 genes</b> significantly differentially expressed (padj < 0.05, |LFC| > 1)")),
              tags$li(HTML("Strong suppression of <b>G2M checkpoint, E2F targets, and MYC targets</b> — confirming antiproliferative effect")),
              tags$li(HTML("Activation of <b>p53 pathway</b> — pro-apoptotic signaling engaged")),
              tags$li(HTML("Suppression of <b>mTORC1 signaling</b> — downstream EGFR pathway blockade confirmed")),
              tags$li(HTML("Activation of <b>IL-6/JAK/STAT3</b> — candidate early adaptive resistance mechanism")),
              tags$li(HTML("Activation of <b>interferon alpha/gamma response</b> — potential immunotherapy combination relevance"))
            ),
            br(),
            h4("Analytical Workflow"),
            p("Data ingestion → Metadata curation → DESeq2 differential expression →
              LFC shrinkage (apeglm) → Gene annotation (org.Hs.eg.db) →
              GSEA (Hallmark + KEGG) → ORA → Visualization")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: DE GENE EXPLORER
      # -----------------------------------------------------------------------
      tabItem(tabName = "volcano",
        fluidRow(
          box(width = 3, title = "Filters", status = "primary", solidHeader = TRUE,
            sliderInput("lfc_cut",  "Min |Log2 FC|", min = 0, max = 5, value = 1, step = 0.25),
            sliderInput("padj_cut", "Max adj. p-value", min = 0.001, max = 0.2, value = 0.05, step = 0.005),
            sliderInput("label_n", "Genes to label", min = 0, max = 40, value = 15, step = 5),
            hr(),
            h5("DE Gene Summary"),
            verbatimTextOutput("de_summary")
          ),
          box(width = 9, title = "Volcano Plot", status = "primary", solidHeader = TRUE,
            plotlyOutput("volcano_plot", height = "500px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Differentially Expressed Genes Table",
              status = "primary", solidHeader = TRUE,
            p("Showing genes filtered by the sliders above. Click column headers to sort."),
            DTOutput("de_table")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: PATHWAY ANALYSIS
      # -----------------------------------------------------------------------
      tabItem(tabName = "pathways",
        fluidRow(
          tabBox(width = 12, title = "Pathway Enrichment Results",
            tabPanel("Hallmark GSEA",
              fluidRow(
                column(3,
                  sliderInput("hallmark_padj", "Max adj. p-value",
                              min = 0.001, max = 0.2, value = 0.05, step = 0.005),
                  radioButtons("hallmark_dir", "Show",
                               choices = c("All", "Activated", "Suppressed"), selected = "All")
                ),
                column(9,
                  plotlyOutput("hallmark_plot", height = "550px")
                )
              ),
              hr(),
              DTOutput("hallmark_table")
            ),
            tabPanel("KEGG GSEA",
              fluidRow(
                column(3,
                  sliderInput("kegg_padj", "Max adj. p-value",
                              min = 0.001, max = 0.2, value = 0.05, step = 0.005),
                  sliderInput("kegg_topn", "Top N pathways",
                              min = 5, max = 40, value = 20, step = 5)
                ),
                column(9,
                  plotlyOutput("kegg_plot", height = "550px")
                )
              ),
              hr(),
              DTOutput("kegg_table")
            )
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: HEATMAP VIEWER
      # -----------------------------------------------------------------------
      tabItem(tabName = "heatmap",
        fluidRow(
          box(width = 3, title = "Heatmap Controls", status = "primary", solidHeader = TRUE,
            sliderInput("heatmap_n",   "Number of top genes", min = 10, max = 100, value = 50, step = 10),
            radioButtons("heatmap_dir", "Gene selection",
                         choices = c("Top DE (both)", "Upregulated only", "Downregulated only"),
                         selected = "Top DE (both)"),
            sliderInput("heatmap_lfc", "Min |LFC|", min = 0, max = 4, value = 1, step = 0.5),
            hr(),
            p(style = "font-size:12px; color:grey;",
              "Heatmap shows Z-score scaled normalized counts across samples.")
          ),
          box(width = 9, title = "Top DE Genes Heatmap",
              status = "primary", solidHeader = TRUE,
            plotOutput("heatmap_plot", height = "650px")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: ABOUT
      # -----------------------------------------------------------------------
      tabItem(tabName = "about",
        fluidRow(
          box(width = 12, title = "About This Project", status = "primary", solidHeader = TRUE,
            h4("Project Overview"),
            p(HTML(
              "This dashboard is part of a computational drug discovery portfolio project demonstrating
              RNA-seq analysis skills relevant to pharmaceutical research. The analysis characterizes
              the transcriptomic response to osimertinib — a third-generation EGFR tyrosine kinase
              inhibitor — in EGFR-mutant NCI-H1975 non-small cell lung cancer cells."
            )),
            hr(),
            h4("Dataset"),
            tags$ul(
              tags$li(HTML("<b>GEO Accession:</b> GSE140941")),
              tags$li(HTML("<b>Cell Line:</b> NCI-H1975 (EGFR L858R + T790M)")),
              tags$li(HTML("<b>Treatment:</b> Osimertinib 1 µM, 12 hours vs DMSO vehicle")),
              tags$li(HTML("<b>Replicates:</b> n = 4 per condition")),
              tags$li(HTML("<b>Sequencing:</b> Bulk RNA-seq, aligned to GRCh38.p12"))
            ),
            hr(),
            h4("Analytical Methods"),
            tags$ul(
              tags$li(HTML("<b>DE Analysis:</b> DESeq2 v1.46 with apeglm LFC shrinkage")),
              tags$li(HTML("<b>Thresholds:</b> padj < 0.05, |Log2FC| > 1")),
              tags$li(HTML("<b>GSEA:</b> fgsea v1.32 — Hallmark (MSigDB); gseKEGG (clusterProfiler v4.14)")),
              tags$li(HTML("<b>ORA:</b> enrichKEGG (clusterProfiler)")),
              tags$li(HTML("<b>Gene annotation:</b> org.Hs.eg.db"))
            ),
            hr(),
            h4("Key References"),
            tags$ol(
              tags$li("Love, M. I., Huber, W., & Anders, S. (2014). DESeq2. Genome Biology, 15(12), 550."),
              tags$li("Wu, T., et al. (2021). clusterProfiler 4.0. The Innovation, 2(3), 100141."),
              tags$li("Soria, J.-C., et al. (2018). Osimertinib in EGFR-mutated NSCLC. NEJM, 378(2), 113–125."),
              tags$li("Korotkevich, G., et al. (2016). fgsea. bioRxiv.")
            ),
            hr(),
            h4("Author & Repository"),
            p(HTML(
              "Author: <b>Matt Muslu</b><br>
               GitHub: <a href='https://github.com/amuslu87' target='_blank'>github.com/amuslu87</a><br>
               All analysis code is available in the project repository."
            ))
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # --- Filtered DE data -------------------------------------------------------
  filtered_de <- reactive({
    dat$res_df %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(
        significance = case_when(
          padj < input$padj_cut & log2FoldChange >  input$lfc_cut ~ "Up",
          padj < input$padj_cut & log2FoldChange < -input$lfc_cut ~ "Down",
          TRUE ~ "NS"
        ),
        neg_log10_padj = -log10(padj)
      )
  })

  # --- DE Summary text --------------------------------------------------------
  output$de_summary <- renderText({
    df <- filtered_de()
    paste0(
      "Upregulated:   ", sum(df$significance == "Up"), "\n",
      "Downregulated: ", sum(df$significance == "Down"), "\n",
      "Not significant: ", sum(df$significance == "NS")
    )
  })

  # --- Volcano plot -----------------------------------------------------------
  output$volcano_plot <- renderPlotly({
    make_volcano(filtered_de(), input$padj_cut, input$lfc_cut, input$label_n)
  })

  # --- DE table ---------------------------------------------------------------
  output$de_table <- renderDT({
    filtered_de() %>%
      filter(significance != "NS") %>%
      arrange(padj) %>%
      dplyr::select(Gene = gene, Log2FC = log2FoldChange, AdjPval = padj,
                    BaseMean = baseMean, Status = significance, Description = description) %>%
      mutate(
        Log2FC  = round(Log2FC, 3),
        AdjPval = signif(AdjPval, 3),
        BaseMean = round(BaseMean, 1)
      ) %>%
      datatable(
        filter    = "top",
        rownames  = FALSE,
        options   = list(pageLength = 15, scrollX = TRUE),
        class     = "cell-border stripe"
      ) %>%
      formatStyle("Status",
        color            = styleEqual(c("Up", "Down"), c("#d73027", "#4575b4")),
        fontWeight       = "bold"
      )
  })

  # --- Hallmark plot ----------------------------------------------------------
  output$hallmark_plot <- renderPlotly({
    df <- dat$hallmark
    if (input$hallmark_dir != "All")
      df <- df %>% filter(direction == input$hallmark_dir)
    make_hallmark_plot(df, input$hallmark_padj)
  })

  # --- Hallmark table ---------------------------------------------------------
  output$hallmark_table <- renderDT({
    df <- dat$hallmark
    if (input$hallmark_dir != "All")
      df <- df %>% filter(direction == input$hallmark_dir)
    df %>%
      filter(padj < input$hallmark_padj) %>%
      dplyr::select(Pathway = pathway_clean, NES, `Adj P-val` = padj,
                    Size = size, Direction = direction) %>%
      mutate(NES = round(NES, 3), `Adj P-val` = signif(`Adj P-val`, 3)) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
                class = "cell-border stripe") %>%
      formatStyle("Direction",
        color      = styleEqual(c("Activated", "Suppressed"), c("#d73027", "#4575b4")),
        fontWeight = "bold"
      )
  })

  # --- KEGG plot --------------------------------------------------------------
  output$kegg_plot <- renderPlotly({
    make_kegg_plot(dat$kegg, input$kegg_padj, input$kegg_topn)
  })

  # --- KEGG table -------------------------------------------------------------
  output$kegg_table <- renderDT({
    dat$kegg %>%
      filter(p.adjust < input$kegg_padj) %>%
      dplyr::select(Pathway = Description, NES, `Adj P-val` = p.adjust,
                    Size = setSize, Direction = direction) %>%
      mutate(NES = round(NES, 3), `Adj P-val` = signif(`Adj P-val`, 3)) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE),
                class = "cell-border stripe") %>%
      formatStyle("Direction",
        color      = styleEqual(c("Activated", "Suppressed"), c("#d73027", "#4575b4")),
        fontWeight = "bold"
      )
  })

  # --- Heatmap ----------------------------------------------------------------
  output$heatmap_plot <- renderPlot({
    n        <- input$heatmap_n
    lfc_cut  <- input$heatmap_lfc

    gene_pool <- dat$res_df %>%
      filter(padj < 0.05, abs(log2FoldChange) >= lfc_cut)

    if (input$heatmap_dir == "Upregulated only")
      gene_pool <- gene_pool %>% filter(log2FoldChange > 0)
    else if (input$heatmap_dir == "Downregulated only")
      gene_pool <- gene_pool %>% filter(log2FoldChange < 0)

    top_genes <- gene_pool %>%
      arrange(padj) %>%
      slice_head(n = n) %>%
      pull(gene)

    req(length(top_genes) >= 2)

    mat <- norm_counts[intersect(top_genes, rownames(norm_counts)), ]
    req(nrow(mat) >= 2)

    mat_scaled <- t(scale(t(mat)))

    col_ann <- data.frame(
      Condition = dat$meta$condition,
      row.names = dat$meta$sample_id
    )

    ann_colors <- list(
      Condition = c("control" = "#4575b4", "treated" = "#d73027")
    )

    p <- pheatmap(
      mat_scaled,
      annotation_col    = col_ann,
      annotation_colors = ann_colors,
      color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
      cluster_rows      = TRUE,
      cluster_cols      = TRUE,
      show_rownames     = nrow(mat) <= 60,
      show_colnames     = TRUE,
      fontsize_row      = max(5, 10 - nrow(mat) %/% 10),
      fontsize_col      = 9,
      main              = paste0("Top ", nrow(mat), " DE Genes (Z-score)"),
      border_color      = NA,
      silent            = TRUE
    )
    print(p)
  })

}

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)



