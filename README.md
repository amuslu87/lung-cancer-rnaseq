# Transcriptomic Response to Osimertinib in EGFR-Mutant Lung Cancer

A reproducible bulk RNA-seq analysis characterizing the transcriptional response to osimertinib treatment in NCI-H1975 lung cancer cells, framed in a drug discovery context. Includes differential expression, pathway enrichment with Shiny app.

---

## Background

Osimertinib is a third-generation EGFR tyrosine kinase inhibitor approved for EGFR-mutant non-small cell lung cancer (NSCLC). Understanding the transcriptomic consequences of EGFR inhibition is essential for identifying pharmacodynamic biomarkers, anticipating resistance mechanisms, and designing rational combination therapies.

**Scientific question:** What transcriptomic pathways are altered in EGFR-mutant lung cancer cells after osimertinib treatment, and which response patterns are relevant to mechanism-of-action interpretation, pharmacodynamic signaling, and candidate biomarker discovery?

---

## Dataset

| Field | Value |
|---|---|
| GEO Accession | [GSE140941](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140941) |
| Cell Line | NCI-H1975 (EGFR L858R + T790M) |
| Treatment | Osimertinib 1 µM, 12 hours |
| Control | DMSO vehicle |
| Replicates | n = 4 per condition |
| Sequencing | Bulk RNA-seq, featureCounts, GRCh38.p12 / GENCODE v28 |

---

## Key Findings

- **4,308 genes** significantly differentially expressed (padj < 0.05, |log₂FC| > 1)
- **Cell cycle arrest confirmed:** Strong suppression of G2M checkpoint, E2F targets, MYC targets, and mitotic spindle gene sets
- **On-target pathway engagement:** mTORC1 signaling suppressed — downstream EGFR blockade confirmed
- **Pro-apoptotic signaling:** p53 pathway activated
- **Candidate resistance mechanism:** IL-6/JAK/STAT3 activated at 12h — early adaptive bypass signaling
- **Immunotherapy relevance:** Interferon alpha/gamma response activated — potential combination opportunity with checkpoint inhibitors

---

## Project Structure

```
lung-cancer-rnaseq/
├── data/
│   ├── raw/              # Raw featureCounts files (excluded from repo — see below)
│   ├── processed/        # Merged count matrix, metadata, DESeq2 object
│   └── metadata/         # Sample annotation CSV
├── scripts/
│   ├── 01_load_and_inspect.R    # Data ingestion and QC
│   ├── 02_deseq2_analysis.R     # Differential expression (DESeq2 + apeglm)
│   ├── 03_figures.R             # PCA, volcano, heatmaps
│   └── 04_enrichment.R          # GSEA (Hallmark + KEGG) and ORA
├── results/
│   ├── figures/          # PCA, volcano, heatmap, sample distance plots
│   ├── tables/           # DESeq2 results CSVs
│   ├── enrichment/       # GSEA and ORA result tables
│   └── findings_summary_report.Rmd  # Full findings report (APA style)
├── shiny_app/
│   ├── app.R             # Interactive dashboard
│   └── utils.R           # Data loading and plot helper functions
├── docs/                 # Citations and dataset metadata
└── README.md
```

---

## Analytical Workflow

```
Raw counts (featureCounts)
    ↓
Data loading + QC (script 01)
    ↓
DESeq2 normalization + LFC shrinkage (script 02)
    ↓
Visualization: PCA, volcano, heatmap (script 03)
    ↓
GSEA: Hallmark + KEGG | ORA (script 04)
    ↓
Interactive Shiny dashboard
```

---

## Interactive Dashboard

The Shiny dashboard provides interactive exploration of all results:

- **Overview** — project summary and key statistics
- **DE Gene Explorer** — interactive volcano plot + filterable gene table
- **Pathway Analysis** — Hallmark and KEGG GSEA dot/bar plots + tables
- **Heatmap Viewer** — dynamic heatmap with gene count and direction controls
- **About** — methods, references, reproducibility

To launch locally:

```r
shiny::runApp("shiny_app")
```

---

## How to Reproduce

### 1. Clone the repository

```bash
git clone https://github.com/amuslu87/lung-cancer-rnaseq-portfolio.git
cd lung-cancer-rnaseq-portfolio
```

### 2. Download raw data

Download `GSE140941_RAW.tar` from [GEO GSE140941](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140941) and extract into `data/raw/`.

### 3. Install R dependencies

```r
install.packages(c("tidyverse", "janitor", "pheatmap", "RColorBrewer",
                   "msigdbr", "shiny", "shinydashboard", "DT", "plotly"))

BiocManager::install(c("DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db",
                       "clusterProfiler", "fgsea", "enrichplot", "EnhancedVolcano"))
```

### 4. Run scripts in order

```r
source("scripts/01_load_and_inspect.R")
source("scripts/02_deseq2_analysis.R")
source("scripts/03_figures.R")
source("scripts/04_enrichment.R")
```

### 5. Launch dashboard

```r
shiny::runApp("shiny_app")
```

---

## Tools and Packages

| Tool | Version | Purpose |
|---|---|---|
| R | 4.4.x | Primary analysis language |
| DESeq2 | 1.46 | Differential expression |
| apeglm | — | LFC shrinkage |
| clusterProfiler | 4.14 | GSEA and ORA |
| fgsea | 1.32 | GSEA permutation testing |
| msigdbr | 10.x | MSigDB Hallmark gene sets |
| EnhancedVolcano | — | Volcano plot |
| pheatmap | — | Heatmap visualization |
| Shiny + shinydashboard | — | Interactive dashboard |
| plotly | — | Interactive figures |
| org.Hs.eg.db | — | Human gene annotation |

---

## Clinical Relevance

This analysis supports the following translational insights:

1. **Pharmacodynamic biomarkers:** E2F target suppression and p21 upregulation as candidate response markers
2. **Resistance prevention:** IL-6/JAK/STAT3 activation at 12h suggests early co-targeting with JAK inhibitors (e.g., ruxolitinib) may delay acquired resistance
3. **Combination strategy:** Interferon pathway activation supports rationale for osimertinib + anti-PD-L1 combinations
4. **Dual pathway suppression:** mTORC1 co-suppression supports osimertinib + mTOR inhibitor synergy studies

---

## Future Directions

- Time-course RNA-seq to characterize kinetics of adaptive resistance
- Single-cell RNA-seq to identify drug-tolerant persister cell subpopulations
- Osimertinib + ruxolitinib combination validation (JAK/STAT3 co-inhibition)
- Cross-validation in patient biopsy datasets (FLAURA / ADAURA trial cohorts)
- siRNA/ASO co-targeting of MYC in combination with osimertinib

---

## Author

**Matt Muslu**  
GitHub: [github.com/amuslu87](https://github.com/amuslu87)

*This project is part of a computational drug discovery portfolio demonstrating RNA-seq analysis, pathway interpretation, and reproducible workflow development in the context of targeted cancer therapy.*
