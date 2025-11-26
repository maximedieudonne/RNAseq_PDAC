# Données utilisées

## 1. TCGA-PAAD (bulk RNA-seq + clinique)

Source : The Cancer Genome Atlas (TCGA), projet **TCGA-PAAD (Pancreatic Adenocarcinoma)** via le portail GDC.

### Étapes (version R via TCGAbiolinks)

Dans R :

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# Expression
query.exp <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(query.exp)
exp.data <- GDCprepare(query.exp)

# Données cliniques
query.cli <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Clinical"
)
GDCdownload(query.cli)
cli.data <- GDCprepare_clinic(query.cli, clinical.info = "patient")
