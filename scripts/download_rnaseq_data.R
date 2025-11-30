library(TCGAbiolinks)

external_path <- "/mnt/e/TCGA_PDAC"   
dir.create(external_path, recursive = TRUE, showWarnings = FALSE)

# 1) Préparer la requête expression
query.exp <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 2) Télécharger les fichiers bruts sur le disque externe, dans un sous-dossier "GDC_raw"
gdc_raw_dir <- file.path(external_path, "GDC_raw")
dir.create(gdc_raw_dir, recursive = TRUE, showWarnings = FALSE)

GDCdownload(query.exp, directory = gdc_raw_dir)

# 3) Préparer l’objet SummarizedExperiment à partir des fichiers téléchargés
exp.data <- GDCprepare(query.exp, directory = gdc_raw_dir)


