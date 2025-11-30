#!/usr/bin/env Rscript

## ============================================================
## download_expression_data.R
## Téléchargement et sauvegarde des données d'expression TCGA-PAAD
## (Transcriptome Profiling, STAR - Counts)
## vers un disque externe sous WSL (ex : /mnt/e/TCGA_PDAC)
## ============================================================

# -----------------------------
# 0. Chargement / installation
# -----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment")
}

library(TCGAbiolinks)
library(SummarizedExperiment)

# -----------------------------
# 1. Chemin du disque externe
# -----------------------------

# À adapter si besoin (par ex. "/mnt/d/TCGA_PDAC")
external_path <- "/mnt/e/TCGA_PDAC"

dir.create(external_path, recursive = TRUE, showWarnings = FALSE)
cat("Dossier de sortie :", external_path, "\n")

# Dossier où GDCdownload écrira les fichiers bruts
gdc_raw_dir <- file.path(external_path, "GDC_raw_expression")
dir.create(gdc_raw_dir, recursive = TRUE, showWarnings = FALSE)
cat("Dossier de téléchargement brut :", gdc_raw_dir, "\n")

# -----------------------------
# 2. Requête expression
# -----------------------------

cat("Préparation requête expression TCGA-PAAD (STAR - Counts)...\n")

query.exp <- GDCquery(
  project       = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# -----------------------------
# 3. Téléchargement
# -----------------------------

cat("Téléchargement des fichiers d'expression dans :", gdc_raw_dir, "\n")

GDCdownload(
  query     = query.exp,
  directory = gdc_raw_dir
)

cat("Téléchargement terminé.\n")

# -----------------------------
# 4. Préparation SummarizedExperiment
# -----------------------------

cat("Préparation de l'objet SummarizedExperiment (GDCprepare)...\n")

exp.data <- GDCprepare(
  query     = query.exp,
  directory = gdc_raw_dir
)

cat("Préparation terminée.\n")
cat("Classes :", class(exp.data), "\n")
cat("Assays disponibles :", paste(assayNames(exp.data), collapse = ", "), "\n")

# -----------------------------
# 5. Sauvegarde de l'objet complet
# -----------------------------

se_file <- file.path(external_path, "tcga_paap_star_counts_SE.rds")
saveRDS(exp.data, se_file)
cat("Objet SummarizedExperiment sauvegardé sous :", se_file, "\n")

# -----------------------------
# 6. Extraction et sauvegarde de la matrice de comptes
# -----------------------------

cat("Extraction de la matrice de comptes 'unstranded'...\n")

assay_name <- "unstranded"
if (!assay_name %in% assayNames(exp.data)) {
  stop("L'assay '", assay_name, "' n'est pas disponible dans exp.data. Assays dispo : ",
       paste(assayNames(exp.data), collapse = ", "))
}

counts_mat <- assay(exp.data, assay_name)

counts_rds_file <- file.path(external_path, "tcga_paap_counts_unstranded.rds")
counts_csv_file <- file.path(external_path, "tcga_paap_counts_unstranded.csv")

saveRDS(counts_mat, counts_rds_file)
write.csv(counts_mat, counts_csv_file)

cat("Matrice de comptes sauvegardée :\n")
cat("  - RDS :", counts_rds_file, "\n")
cat("  - CSV :", counts_csv_file, "\n")

# -----------------------------
# 7. Sauvegarde des annotations
# -----------------------------

cat("Sauvegarde des annotations gènes (rowData) et échantillons (colData)...\n")

gene_annot   <- as.data.frame(rowData(exp.data))
sample_annot <- as.data.frame(colData(exp.data))

gene_annot_file   <- file.path(external_path, "tcga_paap_gene_annotations.csv")
sample_annot_file <- file.path(external_path, "tcga_paap_sample_annotations.csv")

## 7.1 Export des annotations gènes
## (rowData est généralement "propre", mais on sécurise au cas où)

is_list_gene <- sapply(gene_annot, is.list)

if (any(is_list_gene)) {
  cat("Colonnes de type 'list' détectées dans gene_annot, conversion en texte pour CSV...\n")
  gene_annot[is_list_gene] <- lapply(
    gene_annot[is_list_gene],
    function(col) {
      sapply(col, function(x) {
        if (length(x) == 0 || all(is.na(x))) {
          NA_character_
        } else {
          paste(as.character(x), collapse = ";")
        }
      })
    }
  )
}

write.csv(gene_annot, gene_annot_file, row.names = TRUE)

## 7.2 Export des annotations échantillons (colData)
## Ici, on sait qu'il y a des colonnes list → aplatissement obligatoire

is_list_sample <- sapply(sample_annot, is.list)

if (any(is_list_sample)) {
  cat("Colonnes de type 'list' détectées dans sample_annot, conversion en texte pour CSV...\n")
  sample_annot[is_list_sample] <- lapply(
    sample_annot[is_list_sample],
    function(col) {
      sapply(col, function(x) {
        if (length(x) == 0 || all(is.na(x))) {
          NA_character_
        } else {
          paste(as.character(x), collapse = ";")
        }
      })
    }
  )
}

write.csv(sample_annot, sample_annot_file, row.names = TRUE)

cat("Annotations sauvegardées :\n")
cat("  - Gènes   :", gene_annot_file, "\n")
cat("  - Samples :", sample_annot_file, "\n")

# -----------------------------
# 8. Aperçu rapide
# -----------------------------

cat("\nDimensions de la matrice de comptes (gènes x échantillons) :\n")
print(dim(counts_mat))

cat("\nAperçu des 5 premiers gènes (5 colonnes) :\n")
print(counts_mat[1:5, 1:5])

