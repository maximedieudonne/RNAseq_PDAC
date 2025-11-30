#!/usr/bin/env Rscript

## ============================================================
## download_clinical_data.R
## Récupération et sauvegarde des données cliniques TCGA-PAAD
## vers un disque externe sous WSL (ex : /mnt/e/TCGA_PDAC)
## Utilise GDCquery_clinic() (pas de XML, pas de GDCprepare_clinic)
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

library(TCGAbiolinks)

# -----------------------------
# 1. Chemin du disque externe
# -----------------------------

# À adapter si besoin (par ex. "/mnt/d/TCGA_PDAC")
external_path <- "/mnt/e/TCGA_PDAC"

dir.create(external_path, recursive = TRUE, showWarnings = FALSE)
cat("Dossier de sortie :", external_path, "\n")

# -----------------------------
# 2. Récupération clinique
# -----------------------------
# On interroge directement l'API GDC pour obtenir un data.frame
# clinique (niveau patient) sans manipuler de fichiers sur disque.

cat("Récupération des données cliniques TCGA-PAAD via GDCquery_clinic()...\n")

cli.data <- GDCquery_clinic(
  project = "TCGA-PAAD",
  type    = "clinical"
)

cat("Récupération terminée. Nombre de patients :", nrow(cli.data), "\n")

# -----------------------------
# 3. Sauvegarde sur disque externe
# -----------------------------

rds_file <- file.path(external_path, "tcga_paap_clinical_df.rds")
csv_file <- file.path(external_path, "tcga_paap_clinical_df.csv")

# 3.1 Sauvegarde RDS (conserve les colonnes liste, idéal pour R)
saveRDS(cli.data, rds_file)

# 3.2 Création d'une version "aplatie" pour export CSV
cli.data.flat <- cli.data

# repérer les colonnes de type list
is_list_col <- sapply(cli.data.flat, is.list)

if (any(is_list_col)) {
  cat("Colonnes de type 'list' détectées, conversion en texte pour CSV...\n")
  
  cli.data.flat[is_list_col] <- lapply(
    cli.data.flat[is_list_col],
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
} else {
  cat("Aucune colonne de type 'list' trouvée, export CSV direct.\n")
}

# 3.3 Sauvegarde CSV à partir de la version aplatie
write.csv(cli.data.flat, csv_file, row.names = FALSE)

cat("Données cliniques sauvegardées :\n")
cat("  - RDS : ", rds_file, "\n")
cat("  - CSV : ", csv_file, "\n")

# -----------------------------
# 4. Aperçu rapide
# -----------------------------

cat("\nAperçu des premières colonnes (version RDS originale) :\n")
print(head(cli.data[, 1:min(10, ncol(cli.data))]))
