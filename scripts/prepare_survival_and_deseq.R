#!/usr/bin/env Rscript

## ============================================================
## prepare_survival_and_deseq.R
##
## - Recharge l'expression (SummarizedExperiment) et la clinique
## - Construit les variables de survie (os_time, os_event)
## - Crée un objet DESeq2 prêt pour l'analyse différentielle
##   (design ~ risk_group basé sur la médiane d'OS)
##
## Fichiers attendus dans /mnt/e/TCGA_PDAC :
##   - tcga_paap_star_counts_SE.rds
##   - tcga_paap_clinical_df.rds
## ============================================================

# -----------------------------
# 0. Chargement / installation
# -----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment")
}

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(DESeq2)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)

# -----------------------------
# 1. Chemin des fichiers
# -----------------------------

external_path <- "/mnt/e/TCGA_PDAC"   # adapte si besoin

se_file   <- file.path(external_path, "tcga_paap_star_counts_SE.rds")
clin_file <- file.path(external_path, "tcga_paap_clinical_df.rds")

cat("Lecture des données depuis :", external_path, "\n")
cat("  - Expression SE :", se_file, "\n")
cat("  - Clinique      :", clin_file, "\n")

# -----------------------------
# 2. Chargement des données
# -----------------------------

if (!file.exists(se_file)) {
  stop("Fichier SummarizedExperiment introuvable : ", se_file)
}
if (!file.exists(clin_file)) {
  stop("Fichier clinique introuvable : ", clin_file)
}

exp.data <- readRDS(se_file)
cli.data <- readRDS(clin_file)

cat("Objet expression (SummarizedExperiment) chargé.\n")
cat("Dimensions assays :", dim(assay(exp.data)), "\n")
cat("Nombre de patients (clinique) :", nrow(cli.data), "\n")

# -----------------------------
# 3. Construction des variables de survie
# -----------------------------
# On suppose les colonnes suivantes dans cli.data :
# - "bcr_patient_barcode"
# - "vital_status" ("Dead"/"Alive")
# - "days_to_death"
# - "days_to_last_followup"
# (si besoin, adapter ici en fonction des noms réels)

stopifnot("bcr_patient_barcode" %in% colnames(cli.data))

# On convertit les colonnes de temps en numérique (elles sont souvent char)
time_vars <- c("days_to_death", "days_to_last_followup")
for (v in time_vars) {
  if (v %in% colnames(cli.data)) {
    cli.data[[v]] <- suppressWarnings(as.numeric(as.character(cli.data[[v]])))
  }
}

# Temps d'OS : priorité à days_to_death, sinon days_to_last_followup
cli.data <- cli.data %>%
  mutate(
    os_time = dplyr::case_when(
      !is.na(days_to_death) ~ days_to_death,
      !is.na(days_to_last_followup) ~ days_to_last_followup,
      TRUE ~ NA_real_
    ),
    os_event = ifelse(vital_status == "Dead", 1L, 0L)
  )

cat("Variables de survie construites (os_time, os_event).\n")

# On garde un sous-ensemble propre (optionnel)
clin_surv <- cli.data %>%
  select(bcr_patient_barcode, vital_status, days_to_death,
         days_to_last_followup, os_time, os_event, everything())

# -----------------------------
# 4. Harmonisation des IDs avec l'expression
# -----------------------------
# Dans exp.data (SummarizedExperiment), colData contient des barcodes complets
# du style "TCGA-AB-1234-01A-01R-..." ; le patient = 12 premiers caractères.

coldata <- as.data.frame(colData(exp.data))

# nom de la colonne qui contient le barcode complet
# TCGAbiolinks ajoute souvent "barcode" ou "cases"
barcode_col <- NULL
if ("barcode" %in% colnames(coldata)) {
  barcode_col <- "barcode"
} else if ("cases" %in% colnames(coldata)) {
  barcode_col <- "cases"
} else {
  # Si nécessaire, afficher les noms pour debug
  cat("Colonnes disponibles dans colData(exp.data) :\n")
  print(colnames(coldata))
  stop("Impossible de trouver la colonne 'barcode' ou 'cases' dans colData(exp.data).")
}

coldata$sample_barcode  <- as.character(coldata[[barcode_col]])
coldata$patient_barcode <- substr(coldata$sample_barcode, 1, 12)

# On renomme côté clinique pour merger
clin_surv <- clin_surv %>%
  rename(patient_barcode = bcr_patient_barcode)

# Merge colData + clinique
merged_coldata <- coldata %>%
  left_join(clin_surv, by = "patient_barcode")

cat("Fusion expression / clinique effectuée.\n")
cat("Nombre d'échantillons avec os_time non NA : ",
    sum(!is.na(merged_coldata$os_time)), "\n")

# On remet ceci dans colData(exp.data)
colData(exp.data) <- S4Vectors::DataFrame(merged_coldata)

# -----------------------------
# 5. Création de l'objet DESeq2
# -----------------------------
# On choisit l'assay "unstranded" comme matrice de comptes.
# On filtre les gènes avec peu de counts, puis on définit une variable
# de groupe de risque (high/low OS) basée sur la médiane d'os_time.

assay_name <- "unstranded"
if (!assay_name %in% assayNames(exp.data)) {
  stop("L'assay '", assay_name, "' n'existe pas dans exp.data. Assays dispo : ",
       paste(assayNames(exp.data), collapse = ", "))
}

counts <- assay(exp.data, assay_name)

# Filtrage simple : garder les gènes exprimés (>= 10 counts dans au moins 10 échantillons)
keep_genes <- rowSums(counts >= 10) >= 10
counts_filt <- counts[keep_genes, ]

cat("Filtrage des gènes : ", sum(keep_genes), " gènes conservés.\n")

# On met à jour l'objet SummarizedExperiment avec ces gènes filtrés
exp_filt <- exp.data[keep_genes, ]
coldata_filt <- as.data.frame(colData(exp_filt))

# On ne garde que les échantillons avec os_time non NA pour DESeq
non_na_os <- !is.na(coldata_filt$os_time)
exp_filt <- exp_filt[, non_na_os]
coldata_filt <- as.data.frame(colData(exp_filt))

cat("Nombre d'échantillons avec OS utilisé pour DESeq : ",
    ncol(exp_filt), "\n")

# Construction de la variable de risk_group (high/low OS)
os_time_vec <- coldata_filt$os_time
median_os   <- median(os_time_vec, na.rm = TRUE)

risk_group <- ifelse(os_time_vec >= median_os, "high", "low")
coldata_filt$risk_group <- factor(risk_group, levels = c("low", "high"))

colData(exp_filt) <- S4Vectors::DataFrame(coldata_filt)

# Création de l'objet DESeqDataSet
dds <- DESeqDataSet(
  se      = exp_filt,
  design  = ~ risk_group
)

cat("Objet DESeqDataSet créé avec design ~ risk_group.\n")
cat("Dimensions dds (gènes x échantillons) :", dim(dds), "\n")

# -----------------------------
# 6. Sauvegarde des résultats
# -----------------------------

dds_file <- file.path(external_path, "dds_tcga_paap_survival.rds")
clin_surv_file <- file.path(external_path, "tcga_paap_clinical_with_survival.csv")

saveRDS(dds, dds_file)
write.csv(coldata_filt, clin_surv_file, row.names = FALSE)

cat("Objet DESeq2 sauvegardé sous : ", dds_file, "\n")
cat("Clinique + survie + risk_group sauvegardé sous : ", clin_surv_file, "\n")

# -----------------------------
# 7. Aperçu
# -----------------------------

cat("\nAperçu du colData(dds) :\n")
print(head(as.data.frame(colData(dds))[, c("patient_barcode", "os_time", "os_event", "risk_group")]))
