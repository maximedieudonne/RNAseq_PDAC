#!/usr/bin/env Rscript

## ============================================================
## prepare_survival_and_deseq.R
##
## - Charge l'expression (SummarizedExperiment) et la clinique
## - Construit os_time / os_event / risk_group
## - Fusionne clinique + colData de l'expression
## - Crée un objet DESeq2 (design ~ risk_group)
## ============================================================

############### 0. Packages ###############################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_bioc <- c("DESeq2", "SummarizedExperiment", "TCGAbiolinks")
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)
}

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(DESeq2)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)

############### 1. Chemins #################################

external_path <- "/mnt/e/TCGA_PDAC"

se_file   <- file.path(external_path, "tcga_paap_star_counts_SE.rds")
clin_file <- file.path(external_path, "tcga_paap_clinical_df.rds")

cat("Lecture des fichiers:\n")
cat("  Expression :", se_file, "\n")
cat("  Clinique   :", clin_file, "\n\n")

if (!file.exists(se_file)) stop("Fichier SE manquant :", se_file)
if (!file.exists(clin_file)) stop("Fichier clinique manquant :", clin_file)

############### 2. Chargement ################################

exp.data <- readRDS(se_file)
cli.data <- readRDS(clin_file)

cat("Expression chargée : ", dim(assay(exp.data))[1], "gènes x", 
    dim(assay(exp.data))[2], "échantillons\n")
cat("Clinique chargée :", nrow(cli.data), "patients\n\n")

############### 3. Construction de os_time et os_event ########

# Colonnes possibles pour les temps
death_col <- "days_to_death"
followup_candidates <- c("days_to_last_follow_up", "days_to_last_known_disease_status")
followup_col <- intersect(followup_candidates, colnames(cli.data))
followup_col <- ifelse(length(followup_col) > 0, followup_col[1], NA)

cat("Colonne décès :", death_col, "\n")
cat("Colonne follow-up :", followup_col, "\n\n")

# Conversion numériques
for (v in c(death_col, followup_col)) {
  if (!is.na(v) && v %in% colnames(cli.data)) {
    cli.data[[v]] <- suppressWarnings(as.numeric(as.character(cli.data[[v]])))
  }
}

# os_time
cli.data$os_time <- NA_real_
if (death_col %in% colnames(cli.data)) {
  cli.data$os_time <- ifelse(!is.na(cli.data[[death_col]]), cli.data[[death_col]], NA)
}
if (!is.na(followup_col)) {
  cli.data$os_time <- ifelse(is.na(cli.data$os_time) & !is.na(cli.data[[followup_col]]),
                             cli.data[[followup_col]],
                             cli.data$os_time)
}

# os_event
cli.data$os_event <- ifelse(tolower(cli.data$vital_status) == "dead", 1L, 0L)

cat("Variables de survie créées : os_time + os_event\n\n")

# Sous-table clinique standardisée
clin_surv <- cli.data %>%
  select(
    bcr_patient_barcode, vital_status, all_of(death_col),
    all_of(followup_col), os_time, os_event, everything()
  )

############### 4. Fusion expression + clinique ###############

# Extraction colData
coldata <- as.data.frame(colData(exp.data))

# Recherche de la colonne contenant le barcode complet
barcode_col <- NULL

if ("barcode" %in% colnames(coldata)) {
  barcode_col <- "barcode"
} else if ("cases" %in% colnames(coldata)) {
  barcode_col <- "cases"
} else {
  cat("Colonnes disponibles dans colData(exp.data) :\n")
  print(colnames(coldata))
  stop("Impossible de trouver la colonne 'barcode' ou 'cases' dans colData(exp.data).")
}

# Nettoyage des ID
coldata$sample_barcode  <- as.character(coldata[[barcode_col]])
coldata$patient_barcode <- substr(coldata$sample_barcode, 1, 12)

# Détection de la colonne patient dans la clinique
patient_id_candidates <- c("bcr_patient_barcode", "submitter_id", "case_id", "patient_id")
patient_cols <- intersect(patient_id_candidates, colnames(clin_surv))

if (length(patient_cols) == 0) {
  print(colnames(clin_surv))
  stop("Aucune colonne patient trouvée dans la clinique.")
}

# Priorité : bcr_patient_barcode
if ("bcr_patient_barcode" %in% patient_cols) {
  patient_col <- "bcr_patient_barcode"
} else {
  # sinon, prendre le premier trouvé
  patient_col <- patient_cols[1]
}

cat("Colonne patient utilisée pour fusion :", patient_col, "\n")

clin_surv <- clin_surv %>%
  rename_with(~ "patient_barcode", all_of(patient_col))


# Fusion
merged_coldata <- coldata %>% left_join(clin_surv, by = "patient_barcode")

cat("Fusion expression + clinique OK.\n")
cat("Échantillons avec données OS :", sum(!is.na(merged_coldata$os_time)), "\n\n")

colData(exp.data) <- S4Vectors::DataFrame(merged_coldata)

############### 5. Génération DESeq2 ###########################

assay_name <- "unstranded"
if (!assay_name %in% assayNames(exp.data))
  stop("Assay ", assay_name, " non trouvé.")

counts <- assay(exp.data, assay_name)

# Filtrage des gènes
keep_genes <- rowSums(counts >= 10) >= 10
exp_filt <- exp.data[keep_genes, ]
cat("Gènes conservés :", sum(keep_genes), "\n")

coldata_filt <- as.data.frame(colData(exp_filt))

# Garder seulement les échantillons avec OS
sel <- !is.na(coldata_filt$os_time)
exp_filt <- exp_filt[, sel]
coldata_filt <- as.data.frame(colData(exp_filt))

cat("Échantillons gardés pour DESeq :", ncol(exp_filt), "\n")

# Construction risk_group
median_os <- median(coldata_filt$os_time, na.rm = TRUE)
coldata_filt$risk_group <- ifelse(coldata_filt$os_time >= median_os, "high", "low")
coldata_filt$risk_group <- factor(coldata_filt$risk_group, levels = c("low", "high"))

colData(exp_filt) <- S4Vectors::DataFrame(coldata_filt)

dds <- DESeqDataSet(exp_filt, design = ~ risk_group)

cat("DESeqDataSet créé : ", dim(dds)[1], "gènes x", dim(dds)[2], "échantillons\n\n")



############### 6. Sauvegarde #################################

dds_file <- file.path(external_path, "dds_tcga_paap_survival.rds")
clin_surv_file <- file.path(external_path, "tcga_paap_clinical_with_survival.csv")

# Sauvegarde de l'objet DESeq2
saveRDS(dds, dds_file)

cat("Conversion des colonnes 'list' dans coldata_filt pour export CSV...\n")

# Détecter les colonnes de type list
is_list_col <- sapply(coldata_filt, is.list)

if (any(is_list_col)) {
  cat("Colonnes list détectées :", paste(names(coldata_filt)[is_list_col], collapse = ", "), "\n")
  
  # Conversion list → chaîne de texte
  coldata_filt[is_list_col] <- lapply(
    coldata_filt[is_list_col],
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

# Export en CSV
write.csv(coldata_filt, clin_surv_file, row.names = FALSE)

cat("Fichiers sauvegardés :\n")
cat("  - dds :", dds_file, "\n")
cat("  - clinique enrichie :", clin_surv_file, "\n\n")


############### 7. Aperçu #####################################

print(head(as.data.frame(colData(dds))[, c("patient_barcode", "os_time", "os_event", "risk_group")]))
