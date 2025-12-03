#!/usr/bin/env Rscript

## ============================================================
## scrna_GSE212966_seurat.R
##
## - Télécharge les données GSE212966 (si besoin)
## - Charge la matrice single-cell dans Seurat
## - QC / filtrage basique
## - Normalisation, PCA, UMAP, clustering
## - Sauvegarde de l’objet Seurat + quelques plots
## ============================================================

## ---------- 0. Packages ----------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(GEOquery)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

## ---------- 1. Chemins ----------

external_path <- "/mnt/e/TCGA_PDAC"  # adapte si besoin
gse_id        <- "GSE212966"

gse_dir       <- file.path(external_path, gse_id)
dir.create(gse_dir, recursive = TRUE, showWarnings = FALSE)

seurat_rds    <- file.path(gse_dir, paste0(gse_id, "_seurat_raw.rds"))
umap_plot_png <- file.path(gse_dir, paste0(gse_id, "_UMAP_clusters.png"))
qc_violin_png <- file.path(gse_dir, paste0(gse_id, "_QC_violin.png"))

cat("Dossier GSE212966 :", gse_dir, "\n")

## ---------- 2. Téléchargement des fichiers GEO (si besoin) ----------

if (length(list.files(gse_dir, recursive = TRUE)) == 0) {
  cat("Aucun fichier trouvé dans", gse_dir, "→ téléchargement depuis GEO...\n")
  getGEOSuppFiles(gse_id, baseDir = gse_dir)
  cat("Téléchargement terminé.\n")
} else {
  cat("Des fichiers existent déjà dans", gse_dir, "→ pas de téléchargement.\n")
}

## ---------- 3. Détection du format (h5 vs 10x mtx) ----------

all_files <- list.files(gse_dir, recursive = TRUE, full.names = TRUE)
cat("Nombre de fichiers trouvés (supplementary) :", length(all_files), "\n")

# On cherche d'abord un .h5 (format 10x récent)
h5_files <- grep("\\.h5$", all_files, value = TRUE)

# Sinon, on cherche un dossier de type 10x (matrix.mtx + barcodes + genes/features)
mtx_files <- grep("matrix.mtx", all_files, ignore.case = TRUE, value = TRUE)

if (length(h5_files) > 0) {
  cat("Fichier(s) .h5 détecté(s) :\n")
  print(h5_files)
  h5_file <- h5_files[1]
  cat("→ On utilise :", h5_file, "\n")
  
  counts <- Read10X_h5(h5_file)
  
} else if (length(mtx_files) > 0) {
  cat("Fichier(s) matrix.mtx détecté(s) :\n")
  print(mtx_files)
  mtx_file <- mtx_files[1]
  
  tenx_dir <- dirname(mtx_file)
  cat("→ On suppose un dossier 10x, on utilise :", tenx_dir, "\n")
  
  counts <- Read10X(data.dir = tenx_dir)
  
} else {
  stop("Impossible de trouver un .h5 ou un matrix.mtx dans les suppléments de GSE212966.")
}

cat("Dimensions de la matrice de comptes (features x cellules) :", dim(counts), "\n")

## ---------- 4. Création de l’objet Seurat ----------

# On part du principe que c'est de l'RNA
sc <- CreateSeuratObject(
  counts = counts,
  project = "GSE212966_PDAC",
  min.cells = 3,   # gene gardé si présent dans ≥3 cellules
  min.features = 200  # cellule gardée si ≥200 gènes détectés
)

cat("Objet Seurat brut créé. Dimensions (cellules x gènes) :\n")
print(dim(sc))

## ---------- 5. QC : nFeature, nCount, % mito ----------

# Les gènes mitochondriaux humains commencent souvent par "MT-"
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# VlnPlot QC
qc_plot <- VlnPlot(
  sc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(qc_violin_png, qc_plot, width = 9, height = 4)
cat("Plot QC (violons) sauvegardé sous :", qc_violin_png, "\n")

# Filtrage basique : paramètres standards, à affiner
sc <- subset(
  sc,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 20
)

cat("Après filtrage QC, dimensions (cellules x gènes) :\n")
print(dim(sc))

## ---------- 6. Normalisation + réduction de dimension ----------

# 6.1 Normalisation log
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 1e4)

# 6.2 Sélection des gènes variables
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# 6.3 Mise à l’échelle
sc <- ScaleData(sc, features = rownames(sc))

# 6.4 PCA
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

cat("PCA réalisée.\n")

## ---------- 7. Clustering + UMAP ----------

# On choisit 1:30 PC, à ajuster après inspection du ElbowPlot
sc <- FindNeighbors(sc, dims = 1:30)
sc <- FindClusters(sc, resolution = 0.5)

# UMAP
sc <- RunUMAP(sc, dims = 1:30)
cat("UMAP + clustering terminés.\n")

## ---------- 8. Plots UMAP ----------

umap_plot <- DimPlot(sc, reduction = "umap", label = TRUE) +
  ggtitle("GSE212966 - Clusters Seurat")

ggsave(umap_plot_png, umap_plot, width = 7, height = 6)
cat("UMAP clusters sauvegardé sous :", umap_plot_png, "\n")

## ---------- 9. Sauvegarde de l’objet Seurat ----------

saveRDS(sc, seurat_rds)
cat("Objet Seurat sauvegardé sous :", seurat_rds, "\n")

## ---------- 10. Aperçu des clusters ----------

cat("\nRésumé des clusters :\n")
print(table(Idents(sc)))

cat("\nAnalyse single-cell (étape 1 : QC + clustering) terminée.\n")
