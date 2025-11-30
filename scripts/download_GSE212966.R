if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

external_path <- "/mnt/e/TCGA_PDAC"
geo_dir <- file.path(external_path, "GSE212966")

dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO("GSE212966", GSEMatrix = FALSE)

getGEOSuppFiles("GSE212966", baseDir = geo_dir)