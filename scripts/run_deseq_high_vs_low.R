#!/usr/bin/env Rscript

## ============================================================
## run_deseq_high_vs_low.R (version améliorée)
##
## - Charge DDS
## - Exécute DESeq2
## - Convertit ENSG → symbol
## - Filtre gènes significatifs (padj < 0.05, |log2FC| > 1)
## - Génère volcano + MA plot annotés
## - Sort une table finale propre et interprétable
## ============================================================

library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------
# 1. Chemin d’entrée/sortie
# -----------------------------

external_path <- "/mnt/e/TCGA_PDAC"

dds_file        <- file.path(external_path, "dds_tcga_paap_survival.rds")
annot_file      <- file.path(external_path, "tcga_paap_gene_annotations.csv")

res_file        <- file.path(external_path, "DESeq_results_annotated.csv")
volcano_file    <- file.path(external_path, "volcano_high_vs_low.png")
maplot_file     <- file.path(external_path, "maplot_high_vs_low.png")

# -----------------------------
# 2. Chargement DDS + annotations
# -----------------------------

cat("Chargement DDS et annotations...\n")
dds <- readRDS(dds_file)
gene_annot <- read.csv(annot_file)

# Nettoyage ID ENSG (suppression du .12)
gene_annot$ensembl_clean <- str_replace(gene_annot$gene_id, "\\..*", "")

# -----------------------------
# 3. Exécution DESeq2
# -----------------------------

cat("Running DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("risk_group", "high", "low"))
res_df <- as.data.frame(res)

# Ajout des ENSG
res_df$ensembl_id <- rownames(res_df)
res_df$ensembl_clean <- str_replace(res_df$ensembl_id, "\\..*", "")

# -----------------------------
# 4. Fusion avec annotations
# -----------------------------

res_annot <- res_df %>%
  left_join(gene_annot %>% select(ensembl_clean, gene_name, gene_type),
            by = "ensembl_clean")

# Tri par padj
res_annot <- res_annot %>% arrange(padj)

# -----------------------------
# 5. Filtrage significatif
# -----------------------------

sig <- res_annot %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("Gènes significatifs (padj < 0.05 & |log2FC| > 1) :", nrow(sig), "\n")

# -----------------------------
# 6. Volcano plot annoté
# -----------------------------

cat("Génération volcano plot...\n")

# Choix des gènes à annoter (les 10 meilleurs par p-value)
top_labels <- sig %>% head(10)

volcano <- ggplot(res_annot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1),
             alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_text(data = top_labels, aes(label = gene_name),
            vjust = -0.5, size = 3.5, color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot • High OS vs Low OS",
    x = "log2 Fold Change (High vs Low OS)",
    y = "-log10 Adjusted p-value",
    color = "Significatif"
  )

ggsave(volcano_file, volcano, width = 7, height = 6)
cat("Volcano sauvegardé :", volcano_file, "\n")

# -----------------------------
# 7. MA Plot
# -----------------------------

cat("MA plot...\n")
png(maplot_file, width = 800, height = 700)
plotMA(res, ylim = c(-5, 5), main = "MA Plot • High vs Low OS")
dev.off()

cat("MA plot sauvegardé :", maplot_file, "\n")

# -----------------------------
# 8. Sauvegarde des résultats
# -----------------------------

write.csv(res_annot, res_file, row.names = FALSE)
cat("Table DESeq annotée sauvegardée :", res_file, "\n")

# -----------------------------
# 9. Top gènes dans console
# -----------------------------

cat("\nTop 10 gènes up-régulés :\n")
print(sig %>% arrange(desc(log2FoldChange)) %>% select(gene_name, log2FoldChange, padj) %>% head(10))

cat("\nTop 10 gènes down-régulés :\n")
print(sig %>% arrange(log2FoldChange) %>% select(gene_name, log2FoldChange, padj) %>% head(10))

cat("\nAnalyse différentielle terminée.\n")
