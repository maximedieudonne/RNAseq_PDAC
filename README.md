# PDAC Multiomics Project

Intégration de données **bulk RNA-seq** (TCGA-PAAD) et **single-cell RNA-seq** (GSE212966) pour :
- caractériser le microenvironnement tumoral du cancer du pancréas (PDAC),
- construire des **signatures géniques dérivées du single-cell**,
- projeter ces signatures sur les tumeurs TCGA,
- et développer des **modèles de machine learning / survie** pour prédire le risque des patients.

Ce projet est conçu comme un **portfolio** montrant des compétences en :
- Bioinformatique (R & Python),
- Analyse RNA-Seq (bulk & single-cell),
- Intégration multi-omique,
- Machine Learning et modélisation prédictive,
- Reproductibilité (scripts et notebooks structurés).

---

## Structure du projet

```text
pdac-multiomics-project/
│
├── README.md
├── environment.yml
├── requirements.txt
├── data/
│   ├── README_data.md
│   └── (fichiers bruts ou préprocessés, ignorés par git)
│
├── notebooks/
│   ├── 01_bulk_tcga_preprocessing.ipynb
│   ├── 02_scrnaseq_gse212966_seurat.ipynb
│   ├── 03_integration_signatures_bulk.ipynb
│   └── 04_survival_ml_models.ipynb
│
└── src/
    ├── r/
    │   ├── bulk_tcga_pipeline.R
    │   └── gsva_signatures.R
    └── python/
        ├── ml_survival_models.py
        └── utils_data_loading.py