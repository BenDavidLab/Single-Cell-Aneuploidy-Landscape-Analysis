# Single-Cell Aneuploidy Landscape Analysis

This repository contains the R code supporting the analyses in:

> **[Pan-cancer analysis of single-cell RNA sequencing data from 304 human tumors sheds light on the ‘aneuploidy paradox’]** — *[Guy Wolf-Dankovich, Tomer Mashiah, Ron Saad, Einav Somech, Haia Khoury, Itay Tirosh, Uri Ben-David]* 

---

## Overview

This study characterizes the transcriptional consequences of aneuploidy across hundreds of single-cell tumor samples, using a cohort of publicly available scRNA-seq datasets integrated with chromosome-arm-level copy number profiles inferred by the ASCETS algorithm.

The analyses are organized into three scripts:

| File | Contents |
|------|----------|
| `01_main_analysis.R` | Core analysis pipeline, including AS calculation, cell-state grouping, differential gene expression, pathway enrichment, ssGSEA meta-analysis, karyotypic heterogeneity analysis, and CN subclone identification |
| `02_recurrent_aneuploidies.R` | Analysis of the transcriptional consequences associated with recurrent chromosome-arm gains and losses |
| `03_figures.R` | Code to reproduce figures from the main analysis and recurrent aneuploidy analyses |
| `04_immune_analysis.R` | Analysis of the relationship between aneuploidy and tumor cell composition, with code to reproduce the corresponding figures|
| `05_TAM_analysis.R` | Analysis of tumor-associated macrophages, with code to reproduce the corresponding figures |
---

## Data Availability

Gene expression data are publicly available (see paper for accession numbers and sources).  
Copy number profiles were derived using the **ASCETS** algorithm.  
TCGA aneuploidy scores and arm-level CNA calls are from [Taylor et al. 2018, *Cancer Cell*](https://doi.org/10.1016/j.ccell.2018.03.007).

The **sample file map** (`sampleFilesMap`) referenced throughout is a CSV with one row per sample, containing paths to the relevant expression and CNA files. Its columns are described in the code comments.

---

## Dependencies

```r
# CRAN
install.packages(c("tidyverse", "dplyr", "tidyr", "reshape2", "patchwork",
                   "cowplot", "ggrepel", "ggpubr", "effectsize", "metafor",
                   "weights", "sandwich", "lmtest", "pheatmap", "FNN",
                   "dendextend", "cluster", "stringr"))

# Bioconductor
BiocManager::install(c("clusterProfiler", "msigdbr", "GSVA",
                       "Seurat", "enrichplot"))

# Other
install.packages(c("doParallel", "foreach", "data.table",
                   "Matrix", "igraph", "icesTAF", "corto"))
```

---

## sampleFilesMap Structure

The `sampleFilesMap` CSV must contain the following columns:

| Column | Description |
|--------|-------------|
| `Sample` | Sample identifier |
| `Study` | Study/paper identifier |
| `Cancer_type` | Cancer type label |
| `Expression_file_type` | 0 = UMI (10X), 1 = TPM (Smart-seq2), 2 = excluded |
| `Cells_path` | Path to cells CSV (columns: `cell_name`, `sample`, `cell_type`) |
| `Cna_mat` | Path to chromosome-arm CNA matrix (TSV; rows = cells, columns = chr arms) |
| `Norm_Expression_path` | Path to normalized (centered) expression matrix |
| `Norm_Without_Center_Expression_path` | Path to normalized (non-centered) expression matrix |
| `Gene_cna_mat` | Path to gene-level CNA matrix (space-separated) |
