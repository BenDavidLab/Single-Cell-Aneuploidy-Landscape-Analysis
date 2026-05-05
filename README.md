# Single-Cell Aneuploidy Landscape Analysis

This repository contains the R code supporting the analyses in:

> **[Pan-cancer analysis of single-cell RNA sequencing data from 304 human tumors sheds light on the ‘aneuploidy paradox’]** — *[Guy Wolf-Dankovich, Tomer Mashiah, Ron Saad, Einav Somech, Haia Khoury, Itay Tirosh, Uri Ben-David]* 

---

## Overview

This study characterizes the transcriptional consequences of aneuploidy across hundreds of single-cell tumor samples, using a cohort of publicly available scRNA-seq datasets integrated with chromosome-arm-level copy number profiles inferred by the ASCETS algorithm.

The analyses are organized into three scripts:

| File | Contents |
|------|----------|
| `01_main_analysis.R` | Core analysis pipeline, including AS calculation, cell-state grouping, differential gene expression, pathway enrichment, ssGSEA meta-analysis, karyotypic heterogeneity analysis, CN subclone identification, cell cycle scoring, QC aggregation, cell-phase analysis |
| `02_recurrent_aneuploidies.R` | Analysis of the transcriptional consequences associated with recurrent chromosome-arm gains and losses: ssGSEA, meta-analysis, IPTW-controlled DGE, SC vs TCGA correlation |
| `03_figures.R` | Code to reproduce figures from the main analysis and recurrent aneuploidy analyses |
| `04_immune_analysis.R` | Analysis of the relationship between aneuploidy and tumor cell composition, with code to reproduce the corresponding figures|
| `05_TAM_analysis.R` | Analysis of tumor-associated macrophages, with code to reproduce the corresponding figures |
---

## Data Availability

Gene expression data are publicly available (see paper for accession numbers and sources).  
Copy number profiles were derived using the **ASCETS** algorithm.  
TCGA aneuploidy scores and arm-level CNA calls are from [Taylor et al. 2018, *Cancer Cell*](https://doi.org/10.1016/j.ccell.2018.03.007).

---
