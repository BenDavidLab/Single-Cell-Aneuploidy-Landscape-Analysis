################################################################################
# 01_main_analysis.R
#
# Core analytical pipeline for the single-cell aneuploidy study.
#
# This script covers:
#   1.  Aneuploidy Score (AS) calculation
#   2.  Cell grouping by AS (High / Mid / Low)
#   3.  Per-sample differential gene expression (Seurat FindMarkers)
#   4.  Continuous regression of expression on AS (lmFit)
#   5.  Pathway enrichment (clusterProfiler::enricher) for both approaches
#   6.  Per-sample ssGSEA (GSVA) + Welch t-test (high-AS vs low-AS cells)
#   7.  Random-effects meta-analysis across samples (metafor::rma)
#   8.  Pseudo-bulk mean expression & linear regression
#   9.  Karyotypic heterogeneity score (AHS) from pairwise CN correlations
#  10.  CN subclone identification (Louvain/Seurat) and merging
#  11.  QC: cell cycle scoring & expression QC metrics
#  12.  Per-Cell QC Metrics & Aneuploidy Score Collection
#  13. AHS-Stratified Meta-Analysis 
#  14. Subclone-Level Karyotypic Heterogeneity & DGE
#  15. Subclone Enrichment Summary Table
#
# NOTE: This code is provided for transparency. The gene expression data are
# publicly available; paths should be updated to your local copies.
################################################################################

## sampleFilesMap Structure

# The `sampleFilesMap` CSV must contain the following columns:
#   
#   | Column | Description |
#   |--------|-------------|
#   | `Sample` | Sample identifier |
#   | `Study` | Study/paper identifier |
#   | `Cancer_type` | Cancer type label |
#   | `Expression_file_type` | 0 = UMI (10X), 1 = TPM (Smart-seq2), 2 = excluded |
#     | `Cells_path` | Path to cells CSV (columns: `cell_name`, `sample`, `cell_type`) |
#     | `Cna_mat` | Path to chromosome-arm CNA matrix (TSV; rows = cells, columns = chr arms) |
#     | `Norm_Expression_path` | Path to normalized (centered) expression matrix |
#     | `Norm_Without_Center_Expression_path` | Path to normalized (non-centered) expression matrix |
#     | `Gene_cna_mat` | Path to gene-level CNA matrix (space-separated) |
    

# ==============================================================================
# 0.  Libraries & shared aesthetics
# ==============================================================================

suppressMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(stats)
  library(rstatix)
  library(msigdbr)
  library(icesTAF)
  library(clusterProfiler)
  library(doParallel)
  library(foreach)
  library(ggpubr)
  library(Seurat)
  library(Matrix)
  library(GSVA)
  library(metafor)
  library(effectsize)
  library(data.table)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
})

# ---- pathway databases to query ----
pathways_titles_list      <- c("BIOCARTA", "KEGG", "REACTOME", "GOBP", "GOCC", "HALLMARK")
pathways_anontations_list <- c("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "GO:BP", "GO:CC", "")
pathways_categories_list  <- c("C2", "C2", "C2", "C5", "C5", "H")

# ==============================================================================
# 1.  Aneuploidy Score (AS) Calculation
# ==============================================================================
# Input:  cna_df  — data.frame, rows = cells, columns = chromosome arms
#                   (values: "AMP", "DEL", "NEUTRAL", "NC")
#         sample_cells — character vector of cell barcodes to include
#
# Output: cna_df with an added 'AS' column

calculate_AS <- function(cna_df, sample_cells) {
  cna_df <- cna_df[rownames(cna_df) %in% sample_cells, ]
  cna_df$AS <- apply(
    cna_df[, !colnames(cna_df) %in% c("Xp", "Xq", "Yp", "Yq")],
    MARGIN = 1,
    FUN = function(x) length(which(x == "DEL" | x == "AMP"))
  )
  return(cna_df)
}


# ==============================================================================
# 2.  Cell Grouping by AS (High / Mid / Low)
# ==============================================================================
# Cells are labelled:
#   highAS  — AS >= 7
#   midAS   — AS in (3, 7)
#   lowAS   — AS <= 3
#
# Input:  sample_info — data.frame with an 'AS' column; rows = cells
# Output: same data.frame with an added 'pred' column
#
# For more stringent classification, apply alternative AS cutoffs.

divide_cells_by_AS <- function(sample_info) {
  sorted_df   <- sample_info[order(sample_info$AS), ]
  high_value  <- 7
  low_value   <- 3

  high_group <- sorted_df[sorted_df$AS >= high_value, ]
  low_group  <- sorted_df[sorted_df$AS <= low_value,  ]

  sorted_df$pred <- "midAS"
  sorted_df[rownames(high_group), "pred"] <- "highAS"
  sorted_df[rownames(low_group),  "pred"] <- "lowAS"
  return(sorted_df)
}


# ==============================================================================
# 3.  Per-Sample Subclone File Preparation
# ==============================================================================
# For each sample, this function:
#   (a) reads copy-number data and computes AS per cell,
#   (b) labels cells as highAS / midAS / lowAS using divide_cells_by_AS(),
#   (c) loads the normalized expression matrix,
#   (d) saves subclones_info.csv and norm_count_data_<cellTypes>.csv
#       to <outputPrefix>/<Study>/<Sample>/
#
# These intermediate files are required by the DGE and ssGSEA functions below.
#
# Inputs:
#   sampleIndex       — integer row index into sampleFilesMap
#   sample, paper     — identifiers
#   cellTypes         — "malignant" (or "all")
#   outputPrefix      — base output path
#   sampleFilesMap    — data.frame, see README.md for column descriptions

createSubclonesByCellASPerSampleFilteredForNorm <- function(
    sampleIndex, sample, paper, cellTypes, outputPrefix, sampleFilesMap
) {
  cells_df <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])

  if (cellTypes == "all") {
    sample_cells <- cells_df[cells_df$sample == sample, ]$cell_name
  } else {
    sample_cells <- cells_df[
      cells_df$sample == sample &
        tolower(cells_df$cell_type) == tolower(cellTypes), ]$cell_name
  }
  sample_cells <- as.character(sample_cells)

  dir.create(file.path(outputPrefix, paper), showWarnings = FALSE)
  dir.create(file.path(outputPrefix, paper, sample), showWarnings = FALSE)
  path <- file.path(outputPrefix, paper, sample)

  out_file <- file.path(path, "subclones_by_cells_AS_filteredForNorm_AS_thres.csv")
  if (file.exists(out_file)) return(invisible(NULL))

  # ---- CNA matrix ----
  cna_df <- read.csv(sampleFilesMap[sampleIndex, "Cna_mat"], sep = "\t") %>%
    remove_rownames() %>%
    column_to_rownames("sample")
  cna_df <- calculate_AS(cna_df, sample_cells)

  sample_info <- data.frame(AS = cna_df$AS, row.names = rownames(cna_df))
  sample_info$sample   <- sample
  sample_info$paper    <- paper
  sample_info$cell_id  <- rownames(sample_info)
  rownames(sample_info) <- gsub("[.:]", "-", rownames(sample_info))

  # ---- Normalized expression matrix ----
  norm_counts <- fread(sampleFilesMap[sampleIndex, "Norm_Expression_path"],
                       data.table = FALSE) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("[.:]", "-", colnames(norm_counts))

  cols <- gsub("^X", "", colnames(norm_counts))
  expr_filt <- norm_counts[, (colnames(norm_counts) %in% rownames(sample_info)) |
                               (cols %in% rownames(sample_info)), drop = FALSE] %>%
    na.omit()

  filt_cols   <- gsub("^X", "", colnames(expr_filt))
  sample_info <- sample_info[(rownames(sample_info) %in% colnames(expr_filt)) |
                               (rownames(sample_info) %in% filt_cols), ] %>%
    na.omit()

  sample_info <- divide_cells_by_AS(sample_info)
  write.csv(sample_info, out_file)
}


# ---- wrapper: run createSubclonesByCellASPerSampleFilteredForNorm in parallel ----
#
# Inputs: same as per-sample function, plus sampleFilesMap
createSubclonesByCellAS <- function(outputPrefix, sampleFilesMap, cellTypes) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr")
  ) %dopar% {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    message("***** ", paper, " ", sample, " *****")
    tryCatch(
      createSubclonesByCellASPerSampleFilteredForNorm(
        sampleIndex, sample, paper, cellTypes, outputPrefix, sampleFilesMap),
      error = function(e) message("Error in ", paper, "_", sample, ": ", e)
    )
  }
}


# ==============================================================================
# 4a.  Per-Sample AS Summary and t-test
# ==============================================================================
# Computes summary statistics and pairwise t-tests on the AS values of each
# cell group (highAS, midAS, lowAS) within a sample.
#
# Input:  path — directory containing subclones_info.csv
# Output: AS_summary.csv and AS_t_test_res.csv in the same directory

runAsTtest <- function(path) {
  colData <- fread(file.path(path, "subclones_info.csv"),
                   data.table = FALSE, header = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  rownames(colData) <- gsub("\\.", "-", rownames(colData))

  summary <- colData %>%
    group_by(pred) %>%
    get_summary_stats(AS, type = "mean_sd") %>%
    as.data.frame()
  fwrite(summary, file.path(path, "AS_summary.csv"), row.names = TRUE)

  pwc <- colData %>% pairwise_t_test(AS ~ pred, p.adjust.method = "BH")
  if (nrow(pwc) > 0) {
    effect <- colData %>% cohens_d(AS ~ pred)
    pwc$cohens_d_effect <- effect[["effsize"]]
    fwrite(pwc, file.path(path, "AS_t_test_res.csv"), row.names = TRUE)
  }
}


# ==============================================================================
# 4b.  Prepare Data Files for DGE Test
# ==============================================================================
# Reads cell-group labels (from subclones_info) and the normalized expression
# matrix, applies cell filtering, and saves both to disk for use by FindMarkers.
#
# NOTE: This function uses a sampleFilesMap global; adjust the argument list
# if you prefer to pass it explicitly.

createFilesPerSampleForTest <- function(
    sampleIndex, sample, paper, cellTypes, subclones_path, outputPrefix,
    sampleFilesMap
) {
  cells_df <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  if (cellTypes == "all") {
    sample_cells <- cells_df[cells_df$sample == sample, ]$cell_name
  } else {
    sample_cells <- cells_df[
      cells_df$sample == sample &
        tolower(cells_df$cell_type) == tolower(cellTypes), ]$cell_name
  }
  sample_cells <- na.omit(as.character(sample_cells))

  subclones_df <- read.csv(subclones_path) %>%
    column_to_rownames("cell_id") %>%
    filter(cell_id %in% sample_cells | rownames(.) %in% sample_cells) %>%
    mutate(pred = gsub("^[^_]*_", "", pred)) %>%
    filter(!pred %in% c("Unresolved", "Normal")) %>%
    drop_na()

  if (length(unique(subclones_df$pred)) < 2) stop("Only 1 subclone — skipping")

  dir.create(file.path(outputPrefix, paper), showWarnings = FALSE)
  dir.create(file.path(outputPrefix, paper, sample), showWarnings = FALSE)
  path <- file.path(outputPrefix, paper, sample)

  cna_df <- read.csv(sampleFilesMap[sampleIndex, "Cna_mat"], sep = "\t") %>%
    remove_rownames() %>%
    column_to_rownames("sample")
  cna_df <- cna_df[rownames(cna_df) %in% rownames(subclones_df), ]
  cna_df <- cna_df[rownames(subclones_df), ]
  cna_df$AS <- apply(
    cna_df[, !colnames(cna_df) %in% c("Xp", "Xq", "Yp", "Yq")],
    MARGIN = 1, FUN = function(x) length(which(x == "DEL" | x == "AMP"))
  )
  cna_df <- na.omit(cna_df)

  sample_info <- data.frame(pred = subclones_df$pred,
                             row.names = rownames(subclones_df))
  sample_info$sample  <- sample
  sample_info$paper   <- paper
  sample_info$AS      <- cna_df$AS
  rownames(sample_info) <- gsub("[.:]", "-", rownames(sample_info))

  norm_counts <- fread(sampleFilesMap[sampleIndex, "Norm_Expression_path"],
                       data.table = FALSE) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("[.:]", "-", colnames(norm_counts))

  cols      <- gsub("^X", "", colnames(norm_counts))
  expr_filt <- norm_counts[, (colnames(norm_counts) %in% rownames(sample_info)) |
                               (cols %in% rownames(sample_info)), drop = FALSE] %>%
    na.omit()
  filt_cols <- gsub("^X", "", colnames(expr_filt))

  sample_info <- sample_info[
    (rownames(sample_info) %in% colnames(expr_filt)) |
      (rownames(sample_info) %in% filt_cols), ] %>%
    na.omit()
  sample_info <- sample_info[
    ifelse(filt_cols != colnames(expr_filt), filt_cols, colnames(expr_filt)), ]
  if (all(filt_cols == rownames(sample_info))) {
    colnames(expr_filt) <- filt_cols
  }

  write.csv(sample_info, file.path(path, "subclones_info.csv"))
  write.csv(expr_filt,   file.path(path, paste0("norm_count_data_", cellTypes, ".csv")))
}


# ==============================================================================
# 5a.  Group-Based DGE: Seurat FindMarkers (highAS vs lowAS)
# ==============================================================================
# For each sample that has ≥30 cells in both the highAS and lowAS groups,
# runs Seurat FindMarkers (default: Wilcoxon test, min.pct = 0.25).
#
# The comparison is always highAS (ident.1) vs lowAS (ident.2).
#
# Inputs:
#   files_path    — base path containing <Study>/<Sample>/ subdirectories
#   cellTypes     — "malignant"
#   sampleFilesMap — sample metadata table
#   minPct        — min.pct argument to FindMarkers
#   test.use      — statistical test (default "wilcox")
#   output_path   — where to write results

runFindMarkersTest <- function(files_path, cellTypes, sampleFilesMap,
                               minPct, test.use, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "rstatix", "Seurat", "Matrix",
                  "msigdbr", "clusterProfiler", "ggpubr", "icesTAF")
  ) %dopar% {
    sample      <- sampleFilesMap[sampleIndex, "Sample"]
    paper       <- sampleFilesMap[sampleIndex, "Study"]
    cancer_type <- sampleFilesMap[sampleIndex, "Cancer_type"]
    if (sampleFilesMap[sampleIndex, "Expression_file_type"] == 2) return(NULL)

    t_test_path <- file.path(files_path, paper, sample, "AS_t_test_res.csv")
    if (!file.exists(t_test_path)) return(NULL)

    as_t_res <- fread(t_test_path, data.table = FALSE)
    tryCatch(
      runSeuratFindMarkersHighVsLow(
        as_t_res, files_path, cellTypes, sampleFilesMap,
        cancer_type, paper, sample, minPct, test.use, output_path
      ),
      error = function(e) message("FindMarkers error: ", paper, "_", sample, " — ", e)
    )
  }
}

# ---- internal: run Seurat FindMarkers for highAS vs lowAS ----
runSeuratFindMarkersHighVsLow <- function(
    as_t_res, files_path, cellTypes, sampleFilesMap,
    cancerType, paper, sample, minPct, test.use, output_path
) {
  min_cells <- 30

  for (i in seq_len(nrow(as_t_res))) {
    group1 <- as_t_res[i, "group1"]
    group2 <- as_t_res[i, "group2"]
    if (as_t_res[i, "p.adj"] >= 0.25) next

    # Determine reference (lower AS) and clone of interest (higher AS)
    ref            <- if (as_t_res[i, "cohens_d_effect"] < 0) group1 else group2
    cloneOfInterest <- if (ref == group1) group2 else group1

    # Only run highAS vs lowAS
    if (!(cloneOfInterest == "highAS" && ref == "lowAS")) next

    norm_counts <- fread(
      file.path(files_path, paper, sample,
                paste0("norm_count_data_", cellTypes, ".csv")),
      data.table = FALSE, header = TRUE
    ) %>%
      remove_rownames() %>%
      column_to_rownames("V1") %>%
      drop_na()
    colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))

    # Keep variable genes
    norm_counts$sd <- apply(norm_counts, 1, sd)
    norm_counts    <- norm_counts[norm_counts$sd != 0, ]
    norm_counts$sd <- NULL

    colData <- fread(
      file.path(files_path, paper, sample, "subclones_info.csv"),
      data.table = FALSE, header = TRUE
    ) %>%
      drop_na()
    colData[["V1"]] <- gsub("\\.", "-", colData[["V1"]])

    groupsColData <- colData[colData$pred %in% c(group1, group2), ] %>%
      remove_rownames() %>%
      column_to_rownames("V1")
    counts_by_group <- table(groupsColData$pred)

    if (!all(c(cloneOfInterest, ref) %in% names(counts_by_group)) ||
        any(counts_by_group < min_cells)) {
      message("Skipping ", paper, " ", sample, " — insufficient cells per group")
      next
    }

    groupsNormCount <- norm_counts[, rownames(groupsColData)]
    groupsNormCount$sd <- apply(groupsNormCount, 1, sd)
    groupsNormCount    <- groupsNormCount[groupsNormCount$sd != 0, ]
    groupsNormCount$sd <- NULL
    groupsNormCount    <- groupsNormCount[, rownames(groupsColData)]

    if (!all(colnames(groupsNormCount) == rownames(groupsColData))) {
      message("Column/row mismatch for ", paper, " ", sample, " — skipping")
      next
    }

    comp <- paste0(cloneOfInterest, "_vs_", ref)
    message("Running FindMarkers for ", nrow(groupsNormCount), " genes in ",
            paper, " ", sample)

    mat <- Matrix(as.matrix(groupsNormCount), sparse = TRUE)
    s   <- CreateSeuratObject(counts = mat, meta.data = groupsColData)
    s   <- SetAssayData(s, slot = "data", new.data = mat)
    s   <- FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE)
    s   <- ScaleData(s, features = rownames(s), verbose = FALSE)
    s   <- RunPCA(s, npcs = 30, verbose = FALSE) %>%
      RunUMAP(dims = 1:30, verbose = FALSE) %>%
      FindNeighbors(dims = 1:30, verbose = FALSE) %>%
      FindClusters(resolution = 0.7, verbose = FALSE)
    Idents(s) <- "pred"

    de <- FindMarkers(s, ident.1 = cloneOfInterest, ident.2 = ref,
                      min.pct = minPct, test.use = test.use,
                      min.cells.group = min_cells) %>%
      as.data.frame()
    de$group1 <- cloneOfInterest
    de$group2 <- ref

    sig <- filter(de, p_val_adj < 0.25) %>% arrange(p_val_adj)

    out_dir <- file.path(output_path, cancerType, paper, sample, comp)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(de,  file.path(out_dir, paste0("findMarkers_test_res_tbl_",
                                          minPct, "_", test.use,
                                          "_groupSize_", min_cells, ".csv")),
           row.names = TRUE)
    if (nrow(sig) > 0) {
      fwrite(sig, file.path(out_dir, paste0("findMarkers_test_sig_res_",
                                            minPct, "_", test.use,
                                            "_groupSize_", min_cells, ".csv")),
             row.names = TRUE)
    }
  }
}


# ==============================================================================
# 5b.  Pathway Enrichment on FindMarkers Results
# ==============================================================================
# For a given sample, reads the FindMarkers output and runs clusterProfiler
# enricher() on the up-regulated (positive) and down-regulated (negative)
# gene sets separately.
#
# Threshold: padj < 0.25 and |log2FC| > 0  (>10 significant genes required)

runSeuratFindMarkersEnricher <- function(files_path, cellTypes, sampleFilesMap,
                                         minPct, test.use, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "msigdbr", "clusterProfiler",
                  "Seurat", "Matrix", "icesTAF")
  ) %dopar% {
    sample      <- sampleFilesMap[sampleIndex, "Sample"]
    paper       <- sampleFilesMap[sampleIndex, "Study"]
    comp        <- "highAS_vs_lowAS"
    min_cells   <- 30
    path        <- file.path(output_path,
                             sampleFilesMap[sampleIndex, "Cancer_type"],
                             paper, sample, comp)
    res_file    <- file.path(path, paste0("findMarkers_test_res_tbl_",
                                          minPct, "_", test.use,
                                          "_groupSize_", min_cells, ".csv"))
    if (!file.exists(res_file)) return(NULL)
    tryCatch(
      {
        enricher_direction(sample, paper, minPct, test.use, min_cells,
                           "positive", path)
        enricher_direction(sample, paper, minPct, test.use, min_cells,
                           "negative", path)
      },
      error = function(e) message("Enricher error: ", paper, "_", sample, " — ", e)
    )
  }
}

# Input:
#   direction — "positive" (up in highAS) or "negative" (down in highAS)
#   path      — directory containing the FindMarkers result CSV
#   res_file  — full path to the result CSV (columns include gene, avg_log2FC,
#               p_val_adj)
enricher_direction <- function(sample, paper, minPct, test.use, min_cells,
                                direction, path,
                                res_file = NULL, gene_col = "gene",
                                fc_col   = "avg_log2FC",
                                padj_col = "p_val_adj") {
  if (is.null(res_file)) {
    res_file <- file.path(path, paste0("findMarkers_test_res_tbl_",
                                       minPct, "_", test.use,
                                       "_groupSize_", min_cells,
                                       "atLeast10SigGenes.csv"))
  }
  res_df <- fread(res_file, data.table = FALSE)
  names(res_df)[1] <- "gene"
  padj_thres <- 0.25

  genes_df <- if (direction == "positive") {
    res_df[res_df[[fc_col]] > 0 & res_df[[padj_col]] < padj_thres, ]
  } else {
    res_df[res_df[[fc_col]] < 0 & res_df[[padj_col]] < padj_thres, ]
  }

  if (nrow(genes_df) <= 10) return(invisible(NULL))

  enrich_dir <- file.path(path, paste0("findMarkers_enricher_", minPct,
                                        "_", direction))
  dir.create(enrich_dir, showWarnings = FALSE)

  for (i in seq_along(pathways_titles_list)) {
    m_t2g <- msigdbr(species = "Homo sapiens",
                     category    = pathways_categories_list[i],
                     subcategory = pathways_anontations_list[i]) %>%
      select(gs_name, gene_symbol)
    tryCatch({
      res <- enricher(genes_df$gene, TERM2GENE = m_t2g,
                      pvalueCutoff = 1, qvalueCutoff = 1)
      if (!is.null(res) && nrow(res@result) > 0) {
        write.csv(res@result,
                  file.path(enrich_dir, paste0(pathways_titles_list[i],
                                               "_enricher_res.csv")))
      }
    }, error = function(e)
      message("Enricher failed for DB ", pathways_titles_list[i], ": ", e))
  }
}


# ==============================================================================
# 6a.  Continuous Regression: lmFit (per sample)
# ==============================================================================
# Fits a linear model:  gene_expression ~ AS
# for each gene in each sample, to identify genes whose expression correlates
# continuously with the per-cell aneuploidy score.
#
# Input:  files_path   — base path with <Study>/<Sample>/ structure
#         sampleFilesMap
#         output_path  — where to write lmFit_Results_Genes_AS.csv

runLmFitTest <- function(files_path, sampleFilesMap, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "rstatix")
  ) %dopar% {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    if (sampleFilesMap[sampleIndex, "Expression_file_type"] == 2) return(NULL)
    tryCatch(
      runLmFit(files_path, sampleFilesMap, sampleFilesMap[sampleIndex, "Cancer_type"],
               paper, sample, output_path),
      error = function(e)
        message("LmFit error: ", paper, "_", sample, " — ", e)
    )
  }
}

runLmFit <- function(files_path, sampleFilesMap, cancerType, paper, sample,
                     output_path) {
  columns <- c("gene", "coefficiants", "p_value")

  norm_counts <- fread(
    file.path(files_path, paper, sample, "norm_count_data_malignant.csv"),
    data.table = FALSE, header = TRUE
  ) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))

  norm_counts$sd <- apply(norm_counts, 1, sd)
  norm_counts    <- norm_counts[norm_counts$sd != 0, ]
  norm_counts$sd <- NULL

  gene_expr           <- t(norm_counts) %>% as.data.frame()
  gene_expr$cellName  <- colnames(norm_counts)
  colnames(gene_expr) <- gsub("-", "_", colnames(gene_expr))

  colData <- fread(
    file.path(files_path, paper, sample, "subclones_info.csv"),
    data.table = FALSE, header = TRUE
  ) %>%
    drop_na()
  names(colData)[1] <- "cellName"
  colData[["cellName"]] <- gsub("\\.", "-", colData[["cellName"]])

  final_data <- merge(colData[, c("cellName", "AS")],
                      gene_expr, by = "cellName")
  res <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) <- columns

  for (i in 3:ncol(final_data)) {
    form   <- as.formula(paste(colnames(final_data)[i], "AS", sep = "~"))
    lm_fit <- lm(form, data = final_data)
    res    <- rbind(res, c(colnames(final_data)[i],
                            lm_fit$coefficients["AS"],
                            summary(lm_fit)$coefficients[, 4]["AS"]))
  }
  colnames(res) <- columns
  res$qval      <- p.adjust(res$p_value, method = "BH")

  out_dir <- file.path(output_path, cancerType, paper, sample)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(res, file.path(out_dir, "lmFit_Results_Genes_AS.csv"), row.names = TRUE)
}


# ==============================================================================
# 6b.  Pathway Enrichment on lmFit Results (per sample)
# ==============================================================================
# Same enricher_direction() helper as in Section 5b, but applied to per-sample
# lmFit coefficients (positive coefficient = gene increases with AS).

runEnricher_lm <- function(sampleFilesMap, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "msigdbr", "clusterProfiler",
                  "Seurat", "Matrix", "icesTAF")
  ) %dopar% {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    path   <- file.path(output_path,
                        sampleFilesMap[sampleIndex, "Cancer_type"],
                        paper, sample)
    lm_res_file <- file.path(path, "lmFit_Results_Genes_AS.csv")
    if (!file.exists(lm_res_file)) return(NULL)

    tryCatch({
      runLmEnricher_direction(paper, sample, "positive", path, lm_res_file)
      runLmEnricher_direction(paper, sample, "negative", path, lm_res_file)
    },
    error = function(e) message("LmEnricher error: ", paper, "_", sample, " — ", e))
  }
}

runLmEnricher_direction <- function(paper, sample, direction, path, res_file) {
  res_df     <- fread(res_file, data.table = FALSE)
  padj_thres <- 0.25
  genes_df   <- if (direction == "positive") {
    res_df[res_df$coefficiants > 0 & res_df$qval < padj_thres, ]
  } else {
    res_df[res_df$coefficiants < 0 & res_df$qval < padj_thres, ]
  }
  if (nrow(genes_df) == 0) return(invisible(NULL))

  enrich_dir <- file.path(path, "lm", paste0("lm_enricher_", direction))
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(pathways_titles_list)) {
    m_t2g <- msigdbr(species = "Homo sapiens",
                     category    = pathways_categories_list[i],
                     subcategory = pathways_anontations_list[i]) %>%
      select(gs_name, gene_symbol)
    tryCatch({
      res <- enricher(genes_df$gene, TERM2GENE = m_t2g,
                      pvalueCutoff = 1, qvalueCutoff = 1)
      if (!is.null(res) && nrow(res@result) > 0) {
        write.csv(res@result,
                  file.path(enrich_dir,
                            paste0(pathways_titles_list[i], "_lm_enricher_res.csv")))
      }
    }, error = function(e)
      message("LmEnricher DB ", pathways_titles_list[i], " failed: ", e))
  }
}


# ==============================================================================
# 7a.  ssGSEA: Per-Cell Pathway Enrichment Scores
# ==============================================================================
# Runs ssGSEA (via GSVA) on the normalized expression matrix of each sample,
# producing a cells × pathways score matrix.
#
# This function is the *single shared ssGSEA entry point*. It is called with
# different inputs for:
#   - The global high-AS vs low-AS analysis (Sec 7b–7c)
#   - The karyotypic heterogeneity stratified analysis (Sec 8e)
#   - The recurrent aneuploidy analysis (02_recurrent_aneuploidies.R)
#
# Input:
#   files_path     — path containing <Study>/<Sample>/norm_count_data_*.csv
#   cellTypes      — "malignant"
#   cancerType, paper, sample — identifiers
#   output_path    — where to write ssGsea_all_<DB>_res.csv
#   pathway_index  — integer index(es) into pathways_titles_list to run
#                    (default: 6, i.e. HALLMARK only)

ssGsea_all_set <- function(files_path, cellTypes, cancerType, paper, sample,
                            output_path, pathway_index = 6L) {
  norm_counts <- fread(
    file.path(files_path, paper, sample,
              paste0("norm_count_data_", cellTypes, ".csv")),
    data.table = FALSE, header = TRUE
  ) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("[.:]", "-", colnames(norm_counts))

  norm_counts$sd <- apply(norm_counts, 1, sd)
  norm_counts    <- norm_counts[norm_counts$sd != 0, ]
  norm_counts$sd <- NULL

  dir.create(file.path(output_path, paper, sample),
             recursive = TRUE, showWarnings = FALSE)

  for (i in pathway_index) {
    out_file <- file.path(output_path, paper, sample,
                          paste0("ssGsea_all_", pathways_titles_list[i], "_res.csv"))
    if (file.exists(out_file)) next

    m_t2g  <- msigdbr(species = "Homo sapiens",
                      category    = pathways_categories_list[i],
                      subcategory = pathways_anontations_list[i]) %>%
      select(gs_name, gene_symbol)
    gene_sets <- split(m_t2g$gene_symbol, m_t2g$gs_name)

    mat <- Matrix(as.matrix(norm_counts), sparse = TRUE)
    res <- GSVA::gsva(mat, gene_sets, verbose = FALSE, method = "ssgsea") %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()

    fwrite(res, out_file, row.names = TRUE)
  }
}

# ---- merge ssGSEA scores with cell metadata ----
createMergedSampleInfo_withSsgsea <- function(files_path, sample, paper,
                                               output_path, pathway_index = 6L) {
  sample_info <- read.csv(file.path(files_path, paper, sample,
                                    "subclones_info.csv")) %>%
    remove_rownames() %>%
    column_to_rownames("X")

  path <- file.path(output_path, paper, sample)
  merged <- data.frame()

  for (i in pathway_index) {
    f <- file.path(path, paste0("ssGsea_all_", pathways_titles_list[i], "_res.csv"))
    if (!file.exists(f)) next
    scores <- read.csv(f) %>%
      remove_rownames() %>%
      column_to_rownames("X")
    # Keep only pathway columns matching the DB prefix
    keep_cols <- grep(paste0(pathways_titles_list[i], "_"),
                      colnames(scores), value = TRUE)
    scores <- scores[, keep_cols, drop = FALSE]

    if (nrow(merged) > 0) {
      merged <- merge(merged, scores, by = 0) %>%
        remove_rownames() %>%
        column_to_rownames("Row.names")
    } else {
      merged <- scores
    }
  }

  out <- merge(sample_info, merged, by.x = 0, by.y = 0)
  write.csv(out, file.path(path,
                            paste0("subclones_info_with__ssgsea_",
                                   pathways_titles_list[max(pathway_index)],
                                   ".csv")),
            row.names = TRUE)
}


# ==============================================================================
# 7b.  ssGSEA t-test: highAS vs lowAS per sample
# ==============================================================================
# For each sample that has ≥30 highAS and ≥30 lowAS cells:
#   - Loads the merged subclone + ssGSEA scores file
#   - Runs Welch's t-test for each pathway (highAS vs lowAS cells)
#   - Computes Cohen's d effect size
#   - Saves results to ssgsea_HALLMARK_tTestRes.csv
#
# Input:
#   file_path      — path to ssGSEA merged files
#   output_path    — where to write t-test results
#   sampleFilesMap

ssGSEA_tTest_perSample <- function(sampleFilesMap, file_path, output_path,
                                    excluded_studies = c()) {
  for (sampleIndex in seq_len(nrow(sampleFilesMap))) {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    if (paper %in% excluded_studies) next

    merged_file <- file.path(file_path, paper, sample,
                              "subclones_info_with__ssgsea_HALLMARK.csv")
    highlow_dir <- # must contain the highAS_vs_lowAS comparison
      file.path(file_path, sampleFilesMap[sampleIndex, "Cancer_type"],
                paper, sample, "highAS_vs_lowAS")

    if (!file.exists(merged_file) || !file.exists(highlow_dir)) next

    colData <- read.csv(merged_file) %>%
      remove_rownames() %>%
      column_to_rownames("Row.names") %>%
      .[, !colnames(.) %in% c("X.1", "X")] %>%
      mutate(sample = as.character(sample))

    pathway_columns <- grep("HALLMARK_", colnames(colData), value = TRUE)
    results <- data.frame(
      pathway = character(), p_value = numeric(),
      effect_size = numeric(), stderr = numeric(),
      mean_difference = numeric(), sampleID = character(),
      stringsAsFactors = FALSE
    )

    for (pw in pathway_columns) {
      low_vals  <- colData[colData$pred == "lowAS",  pw]
      high_vals <- colData[colData$pred == "highAS", pw]
      tt        <- t.test(high_vals, low_vals, var.equal = FALSE)
      cd        <- cohens_d(high_vals, low_vals, pooled_sd = FALSE)
      results   <- rbind(results, data.frame(
        pathway        = pw,
        p_value        = tt$p.value,
        effect_size    = cd$Cohens_d,
        stderr         = tt$stderr,
        mean_difference = diff(tt$estimate),
        sampleID       = paste0(paper, "_", sample),
        stringsAsFactors = FALSE
      ))
    }
    results$qvalues <- p.adjust(results$p_value, method = "BH")

    dir.create(file.path(output_path, paper, sample),
               recursive = TRUE, showWarnings = FALSE)
    write.csv(results, file.path(output_path, paper, sample,
                                  "ssgsea_HALLMARK_tTestRes.csv"))
  }
}


# ==============================================================================
# 7c.  Random-Effects Meta-Analysis Across Samples
# ==============================================================================
# Aggregates per-sample ssGSEA t-test results using a random-effects model
# (REML) via metafor::rma.  Produces one pooled effect size per pathway.
#
# Input:
#   output_path — directory containing per-sample ssgsea_HALLMARK_tTestRes.csv
#                 files (recursive search)
#   result_file — where to write the summary (default: same directory)

ssgsea_metaAnalysis <- function(output_path,
                                 result_file = file.path(output_path,
                                                          "ssgsea_meta_analysis_tTest_HALLMARKresults.csv"),
                                 metadata_path = NULL) {
  files <- list.files(output_path, pattern = "ssgsea_HALLMARK_tTestRes.csv",
                      full.names = TRUE, recursive = TRUE)
  combined <- bind_rows(lapply(files, read.csv), .id = "id") %>%
    mutate(sampleID = gsub("[. -]", "_", sampleID))

  meta_data <- combined %>%
    group_by(pathway) %>%
    summarise(
      EffectSize = list(effect_size),
      SE         = list(stderr),
      SampleIDs  = list(sampleID),
      .groups    = "drop"
    )

  meta_results <- vector("list", nrow(meta_data))
  for (k in seq_len(nrow(meta_data))) {
    es  <- unlist(meta_data$EffectSize[[k]])
    ses <- unlist(meta_data$SE[[k]])
    meta_results[[k]] <- rma(yi = es, sei = ses, method = "REML")
  }
  names(meta_results) <- meta_data$pathway

  summary_results <- data.frame(
    Pathway         = names(meta_results),
    PooledEffectSize = sapply(meta_results, \(x) x$b),
    CI_Lower        = sapply(meta_results, \(x) x$ci.lb),
    CI_Upper        = sapply(meta_results, \(x) x$ci.ub),
    PValue          = sapply(meta_results, \(x) x$pval)
  )
  summary_results$Qvalues <- p.adjust(summary_results$PValue, method = "BH")
  write.csv(summary_results, result_file, row.names = FALSE)
  invisible(summary_results)
}


# ==============================================================================
# 8a.  Pseudo-Bulk Mean Expression per Sample
# ==============================================================================
# Computes the per-sample mean normalized expression (across all cells) for
# downstream sample-level regression.
#
# Input:  sampleIndex, sample, paper, outputPrefix, sampleFilesMap
# Output: mean_norm_counts_without_center_data.csv per sample

createMeanNormWithoutCenterCountPerSample <- function(
    sampleIndex, sample, paper, outputPrefix, sampleFilesMap
) {
  dir.create(file.path(outputPrefix, paper, sample),
             recursive = TRUE, showWarnings = FALSE)
  path <- file.path(outputPrefix, paper, sample)

  norm_counts <- fread(
    sampleFilesMap[sampleIndex, "Norm_Without_Center_Expression_path"],
    data.table = FALSE
  ) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))

  mean_expr           <- data.frame(rowMeans(norm_counts))
  colnames(mean_expr) <- sample
  write.csv(mean_expr, file.path(path,
                                  "mean_norm_counts_without_center_data.csv"))
}


# ==============================================================================
# 8b.  Merge Pseudo-Bulk Expression Across Samples
# ==============================================================================
# Reads per-sample mean expression files and concatenates into a
# genes × samples matrix.

createMergeMeanNormWithoutCenterExpFile <- function(outputPrefix, sampleFilesMap,
                                                     excluded_studies = c()) {
  df <- data.frame()
  for (sampleIndex in seq_len(nrow(sampleFilesMap))) {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    if (paper %in% excluded_studies) next
    path <- file.path(outputPrefix, paper, sample,
                      "mean_norm_counts_without_center_data.csv")
    if (!file.exists(path)) next
    d <- read.csv(path) %>% column_to_rownames("X")
    colnames(d) <- paste0(paper, "/", sample)
    if (nrow(df) == 0) {
      df <- d
    } else {
      df       <- merge(df, d, by = 0, all = TRUE)
      df[is.na(df)] <- 0
      df <- df %>% column_to_rownames("Row.names")
    }
  }
  write.csv(df, file.path(outputPrefix, "MeanNormWithoutCenterExp_withMinExpGenes.csv"))
}


# ==============================================================================
# 8c.  Pseudo-Bulk Linear Regression (sample-level)
# ==============================================================================
# Fits:  gene_expression ~ Malignant_Mean_AS + Cancer_type
# for each gene using the pseudo-bulk mean expression values.
# The Malignant_Mean_AS covariate is the mean AS across malignant cells.
#
# Input:
#   metadata       — data.frame with columns: sample_name (sample ID),
#                    Malignant_Mean_AS, Cancer_type
#   gene_expression — genes × samples matrix (samples as rows after t())
# Output: lmFit_Cancer_type_Results_Genes_*.csv

pseudoBulk_lmFit <- function(metadata, gene_expression, output_file) {
  columns      <- c("gene", "coefficiants", "p_value")
  score_option <- "Malignant_Mean_AS"

  gene_expr_t           <- t(gene_expression) %>% as.data.frame()
  gene_expr_t$name      <- rownames(gene_expr_t)
  colnames(gene_expr_t) <- gsub("-", "_", colnames(gene_expr_t))

  gene_expr_t <- gene_expr_t[rownames(gene_expr_t) %in% metadata$sample_name, ]
  final_data  <- merge(metadata[, c("sample_name", score_option, "Cancer_type")],
                       gene_expr_t, by.x = "sample_name", by.y = "name")

  res <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) <- columns
  j <- 4   # first gene column

  for (i in j:ncol(final_data)) {
    form   <- as.formula(paste0("`",
      colnames(final_data)[i], , "` ~ ",
      paste(score_option, "Cancer_type", sep = "+"),
      sep = "~"
    ))
    lm_fit <- lm(form, data = final_data)
    res    <- rbind(res, c(
      colnames(final_data)[i],
      lm_fit$coefficients[score_option],
      summary(lm_fit)$coefficients[, 4][score_option]
    ))
  }
  colnames(res) <- columns
  res$qval      <- p.adjust(res$p_value, method = "BH")
  res           <- drop_na(res)

  write.table(res, output_file, sep = ",", row.names = FALSE)
  invisible(res)
}


# ==============================================================================
# 8d.  Pseudo-Bulk Enrichment
# ==============================================================================
# Wrapper that reads the pseudo-bulk lmFit output, builds a ranked gene list
# (filtered at q < 0.25), and runs clusterProfiler::enricher() on the
# positively and negatively associated genes.
#
# This re-uses the same enricher_direction() helper from Section 5b/6b.

pseudoBulk_enrichment <- function(lmFit_file, output_dir) {
  res <- read.csv(lmFit_file)
  names(res)[1] <- "gene"
  res$gene      <- gsub("\\.{2}.*", "", res$gene)
  sig_res       <- res[res$qval < 0.25, ]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (direction in c("positive", "negative")) {
    genes <- if (direction == "positive") {
      sig_res[sig_res$coefficiants > 0, "gene"]
    } else {
      sig_res[sig_res$coefficiants < 0, "gene"]
    }
    if (length(genes) == 0) next

    enrich_dir <- file.path(output_dir, paste0("enricher_", direction))
    dir.create(enrich_dir, showWarnings = FALSE)
    for (i in seq_along(pathways_titles_list)) {
      m_t2g <- msigdbr(species = "Homo sapiens",
                       category    = pathways_categories_list[i],
                       subcategory = pathways_anontations_list[i]) %>%
        select(gs_name, gene_symbol)
      tryCatch({
        res_enr <- enricher(genes, TERM2GENE = m_t2g, pvalueCutoff = 1)
        if (!is.null(res_enr) && nrow(res_enr@result) > 0) {
          write.csv(res_enr@result,
                    file.path(enrich_dir,
                              paste0(pathways_titles_list[i], "_enricher_res.csv")))
        }
      }, error = function(e)
        message("Enricher (pseudo-bulk) failed for DB ", pathways_titles_list[i]))
    }
  }
}


# ==============================================================================
# 9.  Karyotypic Heterogeneity Score (AHS)
# ==============================================================================
# The AHS is defined as 1 − mean(pairwise Pearson correlation of CNA profiles
# among malignant cells in a sample).  Higher AHS = more diverse CNA profiles.
#
# calculateCNHeterogeneity_pairwiseCor: computes the cell × cell correlation
#   matrix for a single sample and writes it to disk.
#
# summarize_CNheterogeneity: reads all per-sample correlation matrices and
#   computes AHS statistics (mean) per sample.

calculateCNHeterogeneity_pairwiseCor <- function(
    sampleIndex, sample, paper, outputPrefix, sampleFilesMap
) {
  cells_df     <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  sample_cells <- cells_df[
    cells_df$sample == sample &
      tolower(cells_df$cell_type) == "malignant", "cell_name"] %>%
    gsub("\\.", "-", .)

  gene_cna <- fread(sampleFilesMap[sampleIndex, "Gene_cna_mat"],
                    data.table = FALSE, sep = " ") %>%
    remove_rownames() %>%
    column_to_rownames("V1")
  colnames(gene_cna) <- gsub("^X|\\.|-", c("", "", "-"), colnames(gene_cna))

  gene_cna <- gene_cna[, colnames(gene_cna) %in% sample_cells] %>% na.omit()
  corr_mat <- cor(gene_cna, method = "pearson")

  dir.create(file.path(outputPrefix, paper, sample),
             recursive = TRUE, showWarnings = FALSE)
  write.csv(corr_mat,
            file.path(outputPrefix, paper, sample,
                      "cells_CN_pairwise_correlation_df_malignantCells.csv"))
}

summarize_CNheterogeneity <- function(outputPrefix, metadata) {
  df <- data.frame()
  for (i in seq_len(nrow(metadata))) {
    paper  <- metadata[i, "Paper"]
    sample <- metadata[i, "Sample"]
    f <- file.path(outputPrefix, paper, sample,
                   "cells_CN_pairwise_correlation_df_malignantCells.csv")
    if (!file.exists(f)) next
    corr <- fread(f, data.table = FALSE) %>%
      remove_rownames() %>%
      column_to_rownames("V1")
    vals <- corr[upper.tri(corr, diag = FALSE)]
    df   <- bind_rows(df, data.frame(
      Paper  = paper,  Sample = sample,
      mean = mean(vals),
      stringsAsFactors = FALSE
    ))
  }
  write.csv(df, file.path(outputPrefix,
                            "samples_CNheterogeneity_pairwaiseCor_allGenes_df.csv"),
            row.names = FALSE)
  invisible(df)
}


# ==============================================================================
# 10. CN Subclone Identification (Louvain / Seurat)
# ==============================================================================
# Each sample's gene-level CNA matrix is filtered to the top 67% most variable
# genes, then clustered with Seurat's Louvain algorithm.
# Adjacent clusters with similar arm-level CNA profiles (using a ±0.15 threshold)
# are merged into a single subclone.
#
# Output: sample_clusters_louvain_seurat.csv with columns:
#   cell_id, cluster (original), MergedCluster (post-merge), sample, cell_type

LouvainClustering_seurat <- function(outputPrefix, sampleFilesMap,
                                      chromosome_arm_info) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "Seurat", "FNN", "igraph")
  ) %dopar% {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    message(paper, " ", sample)
    tryCatch(
      louvain_one_sample(sampleIndex, sample, paper, outputPrefix,
                         sampleFilesMap, chromosome_arm_info),
      error = function(e) message("Louvain error: ", paper, "_", sample, " — ", e)
    )
  }
}

louvain_one_sample <- function(sampleIndex, sample, paper, outputPrefix,
                                sampleFilesMap, chromosome_arm_info) {
  cells_df     <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  sample_cells <- cells_df[
    cells_df$sample == sample &
      tolower(cells_df$cell_type) == "malignant", "cell_name"] %>%
    gsub("\\.", "-", .)

  gene_cna <- fread(sampleFilesMap[sampleIndex, "Gene_cna_mat"],
                    data.table = FALSE, sep = " ") %>%
    remove_rownames() %>%
    column_to_rownames("V1")
  colnames(gene_cna) <- gsub("^X", "", gsub("\\.", "-", colnames(gene_cna)))
  gene_cna <- gene_cna[, colnames(gene_cna) %in% sample_cells] %>% na.omit()

  # Keep top 67% of genes by mean absolute CNA signal
  gene_signal  <- rowMeans(abs(gene_cna))
  top_n        <- ceiling(length(gene_signal) * 0.67)
  top_genes    <- names(sort(gene_signal, decreasing = TRUE))[seq_len(top_n)]
  cn_matrix    <- gene_cna[top_genes, ]

  # Seurat clustering
  seu <- CreateSeuratObject(counts = cn_matrix)
  seu <- SetAssayData(seu, assay = "RNA", slot = "data", new.data = cn_matrix)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = rownames(seu),
                npcs = min(ncol(cn_matrix) - 1, 30))
  seu <- FindNeighbors(seu, features = rownames(seu),
                       k.param = min(ncol(cn_matrix) - 1, 15),
                       dims    = seq_len(min(ncol(cn_matrix) - 1, 15)))
  seu <- FindClusters(seu, resolution = 0.5)

  sample_clusters <- data.frame(
    cell_id = colnames(seu),
    cluster = as.numeric(as.character(seu$seurat_clusters)) + 1
  )

  # ---- Arm-level CNA per cluster ----
  cna_matrix <- cn_matrix[rownames(cn_matrix) %in% chromosome_arm_info$gene, ]
  gene_info  <- chromosome_arm_info %>%
    filter(gene %in% rownames(cna_matrix)) %>%
    .[match(rownames(cna_matrix), .$gene), ]

  clusters  <- setNames(sample_clusters$cluster, sample_clusters$cell_id)
  avg_by_arm <- lapply(unique(gene_info$chrom_arm), function(arm) {
    genes   <- gene_info$gene[gene_info$chrom_arm == arm]
    if (length(genes) < 2) return(NULL)
    mat_avg <- colMeans(cna_matrix[genes, , drop = FALSE])
    data.frame(Cell = names(mat_avg), Arm = arm, AvgCNA = mat_avg)
  }) %>%
    bind_rows() %>%
    mutate(Cluster = clusters[Cell])

  arm_means <- avg_by_arm %>%
    group_by(Cluster, Arm) %>%
    summarise(MeanCNA = mean(AvgCNA), .groups = "drop") %>%
    pivot_wider(names_from = Arm, values_from = MeanCNA) %>%
    as.data.frame()
  rownames(arm_means) <- arm_means$Cluster
  arm_means$Cluster   <- NULL

  assign_status <- function(x) ifelse(x > 0.15, "Amp", ifelse(x < -0.15, "Del", "Neu"))
  arm_status    <- as.data.frame(lapply(arm_means, assign_status))
  rownames(arm_status) <- rownames(arm_means)

  # ---- Merge similar clusters ----
  merged_clusters <- rownames(arm_status)
  if (nrow(arm_status) > 1) {
    for (i in seq_len(nrow(arm_status) - 1)) {
      for (j in (i + 1):nrow(arm_status)) {
        c1 <- rownames(arm_status)[i]
        c2 <- rownames(arm_status)[j]
        identical_status <- all(arm_status[c1, ] == arm_status[c2, ])
        max_diff         <- max(abs(arm_means[c1, ] - arm_means[c2, ]), na.rm = TRUE)
        if (identical_status || max_diff < 0.15) {
          merged_clusters[merged_clusters == c2] <- c1
        }
      }
    }
  }
  names(merged_clusters)        <- as.character(seq_along(merged_clusters))
  sample_clusters$MergedCluster <- merged_clusters[as.character(sample_clusters$cluster)]

  cells_df   <- cells_df[cells_df$sample == sample, ]
  cells_df$cell_name <- gsub("\\.", "-", cells_df$cell_name)

  out <- merge(sample_clusters[, c("cell_id", "cluster", "MergedCluster")],
               cells_df[, c("cell_name", "sample", "cell_type")],
               by.x = "cell_id", by.y = "cell_name", all.x = TRUE)

  dir.create(file.path(outputPrefix, paper, sample),
             recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file.path(outputPrefix, paper, sample,
                            "sample_clusters_louvain_seurat.csv"))
}


# ==============================================================================
# 11. Cell Cycle Scoring (QC / Validation)
# ==============================================================================
# Scores cells for S-phase and G2/M-phase activity using Seurat's
# CellCycleScoring().  G1 cells (Phase == "G1") are subsequently used
# as a cell-cycle-controlled subset for the ssGSEA and pseudo-bulk analyses.
#
# Gene signatures follow Tirosh et al. 2016 (included inline below).

cellCycle_G2M <- c(
  "TOP2A","UBE2C","HMGB2","NUSAP1","CENPF","CCNB1","TPX2","CKS2","BIRC5","PRC1",
  "PTTG1","KPNA2","MKI67","CDC20","CDK1","CCNB2","CDKN3","SMC4","NUF2","ARL6IP1",
  "CKAP2","ASPM","PLK1","CKS1B","CCNA2","AURKA","MAD2L1","GTSE1","HMMR","UBE2T",
  "CENPE","CENPA","KIF20B","AURKB","CDCA3","CDCA8","UBE2S","KNSTRN","KIF2C","PBK",
  "TUBA1B","DLGAP5","TACC3","STMN1","DEPDC1","ECT2","CENPW","ZWINT","HIST1H4C","KIF23"
)

cellCycle_G1S <- c(
  "PCNA","RRM2","FEN1","GINS2","TYMS","MCM3","GMNN","HIST1H4C","CLSPN","ATAD2",
  "TK1","KIAA0101","DUT","HELLS","MCM7","UBE2T","MCM4","CENPU","DHFR","ZWINT",
  "ASF1B","MCM5","DNAJC9","RFC4","HMGB2","CDC6","RRM1","ORC6","CDK1","RAD51AP1",
  "RNASEH2A","CHAF1A","CENPK","CDCA5","SLBP","MCM6","TMEM106C","CENPM","MYBL2","E2F1",
  "USP1","DNMT1","PKMYT1","MAD2L1","PSMC3IP","CDCA4","RFC2","CDC45","UHRF1","MCM2"
)

# Input:
#   sampleIndex, sample, paper — identifiers
#   sampleFilesMap             — sample metadata table
# Output: data.frame with S.Score, G2M.Score, Phase, cell_id, AS per cell

score_cell_cycle_one_sample <- function(sampleIndex, sample, paper, sampleFilesMap) {
  cells_df     <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  sample_cells <- cells_df[
    cells_df$sample == sample &
      tolower(cells_df$cell_type) == "malignant", "cell_name"] %>%
    gsub("[.:]", "-", .)

  cna_df <- read.csv(sampleFilesMap[sampleIndex, "Cna_mat"], sep = "\t") %>%
    remove_rownames() %>%
    column_to_rownames("sample")
  cna_df$AS <- apply(
    cna_df[, !colnames(cna_df) %in% c("Xp", "Xq", "Yp", "Yq")],
    MARGIN = 1, FUN = function(x) length(which(x == "DEL" | x == "AMP"))
  )
  rownames(cna_df) <- gsub("[.:]", "-", rownames(cna_df))
  cna_df <- cna_df[rownames(cna_df) %in% sample_cells, ]

  norm_counts <- fread(
    sampleFilesMap[sampleIndex, "Norm_Without_Center_Expression_path"],
    data.table = FALSE
  ) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    as.data.frame()
  colnames(norm_counts) <- gsub("[.:]", "-", gsub("^X", "", colnames(norm_counts)))
  norm_counts <- norm_counts[, colnames(norm_counts) %in% sample_cells] %>% drop_na()
  rownames(norm_counts) <- gsub("_", "-", rownames(norm_counts))

  mat <- as(as.matrix(norm_counts), "dgCMatrix")
  seu <- CreateSeuratObject(counts = mat)
  seu <- SetAssayData(seu, slot = "data", new.data = mat)

  s_genes   <- intersect(cellCycle_G1S, rownames(seu))
  g2m_genes <- intersect(cellCycle_G2M, rownames(seu))
  seu       <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes)

  meta          <- seu@meta.data
  meta$cell_id  <- rownames(meta)
  meta$sample   <- sample
  meta$paper    <- paper
  meta$AS       <- cna_df$AS[match(meta$cell_id, rownames(cna_df))]
  meta
}

# ==============================================================================
# 12.  Per-Cell QC Metrics & Aneuploidy Score Collection
# ==============================================================================
# Iterates over all samples in sampleFilesMap, loading expression matrices
# and CNA calls to compute per-cell QC metrics and Aneuploidy Scores.
#
# Input:  sampleFilesMap 
#
# Output: df — data.frame accumulated across all samples, columns:
#           paper, sample, cell_id, cell_type,
#           expressed_genes_per_cell, library_size, percent_mito, AS

collect_qc_and_AS <- function(sampleFilesMap) {
  df <- data.frame()
  
  for (sampleIndex in seq_len(nrow(sampleFilesMap))) {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    print(paste0("***** ", paper, " ", sample, " *****"))
    
    cells_df <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
    genes_df <- read.csv(sampleFilesMap[sampleIndex, "Genes_path"], header = FALSE)
    
    exp_data <- ReadMtx(
      sampleFilesMap[sampleIndex, "Expression_path"],
      cells          = gsub("csv", "tsv", sampleFilesMap[sampleIndex, "Cells_path"]),
      features       = gsub("txt", "tsv", sampleFilesMap[sampleIndex, "Genes_path"]),
      cell.column    = 1,
      feature.column = 1,
      skip.cell      = 1
    )
    
    sample_cells      <- cells_df[cells_df$sample == sample, ]
    sample_cells_name <- gsub("[.:]", "-", as.character(sample_cells$cell_name))
    
    sample_expression <- as.data.frame(
      exp_data[, colnames(exp_data) %in% sample_cells_name]
    )
    colnames(sample_expression) <- gsub("^X", "", gsub("[.:]", "-", colnames(sample_expression)))
    
    cna_df <- read.csv(sampleFilesMap[sampleIndex, "Cna_mat"], sep = "\t") %>%
      remove_rownames() %>%
      column_to_rownames(var = "sample")
    cna_df <- calculate_AS(cna_df, sample_cells_name)
    
    # QC metrics
    expressed_genes_per_cell <- colSums(sample_expression > 0)
    library_size             <- colSums(sample_expression, na.rm = TRUE)
    mito_genes               <- grep("^MT-", rownames(sample_expression), value = TRUE)
    mito_counts              <- colSums(sample_expression[mito_genes, , drop = FALSE])
    percent_mito             <- (mito_counts / library_size) * 100
    
    sample_cells$cell_name <- gsub("[.:]", "-", sample_cells$cell_name)
    
    sample_info <- data.frame(expressed_genes_per_cell,
                              row.names = names(expressed_genes_per_cell))
    sample_info$sample       <- sample
    sample_info$paper        <- paper
    sample_info$cell_id      <- rownames(sample_info)
    sample_info$library_size <- library_size[match(sample_info$cell_id, names(library_size))]
    sample_info$percent_mito <- percent_mito[match(sample_info$cell_id, names(percent_mito))]
    sample_info$cell_type    <- sample_cells$cell_type[match(sample_info$cell_id, sample_cells$cell_name)]
    sample_info$AS           <- cna_df$AS[match(sample_info$cell_id, rownames(cna_df))]
    
    rownames(sample_info) <- gsub("[.:]", "-", rownames(sample_info))
    
    df <- if (nrow(df) == 0) sample_info else bind_rows(df, sample_info)
  }
  
  df$name <- paste0(df$paper, " ", df$sample)
  return(df)
}

# ==============================================================================
# 13. AHS-Stratified Meta-Analysis 
# ==============================================================================
# Re-runs the ssGSEA meta-analysis separately on:
#   - "homogeneous" samples: bottom tercile of AHS (≤ 33rd percentile)
#   - "heterogeneous" samples: top tercile of AHS (≥ 67th percentile)
# This tests whether the transcriptional consequences of aneuploidy differ
# based on karyotypic heterogeneity.
#
# Input:
#   ssgsea_ttest_path  — directory containing per-sample ssgsea_tTestRes_<DB>.csv
#   ahs_df             — data.frame with columns: name (sample ID), AHS
#   result_file_homo   — output CSV for homogeneous samples
#   result_file_hetro  — output CSV for heterogeneous samples
#   pathway_index      — integer index into pathways_titles_list (default 6 = HALLMARK)

ssgsea_metaAnalysis_AHS_stratified <- function(ssgsea_ttest_path,
                                               ahs_df,
                                               result_file_homo,
                                               result_file_hetro,
                                               pathway_index = 6L) {
  quantiles <- quantile(ahs_df$AHS, probs = c(1/3, 2/3), na.rm = TRUE)
  homo_names  <- ahs_df$name[ahs_df$AHS <= quantiles[1]]
  hetro_names <- ahs_df$name[ahs_df$AHS >= quantiles[2]]
  
  run_meta <- function(sample_names, result_file) {
    files <- list.files(ssgsea_ttest_path,
                        pattern = paste0("ssgsea_tTestRes_",
                                         pathways_titles_list[pathway_index], ".csv"),
                        full.names = TRUE, recursive = TRUE)
    combined <- bind_rows(lapply(files, read.csv), .id = "id") %>%
      mutate(sampleID = gsub("[. -]", "_", sampleID)) %>%
      filter(sampleID %in% gsub("[. -]", "_", sample_names))
    
    meta_data <- combined %>%
      group_by(pathway) %>%
      summarise(EffectSize = list(effect_size), SE = list(stderr), .groups = "drop")
    
    meta_results <- vector("list", nrow(meta_data))
    for (k in seq_len(nrow(meta_data))) {
      es  <- unlist(meta_data$EffectSize[[k]])
      ses <- unlist(meta_data$SE[[k]])
      meta_results[[k]] <- rma(yi = es, sei = ses, method = "REML")
    }
    names(meta_results) <- meta_data$pathway
    
    summary_results <- data.frame(
      Pathway          = names(meta_results),
      PooledEffectSize = sapply(meta_results, \(x) x$b),
      CI_Lower         = sapply(meta_results, \(x) x$ci.lb),
      CI_Upper         = sapply(meta_results, \(x) x$ci.ub),
      PValue           = sapply(meta_results, \(x) x$pval),
      I2               = sapply(meta_results, \(x) x$I2),
      Tau2             = sapply(meta_results, \(x) x$tau2),
      Q                = sapply(meta_results, \(x) x$QE),
      QPValue          = sapply(meta_results, \(x) x$QEp)
    )
    summary_results$Qvalues <- p.adjust(summary_results$PValue, method = "BH")
    write.csv(summary_results, result_file, row.names = FALSE)
    invisible(summary_results)
  }
  
  run_meta(homo_names,  result_file_homo)
  run_meta(hetro_names, result_file_hetro)
}


# ==============================================================================
# 14a. Subclone-Level Karyotypic Heterogeneity
# ==============================================================================
# For each Louvain subclone within a sample, computes:
#   - Mean pairwise Pearson correlation of gene-level CNA profiles among cells
#     in that subclone (AHS = 1 − mean correlation)
#
# This produces a within-subclone heterogeneity measure, separate from the
# sample-level AHS computed in Section 9.
#
# Input:
#   sampleIndex, sample, paper — identifiers
#   louvain_base_path  — base path to sample_clusters_louvain_seurat.csv files
#   outputPrefix       — base path to write clusters_data.csv per sample
#   sampleFilesMap

calculate_subclone_AHS_per_sample <- function(sampleIndex, sample, paper,
                                              louvain_base_path, outputPrefix,
                                              sampleFilesMap) {
  cells_df     <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  sample_cells <- cells_df[
    cells_df$sample == sample &
      tolower(cells_df$cell_type) == "malignant", "cell_name"] %>%
    gsub("\\.", "-", .)
  
  gene_cna_data <- fread(sampleFilesMap[sampleIndex, "Gene_cna_mat"],
                         data.table = FALSE, sep = " ") %>%
    as_tibble() %>%
    rename(gene_name = V1) %>%
    rename_with(~gsub("^X", "", .), starts_with("X")) %>%
    rename_with(~gsub("\\.", "-", .))
  
  gene_cna_data <- gene_cna_data %>%
    select(gene_name, any_of(sample_cells)) %>%
    na.omit()
  
  louvain_path <- file.path(louvain_base_path, paper, sample,
                            "sample_clusters_louvain_seurat.csv")
  if (!file.exists(louvain_path)) {
    message("Louvain file missing for ", paper, " ", sample)
    return(invisible(NULL))
  }
  l_df <- read.csv(louvain_path) %>%
    as_tibble() %>%
    mutate(cell_id = as.character(cell_id))
  
  common_cells <- intersect(setdiff(colnames(gene_cna_data), "gene_name"),
                            l_df$cell_id)
  gene_cna_filt <- gene_cna_data %>% select(gene_name, all_of(common_cells))
  l_df_filt     <- l_df %>%
    filter(cell_id %in% common_cells) %>%
    arrange(match(cell_id, common_cells))
  
  cn_mat <- as.matrix(gene_cna_filt %>% select(-gene_name))
  rownames(cn_mat) <- gene_cna_filt$gene_name
  clustered_cn <- as.data.frame(t(cn_mat)) %>%
    rownames_to_column("cell_id") %>%
    left_join(l_df_filt %>% select(cell_id, MergedCluster), by = "cell_id") %>%
    rename(cluster = MergedCluster) %>%
    as_tibble()
  
  gene_cols    <- setdiff(colnames(clustered_cn), c("cell_id", "cluster"))
  cluster_data <- data.frame(Cluster = character(), MeanCorrelation = numeric(),
                             Mean_SD = numeric(), Total_Var = numeric(),
                             Size = numeric(), stringsAsFactors = FALSE)
  
  for (cl in unique(clustered_cn$cluster)) {
    cl_mat <- as.matrix(clustered_cn %>% filter(cluster == cl) %>% select(all_of(gene_cols)))
    if (nrow(cl_mat) < 2) {
      cluster_data <- rbind(cluster_data,
                            data.frame(Cluster = as.character(cl),
                                       MeanCorrelation = NA, Mean_SD = NA,
                                       Total_Var = NA, Size = nrow(cl_mat)))
      next
    }
    var_cols   <- apply(cl_mat, 2, var, na.rm = TRUE) > 0
    if (sum(var_cols) < 2) {
      mean_cor <- NA
    } else {
      corr_mat <- cor(t(cl_mat[, var_cols]), method = "pearson",
                      use = "pairwise.complete.obs")
      mean_cor <- mean(corr_mat[upper.tri(corr_mat, diag = FALSE)], na.rm = TRUE)
    }
    gene_sd  <- apply(cl_mat, 2, sd,  na.rm = TRUE)
    gene_var <- apply(cl_mat, 2, var, na.rm = TRUE)
    cluster_data <- rbind(cluster_data,
                          data.frame(Cluster = as.character(cl),
                                     MeanCorrelation = mean_cor,
                                     Size = nrow(cl_mat)))
  }
  
  dir.create(file.path(outputPrefix, paper, sample),
             recursive = TRUE, showWarnings = FALSE)
  write.csv(cluster_data, file.path(outputPrefix, paper, sample, "clusters_data.csv"))
  invisible(cluster_data)
}

# ---- aggregate per-sample clusters_data.csv into a single table ----
# Input: outputPrefix, metadata — data.frame with Paper and Sample columns
summarize_subclone_AHS <- function(outputPrefix, metadata) {
  df <- data.frame()
  for (i in seq_len(nrow(metadata))) {
    f <- file.path(outputPrefix, metadata[i, "Paper"], metadata[i, "Sample"],
                   "clusters_data.csv")
    if (!file.exists(f)) next
    cd <- fread(f, data.table = FALSE)
    for (j in seq_len(nrow(cd))) {
      df <- bind_rows(df, data.frame(
        Paper          = metadata[i, "Paper"],
        Sample         = metadata[i, "Sample"],
        Cluster        = cd[j, "Cluster"],
        MeanCorrelation = cd[j, "MeanCorrelation"],
        Size           = cd[j, "Size"]
      ))
    }
  }
  write.csv(df, file.path(outputPrefix, "samples_clusters_data.csv"), row.names = FALSE)
  invisible(df)
}


# ==============================================================================
# 14b. Louvain Subclone One-vs-All DGE
# ==============================================================================
# For each subclone in a sample, runs Seurat FindMarkers comparing that
# subclone against all other cells ("the_rest").  
#
# The enricher_direction() helper is reused for downstream
# pathway enrichment on the resulting gene lists.
#
# Input:
#   files_path    — base path with subclones_info and norm_count_data files
#   cellTypes     — "malignant"
#   sampleFilesMap
#   minPct, test.use — Seurat FindMarkers parameters
#   output_path   — where to write per-subclone results

runFindMarkersTest_oneVsAll <- function(files_path, cellTypes, sampleFilesMap,
                                        minPct, test.use, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "Seurat", "Matrix", "msigdbr",
                  "clusterProfiler", "ggpubr", "icesTAF")
  ) %dopar% {
    sample      <- sampleFilesMap[sampleIndex, "Sample"]
    paper       <- sampleFilesMap[sampleIndex, "Study"]
    cancer_type <- sampleFilesMap[sampleIndex, "Cancer_type"]
    tryCatch(
      run_oneVsAll_one_sample(
        files_path, cellTypes, sampleFilesMap, cancer_type,
        paper, sample, minPct, test.use, output_path
      ),
      error = function(e)
        message("one-vs-all error: ", paper, "_", sample, " — ", e)
    )
  }
}

run_oneVsAll_one_sample <- function(files_path, cellTypes, sampleFilesMap,
                                    cancerType, paper, sample,
                                    minPct, test.use, output_path) {
  min_cells <- 30
  
  norm_counts <- fread(
    file.path(files_path, paper, sample,
              paste0("norm_count_data_", cellTypes, ".csv")),
    data.table = FALSE, header = TRUE
  ) %>%
    remove_rownames() %>% column_to_rownames("V1") %>% drop_na()
  colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))
  norm_counts$sd <- apply(norm_counts, 1, sd)
  norm_counts    <- norm_counts[norm_counts$sd != 0, ]
  norm_counts$sd <- NULL
  
  colData <- fread(
    file.path(files_path, paper, sample, "subclones_info.csv"),
    data.table = FALSE, header = TRUE
  ) %>% drop_na()
  colData[["V1"]] <- gsub("\\.", "-", colData[["V1"]])
  
  subclones <- sort(unique(colData$pred))
  
  for (subclone in subclones) {
    comp <- paste0(subclone, "_vs_the_rest")
    
    colData$new_pred <- ifelse(colData$pred == subclone, subclone, "the_rest")
    groupsColData <- colData %>%
      remove_rownames() %>% column_to_rownames("V1")
    
    counts_by_group <- table(groupsColData$new_pred)
    if (!all(c(subclone, "the_rest") %in% names(counts_by_group)) ||
        any(counts_by_group < min_cells)) next
    
    groupsNormCount <- norm_counts[, rownames(groupsColData)]
    groupsNormCount$sd <- apply(groupsNormCount, 1, sd)
    groupsNormCount    <- groupsNormCount[groupsNormCount$sd != 0, ]
    groupsNormCount$sd <- NULL
    groupsNormCount    <- groupsNormCount[, rownames(groupsColData)]
    
    if (!all(colnames(groupsNormCount) == rownames(groupsColData))) next
    
    mat <- Matrix(as.matrix(groupsNormCount), sparse = TRUE)
    s   <- CreateSeuratObject(counts = mat, meta.data = groupsColData)
    s   <- SetAssayData(s, slot = "data", new.data = mat)
    s   <- FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE)
    s   <- ScaleData(s, features = rownames(s), verbose = FALSE)
    s   <- RunPCA(s, npcs = 30, verbose = FALSE) %>%
      RunUMAP(dims = 1:30, verbose = FALSE) %>%
      FindNeighbors(dims = 1:30, verbose = FALSE) %>%
      FindClusters(resolution = 0.7, verbose = FALSE)
    Idents(s) <- "new_pred"
    
    de <- FindMarkers(s, ident.1 = subclone, ident.2 = "the_rest",
                      min.pct = minPct, test.use = test.use,
                      min.cells.group = min_cells) %>%
      as.data.frame()
    de$group1 <- subclone; de$group2 <- "the_rest"
    
    sig <- filter(de, p_val_adj < 0.25) %>% arrange(p_val_adj)
    
    out_dir <- file.path(output_path, cancerType, paper, sample, comp)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(de,  file.path(out_dir, paste0("findMarkers_test_res_tbl_",
                                          minPct, "_", test.use,
                                          "_groupSize_", min_cells, ".csv")),
           row.names = TRUE)
    if (nrow(sig) > 0) {
      fwrite(sig, file.path(out_dir, paste0("findMarkers_test_sig_res_",
                                            minPct, "_", test.use,
                                            "_groupSize_", min_cells, ".csv")),
             row.names = TRUE)
    }
  }
}

# ---- wrapper to run enricher on all one-vs-all results ----
runEnricher_oneVsAll <- function(files_path, cellTypes, sampleFilesMap,
                                 minPct, test.use, output_path) {
  foreach(
    sampleIndex = seq_len(nrow(sampleFilesMap)),
    .combine = "c",
    .export  = ls(globalenv()),
    .packages = c("tidyverse", "dplyr", "msigdbr", "clusterProfiler")
  ) %dopar% {
    sample <- sampleFilesMap[sampleIndex, "Sample"]
    paper  <- sampleFilesMap[sampleIndex, "Study"]
    colData <- tryCatch(
      fread(file.path(files_path, paper, sample, "subclones_info.csv"),
            data.table = FALSE) %>% drop_na(),
      error = function(e) NULL
    )
    if (is.null(colData)) return(NULL)
    colData[["V1"]] <- gsub("\\.", "-", colData[["V1"]])
    subclones <- sort(unique(colData$pred))
    min_cells <- 30
    
    for (subclone in subclones) {
      comp <- paste0(subclone, "_vs_the_rest")
      path <- file.path(output_path, sampleFilesMap[sampleIndex, "Cancer_type"],
                        paper, sample, comp)
      res_file <- file.path(path, paste0("findMarkers_test_res_tbl_",
                                         minPct, "_", test.use,
                                         "_groupSize_", min_cells, ".csv"))
      if (!file.exists(res_file)) next
      tryCatch({
        enricher_direction(sample, paper, minPct, test.use, min_cells,
                           "positive", path, res_file = res_file)
        enricher_direction(sample, paper, minPct, test.use, min_cells,
                           "negative", path, res_file = res_file)
      }, error = function(e)
        message("Enricher one-vs-all error: ", paper, "_", sample, " — ", e))
    }
  }
}


# ==============================================================================
# 15. Subclone Enrichment Summary Table
# ==============================================================================
# Reads all one-vs-all enricher results and builds a summary table with
# the number of significantly enriched gene sets per subclone per database.
#
# Input:
#   results_base_path — root of the one-vs-all output tree
#   minPct, test.use  — parameters used during DGE
#   excluded_studies  — studies to skip (e.g. low-quality samples)
#
# Output: data.frame with columns Cancer_type, Study, Sample, comp,
#   BIOCARTA, KEGG, REACTOME, GOBP, GOCC, HALLMARK

build_subclone_enrichment_summary <- function(results_base_path,
                                              minPct = 0.25,
                                              test.use = "wilcox",
                                              direction = "positive",
                                              excluded_studies = c()) {
  min_cells <- 30
  file_pattern <- paste0("*_vs_the_rest/findMarkers_test_res_tbl_",
                         minPct, "_", test.use, "_groupSize_", min_cells, ".csv")
  files <- list.files(results_base_path, pattern = "findMarkers_test_res_tbl",
                      full.names = TRUE, recursive = TRUE)
  files <- files[grepl("_vs_the_rest", files)]
  
  files_df <- t(as.data.frame(
    lapply(files, function(x) strsplit(
      gsub(paste0(results_base_path, "/"), "", x), "/")[[1]][1:4])
  ))
  colnames(files_df) <- c("Cancer_type", "Study", "Sample", "comp")
  files_df <- as.data.frame(files_df) %>%
    distinct() %>%
    filter(!Study %in% excluded_studies)
  
  enrich_list <- lapply(seq_len(nrow(files_df)), function(idx) {
    row  <- files_df[idx, ]
    dir  <- file.path(results_base_path, row$Cancer_type, row$Study,
                      row$Sample, row$comp)
    enrich_dir <- file.path(dir, paste0("findMarkers_enricher_",
                                        minPct, "_", direction))
    counts <- sapply(pathways_titles_list, function(db) {
      f <- file.path(enrich_dir, paste0(db, "_enricher_res.csv"))
      if (!file.exists(f)) return(0L)
      pw_df <- read.csv(f) %>%
        mutate(qvalues = ifelse(is.na(qvalue), p.adjust, qvalue))
      sum(pw_df$qvalues < 0.25, na.rm = TRUE)
    })
    c(as.list(row), as.list(counts))
  })
  
  bind_rows(enrich_list)
}
