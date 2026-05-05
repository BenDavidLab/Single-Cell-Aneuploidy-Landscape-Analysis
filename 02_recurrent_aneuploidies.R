################################################################################
# 02_recurrent_aneuploidies.R
#
# Transcriptional consequences of specific recurrent chromosome-arm gains/losses.
#
# This script covers:
#   1.  Identify recurrent CNA events from TCGA arm-level gain/loss tables
#       (threshold: >20% of tumors of a given type)
#   2.  For each (chromosome-arm, cancer-type) pair with enough single-cell
#       samples, prepare expression and CNA data files
#   3.  Run ssGSEA per sample using ssGsea_all_set() from 01_main_analysis.R
#   4.  Compare DEL/AMP vs NEUTRAL cells with Welch's t-test per pathway
#   5.  Run random-effects meta-analysis across samples per (arm, cancer-type)
#   6.  Per-sample IPTW-controlled DGE (weighted t-test) and GSEA
#   7.  Aggregate results per cancer type across all arms
#   8.  SC vs TCGA Correlation per (arm, cancer-type)

#
# NOTE: The ssGsea_all_set() and enricher_direction() functions from
# 01_main_analysis.R are reused here and must be sourced first:
#   source("01_main_analysis.R")
#
# The same logic applies to both chromosome-arm LOSSES and GAINS; only the
# comparison group changes (DEL vs NEUTRAL for losses, AMP vs NEUTRAL for gains).
# The example code below shows losses; gains follow the identical pattern.
################################################################################


# ==============================================================================
# 0.  Source shared utilities
# ==============================================================================
# source("01_main_analysis.R")   # loads all shared functions & theme

suppressMessages({
  library(tidyverse); library(dplyr); library(rstatix)
  library(msigdbr); library(clusterProfiler); library(GSVA)
  library(Seurat); library(Matrix); library(doParallel); library(foreach)
  library(metafor); library(effectsize); library(enrichplot)
  library(weights); library(sjstats)
})


# ==============================================================================
# 1.  Build (chromosome-arm, cancer-type) Pairs from TCGA Frequency Tables
# ==============================================================================
# Input:
#   tcga_freq_file — CSV with rows = cancer types, columns = chromosome arms;
#                    values = fraction of tumors with an alteration
#   frequency_threshold — minimum fraction to call an arm "recurrently altered"
#                         (default 0.2)
#
# Output: named list chr_ratio_to_type
#   keys   — chromosome arm strings (e.g. "18q")
#   values — character vector of cancer-type labels

build_chr_cancer_pairs <- function(tcga_freq_file, frequency_threshold = 0.2) {
  df <- read.csv(tcga_freq_file)
  colnames(df) <- gsub("^X", "", colnames(df))
  df <- df[df$Cancer_name != "", ]

  chr_ratio_to_type <- list()
  for (i in seq(3, ncol(df))) {
    chr    <- colnames(df)[i]
    subset <- df[df[[chr]] > frequency_threshold, "Cancer_name"]
    if (length(subset) > 0) chr_ratio_to_type[[chr]] <- subset
  }
  chr_ratio_to_type
}


# ==============================================================================
# 2.  Prepare Per-Sample Files for a Given (arm, cancer-type) Pair
# ==============================================================================
# For each cell in a sample, reads the per-arm CNA status (AMP / DEL / NEUTRAL)
# from the chromosome-arm CNA matrix, then saves:
#   - subclones_info.csv   (cell_id, pred = CNA status, AS, sample, paper)
#   - norm_count_data_malignant.csv
#
# This follows the same structure as createFilesPerSampleForTest() in
# 01_main_analysis.R, but groups cells by CNA status of a *specific arm*
# rather than by global AS quantile.
#
# Input:
#   sampleIndex, sample, paper — identifiers
#   cellTypes    — "malignant"
#   chr          — chromosome arm string, e.g. "18q"
#   cnv_type     — "Loss" or "Gain"
#   cancer_type  — cancer type label matching sampleFilesMap$Cancer_type
#   outputPrefix — base output path

createSubclonesByChrPerSampleFilteredForNorm <- function(
    sampleIndex, sample, paper, cellTypes, chr, cnv_type, cancer_type,
    outputPrefix, sampleFilesMap
) {
  cells_df <- read.csv(sampleFilesMap[sampleIndex, "Cells_path"])
  sample_cells <- if (cellTypes == "all") {
    cells_df[cells_df$sample == sample, "cell_name"]
  } else {
    cells_df[cells_df$sample == sample &
               tolower(cells_df$cell_type) == tolower(cellTypes), "cell_name"]
  }
  sample_cells <- as.character(sample_cells)

  cna_df <- read.csv(sampleFilesMap[sampleIndex, "Cna_mat"], sep = "\t") %>%
    remove_rownames() %>%
    column_to_rownames("sample")
  cna_df <- cna_df[rownames(cna_df) %in% sample_cells, ]
  colnames(cna_df) <- gsub("^X", "", colnames(cna_df))
  cna_df$AS <- apply(
    cna_df[, !colnames(cna_df) %in% c("Xp", "Xq", "Yp", "Yq")],
    MARGIN = 1, FUN = function(x) length(which(x == "DEL" | x == "AMP"))
  )
  cna_df <- na.omit(cna_df)

  # Extract the CNA call for the specific arm of interest
  chr_num     <- gsub("chr", "", chr)
  sample_info <- data.frame(pred = cna_df[[chr_num]], row.names = rownames(cna_df))
  sample_info$sample  <- sample
  sample_info$paper   <- paper
  sample_info$cell_id <- rownames(sample_info)
  sample_info$AS      <- cna_df$AS
  rownames(sample_info) <- gsub("[._:]", "-", rownames(sample_info))

  # Enforce minimum group sizes
  summary_grp <- sample_info %>%
    group_by(pred) %>%
    get_summary_stats(AS, type = "mean_sd") %>%
    as.data.frame()
  min_cells <- 30

  if (cnv_type == "Loss") {
    n_del <- if ("DEL" %in% summary_grp$pred)
      summary_grp[summary_grp$pred == "DEL", "n"] else 0
    n_wt  <- if ("NEUTRAL" %in% summary_grp$pred)
      summary_grp[summary_grp$pred == "NEUTRAL", "n"] else 0
    if (n_del <= min_cells || n_wt <= min_cells) {
      message("Not enough DEL or NEUTRAL cells for ", paper, " ", sample)
      return(invisible(NULL))
    }
  } else {
    n_amp <- if ("AMP" %in% summary_grp$pred)
      summary_grp[summary_grp$pred == "AMP", "n"] else 0
    n_wt  <- if ("NEUTRAL" %in% summary_grp$pred)
      summary_grp[summary_grp$pred == "NEUTRAL", "n"] else 0
    if (n_amp <= min_cells || n_wt <= min_cells) {
      message("Not enough AMP or NEUTRAL cells for ", paper, " ", sample)
      return(invisible(NULL))
    }
  }

  norm_counts <- fread(sampleFilesMap[sampleIndex, "Norm_Expression_path"],
                       data.table = FALSE) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("[._:]", "-", colnames(norm_counts))

  cols      <- gsub("^X", "", colnames(norm_counts))
  expr_filt <- norm_counts[,
    (colnames(norm_counts) %in% rownames(sample_info)) |
      (cols %in% rownames(sample_info)), drop = FALSE] %>%
    na.omit()
  filt_cols <- gsub("^X", "", colnames(expr_filt))

  sample_info <- sample_info[
    (rownames(sample_info) %in% colnames(expr_filt)) |
      (rownames(sample_info) %in% filt_cols), ] %>%
    na.omit()
  sample_info <- sample_info[
    ifelse(filt_cols != colnames(expr_filt), filt_cols, colnames(expr_filt)), ]
  if (all(filt_cols == rownames(sample_info))) colnames(expr_filt) <- filt_cols

  path <- file.path(outputPrefix, cnv_type, chr, cancer_type, paper, sample)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  write.csv(sample_info, file.path(path, "subclones_info.csv"))
  write.csv(expr_filt,   file.path(path, paste0("norm_count_data_", cellTypes, ".csv")))
  fwrite(summary_grp,    file.path(path, "AS_summary.csv"), row.names = TRUE)
}


# ==============================================================================
# 3.  ssGSEA + t-test for Recurrent Aneuploidies
# ==============================================================================
# Re-uses ssGsea_all_set() and createMergedSampleInfo_withSsgsea() from
# 01_main_analysis.R to compute per-cell pathway scores, then runs Welch's
# t-test comparing DEL/AMP cells vs NEUTRAL cells per pathway.
#
# Input:
#   files_path   — base path containing <Study>/<Sample>/ subdirs
#   paper, sample — identifiers
#   cnv_type     — "Loss" or "Gain"
#   output_path  — where to write ssgsea_tTestRes_<DB>.csv
#   pathway_index — which pathway DBs to run (default: 6 = HALLMARK; 3 = REACTOME)

ssGSEA_tTest_perSample_arm <- function(files_path, paper, sample, cnv_type,
                                        output_path, pathway_index = 6L,
                                        sampleFilesMap) {
  for (i in pathway_index) {
    merged_file <- file.path(files_path, paper, sample,
                              paste0("subclones_info_with__ssgsea_",
                                     pathways_titles_list[i], ".csv"))
    if (!file.exists(merged_file)) next

    colData <- read.csv(merged_file) %>%
      remove_rownames() %>%
      column_to_rownames("Row.names") %>%
      .[, !colnames(.) %in% c("X.1", "X")] %>%
      mutate(sample = as.character(sample))

    pw_cols <- grep(paste0(pathways_titles_list[i], "_"),
                    colnames(colData), value = TRUE)
    results <- data.frame(
      pathway = character(), p_value = numeric(),
      effect_size = numeric(), stderr = numeric(),
      mean_difference = numeric(), sampleID = character(),
      stringsAsFactors = FALSE
    )

    for (pw in pw_cols) {
      cn_vals      <- colData[colData$pred == ifelse(cnv_type == "Loss", "DEL", "AMP"), pw]
      neutral_vals <- colData[colData$pred == "NEUTRAL", pw]
      tt           <- t.test(cn_vals, neutral_vals, var.equal = FALSE)
      cd           <- cohens_d(cn_vals, neutral_vals, pooled_sd = FALSE)
      results      <- rbind(results, data.frame(
        pathway         = pw,
        p_value         = tt$p.value,
        effect_size     = cd$Cohens_d,
        stderr          = tt$stderr,
        mean_difference = diff(tt$estimate),
        sampleID        = paste0(paper, "_", sample),
        stringsAsFactors = FALSE
      ))
    }
    results$qvalues <- p.adjust(results$p_value, method = "BH")
    write.csv(results, file.path(output_path,
                                  paste0("ssgsea_tTestRes_",
                                         pathways_titles_list[i], ".csv")))
  }
}


# ==============================================================================
# 4.  Meta-Analysis per (arm, cancer-type) Pair
# ==============================================================================
# Reads all per-sample ssgsea_tTestRes_<DB>.csv files in a given directory
# and runs a random-effects meta-analysis (REML) per pathway.
#
# This is the same as ssgsea_metaAnalysis() in 01_main_analysis.R but
# operates on a single (arm, cancer-type) subdirectory.

ssgsea_metaAnalysis_arm <- function(output_path, pathway_index = 6L) {
  for (i in pathway_index) {
    files <- list.files(output_path,
                        pattern = paste0("ssgsea_tTestRes_",
                                         pathways_titles_list[i], ".csv"),
                        full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) next

    combined <- bind_rows(lapply(files, read.csv), .id = "id") %>%
      mutate(sampleID = gsub("[. -]", "_", sampleID))

    meta_data <- combined %>%
      group_by(pathway) %>%
      summarise(EffectSize = list(effect_size),
                SE         = list(stderr),
                .groups    = "drop")

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
      PValue           = sapply(meta_results, \(x) x$pval)
    )
    summary_results$Qvalues <- p.adjust(summary_results$PValue, method = "BH")
    write.csv(summary_results,
              file.path(output_path,
                        paste0("ssgsea_meta_analysis_tTest_",
                               pathways_titles_list[i], "results.csv")),
              row.names = FALSE)
  }
}


# ==============================================================================
# 5.  IPTW-Controlled DGE (Single Sample)
# ==============================================================================
# To control for the confounding between global aneuploidy (AS) and specific
# arm-level CNA status, we use Inverse Probability of Treatment Weighting (IPTW):
#   - A logistic model predicts arm CNA status from AS
#   - Stabilized weights are derived from propensity scores
#   - Weighted t-tests (wtd.t.test) compare gene expression between
#     DEL/AMP and NEUTRAL cells
#
# Input:
#   arm         — chromosome arm, e.g. "18q"
#   cnv_type    — "Loss" or "Gain"
#   paper, sample — identifiers
#   outputPrefix — base output path containing subclones_info.csv
#                  and norm_count_data_malignant.csv
#
# Similar logic was used at the sample level using the TCGA dataset.

DGE_iptw_perSample <- function(arm, cnv_type, paper, sample, outputPrefix) {
  cellTypes <- "malignant"

  norm_counts <- fread(
    file.path(outputPrefix, paper, sample,
              paste0("norm_count_data_", cellTypes, ".csv")),
    data.table = FALSE, header = TRUE
  ) %>%
    remove_rownames() %>%
    column_to_rownames("V1") %>%
    drop_na()
  colnames(norm_counts) <- gsub("\\.", "-", colnames(norm_counts))
  norm_counts$sd <- apply(norm_counts, 1, sd)
  norm_counts    <- norm_counts[norm_counts$sd != 0, ]
  norm_counts$sd <- NULL

  colData <- fread(
    file.path(outputPrefix, paper, sample, "subclones_info.csv"),
    data.table = FALSE, header = TRUE
  ) %>%
    drop_na()
  colData[["V1"]] <- gsub("\\.", "-", colData[["V1"]])

  clone_of_interest <- ifelse(cnv_type == "Loss", "DEL", "AMP")
  groupsColData <- colData[colData$pred %in% c(clone_of_interest, "NEUTRAL"), ] %>%
    remove_rownames() %>%
    column_to_rownames("V1")

  groupsNormCount <- norm_counts[, rownames(groupsColData)]
  groupsNormCount$sd <- apply(groupsNormCount, 1, sd)
  groupsNormCount    <- groupsNormCount[groupsNormCount$sd != 0, ]
  groupsNormCount$sd <- NULL
  groupsNormCount    <- groupsNormCount[, rownames(groupsColData)]

  if (!all(colnames(groupsNormCount) == rownames(groupsColData))) {
    message("Column/row mismatch for ", paper, " ", sample, " — skipping IPTW")
    return(invisible(NULL))
  }

  groupsColData$cell_id <- rownames(groupsColData)

  # Propensity model: predict arm CNA status from global AS
  anp_data <- groupsColData %>%
    mutate(arm.cn.change = (pred == clone_of_interest)) %>%
    select(cell_id, AS, arm.cn.change)

  message("Fitting propensity model (AS → arm CNA status)...")
  prop_model  <- glm(arm.cn.change ~ AS, data = anp_data, family = "binomial")
  prop_scores <- predict(prop_model, type = "response")
  iptw        <- ifelse(anp_data$arm.cn.change,
                        1 / prop_scores, 1 / (1 - prop_scores))
  anp_data$w  <- ifelse(anp_data$arm.cn.change,
                         iptw * mean(anp_data$arm.cn.change),
                         iptw * mean(!anp_data$arm.cn.change))

  expr_df      <- as.data.frame(t(groupsNormCount))
  expr_df$name <- rownames(expr_df)
  df           <- merge(anp_data, expr_df, by.x = "cell_id", by.y = "name")

  message("Running weighted t-tests for ", ncol(df) - 4, " genes...")
  effect_df <- foreach(
    gene = colnames(df)[-c(1, 2, 3, 4)],
    .init    = data.frame(),
    .combine = rbind
  ) %do% {
    tryCatch({
      x   <- df[df$arm.cn.change, gene]
      y   <- df[!df$arm.cn.change, gene]
      w_x <- df$w[df$arm.cn.change]
      w_y <- df$w[!df$arm.cn.change]
      p   <- unname(weights::wtd.t.test(x, y, w_x, w_y)$coefficients["p.value"])
      w_sd_x <- sjstats::weighted_sd(x, w_x)
      w_sd_y <- sjstats::weighted_sd(y, w_y)
      n_x <- length(x); n_y <- length(y)
      pooled_sd <- sqrt(((n_x - 1) * w_sd_x^2 + (n_y - 1) * w_sd_y^2) /
                          (n_x + n_y - 2))
      cd <- (weights::wtd.mean(x, w_x) - weights::wtd.mean(y, w_y)) / pooled_sd
      data.frame(name = gene, pval = p, effect.size = cd)
    }, error = function(e) {
      message("Error on gene ", gene, ": ", e)
      NULL
    })
  }

  effect_df$qval <- p.adjust(effect_df$pval, method = "BH")
  if (any(effect_df$qval < 0.25, na.rm = TRUE)) {
    write.table(effect_df,
                file.path(outputPrefix, paper, sample, "iptw_AS_Results_Genes.csv"),
                sep = ",", row.names = FALSE)
  }
  invisible(effect_df)
}


# ==============================================================================
# 6.  GSEA on IPTW-Ranked Gene List (Single Sample)
# ==============================================================================
# Creates an effect-size ranked gene list from the IPTW output and runs
# clusterProfiler::GSEA().  Saves a GSEA plot for each pathway in `paths`.
#
# Input:
#   outputPrefix — path containing iptw_AS_Results_Genes.csv
#   cnv_type, chr — for labelling
#   paper, sample — identifiers
#   paths        — character vector of gene set IDs to plot
#   pathway_index — 3 (REACTOME) or 6 (HALLMARK)
#
# Similar logic was used at the sample level using the TCGA dataset.

createGSEAPlot_iptw <- function(outputPrefix, cnv_type, chr, paper, sample,
                                  paths, pathway_index = 3L) {
  iptw_file <- file.path(outputPrefix, paper, sample, "iptw_AS_Results_Genes.csv")
  if (!file.exists(iptw_file)) return(invisible(NULL))

  res <- read.csv(iptw_file) %>%
    filter(!is.na(effect.size)) %>%
    arrange(desc(effect.size))

  gene_list       <- res$effect.size
  names(gene_list) <- res$name

  i     <- pathway_index
  m_t2g <- msigdbr(species = "Homo sapiens",
                   category    = pathways_categories_list[i],
                   subcategory = pathways_anontations_list[i]) %>%
    select(gs_name, gene_symbol)

  set.seed(1234)
  gsea_res <- tryCatch(
    clusterProfiler::GSEA(gene_list, TERM2GENE = m_t2g, pvalueCutoff = 1,
                          seed = TRUE),
    error = function(e) { message("GSEA failed: ", e); NULL }
  )
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) return(invisible(NULL))

  write.csv(gsea_res@result,
            file.path(outputPrefix, paper, sample,
                      paste0("iptw_AS_control_gsea_res_",
                             pathways_titles_list[i], ".csv")))

  sig <- gsea_res@result[gsea_res@result$qvalue < 0.25, ]
  for (gs_id in paths) {
    if (!gs_id %in% sig$ID) next
    clean_name <- gsub(paste0(pathways_titles_list[i], "_"), "", gs_id)
    p <- enrichplot::gseaplot2(
      gsea_res, geneSetID = gs_id,
      title = paste0(paper, " ", sample, "\n", clean_name)
    )
    p[[1]] <- p[[1]] +
      labs(y = "Running\nEnrichment Score") +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 17))
    p[[3]] <- p[[3]] +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 15))

    out_file <- file.path(outputPrefix, paper, sample,
                           paste0(gs_id, "_iptw_GSEA_plot.png"))
    # Print p to save — caller should use ggsave() or save_publication_png()
    print(p)
  }
}

# ==============================================================================
# 7.   Aggregate results per cancer type across all arms
# ==============================================================================
# For each cancer type, reads the per-arm meta-analysis results and merges
# them into a single pathways × arms matrix (one file per cancer type).
# Only pairs with ≥ 5 samples are included.
#
# Input:
#   outputPrefix      — base results path
#   chr_ratio_to_type — named list from build_chr_cancer_pairs()
#   cnv_type          — "Loss" or "Gain"
#   count_df          — output of build_pair_count_summary()
#   pathway_index     — 3 (REACTOME) or 6 (HALLMARK)
#   min_samples       — minimum sample count per pair (default 5)

aggregate_results_per_cancer_type <- function(outputPrefix, chr_ratio_to_type,
                                              cnv_type, count_df,
                                              pathway_index = 6L,
                                              min_samples = 5L) {
  cancer_types <- unique(unlist(chr_ratio_to_type))
  for (type in cancer_types) {
    df <- data.frame()
    for (chr in names(chr_ratio_to_type)) {
      if (!type %in% chr_ratio_to_type[[chr]]) next
      pair_count <- if (nrow(count_df[count_df$Arm == chr &
                                      count_df$Cancer_type == type, ]) > 0)
        count_df[count_df$Arm == chr & count_df$Cancer_type == type,
                 "Sample_Count"] else 0L
      if (pair_count < min_samples) next
      
      result_file <- file.path(outputPrefix, cnv_type, chr, type,
                               paste0("ssgsea_meta_analysis_tTest_",
                                      pathways_titles_list[pathway_index],
                                      "results.csv"))
      if (!file.exists(result_file)) next
      
      chr_df <- read.csv(result_file) %>%
        column_to_rownames("Pathway") %>%
        filter(Qvalues < 0.25) %>%
        select(PooledEffectSize)
      colnames(chr_df) <- chr
      
      if (ncol(df) == 0) {
        df <- chr_df
      } else {
        df <- merge(df, chr_df, by = 0, all = TRUE) %>%
          remove_rownames() %>% column_to_rownames("Row.names")
      }
    }
    if (nrow(df) > 1) {
      out_dir <- file.path(outputPrefix, cnv_type, "byType", "ssGSEA_metaAnalysis")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(df,
                file.path(out_dir,
                          paste0(type, "_ssgsea_meta_analysis_",
                                 pathways_titles_list[pathway_index],
                                 "_allExistArms_results.csv")),
                row.names = TRUE)
    }
  }
}


# ==============================================================================
# 8.  SC vs TCGA Correlation per (arm, cancer-type)
# ==============================================================================
# For each (arm, cancer-type) pair, correlates the SC meta-analysis pooled
# effect sizes against the TCGA IPTW-controlled GSEA NES values using
# Spearman rank correlation.
#
# Both rankings are restricted to pathways significant at q < 0.25 in TCGA.
# A pair is included only if ≥ 8 overlapping pathways exist after filtering.
#
# Inputs:
#   cancer_types       — character vector of cancer type labels
#   TCGA_types         — named list: cancer_type → TCGA cohort abbreviation
#   tcga_gsea_base     — base path to TCGA IPTW GSEA CSVs
#                        (file: <tcga_gsea_base>/<arm>/<tcga_type>/<db>_IPTW_controlled_GSEA.csv)
#   sc_results_base    — base path to per-cancer-type allExistArms CSVs
#                        (output of aggregate_results_per_cancer_type)
#   cnv_type           — "Loss" or "Gain"
#   pathway_index      — 3 (REACTOME) or 6 (HALLMARK)
#
# Output: data.frame with columns Cancer_type, TCGA_cancer_type, Arm,
#         Correlation, P_val

compute_SC_TCGA_correlation <- function(cancer_types, TCGA_types,
                                        tcga_gsea_base, sc_results_base,
                                        cnv_type, pathway_index = 6L) {
  db_subdir <- if (pathway_index == 6L) "H_Null" else
    if (pathway_index == 3L) "C2_CP_REACTOME" else ""
  
  df_count <- data.frame()
  
  for (type in cancer_types) {
    sc_file <- file.path(sc_results_base, cnv_type, "byType", "ssGSEA_metaAnalysis",
                         paste0(type, "_ssgsea_meta_analysis_",
                                pathways_titles_list[pathway_index],
                                "_allExistArms_results.csv"))
    if (!file.exists(sc_file)) next
    sc_df          <- read.csv(sc_file)
    colnames(sc_df)[1] <- "Pathway"
    colnames(sc_df)    <- gsub("^X", "", colnames(sc_df))
    
    for (tcga_type in TCGA_types[[type]]) {
      for (arm in setdiff(colnames(sc_df), "Pathway")) {
        tcga_file <- file.path(tcga_gsea_base, arm, tcga_type,
                               paste0(db_subdir, "_IPTW_controlled_GSEA.csv"))
        if (!file.exists(tcga_file)) next
        
        tcga_gsea <- read.csv(tcga_file) %>%
          select(ID, NES, qvalue) %>%
          filter(qvalue < 0.25) %>%
          arrange(NES) %>%
          mutate(rank = row_number())
        
        sc_arm <- sc_df[, c("Pathway", arm)] %>%
          drop_na() %>%
          arrange(.data[[arm]]) %>%
          mutate(rank = row_number())
        
        merged <- merge(
          tcga_gsea[, c("NES", "ID", "rank")],
          sc_arm[, c(arm, "Pathway", "rank")],
          by.x = "ID", by.y = "Pathway", all = TRUE
        )
        colnames(merged) <- c("Path", "TCGA_NES", "TCGA_rank",
                              "Tirosh_pooledES", "Tirosh_rank")
        merged$TCGA_rank   <- as.numeric(merged$TCGA_rank)
        merged$Tirosh_rank <- as.numeric(merged$Tirosh_rank)
        
        merged_complete <- drop_na(merged)
        if (nrow(merged_complete) > 7) {
          ct  <- cor.test(merged$TCGA_NES, merged$Tirosh_pooledES,
                          method = "spearman", exact = FALSE)
          cor_val <- round(ct$estimate, 4)
          p_val   <- round(ct$p.value, 4)
        } else {
          cor_val <- 0; p_val <- 1
        }
        df_count <- bind_rows(df_count, data.frame(
          Cancer_type     = type,
          TCGA_cancer_type = tcga_type,
          Arm             = arm,
          Correlation     = cor_val,
          P_val           = p_val
        ))
      }
    }
  }
  df_count
}

# ---- standard cancer type → TCGA abbreviation mapping ----
# Adapt if your cohort uses different labels.
TCGA_types <- list(
  "Brain"       = "GBM",  "Colorectal" = "COAD", "HN"     = "HNSC",
  "Kidney"      = "KIRC", "Lung"       = "LUAD", "Sarcoma" = "SARC",
  "Skin"        = "SKCM", "Breast"     = "BRCA", "Liver"   = "LIHC",
  "Pancreas"    = "PAAD", "Hematologic" = "LAML", "Prostate" = "PRAD",
  "Ovarian"     = "OV"
)
