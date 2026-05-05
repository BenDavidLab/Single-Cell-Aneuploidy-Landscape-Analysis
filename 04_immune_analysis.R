library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(cowplot)
library(ragg)

options(bitmapType = "cairo")
set.seed(42)

# =============================================================================
# THEME & DIMENSIONS
# =============================================================================

graph_theme <- theme_minimal() + theme(
  axis.text        = element_text(size = 17),
  axis.title       = element_text(size = 22),
  plot.title       = element_text(hjust = 0.5, face = "bold", size = 24),
  axis.text.x      = element_text(angle = 0, size = 17),
  axis.text.y      = element_text(size = 17),
  legend.key.size  = unit(1, 'cm'),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title     = element_text(size = 17),
  legend.text      = element_text(size = 15),
  panel.border     = element_blank(),
  axis.line        = element_line(color = "black"),
  strip.text       = element_text(size = 18)
)

# =============================================================================
# TOGGLES & CONSTANTS
# =============================================================================

SOLID_TUMORS_ONLY              <- TRUE
SKIP_FILTERS                   <- FALSE
MIN_TME_CELLS_PER_SAMPLE       <- 0
MIN_MALIGNANT_CELLS_PER_SAMPLE <- 0
MAX_OTHER_PCT                  <- 40
EXCLUDE_STUDIES                <- c()
CELLTYPE_COL                   <- "cell_type_full"
SAMPLE_COL                     <- "match_id"
MALIGNANT_VAL                  <- "Malignant"

# =============================================================================
# CANCER TYPE COLOR PALETTE
# =============================================================================

cancer_cols <- c(
  "Brain"          = "dodgerblue4",
  "Breast"         = "deeppink3",
  "Colorectal"     = "forestgreen",
  "HN"             = "orange",
  "Head and Neck"  = "orange",
  "Hematologic"    = "red3",
  "Kidney"         = "yellow",
  "Liver"          = "bisque3",
  "Lung"           = "purple",
  "Neuroendocrine" = "lightblue",
  "Ovarian"        = "lightpink",
  "Pancreas"       = "lightgreen",
  "Pituitary"      = "darkcyan",
  "Prostate"       = "darkred",
  "Sarcoma"        = "darkorange",
  "Skin"           = "chocolate4"
)

# =============================================================================
# PATHS
# =============================================================================

merged_rds <- ".../merged_seurat.rds"

samples_file <- ".../samples_metadata.csv"

meta_bulk_file <- ".../pseuso_bulk_metadata.csv"

# =============================================================================
# LOAD DATA
# =============================================================================

message("\n=== Loading data ===")
if (!file.exists(merged_rds)) stop("File not found: ", merged_rds)

seu      <- readRDS(merged_rds)
all_meta <- seu@meta.data
rm(seu)
gc()

all_meta <- all_meta %>%
  mutate(
    !!sym(CELLTYPE_COL) := case_when(
      .data[[CELLTYPE_COL]] %in% c("Granulocyte", "Granulocytes", "Neutrophil") ~ "Granulocyte",
      .data[[CELLTYPE_COL]] %in% c("Lymphoid", "Myeloid", "TIL")               ~ "Immune_Other",
      .data[[CELLTYPE_COL]] %in% c("Other", "Unknown", "Unassigned")            ~ "Unassigned",
      .data[[CELLTYPE_COL]] %in% c("Dendritic", "pDC")                          ~ "Dendritic",
      TRUE ~ .data[[CELLTYPE_COL]]
    )
  )

samples_metadata <- read.csv(
  samples_file, sep = ",", stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM", check.names = FALSE
)

if (colnames(samples_metadata)[1] %in% c("", NA_character_)) {
  if (all(grepl("^\\d+$", head(samples_metadata[[1]], 10)))) {
    samples_metadata <- samples_metadata[, -1]
  } else {
    colnames(samples_metadata)[1] <- "row_id"
  }
}

samples_metadata <- samples_metadata %>%
  mutate(
    Paper    = trimws(Paper),
    Sample   = trimws(Sample),
    match_id = gsub("\\.", "_", paste(Paper, Sample, sep = "_"))
  )

meta_bulk         <- NULL
mean_as_available <- FALSE
if (file.exists(meta_bulk_file)) {
  try({
    meta_bulk <- read.csv(meta_bulk_file, row.names = 1) %>%
      mutate(
        Paper         = trimws(Paper),
        Sample        = trimws(Sample),
        underscore_id = gsub("\\.", "_", paste(Paper, Sample, sep = "_"))
      ) %>%
      select(underscore_id, Mean_AS) %>%
      distinct()
    mean_as_available <- TRUE
  }, silent = TRUE)
}

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

if (SOLID_TUMORS_ONLY) {
  base_dir <- ".../TME_Comprehensive_Analysis_Solid_Only/"
} else {
  base_dir <-".../TME_Comprehensive_Analysis/"
}

filter_tag <- if (SKIP_FILTERS) "unfiltered" else "filtered"
output_dir <- paste0(base_dir, "Publication_Figures_", filter_tag, "/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", output_dir)

# =============================================================================
# METADATA & BULK INTEGRATION
# =============================================================================

if (SOLID_TUMORS_ONLY) {
  message("  Filtering out Hematologic cancers...")
  samples_metadata <- samples_metadata %>%
    filter(!grepl("(?i)hematologic", Cancer_type))
}

keep_cols <- intersect(
  c("match_id", "Malignant_Mean_AS", "high_vs_low_comp"),
  colnames(samples_metadata)
)

metadata_to_merge <- samples_metadata %>%
  select(all_of(keep_cols)) %>%
  distinct(match_id, .keep_all = TRUE) %>%
  rename(sample = match_id)

if (mean_as_available && !is.null(meta_bulk)) {
  metadata_to_merge <- metadata_to_merge %>%
    left_join(meta_bulk, by = c("sample" = "underscore_id"))
  mean_as_available <- sum(!is.na(metadata_to_merge$Mean_AS)) >= 10
}

# =============================================================================
# COMPUTE CELL TYPE PERCENTAGES
# =============================================================================

all_cell_types <- unique(all_meta[[CELLTYPE_COL]][all_meta$is_malignant != MALIGNANT_VAL])
all_samples    <- unique(metadata_to_merge$sample)

malignant_counts <- all_meta %>%
  filter(is_malignant == MALIGNANT_VAL) %>%
  group_by(.data[[SAMPLE_COL]]) %>%
  summarise(n_malignant = n(), .groups = "drop") %>%
  rename(sample = all_of(SAMPLE_COL))

nonmal_counts <- all_meta %>%
  filter(is_malignant != MALIGNANT_VAL) %>%
  group_by(.data[[SAMPLE_COL]], .data[[CELLTYPE_COL]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  rename(sample = all_of(SAMPLE_COL), cell_type = all_of(CELLTYPE_COL))

nonmal_totals <- all_meta %>%
  filter(is_malignant != MALIGNANT_VAL) %>%
  group_by(.data[[SAMPLE_COL]]) %>%
  summarise(n_nonmalignant = n(), .groups = "drop") %>%
  rename(sample = all_of(SAMPLE_COL))

sample_summary <- nonmal_totals %>%
  left_join(malignant_counts, by = "sample") %>%
  mutate(
    n_malignant       = ifelse(is.na(n_malignant), 0, n_malignant),
    n_total           = n_nonmalignant + n_malignant,
    percent_malignant = (n_malignant / n_total) * 100
  )

full_grid <- expand.grid(
  sample    = all_samples,
  cell_type = all_cell_types,
  stringsAsFactors = FALSE
)

nonmal_counts_complete <- full_grid %>%
  left_join(nonmal_counts, by = c("sample", "cell_type")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  left_join(sample_summary, by = "sample") %>%
  mutate(
    percent_of_nonmalignant = (n / n_nonmalignant) * 100,
    percent_of_all_cells    = (n / n_total)        * 100
  )

percent_all_cells_wide <- nonmal_counts_complete %>%
  select(sample, cell_type, percent_of_all_cells) %>%
  pivot_wider(
    names_from   = cell_type,
    values_from  = percent_of_all_cells,
    names_prefix = "pct_"
  ) %>%
  left_join(
    sample_summary %>%
      select(sample, n_nonmalignant, n_malignant, n_total, percent_malignant),
    by = "sample"
  )

sample_cancer_type <- all_meta %>%
  select(sample = all_of(SAMPLE_COL), cancer_type) %>%
  distinct()

pct_raw <- percent_all_cells_wide %>%
  left_join(sample_cancer_type, by = "sample") %>%
  left_join(metadata_to_merge,  by = "sample") %>%
  select(
    sample, cancer_type,
    n_nonmalignant, n_malignant, n_total, percent_malignant,
    Malignant_Mean_AS,
    everything()
  )

write.csv(pct_raw,
          paste0(output_dir, "cell_type_percent_of_all_cells.csv"),
          row.names = FALSE)

# =============================================================================
# APPLY FILTERS
# =============================================================================

immune_cell_names <- c(
  "B_cell", "Dendritic", "Granulocyte", "Immune_Other", "Macrophage",
  "Mast", "Monocyte", "NK_cell", "Plasma", "T_cell"
)

immune_cols <- paste0("pct_", immune_cell_names)
immune_cols <- immune_cols[immune_cols %in% colnames(pct_raw)]

data <- pct_raw %>%
  filter(!is.na(Malignant_Mean_AS)) %>%
  rowwise() %>%
  mutate(total_immune_pct = sum(c_across(all_of(immune_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(cancer_type = ifelse(cancer_type == "HN", "Head and Neck", cancer_type))

if (!SKIP_FILTERS) {
  if (length(EXCLUDE_STUDIES) > 0) {
    exclude_pattern <- paste0("^", EXCLUDE_STUDIES, collapse = "|")
    data <- data %>% filter(!grepl(exclude_pattern, sample))
  }
  data <- data %>%
    filter(n_nonmalignant >= MIN_TME_CELLS_PER_SAMPLE) %>%
    filter(n_malignant    >= MIN_MALIGNANT_CELLS_PER_SAMPLE)

  if ("pct_Unassigned" %in% colnames(data)) {
    data <- data %>% filter(pct_Unassigned <= MAX_OTHER_PCT)
  }
}

# =============================================================================
# PLOT 2 — CELL COMPOSITION BY CANCER TYPE
# =============================================================================

message("\n--- Plot 2: full composition ---")

all_pct_cols <- colnames(data)[grep("^pct_", colnames(data))]

all_pct_incl_mal <- c(
  all_pct_cols,
  if ("percent_malignant" %in% colnames(data)) "percent_malignant" else NULL
)

cancer_type_n <- data %>%
  group_by(cancer_type) %>%
  summarise(n_samples = n(), .groups = "drop")

full_long <- data %>%
  select(sample, cancer_type, all_of(all_pct_incl_mal)) %>%
  pivot_longer(cols = all_of(all_pct_incl_mal),
               names_to = "cell_type", values_to = "pct") %>%
  mutate(
    cell_type = gsub("^pct_", "", cell_type),
    cell_type = gsub("^percent_malignant$", "Malignant", cell_type),
    cell_type = case_when(
      cell_type %in% c("Ductal_1", "Ductal_2")                               ~ "Ductal",
      grepl("(?i)non[-_ ]?malignant", cell_type) | cell_type == "Unassigned" ~ "Unassigned",
      TRUE ~ cell_type
    ),
    cell_type = gsub("_", " ", cell_type)
  ) %>%
  group_by(cancer_type, sample, cell_type) %>%
  summarise(pct = sum(pct, na.rm = TRUE), .groups = "drop") %>%
  group_by(cancer_type, cell_type) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

bar_totals <- full_long %>%
  group_by(cancer_type) %>%
  summarise(bar_total = sum(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  left_join(cancer_type_n, by = "cancer_type") %>%
  mutate(label = paste0("n=", n_samples))

# Stacking order: Malignant on top, Unassigned on bottom, rest ranked by abundance
abundance_ranking <- full_long %>%
  group_by(cell_type) %>%
  summarise(overall_mean = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  pull(cell_type)

dynamic_cells       <- setdiff(abundance_ranking, c("Malignant", "Unassigned"))
final_order         <- c("Unassigned", rev(dynamic_cells), "Malignant")
final_order         <- intersect(final_order, unique(full_long$cell_type))
full_long$cell_type <- factor(full_long$cell_type, levels = final_order)

base_palette <- c(
  "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A",
  "#B15928", "#A6CEE3", "#FB9A99", "#FDBF6F", "#CAB2D6",
  "#FFFF33", "#01665E", "#D95F02", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#0000CD", "#DC143C",
  "#006400", "#4B0082", "#00CED1", "#8B4513", "#C71585",
  "#2E8B57", "#FF1493", "#DAA520", "#4169E1", "#8B0000",
  "#00FA9A", "#FF6347", "#9400D3", "#3CB371", "#CD853F",
  "#1E90FF", "#FF4500", "#228B22", "#DA70D6", "#B8860B"
)

color_map_plot2 <- setNames(base_palette[seq_along(dynamic_cells)], dynamic_cells)
color_map_plot2["Malignant"]  <- "#000000"
color_map_plot2["Unassigned"] <- "#CCCCCC"

cancer_order <- full_long %>%
  group_by(cancer_type) %>%
  summarise(total = sum(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(cancer_type)

bar_totals <- bar_totals %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_order))

p_full <- ggplot(full_long, aes(x = reorder(cancer_type, -mean_pct), y = mean_pct, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    data        = bar_totals,
    aes(x = cancer_type, y = bar_total, label = label),
    inherit.aes = FALSE,
    vjust       = -0.4,
    size        = 3.25,
    fontface    = "plain"
  ) +
  scale_fill_manual(values = color_map_plot2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  graph_theme +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    legend.key.size = unit(0.3, 'cm'),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10)
  ) +
  guides(fill = guide_legend(ncol = 2, reverse = TRUE)) +
  labs(x = "Cancer Type", y = "Mean % of All Cells", fill = "Cell Type")

# =============================================================================
# PLOT 3 — PAN-CANCER SCATTER (pseudo-bulk Mean_AS)
# =============================================================================

message("\n--- Plot 3: pan-cancer scatter ---")

cor_total <- cor.test(data$Mean_AS, data$total_immune_pct, method = "spearman")
message("  Pan-cancer rho = ", round(cor_total$estimate, 4),
        "  p = ", signif(cor_total$p.value, 4))

present_cancer_types_p3 <- unique(data$cancer_type)
cancer_cols_p3 <- cancer_cols[names(cancer_cols) %in% present_cancer_types_p3]
missing_p3     <- setdiff(present_cancer_types_p3, names(cancer_cols_p3))
if (length(missing_p3) > 0) {
  cancer_cols_p3 <- c(cancer_cols_p3, setNames(rep("grey50", length(missing_p3)), missing_p3))
}

p_pan <- ggplot(data, aes(x = Mean_AS, y = total_immune_pct)) +
  geom_point(aes(color = cancer_type), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", col = alpha('blue', 0.5), level = 0.95) +
  stat_cor(
    method      = "spearman",
    size        = 5,
    label.y.npc = 0.95,
    label.x.npc = "right",
    hjust       = 1,
    output.type = "text",
    aes(label = paste0("rho = ", round(after_stat(r), 2), "\np ",
                       ifelse(after_stat(p) < 2.2e-16, "< 2.2e-16",
                              paste0("= ", signif(after_stat(p), 3)))))
  ) +
  scale_color_manual(values = cancer_cols_p3) +
  graph_theme +
  labs(
    x     = "Pseudo-Bulk Aneuploidy Score",
    y     = "Immune Cell Infiltration ",
    color = "Cancer Type"
  ) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, 'cm'),
    legend.text     = element_text(size = 9),
    legend.title    = element_text(size = 10)
  )

# =============================================================================
# HELPERS — CORRELATION / LM FOR PLOT 5
# =============================================================================

get_cor_results <- function(columns, run_controlled = FALSE, as_col = "Malignant_Mean_AS") {
  res_df <- data.frame()
  for (col in columns) {
    cell_name <- gsub("_", " ", gsub("pct_", "", col))

    if (run_controlled) {
      try({
        fit   <- lm(as.formula(paste(col, "~", as_col, "+ cancer_type")), data = data)
        coefs <- summary(fit)$coefficients[as_col, ]
        res_df <- rbind(res_df, data.frame(
          cell_type = cell_name,
          val       = coefs["Estimate"],
          p_value   = coefs["Pr(>|t|)"],
          type      = "Controlled (LM Estimate)"
        ))
      }, silent = TRUE)
    } else {
      test <- cor.test(data[[as_col]], data[[col]], method = "spearman")
      res_df <- rbind(res_df, data.frame(
        cell_type = cell_name,
        val       = as.numeric(test$estimate),
        p_value   = test$p.value,
        type      = "spearman Rho"
      ))
    }
  }
  res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
  res_df <- res_df %>% mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    )
  )
  return(res_df)
}

make_bar_plot <- function(df, ytitle, w = 10, h = 6, fname) {
  p <- ggplot(df, aes(x = reorder(cell_type, val), y = val)) +
    geom_bar(stat = "identity", aes(fill = val < 0), width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
      guide  = "none"
    ) +
    graph_theme +
    labs(x = "Cell Type", y = ytitle) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(hjust = 1)
    ) +
    geom_text(
      aes(label = ifelse(sig == "ns", "", sig),
          hjust = ifelse(val < 0, 1.1, -0.1)),
      size = 6
    )
  invisible(p)
}

# =============================================================================
# PLOT 5a — IMMUNE CELLS vs ANEUPLOIDY (cancer-type controlled, malignant-cell AS)
# =============================================================================

message("\n--- Plot 5a: immune barplot - Malignant_Mean_AS ---")

res_immune_mal <- get_cor_results(immune_cols, run_controlled = TRUE, as_col = "Malignant_Mean_AS")

make_bar_plot(
  df     = res_immune_mal,
  ytitle = "Assosiation between Immune Cell Infiltration and\nMalignant-Cell Aneuploidy Score (LM coefficient)",
  fname  = "05a_immune_controlled_barplot_malignant_AS"
)

write.csv(res_immune_mal,
          paste0(output_dir, "05a_immune_lm_malignant_AS.csv"),
          row.names = FALSE)

# =============================================================================
# PLOT 5b — IMMUNE CELLS vs ANEUPLOIDY (cancer-type controlled, pseudo-bulk AS)
# =============================================================================

message("\n--- Plot 5b: immune barplot - pseudo-bulk Mean_AS ---")

if (mean_as_available && "Mean_AS" %in% colnames(data) &&
    sum(!is.na(data$Mean_AS)) >= 10) {

  res_immune_bulk <- get_cor_results(immune_cols, run_controlled = TRUE, as_col = "Mean_AS")

  make_bar_plot(
    df     = res_immune_bulk,
    ytitle = "Assosiation between Immune Cell Infiltration and\nPsuedo-Bulk Aneuploidy Score (LM coefficient)",
    fname  = "05b_immune_controlled_barplot_bulk_AS"
  )

  write.csv(res_immune_bulk,
            paste0(output_dir, "05b_immune_lm_bulk_AS.csv"),
            row.names = FALSE)

} else {
  message("  Skipped: Mean_AS not available or insufficient values")
}

# =============================================================================
# PLOT 6 — PSEUDO-BULK AS vs TUMOR PURITY
# =============================================================================

message("\n--- Plot 6: pseudo-bulk AS vs Purity ---")

if (mean_as_available && "Mean_AS" %in% colnames(data) &&
    sum(!is.na(data$Mean_AS) & !is.na(data$percent_malignant)) >= 10) {

  present_cancer_types_p6 <- unique(
    data$cancer_type[!is.na(data$Mean_AS) & !is.na(data$percent_malignant)]
  )
  cancer_cols_p6 <- cancer_cols[names(cancer_cols) %in% present_cancer_types_p6]
  missing_p6     <- setdiff(present_cancer_types_p6, names(cancer_cols_p6))
  if (length(missing_p6) > 0) {
    cancer_cols_p6 <- c(cancer_cols_p6, setNames(rep("grey50", length(missing_p6)), missing_p6))
  }

  p_purity_bulk <- ggplot(
    data %>% filter(!is.na(Mean_AS), !is.na(percent_malignant)),
    aes(x = percent_malignant, y = Mean_AS)
  ) +
    geom_point(aes(color = cancer_type), alpha = 0.55, size = 2) +
    geom_smooth(method = "lm", col = alpha('blue', 0.5), level = 0.95) +
    stat_cor(
      method      = "spearman",
      size        = 5,
      label.y.npc = 0.95,
      label.x.npc = "right",
      hjust       = 1,
      output.type = "text",
      aes(label = paste0("rho = ", round(after_stat(r), 2), "\np ",
                         ifelse(after_stat(p) < 2.2e-16, "< 2.2e-16",
                                paste0("= ", signif(after_stat(p), 3)))))
    ) +
    scale_color_manual(values = cancer_cols_p6, name = "Cancer Type") +
    graph_theme +
    labs(x = "Tumor Purity", y = "Pseudo-Bulk Aneuploidy Score") +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      legend.position = "right",
      legend.key.size = unit(0.4, 'cm'),
      legend.text     = element_text(size = 9),
      legend.title    = element_text(size = 10),
      axis.title.y    = element_text(hjust = 0.85)
    )
} else {
  message("  Skipped: Mean_AS not available or insufficient values")
}

# =============================================================================
# SUMMARY TABLES
# =============================================================================

message("\n--- Summary tables ---")

pancancer_results <- data.frame()
for (col in immune_cols) {
  test <- cor.test(data$Malignant_Mean_AS, data[[col]], method = "spearman")
  pancancer_results <- rbind(pancancer_results, data.frame(
    cell_type = gsub("pct_", "", col),
    rho       = as.numeric(test$estimate),
    p_value   = test$p.value,
    n         = nrow(data)
  ))
}
pancancer_results$p_adj <- p.adjust(pancancer_results$p_value, method = "BH")
pancancer_results       <- pancancer_results %>% arrange(rho)

write.csv(pancancer_results,
          paste0(output_dir, "pancancer_correlations_malignant_AS.csv"),
          row.names = FALSE)

lm_model <- lm(total_immune_pct ~ Malignant_Mean_AS + cancer_type, data = data)
lm_sum   <- summary(lm_model)
as_coef  <- lm_sum$coefficients["Malignant_Mean_AS", ]

write.csv(
  data.frame(
    term      = "Malignant_Mean_AS",
    estimate  = as_coef["Estimate"],
    std_error = as_coef["Std. Error"],
    t_value   = as_coef["t value"],
    p_value   = as_coef["Pr(>|t|)"],
    r_squared = lm_sum$r.squared,
    n_samples = nrow(data)
  ),
  paste0(output_dir, "lm_results_malignant_AS.csv"),
  row.names = FALSE
)

message("\n", strrep("=", 65))
message("DONE — results written to: ", output_dir)
message(strrep("=", 65), "\n")
