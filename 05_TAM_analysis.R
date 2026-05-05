library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)
library(pheatmap)
library(ragg)

options(bitmapType = "cairo")
options(future.globals.maxSize = 8000 * 1024^2)  # needed for Harmony
set.seed(42)

# =============================================================================
# TOGGLES & SETTINGS
# =============================================================================

SOLID_TUMORS_ONLY              <- TRUE
SKIP_FILTERS                   <- FALSE
MIN_TME_CELLS_PER_SAMPLE       <- 0
MIN_MALIGNANT_CELLS_PER_SAMPLE <- 0
MAX_OTHER_PCT                  <- 40
EXCLUDE_STUDIES                <- c()

# A cluster is labelled "Uncertain" if its winning score <= 0 OR the gap
# between 1st and 2nd score < ANN_MIN_MARGIN.
ANN_MIN_MARGIN <- 0.0

# Labels to merge into a single subtype for proportions / stats / plots.
SUBTYPE_COLLAPSE <- c(
  "M2_like"   = "M2_like_TAM",
  "TAM_C1Q"   = "M2_like_TAM",
  "TAM_SPP1"  = "M2_like_TAM",
  "TAM_ISG15" = "M2_like_TAM",
  "TAM_FN1"   = "M2_like_TAM"
)

# =============================================================================
# COLUMN NAMES
# =============================================================================

CELLTYPE_COL  <- "cell_type_full"
SAMPLE_COL    <- "match_id"
MALIGNANT_VAL <- "Malignant"

# =============================================================================
# PATHS
# =============================================================================

if (SOLID_TUMORS_ONLY) {
  base_dir <- ".../TME_Comprehensive_Analysis_Solid_Only/"
} else {
  base_dir <- ".../TME_Comprehensive_Analysis/"
}

merged_rds <- ".../merged_seurat.rds"

samples_file <- ".../samples_metadata.csv"

filter_tag <- if (SKIP_FILTERS) "unfiltered" else "filtered"
output_dir <- paste0(base_dir, "Macrophage_Subtype_Analysis_", filter_tag, "/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Composition table from fig_5.1 (used for pct_Unassigned filter)
composition_file <- paste0(base_dir, "Publication_Figures_", filter_tag,
                            "/cell_type_percent_of_all_cells.csv")

message("Output directory: ", output_dir)

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
# HELPER FUNCTIONS
# =============================================================================

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
    labs(x = "Macrophage Subtype", y = ytitle) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_text(
      aes(label = ifelse(sig == "ns", "", sig),
          hjust = ifelse(val < 0, 1.1, -0.1)),
      size = 6
    )
  invisible(p)
}

# =============================================================================
# STEP 1 — LOAD DATA
# =============================================================================

message("\n=== Step 1: Loading data ===")

seu <- readRDS(merged_rds)
DefaultAssay(seu) <- "RNA"
message("  Total cells: ", ncol(seu))

samples_metadata <- read.csv(
  samples_file, sep = ",", stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM", check.names = FALSE
)
if (colnames(samples_metadata)[1] %in% c("", NA_character_)) {
  samples_metadata <- samples_metadata[, -1]
}

samples_metadata <- samples_metadata %>%
  mutate(
    Paper    = trimws(Paper),
    Sample   = trimws(Sample),
    match_id = gsub("\\.", "_", paste(Paper, Sample, sep = "_"))
  )

if (SOLID_TUMORS_ONLY) {
  samples_metadata <- samples_metadata %>%
    filter(!grepl("(?i)hematologic", Cancer_type, perl = TRUE))
}

as_per_sample <- samples_metadata %>%
  select(match_id, Malignant_Mean_AS) %>%
  filter(!is.na(Malignant_Mean_AS)) %>%
  distinct(match_id, .keep_all = TRUE)

message("  Samples with Malignant_Mean_AS: ", nrow(as_per_sample))

if (!SKIP_FILTERS && file.exists(composition_file)) {
  pct_table <- read.csv(composition_file, stringsAsFactors = FALSE, check.names = FALSE)
  if ("pct_Unassigned" %in% colnames(pct_table)) {
    keep_samples  <- pct_table %>% filter(pct_Unassigned <= MAX_OTHER_PCT) %>% pull(sample)
    n_before      <- nrow(as_per_sample)
    as_per_sample <- as_per_sample %>% filter(match_id %in% keep_samples)
    message("  Unassigned filter (<=", MAX_OTHER_PCT, "%): kept ",
            nrow(as_per_sample), " of ", n_before, " samples")
  } else {
    warning("composition_file found but has no pct_Unassigned column — skipping filter")
  }
} else if (!SKIP_FILTERS) {
  warning("composition_file not found — skipping Unassigned filter: ", composition_file)
}

# =============================================================================
# STEP 2 — ATTACH AS TO SEURAT, SUBSET TO MACROPHAGES
# =============================================================================

message("\n=== Step 2: Subsetting to macrophages ===")

# match() avoids rownames round-trip
idx <- match(seu@meta.data[[SAMPLE_COL]], as_per_sample$match_id)
seu@meta.data$Malignant_Mean_AS <- as_per_sample$Malignant_Mean_AS[idx]

message("  Cells with Malignant_Mean_AS: ",
        sum(!is.na(seu@meta.data$Malignant_Mean_AS)))

mac <- subset(seu, subset = cell_type_full == "Macrophage")

total_cells_per_sample <- seu@meta.data %>%
  group_by(sample = .data[[SAMPLE_COL]]) %>%
  summarise(total_cells = n(), .groups = "drop")

rm(seu); gc()

message("  Macrophage cells: ", ncol(mac))

keep_cells <- rownames(mac@meta.data)[!is.na(mac@meta.data$Malignant_Mean_AS)]
mac        <- subset(mac, cells = keep_cells)
message("  Macrophages with AS score: ", ncol(mac))

mac@meta.data$cancer_type <- ifelse(
  mac@meta.data$cancer_type == "HN", "Head and Neck", mac@meta.data$cancer_type
)

# =============================================================================
# STEP 3 — HARMONY INTEGRATION + UNSUPERVISED CLUSTERING
#
# group.by = cancer_type removes tissue-of-origin batch effects while
# preserving shared macrophage biology across tissues.
# Clustering precedes marker scoring — labels are assigned post-hoc so
# clusters are not defined by the genes we will test.
# =============================================================================

library(harmony)

CLUSTER_RESOLUTIONS <- c(0.3, 0.5)
DEFAULT_RES         <- 0.3

mac_rds <- paste0(output_dir, "mac_clustered.rds")

message("\n=== Step 3: Harmony integration + clustering ===")

res_col_needed <- paste0("cluster_res", DEFAULT_RES)
cache_ok <- file.exists(mac_rds) &&
  (res_col_needed %in% colnames(readRDS(mac_rds)@meta.data))

if (cache_ok) {
  message("  Found cached object with res=", DEFAULT_RES, " — loading: ", mac_rds)
  mac <- readRDS(mac_rds)
} else {
  if (file.exists(mac_rds))
    message("  Cache missing resolution ", DEFAULT_RES, " — re-clustering")

  mac <- mac %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)

  message("  Running Harmony (group.by.vars = 'cancer_type')...")
  mac <- RunHarmony(
    mac,
    group.by.vars  = "cancer_type",
    reduction      = "pca",
    reduction.save = "harmony",
    verbose        = FALSE
  )

  mac <- FindNeighbors(mac, reduction = "harmony", dims = 1:20, verbose = FALSE)

  for (res in CLUSTER_RESOLUTIONS) {
    mac <- FindClusters(mac, resolution = res, verbose = FALSE)
    col <- paste0("cluster_res", res)
    mac@meta.data[[col]] <- mac@meta.data$seurat_clusters
    message("  Resolution ", res, ": ",
            length(unique(mac@meta.data$seurat_clusters)), " clusters")
  }

  mac <- RunUMAP(mac, reduction = "harmony", dims = 1:20, verbose = FALSE)

  message("  Saving clustered object to: ", mac_rds)
  saveRDS(mac, mac_rds)
}

Idents(mac) <- paste0("cluster_res", DEFAULT_RES)
mac@meta.data$mac_cluster <- mac@meta.data[[paste0("cluster_res", DEFAULT_RES)]]
all_clusters <- levels(mac@meta.data$mac_cluster)
n_clusters   <- length(all_clusters)

message("  Using resolution ", DEFAULT_RES, " (", n_clusters, " clusters) for downstream analysis")
message("  Cell counts per cluster:")
print(table(mac@meta.data$mac_cluster))

# =============================================================================
# STEP 4 — POST-HOC CLUSTER ANNOTATION
#
# Each cluster is assigned the canonical type with the highest mean module
# score. This is annotation evidence, not the clustering basis.
# =============================================================================

message("\n=== Step 4: Post-hoc cluster annotation ===")

canonical_markers <- list(

  M1_like = c(
    "CXCL9",   "CXCL10",  "CXCL11",
    "IL1A",    "IL1B",    "IL6",     "TNF",
    "CD80",    "CD86",    "CD40",
    "IRF5",    "IRF1",    "IDO1",    "CCL5",
    "HLA-DRA", "HLA-DRB1"
  ),

  M2_like = c(
    "CD163",   "MRC1",    "MSR1",    "CSF1R",
    "CD274",   "PDCD1LG2","CD276",   "VTCN1",
    "IL4R",    "IL1RN",   "IL1R2",   "CD200R1",
    "CCL4",    "CCL13",   "CCL17",   "CCL18",   "CCL22",
    "VEGFA",   "CTSA",    "CTSB",    "CTSD",
    "TGFB1",   "MMP9",    "MMP14",   "MMP19",   "FN1",
    "CLEC7A",  "LYVE1",   "IRF4",    "TNFSF12"
  ),

  TAM_C1Q   = c("C1QA",  "C1QB",   "C1QC"),
  TAM_SPP1  = c("SPP1",  "GPNMB",  "FABP5"),
  TAM_ISG15 = c("ISG15", "IFIT1",  "IFITM1", "CXCL10", "CASP1", "TNFSF10"),
  TAM_FN1   = c("FN1")
)

# Deduplicate each list (duplicate genes get double-weighted in AddModuleScore)
canonical_markers <- lapply(canonical_markers,
                             function(g) unique(intersect(g, rownames(mac))))

mac <- AddModuleScore(mac, features = canonical_markers, name = "AnnScore_")
ann_cols <- grep("^AnnScore_", colnames(mac@meta.data), value = TRUE)
ann_cols <- ann_cols[seq_along(canonical_markers)]
names(ann_cols) <- names(canonical_markers)

cluster_scores <- mac@meta.data %>%
  group_by(cluster = mac_cluster) %>%
  summarise(across(all_of(unname(ann_cols)), mean, na.rm = TRUE),
            .groups = "drop")
colnames(cluster_scores)[-1] <- names(canonical_markers)

write.csv(cluster_scores,
          paste0(output_dir, "cluster_annotation_scores.csv"),
          row.names = FALSE)

message("  Cluster annotation scores (review to verify labels):")
print(cluster_scores)

# Label = type with highest mean score. "Uncertain" if winning score <= 0
# or margin between 1st and 2nd score < ANN_MIN_MARGIN.
cluster_labels <- cluster_scores %>%
  pivot_longer(-cluster, names_to = "type", values_to = "score") %>%
  group_by(cluster) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  summarise(
    top_type  = first(type),
    top_score = first(score),
    sec_score = nth(score, 2),
    margin    = top_score - sec_score,
    .groups   = "drop"
  ) %>%
  mutate(label = case_when(
    top_score <= 0           ~ "Uncertain",
    margin < ANN_MIN_MARGIN  ~ "Uncertain",
    TRUE                     ~ top_type
  ))

message("  Cluster annotation (margin threshold = ", ANN_MIN_MARGIN, "):")
print(cluster_labels %>% select(cluster, label, top_score, margin))

message("  Cluster label assignments:")
print(cluster_labels)

label_map <- setNames(cluster_labels$label, as.character(cluster_labels$cluster))
mac@meta.data$mac_label   <- label_map[as.character(mac@meta.data$mac_cluster)]
mac@meta.data$mac_subtype <- mac@meta.data$mac_label

if (length(SUBTYPE_COLLAPSE) > 0) {
  for (st in names(SUBTYPE_COLLAPSE)) {
    mac@meta.data$mac_subtype[mac@meta.data$mac_subtype == st] <- SUBTYPE_COLLAPSE[[st]]
  }
}

all_subtypes <- sort(unique(na.omit(mac@meta.data$mac_subtype)))
n_subtypes   <- length(all_subtypes)
message("  Subtypes present (", n_subtypes, "): ",
        paste(all_subtypes, collapse = ", "))

# =============================================================================
# STEP 5 — BUILD SAMPLE-LEVEL PROPORTION TABLE
# =============================================================================

message("\n=== Step 5: Building sample-level proportion table ===")

# Full sample universe: all samples with an AS score (including 0-macrophage samples)
all_samples_meta <- as_per_sample %>%
  rename(sample = match_id) %>%
  left_join(
    samples_metadata %>%
      mutate(sample = gsub("\\.", "_", paste(Paper, Sample, sep = "_"))) %>%
      distinct(sample, Cancer_type) %>%
      rename(cancer_type = Cancer_type),
    by = "sample"
  )

cell_counts <- mac@meta.data %>%
  filter(!is.na(mac_subtype)) %>%
  group_by(sample = .data[[SAMPLE_COL]], mac_subtype) %>%
  summarise(n_cells = n(), .groups = "drop")

sample_totals <- cell_counts %>%
  group_by(sample) %>%
  summarise(total_mac = sum(n_cells), .groups = "drop")

sample_props <- expand.grid(
  sample      = all_samples_meta$sample,
  mac_subtype = factor(all_subtypes, levels = all_subtypes),
  stringsAsFactors = FALSE
) %>%
  left_join(cell_counts,            by = c("sample", "mac_subtype")) %>%
  left_join(sample_totals,          by = "sample") %>%
  left_join(total_cells_per_sample, by = "sample") %>%
  mutate(
    n_cells     = replace_na(n_cells, 0),
    total_mac   = replace_na(total_mac, 0),
    total_cells = replace_na(total_cells, 0),
    prop        = ifelse(total_cells == 0, 0, n_cells / total_cells)
  ) %>%
  left_join(all_samples_meta, by = "sample") %>%
  filter(!is.na(cancer_type)) %>%
  group_by(cancer_type) %>%
  mutate(
    as_q30 = quantile(Malignant_Mean_AS, 0.30, na.rm = TRUE),
    as_q70 = quantile(Malignant_Mean_AS, 0.70, na.rm = TRUE),
    aneuploidy_group = factor(case_when(
      Malignant_Mean_AS >= as_q70 ~ "High",
      Malignant_Mean_AS <= as_q30 ~ "Low",
      TRUE                        ~ "Intermediate"
    ), levels = c("Low", "Intermediate", "High"))
  ) %>%
  ungroup()

n_samples <- length(unique(sample_props$sample))
message("  Samples in table: ", n_samples)
message("  High aneuploidy: ",
        length(unique(sample_props$sample[sample_props$aneuploidy_group == "High"])))
message("  Intermediate aneuploidy: ",
        length(unique(sample_props$sample[sample_props$aneuploidy_group == "Intermediate"])))
message("  Low aneuploidy: ",
        length(unique(sample_props$sample[sample_props$aneuploidy_group == "Low"])))

sample_props_wide <- sample_props %>%
  mutate(
    mac_subtype = case_when(
      mac_subtype == "M1_like"     ~ "M1",
      mac_subtype == "M2_like_TAM" ~ "M2_Tam",
      mac_subtype == "Uncertain"   ~ "Uncertain",
      TRUE ~ mac_subtype
    )
  ) %>%
  pivot_wider(
    names_from  = mac_subtype,
    values_from = c(n_cells, prop),
    values_fill = list(n_cells = 0, prop = 0)
  ) %>%
  rename(
    n_M1        = n_cells_M1,
    n_M2_Tam    = n_cells_M2_Tam,
    n_Uncertain = n_cells_Uncertain,
    `%M1`       = prop_M1,
    `%M2_Tam`   = prop_M2_Tam,
    `%Uncertain` = prop_Uncertain
  ) %>%
  select(sample, starts_with("n_"), starts_with("%"), everything())

write.csv(sample_props_wide,
          paste0(output_dir, "sample_subtype_proportions.csv"),
          row.names = FALSE)

# =============================================================================
# STEP 6 — STATISTICAL ANALYSIS
# =============================================================================

message("\n=== Step 6: Statistical analysis ===")

# Wilcoxon: High vs Low per subtype
wilcox_res <- sample_props %>%
  filter(aneuploidy_group %in% c("High", "Low")) %>%
  group_by(mac_subtype) %>%
  summarise(
    n_high    = sum(aneuploidy_group == "High"),
    n_low     = sum(aneuploidy_group == "Low"),
    mean_high = mean(prop[aneuploidy_group == "High"], na.rm = TRUE),
    mean_low  = mean(prop[aneuploidy_group == "Low"],  na.rm = TRUE),
    p_value   = tryCatch(
      wilcox.test(prop[aneuploidy_group == "High"], prop[aneuploidy_group == "Low"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    )
  )

write.csv(wilcox_res,
          paste0(output_dir, "wilcoxon_subtype_high_vs_low.csv"),
          row.names = FALSE)

message("  Wilcoxon results:")
print(wilcox_res %>% select(mac_subtype, mean_high, mean_low, p_adj, sig))

# Cancer-type-controlled LM per subtype
lm_res <- lapply(all_subtypes, function(st) {
  df_s <- sample_props %>%
    filter(mac_subtype == st, !is.na(Malignant_Mean_AS), !is.na(cancer_type))
  if (nrow(df_s) < 10) return(NULL)
  tryCatch({
    fit   <- lm(prop ~ Malignant_Mean_AS + cancer_type, data = df_s)
    coefs <- summary(fit)$coefficients["Malignant_Mean_AS", ]
    data.frame(
      subtype  = st,
      estimate = coefs[["Estimate"]],
      p_value  = coefs[["Pr(>|t|)"]]
    )
  }, error = function(e) NULL)
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
      p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
    )
  )

write.csv(lm_res,
          paste0(output_dir, "lm_subtype_cancer_controlled.csv"),
          row.names = FALSE)

message("  LM results:")
print(lm_res %>% select(subtype, estimate, p_adj, sig))

# =============================================================================
# RECODE SUBTYPE DISPLAY NAMES
# =============================================================================

display_names  <- c("M1_like" = "M1-like", "M2_like_TAM" = "M2-like/TAM")
recode_subtype <- function(x) ifelse(x %in% names(display_names), display_names[x], x)

sample_props$mac_subtype  <- recode_subtype(sample_props$mac_subtype)
wilcox_res$mac_subtype    <- recode_subtype(wilcox_res$mac_subtype)
lm_res$subtype            <- recode_subtype(lm_res$subtype)
all_subtypes              <- recode_subtype(all_subtypes)
mac@meta.data$mac_subtype <- recode_subtype(mac@meta.data$mac_subtype)
mac@meta.data$mac_label   <- recode_subtype(mac@meta.data$mac_label)

# =============================================================================
# COLOURS
# =============================================================================

cluster_pal <- setNames(
  colorRampPalette(brewer.pal(min(n_clusters, 12), "Paired"))(n_clusters),
  all_clusters
)

subtype_colors <- setNames(scales::hue_pal()(length(all_subtypes)), all_subtypes)
subtype_colors["Uncertain"] <- "#AAAAAA"
label_colors <- subtype_colors

subtype_colors_034 <- c("M1-like" = "#e2a150", "M2-like/TAM" = "#dfa9df", "Uncertain" = "#AAAAAA")
for (s in setdiff(all_subtypes, names(subtype_colors_034))) {
  subtype_colors_034[s] <- subtype_colors[s]
}

group_colors <- c("High" = "#004D40", "Low" = "#B2DFDB", "Intermediate" = "#26A69A")

# =============================================================================
# PLOT 3a — STACKED BAR (High vs Low, mean % of all cells)
# =============================================================================

message("\n--- Plot 3a: stacked bar ---")

stack_df <- sample_props %>%
  filter(aneuploidy_group %in% c("High", "Low"), total_cells > 0) %>%
  group_by(aneuploidy_group, mac_subtype) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

p3a <- ggplot(stack_df,
              aes(x = aneuploidy_group, y = mean_prop, fill = mac_subtype)) +
  geom_bar(stat = "identity", width = 0.55, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = subtype_colors_034, name = "Macrophage Subtype") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = NULL, x = "Aneuploidy Group", y = "Mean % of All Cells") +
  graph_theme

# =============================================================================
# PLOT 3b — VIOLIN + BOXPLOT (High vs Low per subtype)
# =============================================================================

message("\n--- Plot 3b: violin+boxplot ---")

p3b <- ggplot(
  sample_props %>% filter(aneuploidy_group %in% c("High", "Low"),
                           total_cells > 0, mac_subtype != "Uncertain"),
  aes(x = aneuploidy_group, y = prop, fill = aneuploidy_group)
) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_beeswarm(aes(color = aneuploidy_group), size = 1.2, alpha = 0.5, cex = 2) +
  stat_compare_means(
    method      = "wilcox.test",
    comparisons = list(c("Low", "High")),
    label       = "p.format",
    method.args = list(alternative = "two.sided"),
    size        = 5
  ) +
  facet_wrap(~ mac_subtype, scales = "free_y") +
  scale_fill_manual(
    values = group_colors,
    breaks = c("High", "Low"),
    labels = c("High AS", "Low AS"),
    name   = "Sample AS"
  ) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_x_discrete(labels = c("Low" = "Low AS", "High" = "High AS")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = NULL, x = "Aneuploidy Group", y = "Mean % of All Cells") +
  graph_theme +
  theme(legend.position = "right")

# =============================================================================
# PLOT 4 — SCATTER (proportion vs continuous AS, per subtype)
# =============================================================================

message("\n--- Plot 4: scatter ---")

cor_labels_04 <- sample_props %>%
  filter(total_cells > 0, mac_subtype != "Uncertain", !is.na(Malignant_Mean_AS)) %>%
  group_by(mac_subtype) %>%
  summarise(
    r = cor(Malignant_Mean_AS, prop, use = "complete.obs"),
    p = cor.test(Malignant_Mean_AS, prop)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("rho = ", round(r, 2), "\np ",
                         ifelse(p < 2.2e-16, "< 2.2e-16", paste0("= ", signif(p, 3)))))

p4 <- ggplot(
  sample_props %>% filter(total_cells > 0, mac_subtype != "Uncertain"),
  aes(x = Malignant_Mean_AS, y = prop, color = mac_subtype)
) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", col = alpha("blue", 0.5), level = 0.95) +
  geom_text(data = cor_labels_04,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.5,
            inherit.aes = FALSE, size = 5) +
  facet_wrap(~ mac_subtype, scales = "free_y") +
  scale_color_manual(values = subtype_colors_034, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = NULL, x = "Malignant-Cell Aneuploidy Score", y = "Mean % of All Cells") +
  graph_theme

# =============================================================================
# PLOT 5 — LM BARPLOT (cancer-type controlled effect per subtype)
# =============================================================================

message("\n--- Plot 5: LM barplot ---")

p5 <- ggplot(lm_res %>% filter(subtype != "Uncertain"),
             aes(x = reorder(subtype, estimate), y = estimate,
                 fill = estimate > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(sig == "ns", "", sig),
                hjust = ifelse(estimate < 0, 1.1, -0.1)),
            size = 6) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                    guide  = "none") +
  labs(
    title = NULL,
    x     = "Macrophage Subtype",
    y     = "Assosiation between Macrophages Subtype Infiltration\nand Malignant-Cell Aneuploidy Score (LM coefficient)"
  ) +
  graph_theme +
  theme(axis.title.x = element_text(hjust = 1, size = 20))

# =============================================================================
# PLOT 6 — ALL IMMUNE CELLS + MACROPHAGE SUBTYPES (Macrophage expanded)
# Same as fig_5.1 Plot 5a but replaces the broad "Macrophage" bar
# with M1-like / M2-like/TAM subtypes.
# =============================================================================

message("\n--- Plot 6: immune + macrophage subtypes LM barplot ---")

immune_cell_names_06 <- c(
  "B_cell", "Dendritic", "Granulocyte", "Immune_Other",
  "Mast", "Monocyte", "NK_cell", "Plasma", "T_cell"
)

if (file.exists(composition_file)) {
  pct_data <- read.csv(composition_file, stringsAsFactors = FALSE,
                       check.names = FALSE) %>%
    filter(!is.na(Malignant_Mean_AS)) %>%
    mutate(cancer_type = ifelse(cancer_type == "HN", "Head and Neck", cancer_type))

  if (!SKIP_FILTERS) {
    if (length(EXCLUDE_STUDIES) > 0) {
      exclude_pattern <- paste0("^", EXCLUDE_STUDIES, collapse = "|")
      pct_data <- pct_data %>% filter(!grepl(exclude_pattern, sample))
    }
    pct_data <- pct_data %>%
      filter(n_nonmalignant >= MIN_TME_CELLS_PER_SAMPLE,
             n_malignant    >= MIN_MALIGNANT_CELLS_PER_SAMPLE)
    if ("pct_Unassigned" %in% colnames(pct_data))
      pct_data <- pct_data %>% filter(pct_Unassigned <= MAX_OTHER_PCT)
  }

  immune_cols_06 <- paste0("pct_", immune_cell_names_06)
  immune_cols_06 <- immune_cols_06[immune_cols_06 %in% colnames(pct_data)]

  lm_immune_06 <- lapply(immune_cols_06, function(col) {
    cell_name <- gsub("_", " ", gsub("^pct_", "", col))
    tryCatch({
      fit   <- lm(as.formula(paste(col, "~ Malignant_Mean_AS + cancer_type")),
                  data = pct_data)
      coefs <- summary(fit)$coefficients["Malignant_Mean_AS", ]
      data.frame(cell_type = cell_name,
                 val       = coefs[["Estimate"]],
                 p_value   = coefs[["Pr(>|t|)"]])
    }, error = function(e) NULL)
  }) %>% bind_rows()

  lm_mac_sub_06 <- lm_res %>%
    filter(subtype != "Uncertain") %>%
    select(cell_type = subtype, val = estimate, p_value)

  lm_combined_06 <- bind_rows(lm_immune_06, lm_mac_sub_06) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      sig   = case_when(
        p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
        p_adj < 0.05  ~ "*",   TRUE          ~ "ns"
      )
    )

  write.csv(lm_combined_06,
            paste0(output_dir, "06_immune_mac_subtypes_lm_malignant_AS.csv"),
            row.names = FALSE)

  make_bar_plot(
    df     = lm_combined_06,
    ytitle = "Association between Immune Cell Infiltration\nand Malignant-Cell Aneuploidy Score (LM coefficient)",
    fname  = "06_immune_mac_subtypes_lm_barplot"
  )
} else {
  message("  Skipped: composition file not found at ", composition_file)
}

# =============================================================================
# DONE
# =============================================================================

message("\n", strrep("=", 65))
message("DONE — all files written to: ", output_dir)
message(strrep("=", 65))
