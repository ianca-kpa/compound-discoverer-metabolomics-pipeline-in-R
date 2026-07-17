# ==========================================================
# SETTINGS FILE
# Copy this file to:
# config/settings.R
#
# Then edit the paths below before running the pipeline.
# ==========================================================

# Input files
cd_file_path <- "data/Compounds_POS_example.xlsx"
cd_sheet <- 1                                  # options: sheet index or sheet name

metadata_path <- "data/metadata_example.xlsx"
metadata_sheet <- 1                            # options: sheet index or sheet name

injection_order_path <- ""                     # optional: path to an injection-order file; leave "" to use metadata/sample order

reference_path <- "data/reference_example.xlsx"
reference_sheet <- 1                           # options: sheet index or sheet name
use_reference_file <- TRUE                     # options: TRUE/FALSE (when FALSE, reference file is ignored)

# Optional manual reference-column names (leave "" for auto-detect)
reference_col_metabolite <- ""                 # e.g. "metabolite"
reference_col_ref_ion <- ""                    # e.g. "Ref ion"
reference_col_mz <- ""                         # e.g. "m/z"
reference_col_rt <- ""                         # e.g. "RT C18 25min 0.3mL"

# Comparison group labels used by comparisons, heatmaps, and stats
comparison_group_control <- "WT"
comparison_group_treatment <- "TG"
model_allowed_groups_by_model <- NULL          # optional: named vector, e.g. c("ModelA" = "WT,TG", "ModelB" = "WT,TG")

# Output directory (root)
# This is only the folder name where this run will be written. The pipeline does
# not infer the biological meaning of the run from this name.
# Recommended suffixes:
# - FINAL_QCRSC_WEIGHT_LOG2 for the final biological analysis.
# - MA_COMPATIBLE_WEIGHT_LOG2_PARETO for runs configured to compare with
#   MetaboAnalyst-compatible weight/log2/pareto workflows.
# Choose the name manually before each run to keep outputs easy to compare.
output_dir <- "output"

# Preprocessing scenarios after QC-RSD/IQR filtering:
# Scenario 1: normalization_mode <- "QC_RSC"
# Scenario 2: normalization_mode <- "qc_loess"
# Scenario 3: normalization_mode <- "CYCLIC_LOESS"  # limma cyclic LOESS
# Scenario 4: normalization_mode <- "pqn_sample"    # PQN without QC, sample median reference
# Scenario 5: normalization_mode <- "none"
# Weight-only is written internally as normalization_mode <- "weight" when
# use_weight_normalization <- TRUE and normalization_mode <- "none".

# Weight normalization
use_weight_normalization <- FALSE              # options: TRUE/FALSE
stop_on_invalid_weight <- TRUE                 # options: TRUE/FALSE
invalid_weight_to_NA <- TRUE                   # options: TRUE/FALSE

# Normalization scenario
normalization_mode <- "none"                 # options: "none", "qc_loess", "cyclic_loess", "qcrsc", "pqn_qc", "pqn_sample"; "weight" is used internally for weight-only
make_qc_diagnostics <- FALSE                   # options: TRUE/FALSE. TRUE runs optional QC/normalization audit plots; FALSE is faster for routine runs.
apply_qcrsc_spectral_cleaning <- TRUE          # options: TRUE/FALSE
qc_loess_span <- 0.75                          # options: numeric between 0 and 1
qc_loess_min_qc_points <- 4                    # options: integer >= 3
pqn_min_qc_points <- 3                         # options: integer >= 2

# Filters
missing_exclusion_max_fraction <- 0.50         # options: 0..1 (set >=1 to disable)
presence_filter_min_fraction <- 0.00           # options: 0..1
impute_half_min <- TRUE                        # options: TRUE/FALSE

# QC RSD thresholds (variant creation)
rsd_thresholds <- c(20)                        # options: c(10,15,20,30,...) etc.
rsd_filter_type <- "RSD"                   # options: "QC_RSD" = 100 * SD/mean; "RSD" = SD/mean
active_variant <- "none"                       # options: "none" = no RSD variant; or "QC_RSD<threshold>" / "RSD<threshold>" from rsd_thresholds

# Low-variance filter
low_variance_filter_method <- "none"            # options: "none" or "iqr"
low_variance_filter_fraction <- 0.20           # options: 0..1
low_variance_filter_rounding <- "ceiling"      # options: "floor", "ceiling", "round"
# Use "ceiling" to mimic MetaboAnalyst IQR filtering.

# Transformation
log2_offset <- 1 # options: 0, 1, 0.5 ... (avoid 0 if you may have zeros)

# Statistical thresholds
p_value_cutoff <- 0.05                         # options: numeric between 0 and 1 (significance threshold for p-value based stats)
fdr_cutoff <- 0.05                             # options: numeric between 0 and 1 (significance threshold for FDR based stats)
fc_cutoff_log2 <- 0                            # options: numeric >= 0 (log2 fold change cutoff for volcano plot labeling and significant heatmap filtering; set to 0 to disable fold change cutoff)   
alpha_sig <- p_value_cutoff                    # compatibility alias used by significant heatmap code

# Statistical test configuration
# options for statistical_test_type: "student", "welch", "wilcoxon", "limma"
statistical_test_type <- "student"             # options: "student" (Student's t-test), "welch" (Welch's t-test), "wilcoxon" (Wilcoxon rank-sum test), "limma" (unpaired moderated t-test using empirical Bayes across all features)
test_is_paired <- FALSE                        # options: TRUE/FALSE (if TRUE, uses paired test; if FALSE, uses unpaired test)
pvalue_correction_method <- "FDR"              # options: "raw", "FDR", "Bonferroni", "Holm", "Hochberg", "Hommel", "BY" (method for p-value adjustment; "raw" = no correction)

# Comparison mode:
# "pairwise" = current two-group workflow
# "multigroup" = five primary pairwise comparisons plus exploratory global multi-group test and selected follow-up pairs
# "both" = current two-group workflow plus multi-group analysis
comparison_mode <- "pairwise"                  # options: "pairwise", "multigroup", "both"
multigroup_groups <- character(0)              # options: character vector, e.g. c("pre", "post", "recovery"); empty = auto-detect biological groups per model
multigroup_test <- "kruskal"                   # options: "kruskal", "anova", "welch_anova"
multigroup_pairwise_mode <- "selected"         # recommended: "selected"; options: "none", "selected", or explicit legacy "all"
multigroup_pairwise_pairs <- NULL              # options: NULL or character vector, e.g. c("pre vs post", "pre vs recovery")

# Known-only filter
use_only_known <- TRUE                         # options: TRUE/FALSE

# Duplicate metabolite handling
duplicate_name_strategy <- "collapse_best_qc_rsd"   # options: "reference_or_best_qc_rsd"; "keep_separate"; "collapse_mean"; "collapse_sum"; "collapse_best_qc_rsd"

# Duplicate rounding
dup_mz_digits <- 4                             # options: integer >= 0 (number of decimal places to round m/z values to when determining duplicates)
dup_rt_digits <- 2                             # options: integer >= 0 (number of decimal places to round RT values to when determining duplicates)

# Metrics to run
run_metrics <- "FDR_and_p_value"               # options: any subset of c("FDR", "p_value", "FDR_and_p_value")
heatmap_rank_metrics <- "FDR_and_p_value"      # options: any subset of c("FDR", "p_value", "FDR_and_p_value")

# Exports
export_metaboanalyst_ready <- TRUE             # options: TRUE/FALSE (exports a global raw table for MetaboAnalyst with model_group metadata)
export_metaboanalyst_duplicate_only <- TRUE    # options: TRUE/FALSE (exports MetaboAnalyst-style data treated only by duplicate handling)
save_stats_excel_per_model <- TRUE             # options: TRUE/FALSE (saves an Excel file with all the stats for each model in a separate sheet)
save_sig_metabolites_txt_per_model <- TRUE     # options: TRUE/FALSE (saves a plain .txt list with one significant metabolite per line for each comparison)
make_volcano_plots <- TRUE

# PCA / Heatmaps
pca_scaling <- "pareto"                        # options: "none","pareto","autoscale"
pca_label_samples <- TRUE                      # options: TRUE/FALSE

make_heatmap_by_model <- TRUE                  # options: TRUE/FALSE
make_heatmap_by_model_sex <- TRUE              # options: TRUE/FALSE
heatmap_top_n <- 50

heatmap_scale_method <- "zscore"               # options: "none","zscore","pareto"
heatmap_order_samples_by_group <- TRUE         # options: TRUE/FALSE (if TRUE, samples will be ordered by group in the heatmaps; if FALSE, original order from the input data will be kept)
heatmap_cluster_distance <- "euclidean"        # options: "euclidean", "manhattan"
heatmap_cluster_method <- "ward.D2"            # options: "ward.D2", "complete", "average"

heatmap_palette_n <- 101                
heatmap_breaks_symmetric <- TRUE
heatmap_breaks_limit <- 5                      # options: numeric > 0 (max absolute value for the heatmap color breaks when heatmap_breaks_symmetric = TRUE)

# Significant heatmaps
make_sig_heatmap_by_model <- FALSE             # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated showing only the significant features for each comparison)
make_sig_heatmap_by_model_sex <- FALSE         # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated showing only the significant features for each comparison)
make_sig_heatmap_FvsM_within_group <- FALSE    # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated for FvsM within each group)
sig_heatmap_max_features <- 70                 # options: integer > 0 (max number of features to show in the significant heatmaps; set to a high number to include all significant features)   
sig_heatmap_require_fc_cutoff <- TRUE          # options: TRUE/FALSE (if TRUE, only features that pass both the significance threshold and the fold-change cutoff will be included in the significant heatmaps; if FALSE, all features that pass the significance threshold will be included regardless of fold-change)

# =============================================================================
# Volcano plot settings
# =============================================================================

# Main volcano style (only "classic" is supported)
# "classic" = categorical publication-like volcano
volcano_style <-  "classic"                    # options: "classic"

# Automatic axis scaling per plot
volcano_auto_axis <- TRUE                      # options: TRUE/FALSE

# Axis expansion factor when volcano_auto_axis = TRUE
volcano_axis_expand_mult <- 0.08               # options: numeric > 0

# Labels
volcano_add_labels <- TRUE                     # options: TRUE/FALSE
volcano_label_number <- Inf                    # options: integer >= 0 or Inf (use Inf to label all significant points)

# Use NULL for automatic labels.
# Example: volcano_custom_labels <- c("L-Tryptophan", "Corticosterone")
volcano_custom_labels <- NULL                  # options: NULL or character vector

# Threshold lines
volcano_add_cutoff_lines <- TRUE               # options: TRUE/FALSE

# -------------------------
# Classic style settings
# -------------------------
volcano_classic_point_size <- 2.5              # options: numeric > 0
volcano_classic_point_shape <- 21              # options: 21, 16, etc.

# IMPORTANT ORDER:
# c("Down", "Normal", "Up")
volcano_classic_fills <- c("#3B82F6", "#BDBDBD", "#EF4444")
volcano_classic_colors <- c("#1D4ED8", "#7A7A7A", "#B91C1C")
volcano_classic_legend_title <- "Regulation"

# Text sanitation
sanitize_names_for_exports <- FALSE              # options: TRUE/FALSE
sanitize_mode <- "greek_latin_ascii"             # options: "greek_latin_ascii" or "ascii_translit"
strip_stereo_prefixes_for_names <- TRUE          # options: TRUE/FALSE

# Output control
output_level <- "standard"                    # options: "minimal", "standard", "full_debug"
