# ==========================================================
# SETTINGS FILE
# Copy this file to:
# pipeline/config/settings.R
#
# Then edit the paths below before running the pipeline.
# ==========================================================

# Input files
cd_file_path <- "data/Compounds_POS_example.xlsx"
cd_sheet <- 1                                  # options: sheet index or sheet name

metadata_path <- "data/metadata_example.xlsx"
metadata_sheet <- 1                            # options: sheet index or sheet name

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

# Output directory (root)
output_dir <- "output"

# Weight normalization
# If enabled, sample intensities are divided by the sample-weight column
# detected or manually mapped from metadata before the main normalization step.
use_weight_normalization <- FALSE              # options: TRUE/FALSE
stop_on_invalid_weight <- TRUE                 # options: TRUE/FALSE
invalid_weight_to_NA <- TRUE                   # options: TRUE/FALSE

# Main normalization after optional weight normalization
normalization_mode <- "none"                    # options: "none", "PQN", or "QC_LOESS" ("LOESS" is accepted as a legacy alias for "QC_LOESS")
loess_min_qc_points <- 5                       # minimum valid QC points per feature for QC-LOESS
QC_LOESS_span <- 0.75                          # LOESS smoothing span for QC-LOESS
injection_order_path <- ""                     # optional separate sample/run-order table for QC-LOESS

# Filters
missing_exclusion_max_fraction <- 0.50         # options: 0..1 (set >=1 to disable)
presence_filter_min_fraction <- 0.00           # options: 0..1
impute_half_min <- TRUE                        # options: TRUE/FALSE

# RSD thresholds (variant creation)
# "qc_rsd" creates variants such as QC_RSD20 from QC sample RSD.
# "rsd" creates variants such as RSD20 from biological/sample RSD.
# "none" skips RSD-based variant filtering.
rsd_filter_metric <- "none"                    # options: "none", "qc_rsd", "rsd"
rsd_thresholds <- c(20)                        # options: c(10,15,20,30,...) etc.
active_variant <- "none"                       # options: "none", paste0("QC_RSD", rsd_thresholds), or paste0("RSD", rsd_thresholds)

# Low-variance filter
low_variance_filter_method <- "none"           # options: "none" or "iqr"
low_variance_filter_fraction <- 0.20           # options: 0..1

# Transformation
log2_offset <- 1 # options: 0, 1, 0.5 ... (avoid 0 if you may have zeros)

# Statistical thresholds
p_value_cutoff <- 0.05                         # options: numeric between 0 and 1 (significance threshold for p-value based stats)
fdr_cutoff <- 0.05                             # options: numeric between 0 and 1 (significance threshold for FDR based stats)
fc_cutoff_log2 <- 0                            # options: numeric >= 0 (log2 fold change cutoff for volcano plot labeling and significant heatmap filtering; set to 0 to disable fold change cutoff)   
alpha_sig <- p_value_cutoff                    # compatibility alias used by significant heatmap code

# Statistical test configuration
# options for statistical_test_type: "student", "welch", "wilcoxon", "limma"
statistical_test_type <- "student"             # options: "student" (Student's t-test), "welch" (Welch's t-test), "wilcoxon" (Wilcoxon rank-sum test), "limma" (Moderated t-test using empirical Bayes)
test_is_paired <- FALSE                        # options: TRUE/FALSE (if TRUE, uses paired test; if FALSE, uses unpaired test)
pvalue_correction_method <- "FDR"              # options: "raw", "FDR", "Bonferroni", "Holm", "Hochberg", "Hommel", "BY" (method for p-value adjustment; "raw" = no correction)

# Known-only filter
use_only_known <- TRUE                         # options: TRUE/FALSE

# Duplicate metabolite handling
# reference_or_best_qc_rsd uses the reference table to choose the closest RT
# match for duplicate named metabolites, then falls back to best QC RSD.
duplicate_name_strategy <- "collapse_best_qc_rsd"   # options: "reference_or_best_qc_rsd"; "keep_separate"; "collapse_mean"; "collapse_sum"; "collapse_best_qc_rsd"

# Duplicate rounding
dup_mz_digits <- 4                             # options: integer >= 0 (number of decimal places to round m/z values to when determining duplicates)
dup_rt_digits <- 2                             # options: integer >= 0 (number of decimal places to round RT values to when determining duplicates)

# Metrics to run
run_metrics <- "FDR_and_p_value"               # options: any subset of c("FDR", "p_value", "FDR_and_p_value")
heatmap_rank_metrics <- "FDR_and_p_value"      # options: any subset of c("FDR", "p_value", "FDR_and_p_value")

# Exports
export_metaboanalyst_ready <- TRUE             # options: TRUE/FALSE (exports a table with the exact format required for MetaboAnalyst enrichment analysis)
save_stats_excel_per_model <- TRUE             # options: TRUE/FALSE (saves an Excel file with all the stats for each model in a separate sheet)
save_sig_metabolites_txt_per_model <- TRUE     # options: TRUE/FALSE (saves a plain .txt list with one significant metabolite per line for each comparison)
make_volcano_plots <- TRUE

# PCA / Heatmaps
pca_scaling <- "pareto"                        # options: "none","pareto","autoscale"

make_heatmap_by_model <- TRUE                  # options: TRUE/FALSE
make_heatmap_by_model_sex <- TRUE              # options: TRUE/FALSE
heatmap_top_n <- 80

heatmap_scale_method <- "pareto"               # options: "none","zscore","pareto"
heatmap_order_samples_by_group <- TRUE         # options: TRUE/FALSE (if TRUE, samples will be ordered by group in the heatmaps; if FALSE, original order from the input data will be kept)
heatmap_cluster_distance <- "euclidean"        # options: "euclidean", "manhattan"
heatmap_cluster_method <- "ward.D2"            # options: "ward.D2", "complete", "average"

heatmap_palette_n <- 101                
heatmap_breaks_symmetric <- TRUE
heatmap_breaks_limit <- 5                      # options: numeric > 0 (max absolute value for the heatmap color breaks when heatmap_breaks_symmetric = TRUE)

# Significant heatmaps
make_sig_heatmap_by_model <- FALSE             # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated showing only the significant features for each comparison)
make_sig_heatmap_by_model_sex <- FALSE         # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated showing only the significant features for each comparison)
make_sig_heatmap_FvsM_within_group <- FALSE    # options: TRUE/FALSE (if TRUE, additional heatmaps will be generated showing only the significant features for the comparison)
sig_heatmap_max_features <- 70                 # options: integer > 0 (max number of features to show in the significant heatmaps; set to a high number to include all significant features)   
sig_heatmap_require_fc_cutoff <- TRUE          # options: TRUE/FALSE (if TRUE, only features that pass both the significance threshold and the fold-change cutoff will be included in the significant heatmaps; if FALSE, all features that pass the significance threshold will be included regardless of fold-change)

# =============================================================================
# Volcano plot settings
# =============================================================================

# Main volcano style:
# "classic" = categorical publication-like volcano
# "gradual" = continuous gradient volcano
# "both"    = export both styles
volcano_style <-  "classic"                    # options: "classic", "gradual", "both"

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

# -------------------------
# Gradual style settings
# -------------------------
volcano_gradual_point_shape <- 21
volcano_gradual_point_size_range <- c(1.5, 6)

# Default internal palettes
volcano_gradual_fills <- c("#39489F", "#39BBEC", "#F9ED36", "#F38466", "#B81F25")
volcano_gradual_colors <- c("#17194E", "#68BFE7", "#F9ED36", "#A22F27", "#211F1F")

# If TRUE, use RColorBrewer instead of the vectors above
volcano_gradual_use_RColorBrewer <- FALSE        # options: TRUE/FALSE
volcano_gradual_brewer_palette <- "RdYlBu"       # options: any RColorBrewer palette name (e.g. "RdYlBu", "Spectral", "RdBu", "PiYG", etc.)
volcano_gradual_brewer_n <- 5                    # options: integer >= 3 (number of colors to use from the palette)
volcano_gradual_reverse_brewer <- TRUE           # options: TRUE/FALSE (reverse the order of the colors from the RColorBrewer palette)

volcano_gradual_legend_title <- "Significance"   # options: character string = Legend title
volcano_gradual_legend_breaks <- c(1, 2, 3, 4, 5)
volcano_gradual_legend_limits <- c(0, 5)         # options: numeric vector of length 2 (min and max for the legend gradient)

# Text sanitation
sanitize_names_for_exports <- TRUE               # options: TRUE/FALSE
sanitize_mode <- "greek_latin_ascii"             # options: "none", "greek_latin_ascii" or "ascii_translit"

# Output control
minimal_output <- FALSE                          # options: TRUE/FALSE. If TRUE, selected plots/statistics are kept while selected intermediate global exports are skipped.
