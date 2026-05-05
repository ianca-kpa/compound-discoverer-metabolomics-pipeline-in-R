# =============================================================================
# 01_validation.R
# Settings validation
# =============================================================================

validate_settings <- function() {
  # ---------------------------------------------------------------------------
  # PCA
  # ---------------------------------------------------------------------------
  stopifnot(pca_scaling %in% c("none", "pareto", "autoscale"))

  # ---------------------------------------------------------------------------
  # Filters
  # ---------------------------------------------------------------------------
  stopifnot(low_variance_filter_method %in% c("none", "iqr"))

  # ---------------------------------------------------------------------------
  # Duplicate handling
  # ---------------------------------------------------------------------------
  stopifnot(
    duplicate_name_strategy %in% c(
      "reference_or_best_qc_rsd",
      "keep_separate",
      "collapse_mean",
      "collapse_sum",
      "collapse_best_qc_rsd"))

  # ---------------------------------------------------------------------------
  # Reference file usage / overrides
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(use_reference_file), length(use_reference_file) == 1)

  if (isTRUE(use_reference_file)) {
    stopifnot(is.character(comparison_path), length(comparison_path) == 1)
    stopifnot(nzchar(trimws(comparison_path)))
  }

  stopifnot(is.character(reference_col_metabolite), length(reference_col_metabolite) == 1)
  stopifnot(is.character(reference_col_ref_ion), length(reference_col_ref_ion) == 1)
  stopifnot(is.character(reference_col_mz), length(reference_col_mz) == 1)
  stopifnot(is.character(reference_col_rt), length(reference_col_rt) == 1)

  # ---------------------------------------------------------------------------
  # Name sanitation
  # ---------------------------------------------------------------------------
  stopifnot(sanitize_mode %in% c("greek_latin_ascii", "ascii_translit"))

  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------
  stopifnot(all(run_metrics %in% c("FDR", "p_value", "FDR_and_p_value")))
  if (exists("heatmap_rank_metrics", inherits = TRUE)) {
    stopifnot(all(heatmap_rank_metrics %in% c("FDR", "p_value", "FDR_and_p_value")))
  }

  # ---------------------------------------------------------------------------
  # Heatmap
  # ---------------------------------------------------------------------------
  stopifnot(heatmap_scale_method %in% c("none", "zscore", "pareto"))
  stopifnot(heatmap_cluster_distance %in% c("euclidean", "manhattan"))
  stopifnot(heatmap_cluster_method %in% c("ward.D2", "complete", "average"))

  # ---------------------------------------------------------------------------
  # Volcano main style
  # ---------------------------------------------------------------------------
  stopifnot(volcano_style %in% c("classic", "gradual", "both"))

  # ---------------------------------------------------------------------------
  # Volcano classic legend
  # ---------------------------------------------------------------------------
  stopifnot(is.character(volcano_classic_legend_title))
  stopifnot(length(volcano_classic_legend_title) == 1)

  # classic colors must be exactly 3
  stopifnot(length(volcano_classic_fills) == 3)
  stopifnot(length(volcano_classic_colors) == 3)

  # ---------------------------------------------------------------------------
  # Volcano gradual legend
  # ---------------------------------------------------------------------------
  stopifnot(is.character(volcano_gradual_legend_title))
  stopifnot(length(volcano_gradual_legend_title) == 1)

  stopifnot(is.numeric(volcano_gradual_legend_breaks))
  stopifnot(length(volcano_gradual_legend_breaks) >= 2)

  stopifnot(is.numeric(volcano_gradual_legend_limits))
  stopifnot(length(volcano_gradual_legend_limits) == 2)

  # ---------------------------------------------------------------------------
  # Volcano gradual palette
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_gradual_use_RColorBrewer))

  if (!volcano_gradual_use_RColorBrewer) {
    stopifnot(length(volcano_gradual_fills) >= 3)
    stopifnot(length(volcano_gradual_colors) >= 3)
  }

  # ---------------------------------------------------------------------------
  # Volcano labels
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_add_labels))
  stopifnot(is.numeric(volcano_label_number))
  stopifnot(volcano_label_number >= 0)

  # ---------------------------------------------------------------------------
  # Axis scaling
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_auto_axis))
  stopifnot(is.numeric(volcano_axis_expand_mult))
  stopifnot(volcano_axis_expand_mult > 0)
}
