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
  if (!exists("rsd_thresholds", inherits = TRUE)) {
    rsd_thresholds <<- c(20)
  }
  if (!exists("active_variant", inherits = TRUE)) {
    active_variant <<- "none"
  }
  if (!exists("rsd_filter_metric", inherits = TRUE)) {
    if (exists("active_variant", inherits = TRUE) && grepl("^RSD", as.character(active_variant), ignore.case = TRUE)) {
      rsd_filter_metric <<- "rsd"
    } else if (exists("active_variant", inherits = TRUE) && grepl("^QC_RSD", as.character(active_variant), ignore.case = TRUE)) {
      rsd_filter_metric <<- "qc_rsd"
    } else {
      rsd_filter_metric <<- "none"
    }
  }
  rsd_filter_metric <- tolower(as.character(rsd_filter_metric))
  stopifnot(rsd_filter_metric %in% c("none", "qc_rsd", "rsd"))
  stopifnot(is.numeric(rsd_thresholds), length(rsd_thresholds) >= 1)
  stopifnot(is.character(active_variant), length(active_variant) == 1)

  # ---------------------------------------------------------------------------
  # Normalization
  # ---------------------------------------------------------------------------
  if (!exists("normalization_mode", inherits = TRUE)) {
    normalization_mode <<- "PQN"
  }
  normalization_mode <- toupper(as.character(normalization_mode))
  if (normalization_mode == "LOESS") normalization_mode <- "QC_LOESS"
  stopifnot(normalization_mode %in% c("NONE", "PQN", "QC_LOESS"))

  if (!exists("loess_min_qc_points", inherits = TRUE)) {
    loess_min_qc_points <<- 5
  }
  if (!exists("QC_LOESS_span", inherits = TRUE)) {
    if (exists("loess_span", inherits = TRUE)) {
      QC_LOESS_span <<- loess_span
    } else {
      QC_LOESS_span <<- 0.75
    }
  }
  if (!exists("injection_order_path", inherits = TRUE)) {
    injection_order_path <<- ""
  }
  stopifnot(is.numeric(loess_min_qc_points), length(loess_min_qc_points) == 1)
  stopifnot(loess_min_qc_points >= 5)
  stopifnot(is.numeric(QC_LOESS_span), length(QC_LOESS_span) == 1)
  stopifnot(QC_LOESS_span > 0, QC_LOESS_span <= 1)
  stopifnot(is.character(injection_order_path), length(injection_order_path) == 1)

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
    stopifnot(is.character(reference_path), length(reference_path) == 1)
    stopifnot(nzchar(trimws(reference_path)))
  }

  stopifnot(is.character(reference_col_metabolite), length(reference_col_metabolite) == 1)
  stopifnot(is.character(reference_col_ref_ion), length(reference_col_ref_ion) == 1)
  stopifnot(is.character(reference_col_mz), length(reference_col_mz) == 1)
  stopifnot(is.character(reference_col_rt), length(reference_col_rt) == 1)

  # ---------------------------------------------------------------------------
  # Name sanitation
  # ---------------------------------------------------------------------------
  stopifnot(sanitize_mode %in% c("none", "greek_latin_ascii", "ascii_translit"))

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
  if (!exists("volcano_style", inherits = TRUE)) {
    volcano_style <<- "classic"
  }
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
