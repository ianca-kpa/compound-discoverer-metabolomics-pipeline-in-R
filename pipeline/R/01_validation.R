# =============================================================================
# 01_validation.R
# Settings validation
# =============================================================================

apply_missing_setting_defaults <- function() {
  if (!exists("output_level", inherits = TRUE)) {
    legacy_minimal <- isTRUE(get0("minimal_output", ifnotfound = FALSE, inherits = TRUE))
    assign("output_level", if (legacy_minimal) "minimal" else "standard", envir = .GlobalEnv)
  }

  defaults <- list(
    cd_sheet = 1,
    metadata_sheet = 1,
    reference_sheet = 1,
    reference_path = "",
    reference_col_metabolite = "",
    reference_col_ref_ion = "",
    reference_col_mz = "",
    reference_col_rt = "",
    comparison_group_control = "WT",
    comparison_group_treatment = "TG",
    use_weight_normalization = FALSE,
    use_reference_file = FALSE,
    stop_on_invalid_weight = TRUE,
    invalid_weight_to_NA = TRUE,
    p_value_cutoff = 0.05,
    fdr_cutoff = 0.05,
    fc_cutoff_log2 = 0,
    run_metrics = "FDR_and_p_value",
    heatmap_rank_metrics = "FDR_and_p_value",
    statistical_test_type = "student",
    test_is_paired = FALSE,
    pvalue_correction_method = "FDR",
    comparison_mode = "pairwise",
    multigroup_groups = character(0),
    multigroup_test = "kruskal",
    multigroup_pairwise_mode = "selected",
    multigroup_pairwise_pairs = NULL,
    normalization_mode = "none",
    make_qc_diagnostics = FALSE,
    qc_loess_span = 0.75,
    qc_loess_min_qc_points = 4,
    pqn_min_qc_points = 3,
    rsd_filter_type = "QC_RSD",
    missing_exclusion_max_fraction = 0.5,
    presence_filter_min_fraction = 0,
    impute_half_min = TRUE,
    rsd_thresholds = c(20),
    active_variant = "none",
    low_variance_filter_method = "none",
    low_variance_filter_fraction = 0.2,
    low_variance_filter_rounding = "ceiling",
    log2_offset = 1,
    pca_scaling = "pareto",
    pca_label_samples = TRUE,
    ellipse_positive = TRUE,
    use_only_known = TRUE,
    duplicate_name_strategy = "collapse_best_qc_rsd",
    dup_mz_digits = 4,
    dup_rt_digits = 2,
    export_metaboanalyst_ready = TRUE,
    export_metaboanalyst_duplicate_only = TRUE,
    make_heatmap_by_model = TRUE,
    make_heatmap_by_model_sex = TRUE,
    heatmap_top_n = 80,
    heatmap_scale_method = "zscore",
    heatmap_order_samples_by_group = TRUE,
    heatmap_cluster_distance = "euclidean",
    heatmap_cluster_method = "ward.D2",
    heatmap_palette_n = 101,
    heatmap_breaks_symmetric = TRUE,
    heatmap_breaks_limit = 5,
    make_sig_heatmap_by_model = FALSE,
    make_sig_heatmap_by_model_sex = FALSE,
    make_sig_heatmap_FvsM_within_group = FALSE,
    sig_heatmap_max_features = 70,
    sig_heatmap_require_fc_cutoff = TRUE,
    save_stats_excel_per_model = TRUE,
    save_sig_metabolites_txt_per_model = TRUE,
    make_volcano_plots = TRUE,
    volcano_style = "classic",
    volcano_auto_axis = TRUE,
    volcano_axis_expand_mult = 0.08,
    volcano_add_labels = TRUE,
    volcano_label_number = Inf,
    volcano_custom_labels = NULL,
    volcano_add_cutoff_lines = TRUE,
    volcano_classic_point_size = 2.5,
    volcano_classic_point_shape = 21,
    volcano_classic_fills = c("#3B82F6", "#BDBDBD", "#EF4444"),
    volcano_classic_colors = c("#1D4ED8", "#7A7A7A", "#B91C1C"),
    volcano_classic_legend_title = "Regulation",
    sanitize_names_for_exports = FALSE,
    sanitize_mode = "greek_latin_ascii",
    strip_stereo_prefixes_for_names = TRUE,
    output_level = "standard",
    minimal_output = FALSE
  )
  for (key in names(defaults)) {
    if (!exists(key, inherits = TRUE)) {
      assign(key, defaults[[key]], envir = .GlobalEnv)
    }
  }

  if (!exists("alpha_sig", inherits = TRUE)) {
    assign("alpha_sig", get0("fdr_cutoff", ifnotfound = 0.05, inherits = TRUE), envir = .GlobalEnv)
  }

  invisible(TRUE)
}

normalize_rsd_filter_type <- function(rsd_filter_type) {
  rsd_filter_type_local <- toupper(trimws(as.character(rsd_filter_type)[1]))

  if (identical(rsd_filter_type_local, "QC_RSD")) {
    return("QC_RSD")
  }

  if (identical(rsd_filter_type_local, "RSD")) {
    return("RSD")
  }

  stop("rsd_filter_type must be 'QC_RSD' or 'RSD'.")
}

normalize_active_variant <- function(active_variant, rsd_thresholds, rsd_filter_type = "QC_RSD") {
  active_variant_local <- trimws(as.character(active_variant)[1])
  rsd_variant_prefix <- normalize_rsd_filter_type(rsd_filter_type)

  normalize_rsd_variant_name <- function(value) {
    value_local <- toupper(trimws(as.character(value)[1]))

    if (identical(value_local, rsd_variant_prefix)) {
      return(value_local)
    }

    if (identical(value_local, "QC_RSD") || identical(value_local, "RSD")) {
      return(rsd_variant_prefix)
    }

    value_local
  }

  if (!nzchar(active_variant_local) || identical(tolower(active_variant_local), "base")) {
    return("BASE")
  }

  if (identical(tolower(active_variant_local), "none")) {
    return("none")
  }

  active_variant_local <- normalize_rsd_variant_name(active_variant_local)

  if (identical(active_variant_local, rsd_variant_prefix)) {
    if (length(rsd_thresholds) == 1 && !is.na(rsd_thresholds[1])) {
      return(paste0(rsd_variant_prefix, rsd_thresholds[1]))
    }

    stop(
      "active_variant = '", rsd_variant_prefix, "' is ambiguous when multiple rsd_thresholds are configured. ",
      "Use one of: ", paste0(rsd_variant_prefix, rsd_thresholds, collapse = ", ")
    )
  }

  if (grepl("^(QC_RSD|RSD)[0-9]+$", toupper(active_variant_local))) {
    active_variant_threshold <- suppressWarnings(as.numeric(sub("^(QC_RSD|RSD)", "", toupper(active_variant_local))))

    if (length(rsd_thresholds) > 0 && !is.na(active_variant_threshold) &&
        !active_variant_threshold %in% rsd_thresholds) {
      stop(
        "active_variant = '", active_variant_local, "' is not available for the configured rsd_thresholds. ",
        "Use one of: ", paste0(rsd_variant_prefix, rsd_thresholds, collapse = ", ")
      )
    }

    return(paste0(rsd_variant_prefix, active_variant_threshold))
  }

  stop(
    "active_variant must be 'none', 'base', '", rsd_variant_prefix, "', or one of: ",
    paste0(rsd_variant_prefix, rsd_thresholds, collapse = ", ")
  )
}

validate_settings <- function() {
  apply_missing_setting_defaults()

  stopifnot(normalize_output_level(output_level) %in% c("minimal", "standard", "full_debug"))

  # ---------------------------------------------------------------------------
  # Normalization
  # ---------------------------------------------------------------------------
  normalization_mode_local <- normalize_normalization_mode(normalization_mode, default = NULL)
  stopifnot(normalization_mode_local %in% c("none", "weight", "qc_loess", "cyclic_loess", "qcrsc", "pqn_qc", "pqn_sample"))
  stopifnot(is.logical(make_qc_diagnostics), length(make_qc_diagnostics) == 1)
  qcrsc_spectral_cleaning_setting <- if (exists(".apply_qcrsc_spectral_cleaning_setting", inherits = TRUE)) {
    get(".apply_qcrsc_spectral_cleaning_setting", inherits = TRUE)
  } else if (exists("apply_qcrsc_spectral_cleaning", inherits = TRUE) &&
             !is.function(get("apply_qcrsc_spectral_cleaning", inherits = TRUE))) {
    get("apply_qcrsc_spectral_cleaning", inherits = TRUE)
  } else {
    TRUE
  }
  stopifnot(is.logical(qcrsc_spectral_cleaning_setting), length(qcrsc_spectral_cleaning_setting) == 1)
  stopifnot(is.numeric(qc_loess_span), qc_loess_span > 0, qc_loess_span <= 1)
  stopifnot(is.numeric(qc_loess_min_qc_points), qc_loess_min_qc_points >= 3)
  stopifnot(is.numeric(pqn_min_qc_points), pqn_min_qc_points >= 2)

  # ---------------------------------------------------------------------------
  # PCA
  # ---------------------------------------------------------------------------
  stopifnot(pca_scaling %in% c("none", "pareto", "autoscale"))

  # ---------------------------------------------------------------------------
  # Filters
  # ---------------------------------------------------------------------------
  stopifnot(low_variance_filter_method %in% c("none", "iqr"))
  stopifnot(low_variance_filter_rounding %in% c("floor", "ceiling", "round"))
  stopifnot(normalize_rsd_filter_type(rsd_filter_type) %in% c("QC_RSD", "RSD"))

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
  # Active variant
  # ---------------------------------------------------------------------------
  if (exists("active_variant", inherits = TRUE)) {
    normalize_active_variant(active_variant, rsd_thresholds, rsd_filter_type)
  }

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
  stopifnot(sanitize_mode %in% c("greek_latin_ascii", "ascii_translit"))

  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------
  stopifnot(all(run_metrics %in% c("FDR", "p_value", "FDR_and_p_value")))
  if (exists("heatmap_rank_metrics", inherits = TRUE)) {
    stopifnot(all(heatmap_rank_metrics %in% c("FDR", "p_value", "FDR_and_p_value")))
  }
  stopifnot(statistical_test_type %in% c("student", "welch", "wilcoxon", "limma"))
  stopifnot(is.logical(test_is_paired), length(test_is_paired) == 1, !is.na(test_is_paired))
  if (identical(statistical_test_type, "limma") && isTRUE(test_is_paired)) {
    stop(
      "statistical_test_type = 'limma' cannot be combined with test_is_paired = TRUE yet. ",
      "Paired limma requires explicit pair identifiers in metadata."
    )
  }
  stopifnot(comparison_mode %in% c("pairwise", "multigroup", "both"))
  stopifnot(multigroup_test %in% c("kruskal", "anova", "welch_anova"))
  stopifnot(multigroup_pairwise_mode %in% c("none", "all", "selected"))

  # ---------------------------------------------------------------------------
  # Heatmap
  # ---------------------------------------------------------------------------
  stopifnot(heatmap_scale_method %in% c("none", "zscore", "pareto"))
  stopifnot(heatmap_cluster_distance %in% c("euclidean", "manhattan"))
  stopifnot(heatmap_cluster_method %in% c("ward.D2", "complete", "average"))

  # ---------------------------------------------------------------------------
  # Volcano main style
  # ---------------------------------------------------------------------------
  stopifnot(volcano_style %in% c("classic"))

  # ---------------------------------------------------------------------------
  # Volcano classic legend
  # ---------------------------------------------------------------------------
  stopifnot(is.character(volcano_classic_legend_title))
  stopifnot(length(volcano_classic_legend_title) == 1)

  # classic colors must be exactly 3
  stopifnot(length(volcano_classic_fills) == 3)
  stopifnot(length(volcano_classic_colors) == 3)

  

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
