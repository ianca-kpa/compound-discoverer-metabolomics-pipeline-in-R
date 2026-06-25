# =============================================================================
# 12_main_pipeline.R
# Main pipeline runner
# =============================================================================

run_untargeted_pipeline <- function(debug_mode = FALSE) {
  pipeline_t0 <- Sys.time()
  runtime_profile_written <- FALSE
  runtime_profile <- NULL
  on.exit({
    if (exists("export_debug_outputs", inherits = FALSE) && isTRUE(export_debug_outputs) && !isTRUE(runtime_profile_written) &&
        exists("runtime_profile_exists", mode = "function") &&
        runtime_profile_exists() &&
        exists("runtime_profile_write", mode = "function")) {
      try(
        runtime_profile_write(file.path(output_dir, "global", "audits_global", "runtime_profile.csv")),
        silent = TRUE
      )
    }
  }, add = TRUE)

  legacy_minimal_output <- isTRUE(get0("minimal_output", ifnotfound = FALSE, inherits = TRUE))
  output_level_local <- normalize_output_level(
    get0("output_level", ifnotfound = NULL, inherits = TRUE),
    legacy_minimal = legacy_minimal_output
  )
  if (isTRUE(debug_mode)) output_level_local <- "full_debug"
  output_flags <- derive_output_flags(output_level_local)
  export_debug_outputs <- output_flags$export_debug_outputs
  export_intermediate_tables <- output_flags$export_intermediate_tables
  export_all_plots <- output_flags$export_all_plots
  export_normalization_audit <- output_flags$export_normalization_audit
  export_qc_summary <- output_flags$export_qc_summary
  export_multigroup_outputs <- output_flags$export_multigroup_outputs
  export_all_pairwise_multigroup <- output_flags$export_all_pairwise_multigroup
  debug_mode <- export_debug_outputs
  minimal_output_enabled <- identical(output_level_local, "minimal")
  comparison_group_control <- get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE)
  comparison_group_treatment <- get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)
  active_variant_requested <- as.character(get0("active_variant", ifnotfound = "none", inherits = TRUE))[1]
  active_variant_requested <- trimws(active_variant_requested)
  if (!nzchar(active_variant_requested) || identical(tolower(active_variant_requested), "base")) {
    active_variant_requested <- "none"
  }
  if (identical(tolower(active_variant_requested), "none")) {
    active_variant_requested <- "BASE"
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  runtime_profile_reset(output_dir)
  log_path <- file.path(output_dir, "PIPELINE_LOG.txt")
  output_readme_path <- file.path(output_dir, "README.txt")
  writeLines(
    c(
      "Untargeted metabolomics pipeline outputs",
      "",
      paste0("Output level: ", output_level_local),
      "Calculations are identical across output levels; only exported files differ.",
      "Minimal: final statistics Excel, principal/combined PCA, primary pairwise volcano plots, README, and pipeline log.",
      "Standard: Minimal plus main heatmaps, summarized QC/normalization audits, and enabled multi-group global/PCA/top-heatmap outputs.",
      "Full / Debug: all final, intermediate, technical, WITH_QC, BIO_ONLY, matrix, audit, plot, CSV, and XLSX artifacts."
    ),
    con = output_readme_path,
    useBytes = TRUE
  )

  with_console_capture_to_file(log_path, {
    message(console_rule())
    message("============================================================")
    message("PIPELINE LOG — Untargeted Metabolomics - Data Processing")
    step_info("Date/time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    message(console_rule())
    message("[INPUTS]")
    step_info("cd_file_path: ", cd_file_path)
    step_info("cd_sheet: ", cd_sheet)
    step_info("metadata_path: ", metadata_path)
    step_info("metadata_sheet: ", metadata_sheet)
    step_info("reference_path: ", if (exists("reference_path", inherits = TRUE) && !is.null(reference_path) && nzchar(as.character(reference_path))) reference_path else "not provided")
    step_info("output_dir: ", output_dir)
    step_info("output_level: ", output_level_local)
    step_info("use_only_known: ", use_only_known)
    step_info("duplicate_name_strategy: ", duplicate_name_strategy)
    step_info("active_variant: ", active_variant_requested)
    # -------------------------------------------------------------------------
    # STEP 02 — Reading input files
    # -------------------------------------------------------------------------
    t_step <- step_start(2, 15, "Reading input files")

    # message("Step 3: Reading input files...")
    cd_raw <- read_any_table(cd_file_path, cd_sheet)
    metadata <- read_any_table(metadata_path, metadata_sheet)
    reference_tbl <- NULL
    use_reference_file_local <- if (exists("use_reference_file", inherits = TRUE)) {
      isTRUE(use_reference_file)
    } else {
      exists("reference_path", inherits = TRUE) &&
        !is.null(reference_path) &&
        nzchar(as.character(reference_path))
    }

    reference_col_overrides <- list(
      metabolite = get_optional_setting_value("reference_col_metabolite", env = .GlobalEnv),
      ref_ion = get_optional_setting_value("reference_col_ref_ion", env = .GlobalEnv),
      mz = get_optional_setting_value("reference_col_mz", env = .GlobalEnv),
      rt = get_optional_setting_value("reference_col_rt", env = .GlobalEnv)
    )

    if (isTRUE(use_reference_file_local) &&
        exists("reference_path", inherits = TRUE) &&
        !is.null(reference_path) &&
        nzchar(as.character(reference_path))) {
      reference_sheet_local <- if (exists("reference_sheet", inherits = TRUE)) reference_sheet else 1
      reference_tbl <- read_any_table(reference_path, reference_sheet_local)
      step_info("Reference table: rows=", nrow(reference_tbl), " cols=", ncol(reference_tbl))
      step_info(
        "Reference column overrides (metabolite/ref_ion/mz/rt): ",
        paste(vapply(reference_col_overrides, function(x) if (is.null(x)) "auto" else x, character(1)), collapse = " / ")
      )
    } else if (!isTRUE(use_reference_file_local)) {
      step_info("Reference table: disabled by use_reference_file = FALSE")
    } else {
      step_info("Reference table: not provided (reference_path empty)")
    }

    step_info("Compound Discoverer table: rows=", nrow(cd_raw), " cols=", ncol(cd_raw))
    step_info("Metadata table: rows=", nrow(metadata), " cols=", ncol(metadata))
    step_ok("Reading input files", t_step)

    # -------------------------------------------------------------------------
    # STEP 03 — Standardizing metadata
    # -------------------------------------------------------------------------
    t_step <- step_start(3, 15, "Standardizing metadata")

    metadata <- clean_metadata(metadata)
    metadata_aligned <- metadata

    step_info("Metadata after cleaning: rows=", nrow(metadata), " cols=", ncol(metadata))
    step_ok("Standardizing metadata", t_step)

    # -------------------------------------------------------------------------
    # STEP 04 — Creating per-model output directory structure
    # -------------------------------------------------------------------------
    t_step <- step_start(4, 15, "Creating per-model output directory structure")

    models_detected <- get_models_from_metadata(metadata_aligned)
    step_info("Models detected: ", paste(models_detected, collapse = ", "))

    paths <- setup_output_dirs(output_dir, model_names = models_detected)
    runtime_profile_set_output_dir(output_dir)
    filter_summary <- tibble()
    summary_csv <- file.path(paths$global$audits, "filter_summary.csv")
    filter_summary_csv <- if (debug_mode) summary_csv else NULL

    step_ok("Creating per-model output directory structure", t_step)

    # -------------------------------------------------------------------------
    # STEP 05 — Building feature table
    # -------------------------------------------------------------------------
    t_step <- step_start(5, 15, "Building feature table")

    feature_tbl <- build_feature_table(
      cd_raw,
      sanitize_names_for_exports = sanitize_names_for_exports,
      sanitize_mode = sanitize_mode,
      mz_digits = dup_mz_digits,
      rt_digits = dup_rt_digits,
      remove_greek_letters_from_names = isTRUE(get0("strip_stereo_prefixes_for_names", ifnotfound = TRUE, inherits = TRUE)),
      remove_isomer_descriptors_from_names = isTRUE(get0("strip_stereo_prefixes_for_names", ifnotfound = TRUE, inherits = TRUE))
    )

    if (debug_mode) {
      out_feat_map <- file.path(paths$global$exports, "02_featureID_to_display_name_map.csv")
      write_csv_safe(
        feature_tbl %>% select(featureID, display_name, Name, Name_canon, dplyr::any_of("Ref ion"), Formula, mz, RT),
        out_feat_map
      )
      log_written_object(log_path, out_feat_map, feature_tbl, note = "feature table")
    }

    step_info("Feature table created: rows=", nrow(feature_tbl), " cols=", ncol(feature_tbl))
    step_ok("Building feature table", t_step)

    # -------------------------------------------------------------------------
    # STEP 06 — Building raw assay
    # -------------------------------------------------------------------------
    t_step <- step_start(6, 15, "Building raw assay")

    # Determine if QC samples are required based on user settings and chosen steps.
    raw_active_variant <- get0("active_variant", ifnotfound = NA, inherits = TRUE)
    normalization_mode_setting <- normalize_normalization_mode(
      get0("normalization_mode", ifnotfound = "qcrsc", inherits = TRUE)
    )
    rsd_thr <- get0("rsd_thresholds", ifnotfound = numeric(0), inherits = TRUE)
    needs_qc_for_normalization <- normalization_mode_setting %in% c("qc_loess", "qcrsc", "pqn_qc")
    requires_injection_order_for_drift <- normalization_mode_setting %in% c("qc_loess", "qcrsc")
    needs_qc_for_diagnostics <- isTRUE(get0("make_qc_diagnostics", ifnotfound = FALSE, inherits = TRUE))
    needs_qc_for_variants <- length(rsd_thr) > 0
    needs_qc_for_duplicates <- duplicate_name_strategy %in% c("reference_or_best_qc_rsd", "collapse_best_qc_rsd")
    require_qc <- needs_qc_for_normalization || needs_qc_for_diagnostics || needs_qc_for_variants || needs_qc_for_duplicates

    assay_bundle <- build_assay_from_cd(
      cd_raw,
      feature_tbl,
      metadata_aligned,
      paths,
      require_qc = require_qc,
      require_injection_order = requires_injection_order_for_drift,
      export_intermediate_tables = export_intermediate_tables
    )

    if (isTRUE(export_debug_outputs)) {
      step_info("DEBUG — assay_obj class: ", paste(class(assay_bundle), collapse = ", "))
      step_info("DEBUG — assay_obj names: ", paste(names(assay_bundle), collapse = ", "))
    }
    assay_num_raw <- assay_bundle$assay_num_raw
    metadata_aligned <- assay_bundle$metadata_aligned
    qc_idx <- assay_bundle$qc_idx
    sample_idx <- assay_bundle$sample_idx
    drift_injection_order <- assay_bundle$drift_injection_order
    rsd_filter_type <- normalize_rsd_filter_type(get0("rsd_filter_type", ifnotfound = "QC_RSD", inherits = TRUE))

    step_info("Raw assay created: rows=", nrow(assay_num_raw), " cols=", ncol(assay_num_raw))
    step_info("Sample index length: ", length(sample_idx))
    step_info("QC index length: ", length(qc_idx))
    step_info("Injection order source: ", assay_bundle$drift_injection_order_source)

    step_ok("Building raw assay", t_step)

    # -------------------------------------------------------------------------
    # STEP 07 — Selecting preprocessing scenario
    # -------------------------------------------------------------------------
    t_step <- step_start(7, 15, "Selecting preprocessing scenario")

    normalization_mode_local <- normalize_normalization_mode(
      get0("normalization_mode", ifnotfound = "qcrsc", inherits = TRUE)
    )
    if (!normalization_mode_local %in% c("none", "weight", "qc_loess", "cyclic_loess", "qcrsc", "pqn_qc", "pqn_sample")) {
      stop("Unsupported normalization_mode: ", normalization_mode_local, ". Use 'none', 'weight', 'qc_loess', 'cyclic_loess', 'qcrsc', 'pqn_qc', or 'pqn_sample'.")
    }

    drift_correction_bundle <- NULL
    pqn_bundle <- NULL
    assay_num_weight <- NULL
    assay_num_normalized_technical <- NULL
    assay_num_prefilter <- assay_num_raw

    step_info("Normalization mode selected: ", normalization_mode_local)
    step_info("Feature filters will run before normalization.")
    if (identical(normalization_mode_local, "qc_loess")) {
      step_info("Scenario: QC-RSD/IQR -> QC-LOESS correction with QCs -> remove QC -> weight normalization on biological samples -> log2.")
    } else if (identical(normalization_mode_local, "cyclic_loess")) {
      step_info("Scenario: filters -> cyclic LOESS normalization matching the tutorial/limma method -> remove QC -> log2.")
      if (isTRUE(use_weight_normalization)) {
        step_info("Weight normalization is ignored for normalization_mode = 'cyclic_loess' to avoid double scale normalization.")
      }
    } else if (identical(normalization_mode_local, "qcrsc")) {
      step_info("Scenario: QC-RSD/IQR -> QC-RSC robust spline correction with QCs -> remove QC -> weight normalization on biological samples -> log2.")
    } else if (identical(normalization_mode_local, "weight")) {
      step_info("Scenario: filters -> weight normalization only -> log2.")
    } else if (identical(normalization_mode_local, "pqn_qc")) {
      step_info("Scenario: QC-RSD/IQR -> PQN using pooled QC reference -> log2.")
      if (isTRUE(use_weight_normalization)) {
        step_info("Weight normalization is ignored for normalization_mode = 'pqn_qc'.")
      }
    } else if (identical(normalization_mode_local, "pqn_sample")) {
      step_info("Scenario: filters -> PQN using biological sample median reference -> remove QC -> log2.")
      if (isTRUE(use_weight_normalization)) {
        step_info("Weight normalization is ignored for normalization_mode = 'pqn_sample'.")
      }
    } else {
      step_info("Scenario: QC-RSD/IQR -> no QC-based normalization -> log2.")
    }
    step_ok("Selecting preprocessing scenario", t_step)

    # -------------------------------------------------------------------------
    # Initialize filter summary table
    # -------------------------------------------------------------------------
    filter_summary <- tibble::tibble(
      step = character(),
      n_features_before = integer(),
      n_features_after = integer(),
      n_removed = integer(),
      pct_removed = numeric()
    )

    # -------------------------------------------------------------------------
    # STEP 09 — Applying missing-value exclusion
    # -------------------------------------------------------------------------
    t_step <- step_start(9, 15, "Applying missing-value exclusion")

    n0 <- ncol(assay_num_prefilter)

    miss_bundle <- filter_missing_exclusion(
      assay_num_prefilter,
      feature_tbl,
      sample_idx,
      max_missing_fraction = missing_exclusion_max_fraction,
      audit_path = if (debug_mode) file.path(paths$global$audits, "missing_exclusion_audit.csv") else NULL
    )

    assay_work <- miss_bundle$assay
    feature_tbl_work <- miss_bundle$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "missing_exclusion",
      n0,
      ncol(assay_work),
      out_csv = filter_summary_csv
    )

    step_info("Features before filter: ", n0)
    step_info("Features after filter: ", ncol(assay_work))
    step_ok("Applying missing-value exclusion", t_step)

    # -------------------------------------------------------------------------
    # STEP 10 — Applying presence filter and imputation
    # -------------------------------------------------------------------------
    t_step <- step_start(10, 15, "Applying presence filter and imputation")

    n1 <- ncol(assay_work)

    miss2 <- presence_filter_and_impute(
      assay_work,
      feature_tbl_work,
      sample_idx,
      min_fraction = presence_filter_min_fraction,
      impute_half_min = impute_half_min,
      audit_path = if (debug_mode) file.path(paths$global$audits, "presence_filter_audit.csv") else NULL
    )

    assay_work <- miss2$assay
    feature_tbl_work <- miss2$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "presence_filter",
      n1,
      ncol(assay_work),
      out_csv = filter_summary_csv
    )

    step_info("Presence filter min fraction: ", presence_filter_min_fraction)
    step_info("Features before filter: ", n1)
    step_info("Features after filter: ", ncol(assay_work))
    step_ok("Applying presence filter and imputation", t_step)

    # -------------------------------------------------------------------------
    # STEP 11 — Applying feature-level filters, normalization, and duplicate handling
    # -------------------------------------------------------------------------
    t_step <- step_start(11, 15, "Applying feature-level filters, normalization, and duplicate handling")

    n2 <- ncol(assay_work)

    known <- filter_known(
      assay_work,
      feature_tbl_work,
      use_only_known,
      audit_path = if (debug_mode) file.path(paths$global$audits, "known_filter_audit.csv") else NULL
    )

    assay_work <- known$assay
    feature_tbl_work <- known$feature

    # -------------------------------------------------------------------------
    # Duplicate handling: apply duplicate collapsing right after Known-only filter
    # -------------------------------------------------------------------------
    step_info("Applying duplicate name handling (post-known filter)...")

    qc_rsd_for_duplicates <- NULL
    if (duplicate_name_strategy %in% c("reference_or_best_qc_rsd", "collapse_best_qc_rsd")) {
      qc_rsd_for_duplicates <- calc_qc_rsd(assay_work, qc_idx, rsd_filter_type = rsd_filter_type)
    }

    dup_out_pre <- collapse_duplicate_names(
      mat = assay_work,
      feature_tbl = feature_tbl_work,
      strategy = duplicate_name_strategy,
      reference_tbl = reference_tbl,
      reference_col_overrides = reference_col_overrides,
      sanitize_mode = sanitize_mode,
      qc_rsd = qc_rsd_for_duplicates,
      audit_path = if (isTRUE(export_debug_outputs)) {
        file.path(
          paths$global$audits,
          paste0("duplicate_name_audit_", duplicate_name_strategy, "_post_Known.csv")
        )
      } else NULL
    )

    assay_work <- dup_out_pre$mat
    feature_tbl_work <- dup_out_pre$feature

    metaboanalyst_export_mat <- assay_work
    metaboanalyst_export_feature <- feature_tbl_work
    metaboanalyst_duplicate_only_mat <- assay_work
    metaboanalyst_duplicate_only_feature <- feature_tbl_work

    filter_summary <- append_filter_summary(
      filter_summary,
      "known_filter",
      n2,
      ncol(assay_work),
      out_csv = filter_summary_csv
    )

    rsd_variant_prefix <- rsd_filter_type
    variants <- list()
    variants$BASE <- list(mat = assay_work, feature = feature_tbl_work)

    if (length(qc_idx) >= 2 && length(rsd_thresholds) > 0) {
      qc_rsd_all <- calc_qc_rsd(assay_work, qc_idx, rsd_filter_type = rsd_filter_type)

      df_qc_rsd_pre <- tibble(
        featureID = names(qc_rsd_all),
        qc_rsd = as.numeric(qc_rsd_all)
      ) %>% arrange(qc_rsd)

      if (debug_mode) {
        out_qc_rsd_pre <- file.path(paths$global$audits, paste0(tolower(rsd_variant_prefix), "_values_pre_variants.csv"))
        write_csv_safe(df_qc_rsd_pre, out_qc_rsd_pre)
        log_written_object(log_path, out_qc_rsd_pre, df_qc_rsd_pre, note = paste0(rsd_variant_prefix, " values (pre-variants)"))
      }

      for (thr in rsd_thresholds) {
        keep <- names(qc_rsd_all)[!is.na(qc_rsd_all) & qc_rsd_all <= thr]
        variants[[paste0(rsd_variant_prefix, thr)]] <- list(
          mat = assay_work[, keep, drop = FALSE],
          feature = feature_tbl_work %>% filter(featureID %in% keep)
        )
      }

      step_info("Variants created based on ", rsd_variant_prefix, " thresholds: ", paste(names(variants), collapse = ", "))

      active_variant_effective <- normalize_active_variant(
        active_variant_requested,
        get0("rsd_thresholds", ifnotfound = numeric(0), inherits = TRUE),
        rsd_filter_type
      )

      if (!active_variant_effective %in% names(variants)) {
        stop(
          "active_variant not found: ", active_variant_effective,
          "\nAvailable: ", paste(names(variants), collapse = ", ")
        )
      }
    } else {
      if (isTRUE(get0("make_qc_diagnostics", ifnotfound = FALSE, inherits = TRUE))) {
        step_info("QC diagnostics skipped because fewer than 2 QC samples are available.")
      }
      if (length(rsd_thresholds) > 0) {
        step_info("QC-based RSD variants skipped because fewer than 2 QC samples are available.")
      }
      active_variant_effective <- "BASE"
    }

    variants$ACTIVE <- variants[[active_variant_effective]]
    step_info("Active variant selected: ", active_variant_effective, " | n_features=", ncol(variants$ACTIVE$mat))

    if (debug_mode) {
      out_post_rsd_mat <- file.path(
        paths$global$exports,
        paste0("10_MATRIX_post_", active_variant_effective, "_postRSD_preLowVar_preDup_ALL.csv")
      )
      dir.create(dirname(out_post_rsd_mat), recursive = TRUE, showWarnings = FALSE)
      write.csv(variants$ACTIVE$mat, file = out_post_rsd_mat, row.names = TRUE)
      log_written_object(
        log_path,
        out_post_rsd_mat,
        variants$ACTIVE$mat,
        note = "post-RSD, pre-lowVar, pre-dup"
      )

      df_post_rsd_named <- as.data.frame(variants$ACTIVE$mat, check.names = FALSE) %>%
        tibble::rownames_to_column("sample")

      df_post_rsd_named <- rename_feature_cols(
        df_post_rsd_named,
        variants$ACTIVE$feature,
        sample_col = "sample"
      )

      out_post_rsd_named <- file.path(
        paths$global$exports,
        paste0("10_TABLE_post_", active_variant_effective, "_postRSD_preLowVar_preDup_ALL_NAMED.csv")
      )
      write_csv_safe(df_post_rsd_named, out_post_rsd_named)
      log_written_object(
        log_path,
        out_post_rsd_named,
        df_post_rsd_named,
        note = "post-RSD named table"
      )
    }

    step_info("Applying feature-level filters, normalization, and duplicate handling", t_step)

    n3 <- ncol(variants$ACTIVE$mat)

    lv <- filter_low_variance_deterministic(
      variants$ACTIVE$mat,
      variants$ACTIVE$feature,
      method = low_variance_filter_method,
      frac = low_variance_filter_fraction,
      rounding = low_variance_filter_rounding,
      sample_idx = sample_idx,
      audit_path = if (debug_mode) file.path(paths$global$audits, "low_variance_iqr_audit_ACTIVE.csv") else NULL
    )

    variants$ACTIVE$mat <- lv$assay
    variants$ACTIVE$feature <- lv$feature
    metaboanalyst_compatible_iqr_mat <- variants$ACTIVE$mat
    metaboanalyst_compatible_iqr_feature <- variants$ACTIVE$feature

    step_info(
      "IQR low-variance filter: fraction=",
      sprintf("%.2f", lv$fraction),
      " | rounding=", lv$rounding,
      " | removed=", lv$n_removed,
      " | kept=", lv$n_after
    )

    filter_summary <- append_filter_summary(
      filter_summary,
      "low_variance_iqr_ACTIVE",
      n3,
      ncol(variants$ACTIVE$mat),
      out_csv = filter_summary_csv
    )

    # Duplicate handling was already applied post-known filter.

    mat_active_filtered <- variants$ACTIVE$mat
    mat_active_normalized <- mat_active_filtered
    # Preserve the pre-normalization matrix for global pre/post comparisons
    mat_active_base <- mat_active_filtered

    run_loess_comparison_diagnostics <- isTRUE(export_all_plots) &&
      isTRUE(make_qc_diagnostics) &&
      normalization_mode_local %in% c("qc_loess", "qcrsc", "cyclic_loess")

    if (isTRUE(run_loess_comparison_diagnostics)) {
      profile_expr("Normalization comparison diagnostics", {
        normalization_compare_dir <- file.path(paths$global$root, "plots_global", "normalization")
        normalization_compare_png <- file.path(normalization_compare_dir, "normalization_comparison_boxplots.png")
        tryCatch(
          plot_normalization_comparison(
            assay_num_base = assay_num_raw,
            metadata_aligned = metadata_aligned,
            sample_idx = sample_idx,
            qc_idx = qc_idx,
            out_png = normalization_compare_png,
            qc_loess_span = qc_loess_span,
            qc_loess_min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order,
            log2_offset = log2_offset
          ),
          error = function(e) {
            message("  - Normalization comparison plot skipped: ", conditionMessage(e))
            NULL
          }
        )

        tutorial_style_png <- file.path(normalization_compare_dir, "normalization_comparison_tutorial_style_raw_boxplots.png")
        tryCatch(
          plot_tutorial_style_normalization_comparison(
            assay_num_raw = assay_num_raw,
            metadata_aligned = metadata_aligned,
            sample_idx = sample_idx,
            qc_idx = qc_idx,
            out_png = tutorial_style_png,
            qc_loess_span = qc_loess_span,
            qc_loess_min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order,
            log2_offset = log2_offset,
            audit_path = file.path(paths$global$audits, "normalization_comparison_tutorial_style_raw_audit.csv"),
            sample_order_path = file.path(paths$global$audits, "normalization_comparison_tutorial_style_sample_order_used.csv")
          ),
          error = function(e) {
            message("  - Tutorial-style normalization comparison plot skipped: ", conditionMessage(e))
            NULL
          }
        )

        qc_loess_weight_png <- file.path(normalization_compare_dir, "qc_loess_weight_comparison.png")
        tryCatch(
          plot_qc_loess_weight_comparison(
            assay_num_base = assay_num_raw,
            metadata_aligned = metadata_aligned,
            sample_idx = sample_idx,
            qc_idx = qc_idx,
            out_png = qc_loess_weight_png,
            qc_loess_span = qc_loess_span,
            qc_loess_min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order,
            log2_offset = log2_offset
          ),
          error = function(e) {
            message("  - QC-LOESS weight comparison plot skipped: ", conditionMessage(e))
            NULL
          }
        )

        qc_qcrsc_weight_png <- file.path(normalization_compare_dir, "qc_qcrsc_weight_comparison.png")
        tryCatch(
          plot_qc_qcrsc_weight_comparison(
            assay_num_base = assay_num_raw,
            metadata_aligned = metadata_aligned,
            sample_idx = sample_idx,
            qc_idx = qc_idx,
            out_png = qc_qcrsc_weight_png,
            min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order,
            log2_offset = log2_offset
          ),
          error = function(e) {
            message("  - QC-RSC weight comparison plot skipped: ", conditionMessage(e))
            NULL
          }
        )
      }, category = "diagnostic")
    } else {
      if (!isTRUE(export_all_plots)) {
        step_info("Technical normalization comparison plots require output_level = 'full_debug'.")
      } else if (!isTRUE(make_qc_diagnostics)) {
        step_info("Normalization comparison plots skipped because make_qc_diagnostics is FALSE.")
      } else if (!normalization_mode_local %in% c("qc_loess", "qcrsc", "cyclic_loess")) {
        step_info("LOESS/QC-RSC normalization comparison plots skipped for normalization_mode = '", normalization_mode_local, "'.")
      } else {
        step_info("Normalization comparison plots skipped because make_qc_diagnostics is FALSE.")
      }
    }

    if (normalization_mode_local %in% c("qc_loess", "qcrsc")) {
      correction_label <- if (identical(normalization_mode_local, "qcrsc")) "QC_RSC" else "QC_LOESS"
      correction_label_pretty <- if (identical(normalization_mode_local, "qcrsc")) "QC-RSC" else "QC-LOESS"
      if (identical(normalization_mode_local, "qcrsc")) {
        drift_correction_bundle <- profile_expr("QC-RSC normalization", {
          normalize_qcrsc_qc_ref(
            mat_active_filtered,
            qc_idx,
            min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order,
            batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL
          )
        }, category = "normalization")
        mat_active_normalized <- drift_correction_bundle$assay_num_qcrsc
        correction_tbl <- drift_correction_bundle$qcrsc_tbl

        qcrsc_spectral_cleaning_enabled <- isTRUE(get0(
          ".apply_qcrsc_spectral_cleaning_setting",
          ifnotfound = TRUE,
          inherits = TRUE
        ))

        if (isTRUE(qcrsc_spectral_cleaning_enabled)) {
          qcrsc_cleaning <- profile_expr("QC-RSC spectral cleaning", {
            apply_qcrsc_spectral_cleaning(
              mat_active_normalized,
              metadata_aligned = metadata_aligned,
              qc_idx = qc_idx,
              batch = if ("batch" %in% names(metadata_aligned)) metadata_aligned$batch else NULL,
              kruskal_p_cutoff = 1e-4,
              wilcoxon_p_cutoff = 1e-14,
              max_qc_rsd = 20
            )
          }, category = "normalization")

          before_cleaning_features <- ncol(mat_active_normalized)
          mat_active_normalized <- qcrsc_cleaning$assay
          keep_ids <- qcrsc_cleaning$keep_feature_ids

          mat_active_filtered <- mat_active_filtered[, keep_ids, drop = FALSE]
          mat_active_base <- mat_active_base[, keep_ids, drop = FALSE]
          feature_tbl_work <- feature_tbl_work %>% dplyr::filter(featureID %in% keep_ids)

          variants$ACTIVE$mat <- variants$ACTIVE$mat[, keep_ids, drop = FALSE]
          variants$ACTIVE$feature <- variants$ACTIVE$feature %>% dplyr::filter(featureID %in% keep_ids)

          metaboanalyst_export_mat <- metaboanalyst_export_mat[, keep_ids, drop = FALSE]
          metaboanalyst_export_feature <- metaboanalyst_export_feature %>% dplyr::filter(featureID %in% keep_ids)
          metaboanalyst_duplicate_only_mat <- metaboanalyst_duplicate_only_mat[, keep_ids, drop = FALSE]
          metaboanalyst_duplicate_only_feature <- metaboanalyst_duplicate_only_feature %>% dplyr::filter(featureID %in% keep_ids)

          filter_summary <- append_filter_summary(
            filter_summary,
            "qcrsc_spectral_cleaning",
            before_cleaning_features,
            ncol(mat_active_normalized),
            out_csv = filter_summary_csv
          )

          if (debug_mode) {
            out_qcrsc_cleaning <- file.path(paths$global$audits, "qcrsc_spectral_cleaning_audit.csv")
            write_csv_safe(qcrsc_cleaning$audit, out_qcrsc_cleaning)
            log_written_object(log_path, out_qcrsc_cleaning, qcrsc_cleaning$audit, note = "QC-RSC spectral cleaning audit")
          }

          step_info(
            "QC-RSC spectral cleaning removed ",
            before_cleaning_features - ncol(qcrsc_cleaning$assay),
            " features (Kruskal-Wallis, Wilcoxon, QC RSD)."
          )
        } else {
          step_info("QC-RSC spectral cleaning disabled by apply_qcrsc_spectral_cleaning <- FALSE.")
        }
      } else {
        drift_correction_bundle <- profile_expr("QC-LOESS normalization", {
          normalize_qc_loess_ref(
            mat_active_filtered,
            qc_idx,
            qc_loess_span = qc_loess_span,
            min_qc_points = qc_loess_min_qc_points,
            injection_order = drift_injection_order
          )
        }, category = "normalization")
        mat_active_normalized <- drift_correction_bundle$assay_num_qc_loess
        correction_tbl <- drift_correction_bundle$qc_loess_tbl
      }

      if (debug_mode) {
        out_correction_audit <- file.path(paths$global$exports, paste0("11_", tolower(correction_label), "_qc_correction_postFilters.csv"))
        write_csv_safe(correction_tbl, out_correction_audit)
        log_written_object(log_path, out_correction_audit, correction_tbl, note = paste(correction_label_pretty, "QC correction table after filters"))

        out_loess_mat <- file.path(paths$global$exports, paste0("11_MATRIX_ACTIVE_postFilters_post", correction_label, "_WITH_QC.csv"))
        write.csv(mat_active_normalized, file = out_loess_mat, row.names = TRUE)
        log_written_object(log_path, out_loess_mat, mat_active_normalized, note = paste("technical ACTIVE matrix after", correction_label_pretty, "with QC"))
      }

      if (isTRUE(export_qc_summary)) {
        profile_expr("QC correction audit", {
          qc_pca_dir <- file.path(paths$global$root, "plots_global", "pca")
          if (isTRUE(export_all_plots)) {
            plot_qc_loess_audit_pca(
              assay_num_pre = mat_active_filtered,
              assay_num_post = mat_active_normalized,
              metadata_aligned = metadata_aligned,
              qc_idx = qc_idx,
              out_dir = qc_pca_dir,
              pca_scaling = pca_scaling,
              log2_offset = log2_offset,
              pca_label_samples = FALSE,
              injection_order = drift_injection_order,
              correction_label = correction_label,
              log_path = log_path
            )
          }
          plot_qc_loess_audit_metrics(
            assay_num_pre = mat_active_filtered,
            assay_num_post = mat_active_normalized,
            metadata_aligned = metadata_aligned,
            qc_idx = qc_idx,
            out_dir = qc_pca_dir,
            audit_dir = paths$global$audits,
            log2_offset = log2_offset,
            injection_order = drift_injection_order,
            correction_label = correction_label,
            log_path = log_path,
            export_plots = export_all_plots,
            export_detailed_table = export_debug_outputs
          )
        }, category = "diagnostic")
      } else {
        step_info("QC summary skipped at output_level = 'minimal'.")
      }

      assay_num_weight <- NULL
      step_info(correction_label_pretty, " QC correction applied before weight normalization.")
    } else if (identical(normalization_mode_local, "weight")) {
      mat_active_normalized <- profile_expr("Weight normalization", {
        normalize_by_weight(
          mat_active_filtered,
          metadata_aligned,
          sample_idx,
          use_weight_normalization = TRUE,
          stop_on_invalid_weight   = stop_on_invalid_weight,
          invalid_weight_to_NA     = invalid_weight_to_NA
        )
      }, category = "normalization")
      assay_num_weight <- mat_active_normalized
      metaboanalyst_export_mat <- normalize_by_weight(
        metaboanalyst_export_mat,
        metadata_aligned,
        sample_idx,
        use_weight_normalization = TRUE,
        stop_on_invalid_weight   = stop_on_invalid_weight,
        invalid_weight_to_NA     = invalid_weight_to_NA
      )

      if (debug_mode) {
        out_weight <- file.path(paths$global$exports, "11_MATRIX_ACTIVE_postFilters_postWeight_ONLY_ALL.csv")
        write.csv(assay_num_weight, file = out_weight, row.names = TRUE)
        log_written_object(log_path, out_weight, assay_num_weight, note = "ACTIVE matrix after filters and weight normalization only")
      }

      step_info("Weight normalization only applied; QC-based normalization skipped.")
    } else if (identical(normalization_mode_local, "pqn_qc")) {
      pqn_bundle <- profile_expr("PQN QC normalization", {
        normalize_pqn_qc_ref(
          mat_active_filtered,
          qc_idx,
          min_qc_points = pqn_min_qc_points
        )
      }, category = "normalization")
      mat_active_normalized <- pqn_bundle$assay_num_pqn

      if (debug_mode) {
        out_pqn <- file.path(paths$global$exports, "11_pqn_qc_normalization_factors_postFilters.csv")
        write_csv_safe(pqn_bundle$pqn_tbl, out_pqn)
        log_written_object(log_path, out_pqn, pqn_bundle$pqn_tbl, note = "PQN QC normalization factors after filters")

        out_pqn_ref <- file.path(paths$global$exports, "11_pqn_qc_reference_spectrum_postFilters.csv")
        write_csv_safe(pqn_bundle$pqn_reference_tbl, out_pqn_ref)
        log_written_object(log_path, out_pqn_ref, pqn_bundle$pqn_reference_tbl, note = "PQN QC reference spectrum after filters")
      }

      assay_num_weight <- mat_active_filtered
      step_info("PQN QC normalization applied after QC-RSD/IQR filters.")
    } else if (identical(normalization_mode_local, "pqn_sample")) {
      pqn_bundle <- profile_expr("PQN sample normalization", {
        normalize_pqn_sample_ref(
          mat_active_filtered,
          sample_idx,
          min_sample_points = 2
        )
      }, category = "normalization")
      mat_active_normalized <- pqn_bundle$assay_num_pqn

      if (debug_mode) {
        out_pqn <- file.path(paths$global$exports, "11_pqn_sample_normalization_factors_postFilters.csv")
        write_csv_safe(pqn_bundle$pqn_tbl, out_pqn)
        log_written_object(log_path, out_pqn, pqn_bundle$pqn_tbl, note = "PQN sample normalization factors after filters")

        out_pqn_ref <- file.path(paths$global$exports, "11_pqn_sample_reference_spectrum_postFilters.csv")
        write_csv_safe(pqn_bundle$pqn_reference_tbl, out_pqn_ref)
        log_written_object(log_path, out_pqn_ref, pqn_bundle$pqn_reference_tbl, note = "PQN sample reference spectrum after filters")
      }

      assay_num_weight <- mat_active_filtered
      step_info("PQN sample normalization applied after filters without requiring QC samples.")
    } else if (identical(normalization_mode_local, "cyclic_loess")) {
      mat_active_normalized <- profile_expr("Cyclic LOESS normalization", {
        normalize_cyclic_loess_ref(
          mat_active_filtered,
          log2_offset = log2_offset,
          output_scale = "linear"
        )
      }, category = "normalization")

      if (debug_mode) {
        out_cyclic_loess_mat <- file.path(paths$global$exports, "11_MATRIX_ACTIVE_postFilters_postCyclicLOESS_WITH_QC.csv")
        write.csv(mat_active_normalized, file = out_cyclic_loess_mat, row.names = TRUE)
        log_written_object(log_path, out_cyclic_loess_mat, mat_active_normalized, note = "technical ACTIVE matrix after cyclic LOESS with QC")
      }

      assay_num_weight <- mat_active_normalized
      step_info("Cyclic LOESS normalization applied with limma::normalizeCyclicLoess.")
    } else {
      assay_num_weight <- mat_active_filtered
      step_info("No QC-based normalization applied after filters.")
    }

    mat_active_base_technical <- mat_active_base
    mat_active_normalized_technical <- mat_active_normalized

    metadata_biological_final <- metadata_aligned %>%
      dplyr::filter(type == "Sample", sample %in% rownames(mat_active_normalized_technical)) %>%
      dplyr::distinct(sample, .keep_all = TRUE)

    final_biological_samples <- metadata_biological_final$sample
    if (length(final_biological_samples) < 2) {
      stop("Final biological matrix has fewer than 2 biological samples after QC removal.")
    }

    mat_active_base_bio <- mat_active_base_technical[final_biological_samples, , drop = FALSE]
    mat_active_normalized_bio <- mat_active_normalized_technical[final_biological_samples, , drop = FALSE]

    qc_samples_in_final <- metadata_aligned %>%
      dplyr::filter(type == "QC", sample %in% rownames(mat_active_normalized_bio)) %>%
      dplyr::pull(sample)

    if (length(qc_samples_in_final) > 0) {
      stop(
        "Final biological matrix contains QC samples: ",
        paste(qc_samples_in_final, collapse = ", "),
        ". QC samples must be removed before final PCA/statistics/volcano/heatmap."
      )
    }

    build_weight_normalization_audit <- function(mat_before, mat_after, metadata_final, applied, stage) {
      tibble::tibble(
        sample = rownames(mat_after),
        type = metadata_final$type[match(rownames(mat_after), metadata_final$sample)],
        weight = metadata_final$weight[match(rownames(mat_after), metadata_final$sample)],
        stage = stage,
        weight_normalization_applied = isTRUE(applied),
        n_features = ncol(mat_after),
        n_missing_before = rowSums(is.na(mat_before)),
        n_missing_after = rowSums(is.na(mat_after)),
        median_intensity_before = apply(mat_before, 1, stats::median, na.rm = TRUE),
        median_intensity_after = apply(mat_after, 1, stats::median, na.rm = TRUE)
      )
    }

    if (normalization_mode_local %in% c("qc_loess", "qcrsc")) {
      correction_label <- if (identical(normalization_mode_local, "qcrsc")) "QC_RSC" else "QC_LOESS"
      correction_label_pretty <- if (identical(normalization_mode_local, "qcrsc")) "QC-RSC" else "QC-LOESS"
      mat_active_post_loess_bio_pre_weight <- mat_active_normalized_bio
      mat_active_normalized_bio <- normalize_by_weight(
        mat_active_normalized_bio,
        metadata_biological_final,
        seq_len(nrow(metadata_biological_final)),
        use_weight_normalization = use_weight_normalization,
        stop_on_invalid_weight   = stop_on_invalid_weight,
        invalid_weight_to_NA     = invalid_weight_to_NA
      )
      assay_num_weight <- mat_active_normalized_bio

      weight_audit <- build_weight_normalization_audit(
        mat_before = mat_active_post_loess_bio_pre_weight,
        mat_after = mat_active_normalized_bio,
        metadata_final = metadata_biological_final,
        applied = use_weight_normalization,
        stage = paste0("post_", correction_label, "_biological_no_QC")
      )

      if (debug_mode) {
        out_weight_audit <- file.path(paths$global$audits, "weight_normalization_audit_biological_samples.csv")
        write_csv_safe(weight_audit, out_weight_audit)
        log_written_object(log_path, out_weight_audit, weight_audit, note = paste("weight normalization audit on biological samples after", correction_label_pretty))

        out_weight <- file.path(paths$global$exports, paste0("12_MATRIX_ACTIVE_post", correction_label, "_postWeight_BIOLOGICAL_NO_QC.csv"))
        write.csv(mat_active_normalized_bio, file = out_weight, row.names = TRUE)
        log_written_object(log_path, out_weight, mat_active_normalized_bio, note = paste("final biological ACTIVE matrix after", correction_label_pretty, "and weight normalization, no QC"))
      }

      step_info("Weight normalization after ", correction_label_pretty, " applied only to biological samples: ", use_weight_normalization)
    } else if (identical(normalization_mode_local, "weight")) {
      weight_audit <- build_weight_normalization_audit(
        mat_before = mat_active_base_bio,
        mat_after = mat_active_normalized_bio,
        metadata_final = metadata_biological_final,
        applied = TRUE,
        stage = "weight_only_biological_no_QC"
      )

      if (debug_mode) {
        out_weight_audit <- file.path(paths$global$audits, "weight_normalization_audit_biological_samples.csv")
        write_csv_safe(weight_audit, out_weight_audit)
        log_written_object(log_path, out_weight_audit, weight_audit, note = "weight normalization audit on biological samples")
      }
    }

    assay_num_normalized_technical <- mat_active_normalized_technical
    variants$ACTIVE$mat_technical_with_qc <- mat_active_normalized_technical
    variants$ACTIVE$mat <- mat_active_normalized_bio

    step_info("ACTIVE variant selected: ", active_variant_effective)
    step_info("Technical ACTIVE matrix with QC after normalization: ", fmt_dims(mat_active_normalized_technical))
    step_info("Final biological ACTIVE matrix without QC: ", fmt_dims(variants$ACTIVE$mat))
    step_info("Duplicate strategy: ", duplicate_name_strategy)

    if (isTRUE(export_normalization_audit)) {
      normalization_summary <- tibble::tibble(
        output_level = output_level_local,
        normalization_mode = normalization_mode_local,
        weight_normalization = isTRUE(use_weight_normalization),
        technical_samples_with_qc = nrow(mat_active_normalized_technical),
        biological_samples_no_qc = nrow(variants$ACTIVE$mat),
        features_final = ncol(variants$ACTIVE$mat),
        qc_samples_available = length(qc_idx),
        injection_order_source = as.character(assay_bundle$drift_injection_order_source)[1]
      )
      normalization_summary_path <- file.path(paths$global$audits, "normalization_summary.csv")
      write_csv_safe(normalization_summary, normalization_summary_path)
      log_written_object(log_path, normalization_summary_path, normalization_summary, note = "summarized normalization audit")
    }

    if (isTRUE(export_qc_summary)) {
      qc_summary <- tibble::tibble(
        output_level = output_level_local,
        qc_samples = length(qc_idx),
        biological_samples = length(sample_idx),
        features_before_normalization = ncol(mat_active_filtered),
        features_final = ncol(variants$ACTIVE$mat),
        qc_diagnostics_requested = isTRUE(make_qc_diagnostics),
        normalization_mode = normalization_mode_local
      )
      qc_summary_path <- file.path(paths$global$audits, "qc_summary.csv")
      write_csv_safe(qc_summary, qc_summary_path)
      log_written_object(log_path, qc_summary_path, qc_summary, note = "QC summary")
    }

    step_ok("Applying feature-level filters, normalization, and duplicate handling", t_step)

    mat_log2_base_technical <- log2_transform(mat_active_base_technical, log2_offset)
    rownames(mat_log2_base_technical) <- rownames(mat_active_base_technical)

    mat_log2_normalized_technical <- log2_transform(mat_active_normalized_technical, log2_offset)
    rownames(mat_log2_normalized_technical) <- rownames(mat_active_normalized_technical)

    mat_log2_base <- log2_transform(mat_active_base_bio, log2_offset)
    rownames(mat_log2_base) <- rownames(mat_active_base_bio)

    mat_log2_normalized <- log2_transform(mat_active_normalized_bio, log2_offset)
    rownames(mat_log2_normalized) <- rownames(mat_active_normalized_bio)

    mat_prelog_final <- mat_active_normalized_bio
    mat_log2_final <- mat_log2_normalized

    # -------------------------------------------------------------------------
    # STEP 12 — Exporting raw data for MetaboAnalyst and final biological matrix
    # -------------------------------------------------------------------------
    t_step <- step_start(12, 15, "Exporting raw data for MetaboAnalyst and final biological matrix")

    metaboanalyst_base <- metaboanalyst_export_mat
    rownames(metaboanalyst_base) <- rownames(metaboanalyst_export_mat)

    base_raw_df <- as.data.frame(metaboanalyst_base) %>%
      tibble::rownames_to_column("sample")

    if (isTRUE(export_intermediate_tables)) {
      out_final_biological_matrix <- file.path(paths$global$exports, "12_MATRIX_ACTIVE_FINAL_BIOLOGICAL_NO_QC_LOG2.csv")
      write.csv(mat_log2_final, file = out_final_biological_matrix, row.names = TRUE)
      log_written_object(log_path, out_final_biological_matrix, mat_log2_final, note = "final biological log2 matrix without QC")
    }

    if (isTRUE(export_intermediate_tables)) {
      export_metaboanalyst_global_raw(
        raw_df = base_raw_df,
        metadata_aligned = metadata_aligned,
        feature_tbl = metaboanalyst_export_feature,
        export_dir = paths$global$exports,
        log_path = log_path,
        value_label = if (identical(normalization_mode_local, "weight")) "weight_ONLY" else "raw"
      )
    }

    step_info("MetaboAnalyst export feature set: post missing/presence-impute/known/duplicate; excludes QC-RSD and IQR.")
    if (identical(normalization_mode_local, "weight")) {
      step_info("MetaboAnalyst export values: weight-normalized only.")
    }
    step_info("MetaboAnalyst export enabled: ", export_metaboanalyst_ready)

    if (isTRUE(export_intermediate_tables) &&
        exists("metaboanalyst_duplicate_only_mat", inherits = FALSE)) {
      dup_only_df <- as.data.frame(metaboanalyst_duplicate_only_mat, check.names = FALSE) %>%
        tibble::rownames_to_column("sample")

      export_metaboanalyst_global_raw(
        raw_df = dup_only_df,
        metadata_aligned = metadata_aligned,
        feature_tbl = metaboanalyst_duplicate_only_feature,
        export_dir = paths$global$exports,
        log_path = log_path,
        value_label = "duplicate_ONLY"
      )
    }

    if (isTRUE(export_intermediate_tables) && exists("metaboanalyst_compatible_iqr_mat", inherits = FALSE)) {
      compatible_iqr_df <- as.data.frame(metaboanalyst_compatible_iqr_mat, check.names = FALSE) %>%
        tibble::rownames_to_column("sample")

      export_metaboanalyst_global_raw(
        raw_df = compatible_iqr_df,
        metadata_aligned = metadata_aligned,
        feature_tbl = metaboanalyst_compatible_iqr_feature,
        export_dir = paths$global$exports,
        log_path = log_path,
        value_label = "duplicate_IQR",
        file_prefix = "MA_COMPATIBLE",
        include_with_qc = FALSE
      )
      step_info("MetaboAnalyst compatible IQR export: duplicate handling + IQR only; do not apply IQR again on the site.")
    }
    step_info("MetaboAnalyst duplicate-only export enabled: ", export_metaboanalyst_duplicate_only)
    step_ok("Exporting raw data for MetaboAnalyst and final biological matrix", t_step)

    # -------------------------------------------------------------------------
    # STEP 13 — Generating PCA plots
    # -------------------------------------------------------------------------
    t_step <- step_start(13, 15, "Generating PCA plots")

    step_info("Technical PCA before-normalization matrix with QC: mat_log2_base_technical | ", fmt_dims(mat_log2_base_technical))
    step_info("Technical PCA after-normalization matrix with QC: mat_log2_normalized_technical | ", fmt_dims(mat_log2_normalized_technical))
    step_info("Biological PCA before-normalization matrix without QC: mat_log2_base | ", fmt_dims(mat_log2_base))
    step_info("Biological PCA after-normalization matrix without QC: mat_log2_normalized | ", fmt_dims(mat_log2_normalized))
    step_info("Final biological PCA/statistics/heatmap matrix without QC: mat_log2_final | ", fmt_dims(mat_log2_final))

    pca_technical_pre_post_max_abs_diff <- if (identical(dim(mat_log2_base_technical), dim(mat_log2_normalized_technical)) &&
                                               identical(rownames(mat_log2_base_technical), rownames(mat_log2_normalized_technical)) &&
                                               identical(colnames(mat_log2_base_technical), colnames(mat_log2_normalized_technical))) {
      diff_values <- abs(mat_log2_base_technical - mat_log2_normalized_technical)
      diff_max <- suppressWarnings(max(diff_values, na.rm = TRUE))
      if (is.finite(diff_max)) diff_max else NA_real_
    } else {
      NA_real_
    }
    step_info(
      "Technical PCA mat_pre vs mat_post max(abs(mat_pre - mat_post), na.rm = TRUE): ",
      ifelse(is.na(pca_technical_pre_post_max_abs_diff), "NA", format(pca_technical_pre_post_max_abs_diff, scientific = TRUE, digits = 6)),
      " | same_dim=", identical(dim(mat_log2_base_technical), dim(mat_log2_normalized_technical)),
      " | same_rows=", identical(rownames(mat_log2_base_technical), rownames(mat_log2_normalized_technical)),
      " | same_cols=", identical(colnames(mat_log2_base_technical), colnames(mat_log2_normalized_technical))
    )

    pca_biological_pre_post_max_abs_diff <- if (identical(dim(mat_log2_base), dim(mat_log2_normalized)) &&
                                                identical(rownames(mat_log2_base), rownames(mat_log2_normalized)) &&
                                                identical(colnames(mat_log2_base), colnames(mat_log2_normalized))) {
      diff_values <- abs(mat_log2_base - mat_log2_normalized)
      diff_max <- suppressWarnings(max(diff_values, na.rm = TRUE))
      if (is.finite(diff_max)) diff_max else NA_real_
    } else {
      NA_real_
    }
    step_info(
      "Biological PCA mat_pre vs mat_post max(abs(mat_pre - mat_post), na.rm = TRUE): ",
      ifelse(is.na(pca_biological_pre_post_max_abs_diff), "NA", format(pca_biological_pre_post_max_abs_diff, scientific = TRUE, digits = 6)),
      " | same_dim=", identical(dim(mat_log2_base), dim(mat_log2_normalized)),
      " | same_rows=", identical(rownames(mat_log2_base), rownames(mat_log2_normalized)),
      " | same_cols=", identical(colnames(mat_log2_base), colnames(mat_log2_normalized))
    )

    suppressWarnings({
      profile_expr("Combined biological PCA", {
        if (isTRUE(export_all_plots)) {
          plot_global_pca_exports(
            mat_log2_pre = mat_log2_base_technical,
            mat_log2_post = mat_log2_normalized_technical,
            metadata_aligned = metadata_aligned,
            paths = paths,
            pca_scaling = pca_scaling,
            ellipse_positive = if (exists("ellipse_positive", inherits = TRUE)) ellipse_positive else TRUE,
            log_path = log_path,
            file_tag = "TECHNICAL_WITH_QC",
            title_label = "technical with QC"
          )
        }

        plot_global_pca_exports(
          mat_log2_pre = mat_log2_base,
          mat_log2_post = mat_log2_normalized,
          metadata_aligned = metadata_biological_final,
          paths = paths,
          pca_scaling = pca_scaling,
          ellipse_positive = if (exists("ellipse_positive", inherits = TRUE)) ellipse_positive else TRUE,
          log_path = log_path,
          file_tag = "BIOLOGICAL_NO_QC",
          title_label = "biological no QC",
          stages = if (isTRUE(export_all_plots)) NULL else "post_normalization",
          include_per_model = export_all_plots
        )
      }, category = "plot")

      if (isTRUE(export_all_plots)) {
        profile_expr("Per-model PCA pre/post", {
          plot_pca_pre_post_per_model(
            mat_log2_pre = mat_log2_base,
            mat_log2_post = mat_log2_normalized,
            metadata_biological_final = metadata_biological_final,
            paths = paths,
            pca_scaling = pca_scaling,
            pca_label_samples = pca_label_samples,
            correction_label = if (exists("correction_label", inherits = FALSE)) correction_label else normalization_mode_local,
            log_path = log_path
          )
        }, category = "plot")
      }

      profile_expr("Per-model PCA final", {
        plot_pca_per_model(
          mat_log2_final,
          metadata_biological_final,
          paths,
          pca_scaling = pca_scaling,
          ellipse_positive = if (exists("ellipse_positive", inherits = TRUE)) ellipse_positive else TRUE,
          pca_label_samples = pca_label_samples,
          comparison_names = if (isTRUE(export_all_plots)) NULL else "ALL_TGvsWT",
          include_secondary_pca = export_all_plots
        )
      }, category = "plot")

      if (isTRUE(export_multigroup_outputs) && tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")) {
        profile_expr("Exploratory multi-group PCA", {
          plot_pca_multigroup_per_model(
            mat_log2 = mat_log2_final,
            metadata_aligned = metadata_biological_final,
            paths = paths,
            multigroup_groups = multigroup_groups,
            pca_scaling = pca_scaling,
            ellipse_positive = if (exists("ellipse_positive", inherits = TRUE)) ellipse_positive else TRUE,
            pca_label_samples = pca_label_samples,
            log_path = log_path
          )
        }, category = "plot")
      }

      step_info("PCA scaling: ", pca_scaling)
      step_ok("Generating PCA plots", t_step)

      # -------------------------------------------------------------------------
      # STEP 14 — Preparing statistics-driven heatmaps
      # -------------------------------------------------------------------------
      t_step <- step_start(14, 15, "Preparing statistics-driven heatmaps")

      heatmap_rank_metrics_cfg <- if (exists("heatmap_rank_metrics", inherits = TRUE)) {
        heatmap_rank_metrics
      } else {
        run_metrics
      }

      step_info("Top/significant heatmaps will be generated after statistics.")
      step_info("Top heatmap metrics configured: ", paste(heatmap_rank_metrics_cfg, collapse = ", "))
      step_info("Significant heatmap metrics configured: ", paste(run_metrics, collapse = ", "))
      step_ok("Preparing statistics-driven heatmaps", t_step)

      # -------------------------------------------------------------------------
      # STEP 15 — Generating statistics, volcano plots, stats Excel, and heatmaps
      # -------------------------------------------------------------------------
      t_step <- step_start(15, 15, "Generating statistics, volcano plots, stats Excel, and heatmaps")
      run_metrics_expanded <- normalize_metric_modes(run_metrics)
      heatmap_rank_metrics_expanded <- normalize_metric_modes(heatmap_rank_metrics_cfg)

      step_info("Primary stats (ALL_TGvsWT, F_TGvsWT, M_TGvsWT, TG_FvsM, WT_FvsM) with pairwise volcano plots.")
      if (tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")) {
        step_info("Exploratory MULTIGROUP_GLOBAL enabled: no directional FC, Up/Down classification, or volcano plot.")
      }

      stats_5sets_by_model <- profile_expr("Statistics and volcano plots", {
        run_all_stats_5sets_per_model(
          mat_log2 = mat_log2_final,
          mat_raw = mat_prelog_final,
          metadata_aligned = metadata_biological_final,
          feature_tbl = variants$ACTIVE$feature,
          paths = paths,
          p_value_cutoff = p_value_cutoff,
          fdr_cutoff = fdr_cutoff,
          fc_cutoff_log2 = fc_cutoff_log2,
          run_metrics = run_metrics,
          make_volcano_plots = TRUE,
          volcano_style = volcano_style,
          comparison_configs = COMPARISON_CONFIGS,
          statistical_test_type = statistical_test_type,
          test_is_paired = test_is_paired,
          pvalue_correction_method = pvalue_correction_method,
          comparison_mode = comparison_mode,
          multigroup_groups = multigroup_groups,
          multigroup_test = multigroup_test,
          multigroup_pairwise_mode = multigroup_pairwise_mode,
          multigroup_pairwise_pairs = multigroup_pairwise_pairs,
          export_all_pairwise_multigroup = export_all_pairwise_multigroup
        )
      }, category = "stats")

      # -----------------------------------------------------------------------
      # Top heatmaps by model
      # -----------------------------------------------------------------------
      if (!isTRUE(minimal_output_enabled)) {
        step_info("TOP heatmaps per model (loop rank metrics)...")

        profile_expr("Top heatmaps by model", {
          suppressWarnings({
            for (rk in heatmap_rank_metrics_expanded) {
              plot_heatmap_top_ttest_per_model(
                mat_log2_final,
                metadata_biological_final,
                variants$ACTIVE$feature,
                paths,
                stats_5sets_by_model = stats_5sets_by_model,
                top_n = heatmap_top_n,
                rank_by = rk,
                split_by_sex = FALSE,
                order_samples_by_group = heatmap_order_samples_by_group,
                scale_method = heatmap_scale_method,
                comparison_mode = "pairwise",
                multigroup_groups = multigroup_groups
              )
              if (isTRUE(export_multigroup_outputs) && tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")) {
                plot_heatmap_top_ttest_per_model(
                  mat_log2_final,
                  metadata_biological_final,
                  variants$ACTIVE$feature,
                  paths,
                  stats_5sets_by_model = stats_5sets_by_model,
                  top_n = heatmap_top_n,
                  rank_by = rk,
                  split_by_sex = FALSE,
                  order_samples_by_group = heatmap_order_samples_by_group,
                  scale_method = heatmap_scale_method,
                  comparison_mode = "multigroup",
                  multigroup_groups = multigroup_groups
                )
              }
            }
          })
        }, category = "plot")
      }

      # -----------------------------------------------------------------------
      # Top heatmaps by model and sex
      # -----------------------------------------------------------------------
      if (isTRUE(export_all_plots)) {
        step_info("TOP heatmaps per model split by sex (loop rank metrics)...")

        profile_expr("Top heatmaps by model and sex", {
          suppressWarnings({
            for (rk in heatmap_rank_metrics_expanded) {
              plot_heatmap_top_ttest_per_model(
                mat_log2_final,
                metadata_biological_final,
                variants$ACTIVE$feature,
                paths,
                stats_5sets_by_model = stats_5sets_by_model,
                top_n = heatmap_top_n,
                rank_by = rk,
                split_by_sex = TRUE,
                order_samples_by_group = heatmap_order_samples_by_group,
                scale_method = heatmap_scale_method,
                comparison_mode = "pairwise",
                multigroup_groups = multigroup_groups
              )
              if (isTRUE(export_multigroup_outputs) && tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")) {
                plot_heatmap_top_ttest_per_model(
                  mat_log2_final,
                  metadata_biological_final,
                  variants$ACTIVE$feature,
                  paths,
                  stats_5sets_by_model = stats_5sets_by_model,
                  top_n = heatmap_top_n,
                  rank_by = rk,
                  split_by_sex = TRUE,
                  order_samples_by_group = heatmap_order_samples_by_group,
                  scale_method = heatmap_scale_method,
                  comparison_mode = "multigroup",
                  multigroup_groups = multigroup_groups
                )
              }
            }
          })
        }, category = "plot")
      }

      {
        step_info("Exporting stats to Excel (README + 5 tabs per model)...")

        profile_expr("Stats Excel export", {
          export_stats_excel_by_model(
            stats_5sets_by_model,
            paths = paths,
            p_value_cutoff = p_value_cutoff,
            fdr_cutoff = fdr_cutoff,
            fc_cutoff_log2 = fc_cutoff_log2,
            active_variant = active_variant_effective,
            log_path = log_path,
            statistical_test_type = statistical_test_type,
            test_is_paired = test_is_paired,
            pvalue_correction_method = pvalue_correction_method,
            comparison_mode = comparison_mode,
            multigroup_test = multigroup_test,
            multigroup_groups = multigroup_groups,
            multigroup_pairwise_mode = multigroup_pairwise_mode,
            run_metrics = run_metrics,
            include_multigroup_outputs = export_multigroup_outputs,
            output_level = output_level_local
          )
        }, category = "export")
      }

      if (isTRUE(export_debug_outputs)) {
        step_info("Exporting significant metabolites TXT files per model/comparison...")

        profile_expr("Significant metabolites TXT export", {
          export_significant_metabolites_txt_by_model(
            stats_5sets_by_model,
            paths = paths,
            p_value_cutoff = p_value_cutoff,
            fdr_cutoff = fdr_cutoff,
            fc_cutoff_log2 = fc_cutoff_log2,
            active_variant = active_variant_effective,
            log_path = log_path,
            require_fc_cutoff = FALSE
          )
        }, category = "export")
      }

      models <- sort(unique(metadata_biological_final$model))
      sexes <- c("F", "M")

      # -----------------------------------------------------------------------
      # Significant heatmaps by model (ALL sex TG vs WT)
      # -----------------------------------------------------------------------
      if (isTRUE(export_all_plots)) {
        profile_expr("Significant heatmaps by model", {
          step_info("Significant heatmaps (ALL sex ", comparison_group_treatment, "/", comparison_group_control, ") for BOTH metrics...")

          for (met in run_metrics_expanded) {
            # Use the appropriate cutoff for this metric
            current_alpha <- if (met == "FDR") fdr_cutoff else p_value_cutoff
          
            for (m in models) {
              mp <- get_model_paths(paths, m)
              model_groups <- resolve_model_group_values(m)
              meta_m <- metadata_biological_final %>%
                dplyr::filter(
                  type == "Sample",
                  model == m,
                  group %in% c(model_groups$control, model_groups$treatment)
                )

              comparison_label <- paste0(model_groups$treatment, "vs", model_groups$control)
              st <- stats_5sets_by_model[[m]][["ALL_TGvsWT"]]
              if (!is.null(st)) {
                out_png <- file.path(
                  mp$plots$heatmap_significant_all,
                  paste0(
                    "HEATMAP_SIG_ACTIVE_ALL_", comparison_label, "_", met, "_lt_", current_alpha,
                    "_model_", m, "_scale_", heatmap_scale_method, ".png"
                  )
                )

                plot_sig_heatmap_from_stats(
                  mat_log2_final,
                  meta_m,
                  variants$ACTIVE$feature,
                  st,
                  sig_metric = met,
                  alpha_sig = current_alpha,
                  fc_cutoff_log2 = fc_cutoff_log2,
                  require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                  sig_max = sig_heatmap_max_features,
                  scale_method = heatmap_scale_method,
                  order_samples_by_group = heatmap_order_samples_by_group,
                  out_png = out_png,
                  title_main = paste0(
                    "SIG (", comparison_label, ") | model=", m,
                    " | sex=ALL | ", met, "<", current_alpha,
                    if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                    " | scale=", heatmap_scale_method
                  )
                )
              }
            }
          }
        }, category = "plot")
      }

      # -----------------------------------------------------------------------
      # Significant heatmaps by model and sex (TG vs WT within F and within M)
      # -----------------------------------------------------------------------
      if (isTRUE(export_all_plots)) {
        profile_expr("Significant heatmaps by model and sex", {
          step_info("Significant heatmaps (", comparison_group_treatment, "/", comparison_group_control, ") BY sex for BOTH metrics...")

          for (met in run_metrics_expanded) {
            # Use the appropriate cutoff for this metric
            current_alpha <- if (met == "FDR") fdr_cutoff else p_value_cutoff
          
            for (m in models) {
              mp <- get_model_paths(paths, m)
              model_groups <- resolve_model_group_values(m)

              for (sx in sexes) {
                meta_m <- metadata_biological_final %>%
                  dplyr::filter(type == "Sample", model == m, sex == sx, group %in% c(model_groups$control, model_groups$treatment))

                st <- if (sx == "F") {
                  stats_5sets_by_model[[m]][["F_TGvsWT"]]
                } else {
                  stats_5sets_by_model[[m]][["M_TGvsWT"]]
                }

                if (is.null(st)) next

                out_png <- file.path(
                  mp$plots$heatmap_significant_by_sex,
                  paste0(
                    "HEATMAP_SIG_ACTIVE_", sx, "_", model_groups$treatment, "vs", model_groups$control, "_", met, "_lt_", current_alpha,
                    "_model_", m, "_scale_", heatmap_scale_method, ".png"
                  )
                )

                plot_sig_heatmap_from_stats(
                  mat_log2_final,
                  meta_m,
                  variants$ACTIVE$feature,
                  st,
                  sig_metric = met,
                  alpha_sig = current_alpha,
                  fc_cutoff_log2 = fc_cutoff_log2,
                  require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                  sig_max = sig_heatmap_max_features,
                  scale_method = heatmap_scale_method,
                  order_samples_by_group = heatmap_order_samples_by_group,
                  out_png = out_png,
                  title_main = paste0(
                    "SIG (", model_groups$treatment, "/", model_groups$control, ") | model=", m,
                    " | sex=", sx,
                    " | ", met, "<", current_alpha,
                    if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                    " | scale=", heatmap_scale_method
                  )
                )
              }
            }
          }
        }, category = "plot")
      }

      # -----------------------------------------------------------------------
      # Significant heatmaps by model and sex (FvsM within treatment and within control)
      # -----------------------------------------------------------------------
      if (isTRUE(export_all_plots)) {
        profile_expr("Significant heatmaps FvsM within group", {
          step_info("Significant heatmaps (FvsM within ", comparison_group_treatment, " and ", comparison_group_control, ") for BOTH metrics...")

          for (met in run_metrics_expanded) {
            # Use the appropriate cutoff for this metric
            current_alpha <- if (met == "FDR") fdr_cutoff else p_value_cutoff
          
            for (m in models) {
              mp <- get_model_paths(paths, m)
              model_groups <- resolve_model_group_values(m)

              meta_tg <- metadata_biological_final %>%
                dplyr::filter(type == "Sample", model == m, group == model_groups$treatment, sex %in% c("F", "M"))
              st_tg <- stats_5sets_by_model[[m]][["TG_FvsM"]]

              if (!is.null(st_tg)) {
                out_png_tg <- file.path(
                  mp$plots$heatmap_significant_tg_f_vs_tg_m,
                  paste0(
                    "HEATMAP_SIG_ACTIVE_", model_groups$treatment, "_FvsM_", met, "_lt_", current_alpha,
                    "_model_", m, "_scale_", heatmap_scale_method, ".png"
                  )
                )

                plot_sig_heatmap_from_stats(
                  mat_log2_final,
                  meta_tg,
                  variants$ACTIVE$feature,
                  st_tg,
                  sig_metric = met,
                  alpha_sig = current_alpha,
                  fc_cutoff_log2 = fc_cutoff_log2,
                  require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                  sig_max = sig_heatmap_max_features,
                  scale_method = heatmap_scale_method,
                  order_samples_by_group = FALSE,
                  out_png = out_png_tg,
                  title_main = paste0(
                    "SIG (FvsM within ", model_groups$treatment, ") | model=", m,
                    " | ", met, "<", current_alpha,
                    if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                    " | log2FC=log2(F/M) | scale=", heatmap_scale_method
                  )
                )
              }

              meta_wt <- metadata_biological_final %>%
                dplyr::filter(type == "Sample", model == m, group == model_groups$control, sex %in% c("F", "M"))
              st_wt <- stats_5sets_by_model[[m]][["WT_FvsM"]]

              if (!is.null(st_wt)) {
                out_png_wt <- file.path(
                  mp$plots$heatmap_significant_wt_f_vs_wt_m,
                  paste0(
                    "HEATMAP_SIG_ACTIVE_", model_groups$control, "_FvsM_", met, "_lt_", current_alpha,
                    "_model_", m, "_scale_", heatmap_scale_method, ".png"
                  )
                )

                plot_sig_heatmap_from_stats(
                  mat_log2_final,
                  meta_wt,
                  variants$ACTIVE$feature,
                  st_wt,
                  sig_metric = met,
                  alpha_sig = current_alpha,
                  fc_cutoff_log2 = fc_cutoff_log2,
                  require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                  sig_max = sig_heatmap_max_features,
                  scale_method = heatmap_scale_method,
                  order_samples_by_group = FALSE,
                  out_png = out_png_wt,
                  title_main = paste0(
                    "SIG (FvsM within ", model_groups$control, ") | model=", m,
                    " | ", met, "<", current_alpha,
                    if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                    " | log2FC=log2(F/M) | scale=", heatmap_scale_method
                  )
                )
              }
            }
          }
        }, category = "plot")
      }

      step_info("Volcano plots enabled: ", make_volcano_plots)
      step_info("Volcano metrics: ", paste(run_metrics_expanded, collapse = ", "))
      step_info("Top heatmaps metrics: ", paste(heatmap_rank_metrics_expanded, collapse = ", "))
      step_info("Significant heatmaps metrics: ", paste(run_metrics_expanded, collapse = ", "))
      step_info("Stats Excel export enabled: ", save_stats_excel_per_model)
      step_info("Top heatmaps by model enabled: ", make_heatmap_by_model)
      step_info("Top heatmaps by model and sex enabled: ", make_heatmap_by_model_sex)
      step_info("Top heatmaps FvsM")
      step_info("Significant heatmaps by model enabled: ", make_sig_heatmap_by_model)
      step_info("Significant heatmaps by model and sex enabled: ", make_sig_heatmap_by_model_sex)
      step_info("Significant heatmaps FvsM enabled: ", make_sig_heatmap_FvsM_within_group)
      step_ok("Generating statistics, volcano plots, stats Excel, and heatmaps", t_step)
      # -------------------------------------------------------------------------
      # Finalize filter summary table and save
      # -------------------------------------------------------------------------

      if (file.exists(summary_csv)) {
        df_sum <- readr::read_csv(summary_csv, show_col_types = FALSE)
        log_written_object(log_path, summary_csv, df_sum, note = "Filter summary table")
      }

      if (isTRUE(export_debug_outputs)) {
        method_summary_path <- file.path(paths$global$exports, "method_summary.txt")
        tryCatch(
          {
            write_method_summary(
              method_summary_path,
              filter_summary = filter_summary,
              injection_order_source = assay_bundle$drift_injection_order_source
            )
            step_info("Method summary exported: ", method_summary_path)
          },
          error = function(e) {
            step_info("Method summary export skipped: ", conditionMessage(e))
          }
        )
      }

      remove_empty_directories(output_dir)
      if (isTRUE(export_debug_outputs)) {
        runtime_profile <- runtime_profile_write(
          file.path(paths$global$audits, "runtime_profile.csv")
        )
        runtime_profile_written <- TRUE
      }

      message("\n", console_rule())
      message("[PIPELINE DONE] Total runtime: ", fmt_time_sec(pipeline_t0), " sec")
      message(console_rule())
      message("[OUTPUT DIR] ", output_dir)
      message("[LOG FILE] ", log_path)
      message(console_rule())

      return(list(
        cd_raw = cd_raw,
        metadata = metadata,
        metadata_aligned = metadata_aligned,
        models_detected = models_detected,
        paths = paths,
        feature_tbl = feature_tbl,
        assay_bundle = assay_bundle,
        assay_num_raw = assay_num_raw,
        assay_num_weight = assay_num_weight,
        assay_num_normalized_technical = assay_num_normalized_technical,
        variants = variants,
        filter_summary = filter_summary,
        stats_5sets_by_model = stats_5sets_by_model,
        runtime_profile = runtime_profile,
        output_level = output_level_local,
        log_path = log_path
      ))
    })
  })
}
