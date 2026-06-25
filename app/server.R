server <- function(input, output, session) {
  status_message <- reactiveVal("Ready.")
  pipeline_log_text <- reactiveVal("No run executed yet.")
  package_status <- reactiveVal("")
  install_status_text <- reactiveVal("Click the button to check and install required packages.")
  missing_packages_state <- reactiveVal(setdiff(required_packages, rownames(installed.packages())))
  initializing <- reactiveVal(TRUE)
  clearing_inputs <- reactiveVal(FALSE)
  output_level_last_confirmed <- reactiveVal("standard")
  output_level_modal_active <- reactiveVal(FALSE)

  process_state <- reactiveValues(
    proc = NULL,
    log_file = NULL,
    pipeline_log_file = NULL,
    running = FALSE,
    started_at = NULL,
    pid = NULL
  )
  action_confirm <- reactiveValues(
    ask_run = TRUE,
    ask_stop = TRUE,
    pending_run_cfg = NULL
  )
  selected_result_image <- reactiveVal(NULL)
  gallery_refresh_tick <- reactiveVal(0)
  gallery_state <- reactiveValues(dir = NULL, prefix = NULL)
  session_started_at <- Sys.time()
  inputs_cleared_timestamp <- reactiveVal(NULL)
  use_reference_file_last_state <- reactiveVal(FALSE)

  pipeline_modal <- function(title, ..., footer, easyClose = TRUE, size = "m") {
    modalDialog(
      title = title,
      size = size,
      easyClose = easyClose,
      tags$div(
        class = "pipeline-modal-body",
        ...
      ),
      footer = tags$div(
        class = "pipeline-modal-footer",
        footer
      )
    )
  }

  refresh_package_status <- function() {
    missing_pkgs <- find_missing_packages()
    missing_packages_state(missing_pkgs)

    if (length(missing_pkgs) == 0) {
      install_status_text("All required packages are installed.")
    } else {
      install_status_text(paste0(
        "Missing packages detected (",length(missing_pkgs), "): ",
        paste(missing_pkgs, collapse = ", ")
      ))
    }

    invisible(missing_pkgs)
  }

  reset_settings_builder_inputs <- function() {
    if (!exists("settings_form_sections", inherits = TRUE)) {
      return(invisible(FALSE))
    }

    reset_one_setting <- function(spec) {
      key <- as.character(spec$key)[1]
      type <- if (!is.null(spec$type) && length(spec$type) > 0) as.character(spec$type)[1] else "text"
      id <- setting_input_id(key)
      default <- spec$default

      try(
        switch(type,
          checkbox = updateCheckboxInput(session, id, value = isTRUE(default)),
          logical_select = updateSelectInput(session, id, selected = if (isTRUE(default)) "TRUE" else "FALSE"),
          numeric = updateNumericInput(session, id, value = suppressWarnings(as.numeric(default)[1])),
          integer = updateNumericInput(session, id, value = suppressWarnings(as.integer(default)[1])),
          select = updateSelectInput(session, id, selected = as.character(default)[1]),
          multiselect = updateSelectizeInput(session, id, selected = setting_default_vector(default)),
          vector_numeric = updateSelectizeInput(session, id, selected = setting_default_numeric_vector(default)),
          vector_text = updateSelectizeInput(session, id, selected = setting_default_vector(default)),
          detected_multiselect = updateCheckboxGroupInput(session, id, selected = setting_default_vector(default)),
          nullable_vector_text = updateSelectizeInput(session, id, selected = setting_default_vector(default)),
          selectize_text = updateSelectizeInput(session, id, selected = as.character(default)[1]),
          sheet = updateTextInput(session, id, value = as.character(default)[1]),
          updateTextInput(session, id, value = as.character(default)[1])
        ),
        silent = TRUE
      )
    }

    for (section in settings_form_sections) {
      for (spec in section$fields) {
        reset_one_setting(spec)
      }
    }

    invisible(TRUE)
  }

  reset_common_inputs <- function() {
    shinyjs::reset("data_file")
    shinyjs::reset("metadata_file")
    shinyjs::reset("injection_order_file")
    shinyjs::reset("reference_file")
    reset_settings_builder_inputs()
    updateTextInput(session, "output_dir", value = "output")
    updateSelectizeInput(session, "allowed_metadata_groups", selected = character(0))
    updateCheckboxInput(session, "use_reference_file", value = FALSE)
    updateCheckboxInput(session, "use_weight_normalization", value = FALSE)
    updateCheckboxInput(session, "manual_metadata_cols", value = FALSE)
    updateCheckboxInput(session, "show_metadata_column_fields", value = FALSE)
    updateCheckboxInput(session, "manual_reference_cols", value = FALSE)
    output_level_last_confirmed("standard")
    reset_metadata_columns()
    reset_reference_columns()
    use_reference_file_last_state(FALSE)
    selected_result_image(NULL)
    pipeline_log_text("No run executed yet.")
    inputs_cleared_timestamp(Sys.time())
  }

  session$onFlushed(function() {
    initializing(TRUE)
    reset_common_inputs()
    status_message("Ready.")
    shinyjs::runjs("setTimeout(function() { Shiny.setInputValue('init_complete', true); }, 50);")
  })

  observeEvent(input$init_complete,
    {
      initializing(FALSE)
    },
    once = TRUE
  )

  observeEvent(input[[setting_input_id("normalization_mode")]],
    {
      if (isTRUE(initializing())) {
        return()
      }
      scenario <- normalize_normalization_mode(
        input[[setting_input_id("normalization_mode")]],
        default = NULL
      )
      if (is.null(scenario) || !scenario %in% c("qc_loess", "qcrsc")) {
        return()
      }
      correction_label <- if (identical(scenario, "qcrsc")) "QC-RSC" else "QC-LOESS"

      showModal(pipeline_modal(
        title = "Weight normalization",
        tags$p(paste0("Apply weight normalization after ", correction_label, "?")),
        radioButtons(
          "modal_weight_norm_choice",
          label = NULL,
          choices = c("Yes" = "yes", "No" = "no"),
          selected = if (isTRUE(input$use_weight_normalization)) "yes" else "no"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_weight_norm", "Confirm")
        ),
        easyClose = TRUE
      ))
    },
    ignoreInit = TRUE
  )

  observeEvent(input$confirm_weight_norm, {
    choice <- isolate(input$modal_weight_norm_choice)
    updateCheckboxInput(session, "use_weight_normalization", value = identical(choice, "yes"))
    removeModal()
  })

  observeEvent(input[[setting_input_id("output_level")]],
    {
      if (isTRUE(initializing()) || isTRUE(clearing_inputs()) || isTRUE(output_level_modal_active())) {
        return()
      }

      selected_level <- normalize_output_level(
        input[[setting_input_id("output_level")]],
        legacy_minimal = FALSE
      )
      if (identical(selected_level, "standard")) {
        output_level_last_confirmed("standard")
        return()
      }

      output_level_modal_active(TRUE)
      level_label <- switch(selected_level,
        minimal = "Minimal",
        full_debug = "Full / Debug",
        selected_level
      )

      showModal(pipeline_modal(
        title = "Confirm output level",
        size = "m",
        easyClose = FALSE,
        tags$p(paste0("You selected ", level_label, " output.")),
        tags$p("Standard is recommended for routine runs. Confirm that you want to use this output level?"),
        footer = tagList(
          actionButton("cancel_output_level_change", "Cancel"),
          actionButton("confirm_output_level_change", "Confirm", class = "btn-warning")
        )
      ))
    },
    ignoreInit = TRUE
  )

  observeEvent(input$confirm_output_level_change,
    {
      selected_level <- normalize_output_level(
        input[[setting_input_id("output_level")]],
        legacy_minimal = FALSE
      )
      output_level_last_confirmed(selected_level)
      output_level_modal_active(FALSE)
      removeModal()
    },
    ignoreInit = TRUE
  )

  observeEvent(input$cancel_output_level_change,
    {
      previous_level <- normalize_output_level(
        output_level_last_confirmed(),
        legacy_minimal = FALSE
      )
      output_level_modal_active(TRUE)
      updateSelectInput(session, setting_input_id("output_level"), selected = previous_level)
      session$onFlushed(function() {
        output_level_modal_active(FALSE)
      }, once = TRUE)
      removeModal()
    },
    ignoreInit = TRUE
  )

  build_settings_builder_config <- function(current_text) {
    cfg <- current_text

    for (section in settings_form_sections) {
      for (spec in section$fields) {
        input_value <- input[[setting_input_id(spec$key)]]
        if (is.null(input_value)) {
          next
        }
        replacement <- switch(spec$type,
          checkbox = setting_value_logical(input_value),
          logical_select = setting_value_logical(identical(as.character(input_value)[1], "TRUE")),
          numeric = setting_value_numeric(input_value, default = spec$default),
          integer = setting_value_integer(input_value, default = spec$default),
          select = setting_value_text(input_value),
          multiselect = setting_value_vector_text(input_value),
          vector_numeric = setting_value_vector_numeric(input_value),
          vector_text = setting_value_vector_text(input_value),
          detected_multiselect = setting_value_vector_text(input_value),
          nullable_vector_text = setting_value_vector_text(input_value, allow_null = TRUE),
          selectize_text = setting_value_text(input_value),
          sheet = setting_value_sheet(input_value),
          setting_value_text(input_value)
        )

        cfg <- replace_or_append(cfg, spec$key, replacement)
      }
    }

    cfg <- apply_active_variant_config(
      cfg,
      input[[setting_input_id("active_variant")]],
      input[[setting_input_id("rsd_thresholds")]]
    )

    scenario <- gsub("_", "", tolower(safe_trimws(extract_config_value(cfg, "normalization_mode"))))
    weight_selected <- isTRUE(input$use_weight_normalization)
    if (identical(scenario, "none") && isTRUE(weight_selected)) {
      cfg <- replace_or_append(cfg, "normalization_mode", dQuote("weight"))
    }
    cfg <- replace_or_append(cfg, "use_weight_normalization", if (isTRUE(weight_selected)) "TRUE" else "FALSE")

    cfg
  }

  output$settings_builder_ui <- renderUI({
    detected_groups <- get_detected_metadata_groups()
    build_settings_builder_ui(dynamic_choices = list(
      multigroup_groups = detected_groups,
      multigroup_pairwise_pairs = make_multigroup_pair_choices(detected_groups)
    ))
  })

  output$settings_glossary_ui <- renderUI({
    build_settings_glossary_ui()
  })

  if (requireNamespace("shinyFiles", quietly = TRUE)) {
    shinyFiles::shinyDirChoose(
      input,
      "browse_output_dir",
      roots = get_shiny_roots(),
      session = session
    )
  }

  output_level_enabled <- function() {
    normalize_output_level(
      input[[setting_input_id("output_level")]],
      legacy_minimal = FALSE
    )
  }

  # `replace_or_append()` moved to pipeline/R/03_helpers_io_log.R for reuse.


  metadata_column_mapping <- function() {
    list(
      sample = safe_trimws(input$metadata_col_sample),
      weight = safe_trimws(input$metadata_col_weight),
      group = safe_trimws(input$metadata_col_group),
      sex = safe_trimws(input$metadata_col_sex),
      model = safe_trimws(input$metadata_col_model)
    )
  }

  metadata_model_groups_input_id <- function(model_name) {
    paste0("metadata_model_groups_", sanitize_input_id_fragment(model_name))
  }

  metadata_groups_text_for_model <- function(model_name, model_groups, fallback = input$allowed_metadata_groups) {
    detected_groups <- model_groups[[model_name]]
    detected_groups <- trimws(as.character(detected_groups))
    detected_groups <- detected_groups[!is.na(detected_groups) & nzchar(detected_groups)]
    detected_groups <- detected_groups[!is_missing_like(detected_groups)]
    detected_groups <- detected_groups[!(toupper(detected_groups) %in% metadata_qc_group_aliases())]
    detected_groups <- unique(detected_groups)

    if (length(detected_groups) >= 2) {
      fallback_groups <- parse_allowed_groups(fallback)
      fallback_groups <- fallback_groups[nzchar(fallback_groups)]
      if (length(fallback_groups) > 0) {
        detected_norm <- toupper(detected_groups)
        fallback_present <- fallback_groups[toupper(fallback_groups) %in% detected_norm]
        remaining <- detected_groups[!(detected_norm %in% toupper(fallback_present))]
        detected_groups <- unique(c(fallback_present, remaining))
      }

      return(paste(detected_groups, collapse = ", "))
    }

    safe_trimws(fallback)
  }

  get_detected_metadata_models <- function() {
    md_path <- resolve_input_file("metadata")
    detect_metadata_models(read_metadata_for_app(md_path, metadata_column_mapping()))
  }

  get_detected_metadata_groups_by_model <- function() {
    md_path <- resolve_input_file("metadata")
    detect_metadata_groups_by_model(read_metadata_for_app(md_path, metadata_column_mapping()))
  }

  get_detected_metadata_groups <- function() {
    md_path <- resolve_input_file("metadata")
    detect_metadata_groups(read_metadata_for_app(md_path, metadata_column_mapping()))
  }

  output$allowed_metadata_groups_ui <- renderUI({
    detected_groups <- get_detected_metadata_groups()
    current_groups <- isolate(parse_allowed_groups(input$allowed_metadata_groups))
    selected_groups <- current_groups
    if (length(selected_groups) == 0) {
      selected_groups <- detected_groups
    }

    selectizeInput(
      "allowed_metadata_groups",
      "Allowed metadata groups",
      choices = unique(c(selected_groups, detected_groups)),
      selected = selected_groups,
      multiple = TRUE,
      options = list(
        create = TRUE,
        persist = TRUE,
        plugins = list("remove_button", "drag_drop"),
        placeholder = "Choose groups (control first, test second)"
      )
    )
  })

  metadata_allowed_groups_by_model <- function() {
    if (!isTRUE(input$manual_metadata_cols)) {
      return(setNames(character(0), character(0)))
    }

    models <- get_detected_metadata_models()
    if (length(models) == 0) {
      return(setNames(character(0), character(0)))
    }

    alias_values <- vapply(models, function(model_name) {
      safe_trimws(input[[metadata_model_groups_input_id(model_name)]])
    }, character(1), USE.NAMES = FALSE)

    names(alias_values) <- models
    alias_values <- alias_values[nzchar(alias_values)]
    alias_values
  }

  resolve_input_path <- function(uploaded, external_path, kind) {
    if (!is.null(uploaded)) {
      return(file.path("data", basename(uploaded$name)))
    }

    ext <- safe_trimws(external_path)
    if (nzchar(ext)) {
      return(ext)
    }

    status_message(paste(
      "No",
      kind,
      "path was provided. Keep existing config value or provide one."
    ))
    NULL
  }

  save_config_and_inputs <- function(config_text) {
    dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(project_root, "data"), recursive = TRUE, showWarnings = FALSE)
    mapping <- metadata_column_mapping()

    if (!is.null(input$data_file)) {
      file.copy(
        input$data_file$datapath,
        file.path(project_root, "data", basename(input$data_file$name)),
        overwrite = TRUE
      )
    }

    if (!is.null(input$metadata_file)) {
      file.copy(
        input$metadata_file$datapath,
        file.path(project_root, "data", basename(input$metadata_file$name)),
        overwrite = TRUE
      )
      persist_metadata_mapping(
        source_path = input$metadata_file$datapath,
        uploaded_name = input$metadata_file$name,
        mapping = mapping
      )
    }

    if (is.null(input$metadata_file) && has_metadata_mapping(mapping)) {
      metadata_source <- resolve_input_file("metadata", prefer_mapped = FALSE)
      if (!is.null(metadata_source) && nzchar(metadata_source) && file.exists(metadata_source)) {
        persist_metadata_mapping(
          source_path = metadata_source,
          mapping = mapping
        )
      }
    }

    if (!is.null(input$injection_order_file)) {
      file.copy(
        input$injection_order_file$datapath,
        file.path(project_root, "data", "Input Files.xlsx"),
        overwrite = TRUE
      )
    }

    if (!is.null(input$reference_file)) {
      file.copy(
        input$reference_file$datapath,
        file.path(project_root, "data", basename(input$reference_file$name)),
        overwrite = TRUE
      )
    }

    clean_config <- normalize_config_text(config_text)
    writeLines(clean_config, active_config_path, useBytes = TRUE)
    updateTextAreaInput(session, "config_text", value = clean_config)
  }

  build_quick_config <- function(current_text) {
    cfg <- build_settings_builder_config(current_text)
    allowed_groups <- parse_allowed_groups(input$allowed_metadata_groups)
    duplicate_strategy_effective <- safe_trimws(input[[setting_input_id("duplicate_name_strategy")]])
    if (!nzchar(duplicate_strategy_effective)) {
      duplicate_strategy_effective <- "collapse_best_qc_rsd"
    }

    output_dir_quick <- sanitize_output_dir_path(input$output_dir)
    cfg <- replace_or_append(cfg, "output_dir", paste0('"', output_dir_quick, '"'))
    normalization_mode_quick <- input[[setting_input_id("normalization_mode")]]
    if (is.null(normalization_mode_quick) || !nzchar(safe_trimws(normalization_mode_quick))) {
      normalization_mode_quick <- extract_config_value(cfg, "normalization_mode")
    }
    normalization_mode_quick_check <- normalize_normalization_mode(normalization_mode_quick, default = NULL)
    if (is.null(normalization_mode_quick) || !normalization_mode_quick_check %in% c("qc_loess", "cyclic_loess", "qcrsc", "weight", "none", "pqn_qc", "pqn_sample")) {
      normalization_mode_quick <- "qcrsc"
    }

    pca_scaling_quick <- safe_trimws(input[[setting_input_id("pca_scaling")]])
    if (!nzchar(pca_scaling_quick)) {
      pca_scaling_quick <- "pareto"
    }
    heatmap_scale_quick <- safe_trimws(input[[setting_input_id("heatmap_scale_method")]])
    if (!nzchar(heatmap_scale_quick)) {
      heatmap_scale_quick <- "zscore"
    }
    normalization_mode_quick <- normalize_normalization_mode(normalization_mode_quick)
    weight_normalization_effective <- isTRUE(input$use_weight_normalization)
    if (identical(normalization_mode_quick, "weight")) {
      normalization_mode_quick <- "none"
    }
    if (identical(normalization_mode_quick, "none") && isTRUE(weight_normalization_effective)) {
      normalization_mode_quick <- "weight"
    }

    cfg <- replace_or_append(cfg, "normalization_mode", dQuote(normalization_mode_quick))
    cfg <- replace_or_append(cfg, "use_weight_normalization", if (isTRUE(weight_normalization_effective)) "TRUE" else "FALSE")
    cfg <- apply_active_variant_config(
      cfg,
      input[[setting_input_id("active_variant")]],
      input[[setting_input_id("rsd_thresholds")]]
    )
    low_variance_quick <- safe_trimws(input[[setting_input_id("low_variance_filter_method")]])
    if (!nzchar(low_variance_quick) || !low_variance_quick %in% c("iqr", "none")) {
      low_variance_quick <- "none"
    }
    cfg <- replace_or_append(cfg, "low_variance_filter_method", dQuote(low_variance_quick))
    cfg <- replace_or_append(cfg, "pca_scaling", dQuote(pca_scaling_quick))
    cfg <- replace_or_append(cfg, "heatmap_scale_method", dQuote(heatmap_scale_quick))
    cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote(duplicate_strategy_effective))
    cfg <- replace_or_append(cfg, "p_value_cutoff", setting_value_numeric(input[[setting_input_id("p_value_cutoff")]], default = 0.05))
    cfg <- replace_or_append(cfg, "fdr_cutoff", setting_value_numeric(input[[setting_input_id("fdr_cutoff")]], default = 0.05))
    run_metrics_effective <- safe_trimws(input[[setting_input_id("run_metrics")]])
    statistical_test_effective <- safe_trimws(input[[setting_input_id("statistical_test_type")]])
    paired_effective <- identical(as.character(input[[setting_input_id("test_is_paired")]])[1], "TRUE")
    cfg <- replace_or_append(cfg, "run_metrics", dQuote(run_metrics_effective))
    cfg <- replace_or_append(cfg, "heatmap_rank_metrics", dQuote(run_metrics_effective))
    cfg <- replace_or_append(cfg, "statistical_test_type", dQuote(statistical_test_effective))
    cfg <- replace_or_append(cfg, "test_is_paired", if (paired_effective) "TRUE" else "FALSE")
    # Derive p-value correction method from run_metrics selection so that
    # `run_metrics` becomes the primary choice for significance metric.
    # If the user selected any option including "FDR", default correction is "FDR" (BH),
    # otherwise use "raw" (no correction) to match p-value usage.
    derived_pvalue_correction <- if (grepl("FDR", run_metrics_effective, ignore.case = TRUE)) {
      "FDR"
    } else {
      "raw"
    }
    cfg <- replace_or_append(cfg, "pvalue_correction_method", dQuote(derived_pvalue_correction))
    cfg <- replace_or_append(cfg, "output_level", dQuote(output_level_enabled()))
    cfg <- replace_or_append(cfg, "use_reference_file", if (isTRUE(input$use_reference_file)) "TRUE" else "FALSE")
    cfg <- replace_or_append(cfg, "reference_col_metabolite", dQuote(safe_trimws(input$reference_col_metabolite)))
    cfg <- replace_or_append(cfg, "reference_col_ref_ion", dQuote(safe_trimws(input$reference_col_ref_ion)))
    cfg <- replace_or_append(cfg, "reference_col_mz", dQuote(safe_trimws(input$reference_col_mz)))
    cfg <- replace_or_append(cfg, "reference_col_rt", dQuote(safe_trimws(input$reference_col_rt)))

    model_allowed_groups <- metadata_allowed_groups_by_model()
    cfg <- replace_or_append(cfg, "model_allowed_groups_by_model", format_named_character_vector(model_allowed_groups))

    if (length(allowed_groups) >= 2) {
      cfg <- replace_or_append(cfg, "comparison_group_control", dQuote(allowed_groups[1]))
      cfg <- replace_or_append(cfg, "comparison_group_treatment", dQuote(allowed_groups[2]))
    }

    # Derive alpha_sig based on run_metrics selection
    # When both metrics are selected, use FDR (more conservative)
    # Otherwise use the cutoff corresponding to the selected metric
    derived_alpha_sig <- if (grepl("FDR", run_metrics_effective, ignore.case = TRUE)) {
      "fdr_cutoff"
    } else {
      "p_value_cutoff"
    }
    cfg <- replace_or_append(cfg, "alpha_sig", derived_alpha_sig)

    data_path <- resolve_input_path(input$data_file, input$external_data_path, "data")
    metadata_path <- resolve_input_path(input$metadata_file, input$external_metadata_path, "metadata")
    reference_path <- resolve_input_path(input$reference_file, input$external_reference_path, "reference")

    if (!is.null(data_path)) {
      cfg <- replace_or_append(cfg, "cd_file_path", dQuote(data_path))
    }

    if (!is.null(metadata_path)) {
      metadata_path_cfg <- metadata_path
      mapping <- metadata_column_mapping()
      mapped_rel <- metadata_effective_rel_path(
        source_path = metadata_path,
        uploaded_name = if (!is.null(input$metadata_file)) input$metadata_file$name else NULL,
        mapping = mapping
      )
      if (!is.null(mapped_rel)) {
        metadata_path_cfg <- mapped_rel
      }
      cfg <- replace_or_append(cfg, "metadata_path", dQuote(metadata_path_cfg))
    }

    if (isTRUE(input$use_reference_file) && !is.null(reference_path)) {
      cfg <- replace_or_append(cfg, "reference_path", dQuote(reference_path))
    } else if (!isTRUE(input$use_reference_file)) {
      cfg <- replace_or_append(cfg, "reference_path", dQuote(""))
    }

    cfg
  }

  normalization_requires_injection_order <- function(config_text) {
    mode <- safe_trimws(extract_config_value(config_text, "normalization_mode"))
    normalize_normalization_mode(mode, default = NULL) %in% c("qc_loess", "qcrsc")
  }

  has_injection_order_file_for_run <- function() {
    if (!is.null(input$injection_order_file)) {
      return(TRUE)
    }

    candidate_paths <- unique(c(
      file.path(project_root, "data", "Input Files.xlsx"),
      list.files(file.path(project_root, "data"), pattern = "^input_order.*\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
    ))
    candidate_paths <- candidate_paths[file.exists(candidate_paths)]

    any(vapply(candidate_paths, function(path) {
      input_files_ref <- read_input_files_reference(path)
      !is.null(input_files_ref) && nrow(input_files_ref) > 0
    }, logical(1)))
  }

  normalized_output_dir <- function() {
    out <- sanitize_output_dir_path(input$output_dir)
    if (!nzchar(out)) {
      return("output")
    }
    out
  }

  build_expected_output_manifest <- function(output_level = output_level_enabled()) {
    out_dir <- normalized_output_dir()
    cfg <- build_quick_config(current_config_text())
    build_expected_output_manifest_from_config(cfg, out_dir, output_level)
  }

  launch_pipeline <- function() {
    rscript_cmd <- file.path(R.home("bin"), "Rscript")
    if (.Platform$OS.type == "windows") {
      rscript_cmd <- paste0(rscript_cmd, ".exe")
    }

    if (!file.exists(rscript_cmd)) {
      status_message("R script executable not found. Cannot run pipeline from UI.")
      return()
    }

    log_file <- tempfile(pattern = "pipeline_run_", fileext = ".log")
    output_log_file <- file.path(resolve_output_dir_abs(), "PIPELINE_LOG.txt")
    pipeline_log_text(paste0(
      "Starting pipeline at ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "...\nWaiting for Rscript output."
    ))

    if (requireNamespace("processx", quietly = TRUE)) {
      proc <- tryCatch(
        processx::process$new(
          command = rscript_cmd,
          args = c("pipeline/run_pipeline.R"),
          wd = project_root,
          stdout = log_file,
          stderr = log_file,
          cleanup = TRUE
        ),
        error = function(e) e
      )

      if (inherits(proc, "error")) {
        status_message(paste("Failed to start pipeline:", conditionMessage(proc)))
        return()
      }

      process_state$proc <- proc
      process_state$log_file <- log_file
      process_state$pipeline_log_file <- output_log_file
      process_state$running <- TRUE
      process_state$started_at <- Sys.time()
      process_state$pid <- tryCatch(proc$get_pid(), error = function(e) NULL)
      status_message("Pipeline is running in background. Open 'Pipeline Log' tab for live output.")
    } else {
      status_message("Package 'processx' not installed. Running pipeline synchronously without live streaming.")

      exit_code <- tryCatch(
        system2(
          rscript_cmd,
          args = c("pipeline/run_pipeline.R"),
          stdout = log_file,
          stderr = log_file,
          wait = TRUE
        ),
        error = function(e) {
          writeLines(paste("Pipeline execution failed:", conditionMessage(e)), log_file)
          1L
        }
      )

      if (file.exists(output_log_file)) {
        pipeline_log_text(read_log_preview(output_log_file))
      } else if (file.exists(log_file)) {
        pipeline_log_text(read_log_preview(log_file))
      } else {
        pipeline_log_text("No log file generated.")
      }

      if (isTRUE(as.integer(exit_code) == 0L)) {
        status_message("Pipeline run completed successfully.")
      } else {
        status_message(
          paste(
            "Pipeline run failed with exit code",
            as.integer(exit_code),
            ". Check Pipeline Log tab."
          )
        )
      }
    }
  }

  run_pipeline_now <- function(cfg_to_run) {
    updateTextAreaInput(session, "config_text", value = cfg_to_run)
    save_config_and_inputs(cfg_to_run)
    launch_pipeline()
  }

  stop_pipeline_now <- function() {
    if (!isTRUE(process_state$running) || is.null(process_state$proc)) {
      status_message("No running pipeline process to stop.")
      return()
    }

    if (process_state$proc$is_alive()) {
      process_state$proc$kill()
      status_message("Pipeline process stopped by user.")
    } else {
      status_message("Pipeline process is not running.")
    }

    process_state$running <- FALSE
  }

  clear_all_inputs <- function() {
    clearing_inputs(TRUE)
    initializing(TRUE)
    on.exit({
      clearing_inputs(FALSE)
      initializing(FALSE)
    }, add = TRUE)
    if (isTRUE(process_state$running) && !is.null(process_state$proc)) {
      try({
        if (process_state$proc$is_alive()) {
          process_state$proc$kill()
        }
      }, silent = TRUE)
    }
    reset_common_inputs()
    package_status("")
    refresh_package_status()
    updateTextInput(session, "external_data_path", value = "")
    updateTextInput(session, "external_metadata_path", value = "")
    updateTextInput(session, "external_reference_path", value = "")
    cfg <- read_initial_config()
    cfg <- replace_or_append(cfg, "cd_file_path", dQuote(""))
    cfg <- replace_or_append(cfg, "metadata_path", dQuote(""))
    cfg <- replace_or_append(cfg, "reference_path", dQuote(""))
    updateTextAreaInput(session, "config_text", value = cfg)
    gallery_state$dir <- NULL
    gallery_state$prefix <- NULL
    selected_result_image(NULL)
    action_confirm$ask_run <- TRUE
    action_confirm$ask_stop <- TRUE
    action_confirm$pending_run_cfg <- NULL
    process_state$log_file <- NULL
    process_state$pipeline_log_file <- NULL
    process_state$proc <- NULL
    process_state$running <- FALSE
    process_state$started_at <- NULL
    process_state$pid <- NULL
    session_started_at <<- Sys.time()
    inputs_cleared_timestamp(session_started_at)
    bump_results_gallery()
    shinyjs::runjs("window.location.hash = '#data'; setTimeout(function() { window.location.hash = '#data'; }, 100);")
    status_message("App reset to startup defaults. Ready for a new analysis.")
  }

  reset_metadata_columns <- function() {
    updateTextInput(session, "metadata_col_sample", value = "")
    updateTextInput(session, "metadata_col_weight", value = "")
    updateTextInput(session, "metadata_col_group", value = "")
    updateTextInput(session, "metadata_col_sex", value = "")
    updateTextInput(session, "metadata_col_model", value = "")
  }

  reset_reference_columns <- function() {
    updateTextInput(session, "reference_col_metabolite", value = "")
    updateTextInput(session, "reference_col_ref_ion", value = "")
    updateTextInput(session, "reference_col_mz", value = "")
    updateTextInput(session, "reference_col_rt", value = "")
  }

  resolve_output_dir_abs <- function() {
    out <- strip_outer_quotes(input$output_dir)
    if (!nzchar(out)) {
      out <- "output"
    }

    if (is_absolute_path(out)) {
      return(normalizePath(out, winslash = "/", mustWork = FALSE))
    }

    normalizePath(file.path(project_root, out), winslash = "/", mustWork = FALSE)
  }

  # Resolve uploaded, external, mapped, or config-backed input paths in priority order.
  resolve_input_file <- function(
    kind = c("data", "metadata", "reference"),
    prefer_mapped = TRUE,
    allow_config_fallback = TRUE
  ) {
    kind <- match.arg(kind)

    uploaded <- switch(kind,
      data = input$data_file,
      metadata = input$metadata_file,
      reference = input$reference_file
    )

    external <- switch(kind,
      data = safe_trimws(input$external_data_path),
      metadata = safe_trimws(input$external_metadata_path),
      reference = safe_trimws(input$external_reference_path)
    )

    cfg_key <- switch(kind,
      data = "cd_file_path",
      metadata = "metadata_path",
      reference = "reference_path"
    )

    clear_ts <- inputs_cleared_timestamp()

    # Ignore stale uploaded temp files from before the latest Clear/reset.
    is_valid_current_file <- function(file_path) {
      if (is.null(file_path) || !file.exists(file_path)) {
        return(FALSE)
      }

      file_mtime <- file.mtime(file_path)

      if (!is.null(clear_ts)) {
        return(file_mtime >= (clear_ts - 1))
      }

      file_mtime >= (session_started_at - 2)
    }

    if (identical(kind, "metadata") && isTRUE(prefer_mapped)) {
      mapping <- metadata_column_mapping()

      uploaded_name <- if (!is.null(uploaded) && !is.null(uploaded$name)) {
        uploaded$name
      } else {
        NULL
      }

      source_hint <- NULL

      if (!is.null(uploaded) &&
        !is.null(uploaded$datapath) &&
        is_valid_current_file(uploaded$datapath)) {
        source_hint <- uploaded$datapath
      } else if (nzchar(external)) {
        source_hint <- external
      } else if (isTRUE(allow_config_fallback)) {
        source_hint <- extract_config_value(input$config_text, cfg_key)
      }

      mapped_rel <- metadata_effective_rel_path(
        source_path = source_hint,
        uploaded_name = uploaded_name,
        mapping = mapping
      )

      if (!is.null(mapped_rel)) {
        mapped_abs <- file.path(project_root, mapped_rel)

        if (file.exists(mapped_abs)) {
          return(mapped_abs)
        }
      }
    }

    if (!is.null(uploaded) &&
      !is.null(uploaded$datapath) &&
      is_valid_current_file(uploaded$datapath)) {
      return(uploaded$datapath)
    }

    if (nzchar(external)) {
      if (is_absolute_path(external)) {
        return(external)
      }

      return(file.path(project_root, external))
    }

    if (isTRUE(allow_config_fallback)) {
      from_cfg <- extract_config_value(input$config_text, cfg_key)

      if (!is.null(from_cfg) && nzchar(from_cfg)) {
        if (is_absolute_path(from_cfg)) {
          return(from_cfg)
        }

        return(file.path(project_root, from_cfg))
      }
    }

    NULL
  }

  get_result_image_files <- function() {
    out_dir <- resolve_output_dir_abs()

    if (!dir.exists(out_dir)) {
      return(character(0))
    }

    files <- list.files(
      out_dir,
      pattern = "\\.(png|jpg|jpeg)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(files) == 0) {
      return(character(0))
    }

    info <- tryCatch(file.info(files), error = function(e) NULL)
    if (is.null(info) || is.null(info$mtime)) {
      return(character(0))
    }

    mtime <- info$mtime
    if (length(mtime) != length(files)) {
      return(character(0))
    }

    valid_mtime <- !is.na(mtime)
    files <- files[valid_mtime]
    mtime <- mtime[valid_mtime]

    if (length(files) == 0) {
      return(character(0))
    }

    clear_ts <- inputs_cleared_timestamp()
    min_time <- if (!is.null(clear_ts)) clear_ts - 1 else session_started_at - 2
    keep_current <- mtime >= min_time
    files <- files[keep_current]
    mtime <- mtime[keep_current]

    if (length(files) == 0) {
      return(character(0))
    }

    files[order(mtime, tolower(basename(files)), decreasing = TRUE)]
  }

  ensure_results_resource_path <- function(out_dir) {
    out_dir_norm <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

    if (!identical(gallery_state$dir, out_dir_norm) || is.null(gallery_state$prefix)) {
      gallery_state$dir <- out_dir_norm
      gallery_state$prefix <- paste0(
        "results_preview_",
        as.integer(Sys.time()),
        "_",
        sample.int(99999, 1)
      )
      addResourcePath(gallery_state$prefix, out_dir_norm)
    }

    gallery_state$prefix
  }

  build_result_image_src <- function(path) {
    out_dir <- resolve_output_dir_abs()

    if (!dir.exists(out_dir) || !file.exists(path)) {
      return(NULL)
    }

    prefix <- ensure_results_resource_path(out_dir)
    rel <- rel_path_from_output(path, out_dir)

    paste0(prefix, "/", utils::URLencode(rel, reserved = TRUE))
  }

  current_config_text <- function(fallback_text = input$config_text) {
    if (!is.null(fallback_text) && nzchar(trimws(fallback_text))) {
      return(fallback_text)
    }
    read_initial_config()
  }

  observeEvent(TRUE,
    {
      cfg <- input$config_text

      output_level_from_cfg <- extract_config_value(cfg, "output_level")
      legacy_minimal <- config_flag_value(cfg, "minimal_output", default = FALSE)
      output_level_from_cfg <- normalize_output_level(output_level_from_cfg, legacy_minimal = legacy_minimal)
      output_level_last_confirmed(output_level_from_cfg)
      updateSelectInput(session, setting_input_id("output_level"), selected = output_level_from_cfg)
    },
    once = TRUE
  )

  observe({
    refresh_package_status()
  })

  bump_results_gallery <- function() {
    gallery_refresh_tick(gallery_refresh_tick() + 1L)
  }

  observeEvent(input$save_settings_form, {
    cfg <- build_settings_builder_config(current_config_text())
    clean_config <- normalize_config_text(cfg)
    writeLines(clean_config, active_config_path, useBytes = TRUE)
    updateTextAreaInput(session, "config_text", value = clean_config)
    status_message("Settings saved from the form to config/settings.R.")
  })

  observeEvent(input$data_file,
    {
      if (is.null(input$data_file) ||
          is.null(input$data_file$name) ||
          !nzchar(input$data_file$name)) {
        return()
      }

      suggested_output <- make_output_subdir_from_data_file(input$data_file$name)
      updateTextInput(session, "output_dir", value = suggested_output)
      status_message(paste(
        "Output directory suggested from the data file; you can edit it before running:",
        suggested_output
      ))
    },
    ignoreInit = TRUE
  )

  observeEvent(input$output_dir,
    {
      output_dir_value <- sanitize_output_dir_path(input$output_dir)
      cfg <- current_config_text()
      cfg <- replace_or_append(cfg, "output_dir", dQuote(output_dir_value))
      updateTextAreaInput(session, "config_text", value = cfg)
      selected_result_image(NULL)
      bump_results_gallery()
    },
    ignoreInit = TRUE
  )

  observeEvent(input$apply_reference_cols, {
    if (!isTRUE(input$use_reference_file)) {
      status_message("Enable 'Use reference file for duplicate matching' before applying reference column names.")
      return()
    }

    cfg <- input$config_text
    cfg <- replace_or_append(cfg, "use_reference_file", "TRUE")
    cfg <- replace_or_append(cfg, "reference_col_metabolite", dQuote(safe_trimws(input$reference_col_metabolite)))
    cfg <- replace_or_append(cfg, "reference_col_ref_ion", dQuote(safe_trimws(input$reference_col_ref_ion)))
    cfg <- replace_or_append(cfg, "reference_col_mz", dQuote(safe_trimws(input$reference_col_mz)))
    cfg <- replace_or_append(cfg, "reference_col_rt", dQuote(safe_trimws(input$reference_col_rt)))

    updateTextAreaInput(session, "config_text", value = cfg)
    status_message("Reference column names applied in editor. Click 'Save config/settings.R' to persist.")
  })

  observeEvent(input$use_reference_file,
    {
      previous_use_reference_file <- isTRUE(use_reference_file_last_state())
      current_use_reference_file <- isTRUE(input$use_reference_file)
      use_reference_file_last_state(current_use_reference_file)

      if (isTRUE(input$use_reference_file)) {
        updateSelectInput(session, setting_input_id("duplicate_name_strategy"), selected = "reference_or_best_qc_rsd")
        cfg <- input$config_text
        cfg <- replace_or_append(cfg, "use_reference_file", "TRUE")
        cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote("reference_or_best_qc_rsd"))
        updateTextAreaInput(session, "config_text", value = cfg)
        if (!isTRUE(initializing()) && !previous_use_reference_file) {
          status_message("Reference file enabled: duplicate handling was automatically set to 'reference_or_best_qc_rsd'.")
        }
      }

      if (!isTRUE(input$use_reference_file)) {
        updateCheckboxInput(session, "manual_reference_cols", value = FALSE)
        updateSelectInput(session, setting_input_id("duplicate_name_strategy"), selected = "collapse_best_qc_rsd")
        cfg <- input$config_text
        cfg <- replace_or_append(cfg, "use_reference_file", "FALSE")
        cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote("collapse_best_qc_rsd"))
        cfg <- replace_or_append(cfg, "reference_path", dQuote(""))
        updateTextAreaInput(session, "config_text", value = cfg)
        if (!isTRUE(initializing()) && !isTRUE(clearing_inputs()) && previous_use_reference_file) {
          status_message("Reference file disabled: duplicate handling reset to 'collapse_best_qc_rsd'.")
        }
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input[[setting_input_id("duplicate_name_strategy")]], {
    cfg <- input$config_text
    cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote(safe_trimws(input[[setting_input_id("duplicate_name_strategy")]])))
    updateTextAreaInput(session, "config_text", value = cfg)
  }, ignoreInit = TRUE)

  observeEvent(input$manual_reference_cols,
    {
      if (!isTRUE(input$manual_reference_cols)) {
        reset_reference_columns()
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$manual_metadata_cols,
    {
      if (!isTRUE(input$manual_metadata_cols)) {
        updateCheckboxInput(session, "show_metadata_column_fields", value = FALSE)
        reset_metadata_columns()
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$show_metadata_column_fields,
    {
      if (!isTRUE(input$show_metadata_column_fields)) {
        reset_metadata_columns()
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$apply_metadata_cols,
    {
      mapping <- metadata_column_mapping()
      model_groups <- metadata_allowed_groups_by_model()
      allowed_groups <- parse_allowed_groups(input$allowed_metadata_groups)
      has_group_settings <- length(model_groups) > 0 || length(allowed_groups) >= 2

      if (!has_metadata_mapping(mapping) && !isTRUE(has_group_settings)) {
        status_message("No metadata mapping or group settings were provided.")
        return()
      }

      if (has_metadata_mapping(mapping) && isTRUE(has_group_settings)) {
        status_message("Metadata column mapping and group settings captured. They will be applied on save/run.")
      } else if (has_metadata_mapping(mapping)) {
        status_message("Metadata column mapping captured. It will be applied on save/run.")
      } else {
        status_message("Metadata group settings captured. They will be applied on save/run.")
      }
    },
    ignoreInit = TRUE
  )

  output$metadata_model_alias_ui <- renderUI({
    if (!isTRUE(input$manual_metadata_cols)) {
      return(NULL)
    }

    models <- get_detected_metadata_models()
    model_groups <- get_detected_metadata_groups_by_model()

    if (length(models) == 0) {
      return(tags$p(
        class = "small-note",
        "Upload a metadata file with a model column to show one allowed-groups box per detected model."
      ))
    }

    # The global selector uses selection order to mean control -> test. Read it
    # only while initializing the per-model fields; later clicks must not rebuild
    # those fields and silently swap their biological meaning.
    allowed_groups_initial <- isolate(input$allowed_metadata_groups)

    alias_boxes <- lapply(models, function(model_name) {
      detected_groups <- model_groups[[model_name]]
      detected_groups_text <- if (length(detected_groups) == 0) {
        "None detected"
      } else {
        paste(detected_groups, collapse = ", ")
      }
      suggested_groups_text <- metadata_groups_text_for_model(
        model_name,
        model_groups,
        fallback = allowed_groups_initial
      )

      tags$div(
        class = "model-allowed-groups-card",
        tags$h5(paste0("Model: ", model_name)),
        tags$p(
          class = "small-note",
          paste0("Groups detected in metadata: ", detected_groups_text)
        ),
        textInput(
          metadata_model_groups_input_id(model_name),
          label = "Allowed metadata groups (control first, test second)",
          value = suggested_groups_text
        )
      )
    })

    tags$div(
      style = "margin-top: 10px;",
      tags$h5("Allowed metadata groups by model"),
      tags$p(
        class = "small-note",
        "Use one box per model. The first value is treated as control and the second as test."
      ),
      tags$div(
        class = "model-allowed-groups-list",
        do.call(tagList, alias_boxes)
      )
    )
  })

  observeEvent(
    list(
      input$metadata_file,
      input$external_metadata_path,
      input$metadata_col_sample,
      input$metadata_col_weight,
      input$metadata_col_group,
      input$metadata_col_sex,
      input$metadata_col_model,
      input$manual_metadata_cols
    ),
    {
      if (!isTRUE(input$manual_metadata_cols)) {
          return()
      }

      models <- get_detected_metadata_models()
      if (length(models) == 0) {
        return()
      }

      model_groups <- get_detected_metadata_groups_by_model()
      allowed_groups_current <- isolate(input$allowed_metadata_groups)
      for (model_name in models) {
        updateTextInput(
          session,
          metadata_model_groups_input_id(model_name),
          value = metadata_groups_text_for_model(
            model_name,
            model_groups,
            fallback = allowed_groups_current
          )
        )
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$browse_output_dir, {
    if (!requireNamespace("shinyFiles", quietly = TRUE)) {
      return()
    }

    selected <- tryCatch(
      shinyFiles::parseDirPath(get_shiny_roots(), input$browse_output_dir),
      error = function(e) character(0)
    )

    if (length(selected) > 0 && nzchar(selected[1])) {
      updateTextInput(session, "output_dir", value = selected[1])
      status_message("Output directory selected from directory chooser.")
    }
  })

  observeEvent(input$save_config, {
    save_config_and_inputs(current_config_text())
    status_message("Saved config/settings.R and copied uploaded files to data/.")
  })

  observeEvent(input$open_output_dir_gallery,
    {
      out_dir <- resolve_output_dir_abs()

      if (!dir.exists(out_dir)) {
        status_message("Output directory not found yet. Run the pipeline first.")
        return()
      }

      opened <- FALSE
      try(
        {
          if (.Platform$OS.type == "windows") {
            shell.exec(out_dir)
          } else {
            utils::browseURL(out_dir)
          }
          opened <- TRUE
        },
        silent = TRUE
      )

      if (isTRUE(opened)) {
        status_message(paste("Opened output directory:", out_dir))
      } else {
        status_message("Could not open output directory automatically.")
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$refresh_results_gallery,
    {
      bump_results_gallery()
      img_files <- get_result_image_files()
      current <- selected_result_image()
      if (length(img_files) > 0 && !is_valid_image_path(current)) {
        selected_result_image(img_files[1])
      }
      status_message("Results gallery refreshed.")
    },
    ignoreInit = TRUE
  )

  observeEvent(input$install_missing_packages, {
    missing_pkgs <- find_missing_packages()

    if (length(missing_pkgs) == 0) {
      install_status_text("All required packages are already installed. Nothing to install.")
      package_status(paste(
        "All required packages are already installed:",
        paste(required_packages, collapse = ", ")
      ))
      return()
    }

    install_status_text(paste(
      "Installing",
      length(missing_pkgs),
      "package(s). Please wait..."
    ))
    package_status(paste("Installing:", paste(missing_pkgs, collapse = ", "), "..."))

    install_result <- tryCatch(
      {
        utils::install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
        TRUE
      },
      error = function(e) {
        install_status_text(paste("Package installation failed:", conditionMessage(e)))
        package_status(paste("Package installation failed:", conditionMessage(e)))
        FALSE
      }
    )

    if (isTRUE(install_result)) {
      still_missing <- find_missing_packages()
      missing_packages_state(still_missing)

      if (length(still_missing) == 0) {
        install_status_text("Package installation complete. All required packages are now available.")
        package_status(paste(
          "Package installation complete. Installed set:",
          paste(required_packages, collapse = ", ")
        ))
        status_message("Packages installed. Reloading the app to refresh the interface...")
        session$reload()
      } else {
        install_status_text(paste0(
          "Installation finished, but some packages are still missing: ",
          paste(still_missing, collapse = ", ")
        ))
        package_status(paste(
          "Installation finished, but still missing:",
          paste(still_missing, collapse = ", ")
        ))
      }
    }
  })

  observeEvent(input$run_pipeline, {
    updateTabsetPanel(session, "main_tabs", selected = "Pipeline Log")
    pipeline_log_text(paste0(
      "Run requested at ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ".\nValidating inputs before starting Rscript..."
    ))

    if (isTRUE(process_state$running)) {
      status_message("Pipeline is already running.")
      return()
    }

    data_path_for_validation <- resolve_input_file("data")
    metadata_path_for_validation <- resolve_input_file("metadata")
    reference_path_for_validation <- if (isTRUE(input$use_reference_file)) {
      resolve_input_file("reference")
    } else {
      NULL
    }

    missing_inputs <- character(0)
    if (!is_valid_file_path(data_path_for_validation)) {
      missing_inputs <- c(missing_inputs, "Data")
    }
    if (!is_valid_file_path(metadata_path_for_validation)) {
      missing_inputs <- c(missing_inputs, "Metadata")
    }
    if (isTRUE(input$use_reference_file) && !is_valid_file_path(reference_path_for_validation)) {
      missing_inputs <- c(missing_inputs, "Reference")
    }

    if (allowed_groups_missing_comma(input$allowed_metadata_groups)) {
      status_message(
        "Metadata validation failed: select at least two groups (control first, then test)."
      )
      return()
    }

    if (length(missing_inputs) > 0) {
      status_message(paste0(
        "Missing required input(s) before running the pipeline: ",
        paste(missing_inputs, collapse = ", "),
        "."
      ))
      pipeline_log_text(paste0(
        "Pipeline did not start.\nMissing required input(s): ",
        paste(missing_inputs, collapse = ", "),
        "."
      ))
      return()
    }

    md_check <- tryCatch(
      validate_metadata_columns(
        metadata_path_for_validation,
        metadata_mapping = metadata_column_mapping(),
        allowed_groups = parse_allowed_groups(input$allowed_metadata_groups),
        model_allowed_groups_by_model = metadata_allowed_groups_by_model()
      ),
      error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )

    if (!isTRUE(md_check$ok)) {
      status_message(paste("Metadata validation failed:", md_check$message))
      pipeline_log_text(paste("Pipeline did not start.\nMetadata validation failed:", md_check$message))
      return()
    }

    cfg_to_run <- build_quick_config(current_config_text())
    if (isTRUE(normalization_requires_injection_order(cfg_to_run)) && !isTRUE(has_injection_order_file_for_run())) {
      msg <- paste(
        "Injection order file is required for normalization_mode = qcrsc or qc_loess.",
        "Please upload the Injection order file before running the pipeline."
      )
      status_message(msg)
      pipeline_log_text(paste("Pipeline did not start.", msg, sep = "\n"))
      return()
    }

    if (isTRUE(action_confirm$ask_run)) {
      action_confirm$pending_run_cfg <- cfg_to_run
      pipeline_log_text(paste0(
        "Inputs validated at ",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        ".\nWaiting for run confirmation..."
      ))

      showModal(pipeline_modal(
        title = "Confirm run",
        size = "m",
        easyClose = TRUE,
        tags$p("You are about to start the pipeline with the current settings."),
        tags$ul(
          tags$li(paste("Output directory:", safe_trimws(input$output_dir))),
          tags$li(paste("Output level:", output_level_enabled())),
          tags$li(paste("Use reference file:", if (isTRUE(input$use_reference_file)) "Enabled" else "Disabled")),
          tags$li(paste("Analysis mode:", safe_trimws(input$settings_comparison_mode))),
          if (safe_trimws(input$settings_comparison_mode) %in% c("multigroup", "both")) {
            tags$li("MULTIGROUP_GLOBAL is exploratory: no directional FC, Up/Down classification, or volcano plot.")
          }
        ),
        checkboxInput("skip_run_confirm", "Do not ask again in this session", value = FALSE),
        footer = tagList(
          actionButton("cancel_run_confirm", "Cancel"),
          actionButton("confirm_run_now", "Run pipeline", class = "btn-primary")
        )
      ))
      return()
    }

    run_pipeline_now(cfg_to_run)
  })

  observeEvent(input$confirm_run_now,
    {
      removeModal()

      if (isTRUE(input$skip_run_confirm)) {
        action_confirm$ask_run <- FALSE
      }

      cfg_to_run <- action_confirm$pending_run_cfg
      if (is.null(cfg_to_run) || !nzchar(cfg_to_run)) {
        cfg_to_run <- build_quick_config(current_config_text())
      }

      if (isTRUE(normalization_requires_injection_order(cfg_to_run)) && !isTRUE(has_injection_order_file_for_run())) {
        msg <- paste(
          "Injection order file is required for normalization_mode = qcrsc or qc_loess.",
          "Please upload the Injection order file before running the pipeline."
        )
        status_message(msg)
        pipeline_log_text(paste("Pipeline did not start.", msg, sep = "\n"))
        action_confirm$pending_run_cfg <- NULL
        return()
      }

      action_confirm$pending_run_cfg <- NULL
      run_pipeline_now(cfg_to_run)
    },
    ignoreInit = TRUE
  )

  observeEvent(input$cancel_run_confirm,
    {
      removeModal()
      action_confirm$pending_run_cfg <- NULL
      status_message("Pipeline run cancelled.")
    },
    ignoreInit = TRUE
  )

  observeEvent(input$stop_pipeline, {
    if (!isTRUE(process_state$running) || is.null(process_state$proc)) {
      status_message("No running pipeline process to stop.")
      return()
    }

    if (isTRUE(action_confirm$ask_stop)) {
      showModal(pipeline_modal(
        title = "Confirm stop",
        size = "m",
        easyClose = TRUE,
        tags$p("Stopping now may leave partial outputs in the output directory."),
        tags$p("Do you want to stop the running pipeline?"),
        checkboxInput("skip_stop_confirm", "Do not ask again in this session", value = FALSE),
        footer = tagList(
          actionButton("cancel_stop_confirm", "Cancel"),
          actionButton("confirm_stop_now", "Stop pipeline", class = "btn-danger")
        )
      ))
      return()
    }

    stop_pipeline_now()
  })

  observeEvent(input$confirm_stop_now,
    {
      removeModal()

      if (isTRUE(input$skip_stop_confirm)) {
        action_confirm$ask_stop <- FALSE
      }

      stop_pipeline_now()
    },
    ignoreInit = TRUE
  )

  observeEvent(input$cancel_stop_confirm,
    {
      removeModal()
      status_message("Stop action cancelled.")
    },
    ignoreInit = TRUE
  )

  observeEvent(input$clear_all,
    {
      showModal(pipeline_modal(
        title = "Reset app?",
        size = "m",
        easyClose = TRUE,
        tags$p("This will clear uploaded files, reset settings to the app startup defaults, clear the pipeline log, and refresh the results state."),
        tags$p("This action cannot be undone. Proceed?"),
        footer = tagList(
          actionButton("cancel_clear", "Cancel"),
          actionButton("confirm_clear", "Reset app", class = "btn-warning")
        )
      ))
    },
    ignoreInit = TRUE
  )

  observeEvent(input$confirm_clear,
    {
      removeModal()
      clear_all_inputs()
    },
    ignoreInit = TRUE
  )

  observeEvent(input$cancel_clear,
    {
      removeModal()
      status_message("Clear action cancelled.")
    },
    ignoreInit = TRUE
  )

  observeEvent(input$selected_result_image_click, {
    img_path <- input$selected_result_image_click

    if (is_valid_image_path(img_path)) {
      selected_result_image(img_path)
      showModal(modalDialog(
        title = "Figure preview",
        easyClose = TRUE,
        footer = NULL,
        uiOutput("gallery_modal_content")
      ))
      status_message("Result image selected.")
    }
  })

  observeEvent(input$prev_image, {
    imgs <- get_result_image_files()
    if (length(imgs) == 0) {
      return()
    }

    current <- selected_result_image()
    idx <- match(current, imgs)
    if (is.na(idx)) {
      selected_result_image(imgs[1])
      return()
    }

    prev_idx <- if (idx <= 1) length(imgs) else idx - 1
    selected_result_image(imgs[prev_idx])
  })

  observeEvent(input$next_image, {
    imgs <- get_result_image_files()
    if (length(imgs) == 0) {
      return()
    }

    current <- selected_result_image()
    idx <- match(current, imgs)
    if (is.na(idx)) {
      selected_result_image(imgs[1])
      return()
    }

    next_idx <- if (idx >= length(imgs)) 1 else idx + 1
    selected_result_image(imgs[next_idx])
  })

  observeEvent(input$close_image_modal, {
    removeModal()
  })

  observe({
    if (!isTRUE(process_state$running)) {
      return()
    }

    invalidateLater(2000, session)

    live_log <- NULL

    if (!is.null(process_state$log_file) && file.exists(process_state$log_file)) {
      live_log <- process_state$log_file
    } else if (!is.null(process_state$pipeline_log_file) && file.exists(process_state$pipeline_log_file)) {
      live_log <- process_state$pipeline_log_file
    }

    if (!is.null(live_log)) {
      pipeline_log_text(read_log_preview(live_log))
    }

    if (!is.null(process_state$proc) &&
      !process_state$proc$is_alive()) {
      exit_code <- process_state$proc$get_exit_status()
      process_state$running <- FALSE

      if (isTRUE(as.integer(exit_code) == 0L)) {
        status_message("Pipeline run completed successfully.")
      } else {
        status_message(
          paste(
            "Pipeline run failed with exit code",
            as.integer(exit_code),
            ". Check Pipeline Log tab."
          )
        )
      }

      if (!is.null(process_state$pipeline_log_file) && file.exists(process_state$pipeline_log_file)) {
        pipeline_log_text(read_log_preview(process_state$pipeline_log_file))
      } else if (!is.null(process_state$log_file) && file.exists(process_state$log_file)) {
        pipeline_log_text(read_log_preview(process_state$log_file))
      }

      selected_result_image(NULL)
      bump_results_gallery()
    }
  })

  observe({
    pipeline_log_text()
    session$sendCustomMessage("scrollPipelineLog", list())
  })

  output$top_status_banner <- renderUI({
    tags$div(
      style = "margin-bottom:10px; padding:10px; border-radius:8px; background:#ecfeff; border:1px solid #99f6e4; color:#134e4a; display:flex; align-items:center; justify-content:space-between; gap:10px;",
      tags$div(
        style = "flex:1; text-align:center;",
        status_message()
      ),
      tags$div(
        style = "flex-shrink:0;",
        actionButton("clear_all", "Clear all", class = "btn-secondary", style = "padding:6px 12px; font-size:12px; background:#242424; border-color:#242424; color:#fff;", icon = icon("trash-alt"))
      )
    )
  })

  output$allowed_metadata_groups_hint <- renderUI({
    hint <- allowed_groups_hint_text(input$allowed_metadata_groups)
    if (!nzchar(hint)) {
      return(NULL)
    }

    tags$p(
      class = "small-note",
      style = "margin-top:6px; color:#b45309;",
      hint
    )
  })

  output$package_status <- renderText({
    package_status()
  })

  output$install_status_text <- renderText({
    install_status_text()
  })

  output$package_management_ui <- renderUI({
    missing_pkgs <- missing_packages_state()
    missing_packages_note <- if (length(missing_pkgs) > 0) {
      tags$p(
        class = "small-note",
        paste("Missing packages:", paste(missing_pkgs, collapse = ", "))
      )
    }

    tags$div(
      tags$hr(),
      h4("Package Management"),
      tags$p(
        class = "small-note",
        "The pipeline relies on several R packages. Click the button to check if all required packages are installed and install any missing ones."
      ),
      missing_packages_note,
      actionButton("install_missing_packages", "Install missing packages"),
      tags$div(
        class = "package-status-output",
        textOutput("package_status", container = span)
      ),
      tags$hr()
    )
  })

  output$output_dir_browser_ui <- renderUI({
    if (requireNamespace("shinyFiles", quietly = TRUE)) {
      shinyFiles::shinyDirButton(
        "browse_output_dir",
        "Browse...",
        "Select output directory"
      )
    } else {
      tags$p(
        class = "small-note",
        "Install package 'shinyFiles' to enable output folder browsing."
      )
    }
  })

  output$output_dir_status <- renderText({
    out <- sanitize_output_dir_path(input$output_dir)
    abs_path <- if (is_absolute_path(out)) {
      normalizePath(out, winslash = "/", mustWork = FALSE)
    } else {
      normalizePath(file.path(project_root, out), winslash = "/", mustWork = FALSE)
    }

    paste("Output will be saved to:", abs_path)
  })

  output$data_overview <- renderUI({
    # Depend on inputs_cleared_timestamp() so overview is refreshed when inputs are cleared/opened
    inputs_cleared_timestamp()

    md_path <- resolve_input_file(
      "metadata",
      allow_config_fallback = FALSE
    )

    data_path <- resolve_input_file(
      "data",
      allow_config_fallback = FALSE
    )
    allowed_groups <- unique(parse_allowed_groups(input$allowed_metadata_groups))
    allowed_groups_norm <- toupper(allowed_groups)
    metadata_mapping <- metadata_column_mapping()

    summarize_table <- function(path, missing_message, not_found_message, builder) {
      if (is.null(path) || !nzchar(path)) {
        return(list(ok = FALSE, msg = missing_message))
      }

      if (!file.exists(path)) {
        return(list(ok = FALSE, msg = not_found_message))
      }

      tryCatch(builder(path), error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    }

        md_info <- summarize_table(
      md_path,
      "No metadata file selected.",
      "Metadata file not found.",
      function(path) {
        md <- read_metadata_with_mapping(path, metadata_mapping)
        cols <- tolower(names(md))

        sample_n <- if ("sample" %in% cols) {
          sample_idx <- which(cols == "sample")
          if (length(sample_idx) > 0) {
            length(unique(na.omit(md[[sample_idx[1]]])))
          } else {
            NA_integer_
          }
        } else {
          NA_integer_
        }

        groups <- if ("group" %in% cols) {
          group_idx <- which(cols == "group")
          if (length(group_idx) > 0) {
            g <- as.character(md[[group_idx[1]]])
            g <- g[!is.na(g)]
            g <- trimws(g)
            g <- g[!is_missing_like(g)]
            sort(unique(g))
          } else {
            character(0)
          }
        } else {
          character(0)
        }

        invalid_groups <- character(0)

        if ("group" %in% cols && length(allowed_groups_norm) > 0) {
          group_idx <- which(cols == "group")
          group_raw <- trimws(as.character(md[[group_idx[1]]]))

          qc_aliases <- metadata_qc_group_aliases()
          row_is_qc <- toupper(group_raw) %in% qc_aliases

          if ("sample" %in% cols) {
            sample_idx <- which(cols == "sample")
            sample_raw <- trimws(as.character(md[[sample_idx[1]]]))
            row_is_qc <- row_is_qc | grepl("^QC", sample_raw, ignore.case = TRUE)
          }

          if (any(cols %in% c("type", "sample_type"))) {
            type_idx <- which(cols %in% c("type", "sample_type"))[1]
            type_raw <- trimws(as.character(md[[type_idx]]))
            row_is_qc <- row_is_qc | (toupper(type_raw) %in% qc_aliases)
          }

          groups_for_validation <- group_raw
          model_allowed_groups <- metadata_allowed_groups_by_model()

          if (
            length(model_allowed_groups) > 0 &&
              "model" %in% cols &&
              length(allowed_groups) >= 2
          ) {
            model_idx <- which(cols == "model")
            model_raw <- trimws(as.character(md[[model_idx[1]]]))

            valid_rows <- !is.na(groups_for_validation) &
              !is_missing_like(groups_for_validation) &
              !is.na(model_raw) &
              !is_missing_like(model_raw)

            if (any(valid_rows)) {
              groups_for_validation[valid_rows] <- normalize_model_group_pairs(
                groups_for_validation[valid_rows],
                model_raw[valid_rows],
                model_allowed_groups,
                allowed_groups[1],
                allowed_groups[2]
              )
            }
          }

          groups_for_validation <- groups_for_validation[!row_is_qc]
          groups_for_validation <- groups_for_validation[!is.na(groups_for_validation)]
          groups_for_validation <- groups_for_validation[!is_missing_like(groups_for_validation)]

          if (length(groups_for_validation) > 0) {
            groups_norm <- toupper(trimws(groups_for_validation))
            invalid_groups <- sort(unique(groups_for_validation[!(groups_norm %in% allowed_groups_norm)]))
          }
        }

        models <- if ("model" %in% cols) {
          model_idx <- which(cols == "model")
          if (length(model_idx) > 0) {
            sort(unique(na.omit(as.character(md[[model_idx[1]]]))))
          } else {
            character(0)
          }
        } else {
          character(0)
        }

        sexes <- if ("sex" %in% cols) {
          sex_idx <- which(cols == "sex")
          if (length(sex_idx) > 0) {
            table(as.character(md[[sex_idx[1]]]), useNA = "ifany")
          } else {
            NULL
          }
        } else {
          NULL
        }

        list(
          ok = TRUE,
          n_rows = nrow(md),
          n_samples = sample_n,
          groups = groups,
          invalid_groups = invalid_groups,
          models = models,
          sexes = sexes,
          path = path
        )
      }
    )

    data_info <- summarize_table(
      data_path,
      "No data file selected.",
      "Data file not found.",
      function(path) {
        dat <- safe_read_table(path)
        area_cols <- grep("^Area", names(dat), ignore.case = TRUE)

        list(
          ok = TRUE,
          n_rows = nrow(dat),
          n_cols = ncol(dat),
          n_area_cols = length(area_cols),
          path = path
        )
      }
    )

    comp_info <- if (isTRUE(input$use_reference_file)) {
      summarize_table(
        resolve_input_file("reference"),
        "No reference file selected.",
        "Reference file not found.",
        function(path) {
          comp <- safe_read_table(path)

          list(
            ok = TRUE,
            n_rows = nrow(comp),
            n_cols = ncol(comp),
            path = path
          )
        }
      )
    } else {
      list(ok = TRUE, disabled = TRUE)
    }

    tags$div(
      style = "border: 1px solid #dbe4ef; border-radius: 10px; padding: 12px; background: #fff;",
      tags$h5("Data matrix"),
      if (isTRUE(data_info$ok)) {
        tags$ul(
          tags$li(paste("File:", data_info$path)),
          tags$li(paste("Rows:", data_info$n_rows)),
          tags$li(paste("Columns:", data_info$n_cols)),
          tags$li(paste("Area columns:", data_info$n_area_cols))
        )
      } else {
        tags$p(style = "color:#b91c1c;", paste("Data summary unavailable:", data_info$msg))
      },
      tags$h5("Metadata"),
      if (isTRUE(md_info$ok)) {
        tags$ul(
          tags$li(paste("File:", md_info$path)),
          tags$li(paste("Rows:", md_info$n_rows)),
          tags$li(paste("Samples:", ifelse(is.na(md_info$n_samples), "N/A", md_info$n_samples))),
          tags$li(paste(
            "Models:",
            if (length(md_info$models) == 0) "N/A" else paste(md_info$models, collapse = ", ")
          )),
          tags$li(paste(
            "Groups:",
            if (length(md_info$groups) == 0) "N/A" else paste(md_info$groups, collapse = ", ")
          )),
          tags$li(paste(
            "Allowed groups (current setting):",
            if (length(allowed_groups) == 0) "N/A" else paste(allowed_groups, collapse = ", ")
          )),
          tags$li(paste(
            "Sex counts:",
            if (is.null(md_info$sexes)) {
              "N/A"
            } else {
              paste(names(md_info$sexes), md_info$sexes, collapse = " | ")
            }
          )),
          if (length(md_info$invalid_groups) > 0) {
            tags$li(
              style = "color:#b91c1c; font-weight:600;",
              paste(
                "Invalid group values found:",
                paste(md_info$invalid_groups, collapse = ", "),
                "- please review metadata file."
              )
            )
          }
        )
      } else {
        tags$p(style = "color:#b91c1c;", paste("Metadata summary unavailable:", md_info$msg))
      },
      tags$h5("Reference"),
      if (isTRUE(comp_info$ok) && isTRUE(comp_info$disabled)) {
        tags$p(
          class = "small-note",
          "Reference file is disabled in Inputs. Reference summary is not required for this run."
        )
      } else if (isTRUE(comp_info$ok)) {
        tags$ul(
          tags$li(paste("File:", comp_info$path)),
          tags$li(paste("Rows:", comp_info$n_rows)),
          tags$li(paste("Columns:", comp_info$n_cols))
        )
      } else {
        if (identical(comp_info$msg, "No reference file selected.")) {
          tags$div(
            style = "color:#b91c1c;", paste("Reference file not selected. Please upload a reference file in the Inputs panel.")
          )
        } else {
          tags$p(style = "color:#b91c1c;", paste("Reference summary unavailable:", comp_info$msg))
        }
      }
    )
  })

  output$qc_pca_comparison_summary <- renderTable({
    gallery_refresh_tick()

    out_dir <- resolve_output_dir_abs()
    summary_path <- file.path(out_dir, "global", "audits_global", "qc_pca_comparison_summary.csv")

    placeholder <- data.frame(
      metric = "No comparison summary available yet",
      value = "Run the pipeline first",
      unit = "",
      stage = "",
      note = "",
      stringsAsFactors = FALSE
    )

    if (!file.exists(summary_path)) {
      return(placeholder)
    }

    info <- file.info(summary_path)
    clear_ts <- inputs_cleared_timestamp()
    min_time <- if (!is.null(clear_ts)) clear_ts - 1 else session_started_at - 2
    if (is.na(info$mtime) || info$mtime < min_time) {
      return(placeholder)
    }

    tryCatch(
      {
        df <- readr::read_csv(summary_path, show_col_types = FALSE)
        if (nrow(df) == 0) placeholder else as.data.frame(df)
      },
      error = function(e) placeholder
    )
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  output$results_gallery_summary <- renderUI({
    gallery_refresh_tick()

    out_dir <- resolve_output_dir_abs()
    img_files <- get_result_image_files()
    summary_path <- file.path(out_dir, "global", "audits_global", "qc_pca_comparison_summary.csv")
    clear_ts <- inputs_cleared_timestamp()
    min_time <- if (!is.null(clear_ts)) clear_ts - 1 else session_started_at - 2
    summary_current <- FALSE
    if (file.exists(summary_path)) {
      summary_info <- file.info(summary_path)
      summary_current <- !is.na(summary_info$mtime) && summary_info$mtime >= min_time
    }

    latest_file <- ""
    latest_time <- ""
    if (length(img_files) > 0) {
      info <- file.info(img_files)
      latest_idx <- which.max(info$mtime)
      latest_file <- rel_path_from_output(img_files[latest_idx], out_dir)
      latest_time <- format(info$mtime[latest_idx], "%Y-%m-%d %H:%M:%S")
    }

    tags$div(
      class = "results-state-bar",
      tags$div(
        tags$strong("Output directory"),
        tags$span(if (dir.exists(out_dir)) out_dir else paste0(out_dir, " (not found)"))
      ),
      tags$div(
        tags$strong("Figures"),
        tags$span(as.character(length(img_files)))
      ),
      tags$div(
        tags$strong("QC/PCA table"),
        tags$span(if (isTRUE(summary_current)) "available" else "not found for current run")
      ),
      tags$div(
        tags$strong("Latest figure"),
        tags$span(if (nzchar(latest_file)) paste0(latest_file, " | ", latest_time) else "none")
      )
    )
  })

  output$results_gallery <- renderUI({
    gallery_refresh_tick()

    out_dir <- resolve_output_dir_abs()

    if (!dir.exists(out_dir)) {
      return(tags$p("Output directory not found yet. Run the pipeline first."))
    }

    img_files <- get_result_image_files()

    if (length(img_files) == 0) {
      return(tags$p("No result images found yet. Run the pipeline to generate figures."))
    }

    current <- selected_result_image()
    if (!is_valid_image_path(current)) {
      selected_result_image(img_files[1])
    }

    prefix <- ensure_results_resource_path(out_dir)

    cards <- lapply(img_files, function(img_path) {
      rel <- rel_path_from_output(img_path, out_dir)
      src <- paste0(prefix, "/", utils::URLencode(rel, reserved = TRUE))

      current_selected <- selected_result_image()
      is_active <- !is.null(current_selected) &&
        identical(
          normalizePath(img_path, winslash = "/", mustWork = FALSE),
          normalizePath(current_selected, winslash = "/", mustWork = FALSE)
        )

      tags$div(
        class = if (isTRUE(is_active)) "results-gallery-card active" else "results-gallery-card",
        tags$div(class = "results-gallery-name", rel),
        tags$a(
          href = "#",
          class = "results-gallery-thumb-link",
          onclick = sprintf(
            "Shiny.setInputValue('selected_result_image_click', '%s', {priority: 'event'}); return false;",
            gsub("'", "\\\\'", img_path)
          ),
          tags$img(src = src, class = "results-gallery-thumb")
        )
      )
    })

    tags$div(class = "results-gallery-grid", do.call(tagList, cards))
  })

  output$gallery_modal_content <- renderUI({
    img_path <- selected_result_image()
    if (!is_valid_image_path(img_path)) {
      return(tags$p("No figure selected."))
    }

    img_src <- build_result_image_src(img_path)
    out_dir <- resolve_output_dir_abs()
    rel <- rel_path_from_output(img_path, out_dir)
    imgs <- get_result_image_files()
    pos <- match(img_path, imgs)
    pos_txt <- if (is.na(pos)) "" else paste0(pos, " / ", length(imgs))

    tags$div(
      class = "gallery-modal-shell",
      tags$button(
        type = "button",
        class = "gallery-modal-close-x",
        title = "Close",
        onclick = "Shiny.setInputValue('close_image_modal', Date.now(), {priority: 'event'});",
        "\u00d7"
      ),
      tags$div(
        class = "gallery-modal-image-wrap",
        actionButton(
          "prev_image",
          "\u2039",
          class = "gallery-modal-nav gallery-modal-nav-left"
        ),
        tags$img(src = img_src, class = "gallery-modal-image"),
        actionButton(
          "next_image",
          "\u203a",
          class = "gallery-modal-nav gallery-modal-nav-right"
        )
      ),
      tags$p(class = "gallery-modal-caption", paste(rel, pos_txt)),
      tags$div(
        class = "results-download-row",
        downloadButton("download_selected_png", "Download PNG"),
        downloadButton("download_selected_jpeg", "Download JPEG")
      )
    )
  })

  output$selected_result_preview <- renderUI({
    img_path <- selected_result_image()

    if (!is_valid_image_path(img_path)) {
      return(tags$div(
        class = "results-preview-box",
        tags$p("Select an image from the gallery to preview and download it.")
      ))
    }

    img_src <- build_result_image_src(img_path)
    out_dir <- resolve_output_dir_abs()
    rel <- rel_path_from_output(img_path, out_dir)

    tags$div(
      class = "results-preview-box",
      tags$div(class = "results-preview-path", rel),
      tags$img(src = img_src, class = "results-preview-image"),
      tags$div(
        class = "results-download-row",
        downloadButton("download_selected_png", "Download PNG"),
        downloadButton("download_selected_jpeg", "Download JPEG")
      )
    )
  })

  create_image_download_handler <- function(format = "png") {
    list(
      filename = function() {
        img_path <- selected_result_image()

        if (!is_valid_image_path(img_path)) {
          return(paste0("result_image.", format))
        }

        paste0(tools::file_path_sans_ext(basename(img_path)), ".", format)
      },
      content = function(file) {
        img_path <- selected_result_image()

        if (!is_valid_image_path(img_path)) {
          stop("No image selected.")
        }

        ext <- tolower(tools::file_ext(img_path))

        if (format == "png" && ext == "png") {
          file.copy(img_path, file, overwrite = TRUE)
          return()
        } else if (format == "jpeg" && ext %in% c("jpg", "jpeg")) {
          file.copy(img_path, file, overwrite = TRUE)
          return()
        }

        if (!requireNamespace("magick", quietly = TRUE)) {
          stop(paste("Package 'magick' is required to convert images to", toupper(format)))
        }

        img <- magick::image_read(img_path)
        magick::image_write(image = img, path = file, format = format)
      }
    )
  }

  output$download_selected_png <- downloadHandler(
    filename = create_image_download_handler("png")$filename,
    content = create_image_download_handler("png")$content
  )

  output$download_selected_jpeg <- downloadHandler(
    filename = create_image_download_handler("jpeg")$filename,
    content = create_image_download_handler("jpeg")$content
  )

  output$script_content <- renderText({
    selected <- input$script_select

    if (is.null(selected) || !nzchar(selected)) {
      return("No script selected.")
    }

    idx <- which(script_names == selected)

    if (length(idx) == 0) {
      return("Selected script was not found.")
    }

    safe_read_file(script_paths[idx[1]])
  })

  output$pipeline_log <- renderText({
    pipeline_log_text()
  })


  cleanup_app_session <- function() {
    if (!is.null(process_state$proc) && isTRUE(process_state$running)) {
      try(
        {
          if (process_state$proc$is_alive()) {
            process_state$proc$kill()
          }
        },
        silent = TRUE
      )
    }

    process_state$proc <- NULL
    process_state$running <- FALSE
    process_state$pipeline_log_file <- NULL

    if (!is.null(process_state$log_file) && file.exists(process_state$log_file)) {
      try(unlink(process_state$log_file), silent = TRUE)
    }
    process_state$log_file <- NULL

    selected_result_image(NULL)
    if (!is.null(gallery_state$prefix)) {
      try(removeResourcePath(gallery_state$prefix), silent = TRUE)
      gallery_state$prefix <- NULL
      gallery_state$out_dir <- NULL
    }

    invisible(TRUE)
  }

  session$onSessionEnded(function() {
    try(cleanup_app_session(), silent = TRUE)
  })

  onStop(function() {
    try(cleanup_app_session(), silent = TRUE)
  })
}
