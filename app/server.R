server <- function(input, output, session) {
  status_message <- reactiveVal("Ready.")
  pipeline_log_text <- reactiveVal("No run executed yet.")
  package_status <- reactiveVal("")
  install_status_text <- reactiveVal("Click the button to check and install required packages.")
  missing_packages_state <- reactiveVal(setdiff(required_packages, rownames(installed.packages())))
  initializing <- reactiveVal(TRUE)
  clearing_inputs <- reactiveVal(FALSE)

  process_state <- reactiveValues(
    proc = NULL,
    log_file = NULL,
    pipeline_log_file = NULL,
    running = FALSE
  )
  minimal_output_guard <- reactiveValues(updating = FALSE)
  action_confirm <- reactiveValues(
    ask_run = TRUE,
    ask_stop = TRUE,
    pending_run_cfg = NULL
  )
  selected_result_image <- reactiveVal(NULL)
  gallery_state <- reactiveValues(dir = NULL, prefix = NULL)
  session_started_at <- Sys.time()
  inputs_cleared_timestamp <- reactiveVal(NULL)
  use_reference_file_last_state <- reactiveVal(FALSE)

  reset_common_inputs <- function() {
    shinyjs::reset("data_file")
    shinyjs::reset("metadata_file")
    shinyjs::reset("reference_file")
    updateTextInput(session, "output_dir", value = "output")
    updateSelectInput(session, "duplicate_name_strategy", selected = "collapse_best_qc_rsd")
    updateSelectInput(session, "run_metrics", selected = "FDR_and_p_value")
    updateSelectInput(session, "statistical_test_type", selected = "student")
    updateSelectInput(session, "test_is_paired", selected = "FALSE")
    updateTextInput(session, "allowed_metadata_groups", value = "WT, TG")
    updateCheckboxInput(session, "use_reference_file", value = FALSE)
    updateCheckboxInput(session, "use_weight_normalization", value = FALSE)
    updateSelectInput(session, "normalization_mode", selected = "PQN")
    updateNumericInput(session, "loess_min_qc_points", value = 5)
    updateNumericInput(session, "QC_LOESS_span", value = 0.75)
    updateCheckboxInput(session, "minimal_output", value = FALSE)
    updateCheckboxInput(session, "manual_metadata_cols", value = FALSE)
    updateCheckboxInput(session, "show_metadata_column_fields", value = FALSE)
    updateCheckboxInput(session, "manual_reference_cols", value = FALSE)
    reset_metadata_columns()
    reset_reference_columns()
    selected_result_image(NULL)
    pipeline_log_text("No run executed yet.")
    inputs_cleared_timestamp(Sys.time())
  }

  session$onFlushed(function() {
    initializing(TRUE)
    default_cfg <- if (file.exists(active_config_path)) {
      safe_read_file(active_config_path)
    } else {
      "# No default settings file found in config/."
    }
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

  observeEvent(input$metadata_file,
    {
      req(input$metadata_file)
      if (isTRUE(initializing())) {
        return()
      }

      showModal(modalDialog(
        title = "Weight normalization",
        radioButtons(
          "modal_weight_norm_choice",
          "Apply weight normalization to samples?",
          choices = c("Yes" = "yes", "No" = "no"),
          selected = if (isTRUE(setting_display_logical(initial_settings_text, "use_weight_normalization", default = FALSE))) "yes" else "no"
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

  set_minimal_output_status <- function(enabled) {
    if (isTRUE(enabled)) {
      status_message("Minimal output enabled: plots and statistics are preserved; selected intermediate global exports are skipped.")
    } else {
      status_message("Minimal output disabled: full set of outputs (including detailed plots/exports) will be generated.")
    }
  }

  render_builder_control <- function(spec) {
    spec_key <- if (!is.null(spec$key) && length(spec$key) > 0 && nzchar(as.character(spec$key)[1])) {
      as.character(spec$key)[1]
    } else {
      "unknown_key"
    }
    spec_label <- if (!is.null(spec$label) && length(spec$label) > 0 && nzchar(as.character(spec$label)[1])) {
      as.character(spec$label)[1]
    } else {
      spec_key
    }
    spec_type <- if (!is.null(spec$type) && length(spec$type) > 0 && nzchar(as.character(spec$type)[1])) {
      as.character(spec$type)[1]
    } else {
      "text"
    }
    spec_choices <- if (!is.null(spec$choices) && length(spec$choices) > 0) {
      as.character(spec$choices)
    } else {
      character(0)
    }

    input_id <- setting_input_id(spec_key)
    default_text <- initial_settings_text

    to_scalar <- function(x, fallback = "") {
      if (is.null(x) || length(x) == 0 || all(is.na(x))) {
        return(as.character(fallback))
      }
      as.character(x[1])
    }

    value <- switch(spec_type,
      checkbox = setting_display_logical(default_text, spec_key, default = isTRUE(spec$default)),
      logical_select = if (setting_display_logical(default_text, spec_key, default = isTRUE(spec$default))) "TRUE" else "FALSE",
      numeric = setting_display_numeric(default_text, spec_key, default = spec$default),
      integer = setting_display_numeric(default_text, spec_key, default = spec$default),
      multiselect = setting_default_vector(setting_display_value(default_text, spec_key, default = spec$default)),
      vector_numeric = setting_default_numeric_vector(setting_display_value(default_text, spec_key, default = spec$default)),
      vector_text = setting_default_vector(setting_display_value(default_text, spec_key, default = spec$default)),
      nullable_vector_text = setting_default_vector(setting_display_value(default_text, spec_key, default = spec$default)),
      setting_display_value(default_text, spec_key, default = spec$default)
    )

    vector_choices <- unique(c(as.character(value), as.character(spec$default), spec_choices))
    vector_choices <- vector_choices[nzchar(vector_choices)]

    control <- switch(spec_type,
      checkbox = checkboxInput(input_id, spec_label, value = isTRUE(value)),
      logical_select = selectInput(
        input_id,
        spec_label,
        choices = c("TRUE", "FALSE"),
        selected = value
      ),
      numeric = numericInput(
        input_id,
        spec_label,
        value = value,
        min = spec$min,
        max = spec$max,
        step = if (!is.null(spec$step)) spec$step else 0.1
      ),
      integer = numericInput(
        input_id,
        spec_label,
        value = value,
        min = spec$min,
        max = spec$max,
        step = if (!is.null(spec$step)) spec$step else 1
      ),
      select = selectInput(
        input_id,
        spec_label,
        choices = spec_choices,
        selected = {
          value_chr <- trimws(to_scalar(value, fallback = ""))
          if (identical(spec_key, "normalization_mode")) {
            value_norm <- toupper(value_chr)
            if (value_norm == "LOESS") value_norm <- "QC_LOESS"
            value_chr <- if (value_norm == "NONE") "none" else value_norm
          }
          if (length(spec_choices) == 0) {
            ""
          } else if (!nzchar(value_chr) || !(value_chr %in% spec_choices)) {
            as.character(spec_choices[1])
          } else {
            value_chr
          }
        }
      ),
      multiselect = selectizeInput(
        input_id,
        spec_label,
        choices = spec_choices,
        selected = intersect(as.character(value), spec_choices),
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      vector_numeric = selectizeInput(
        input_id,
        spec_label,
        choices = vector_choices,
        selected = if (length(value) == 0) character(0) else value,
        multiple = TRUE,
        options = list(
          create = TRUE,
          persist = TRUE,
          placeholder = if (!is.null(spec$placeholder)) as.character(spec$placeholder)[1] else "Add numeric values and press Enter"
        )
      ),
      vector_text = selectizeInput(
        input_id,
        spec_label,
        choices = vector_choices,
        selected = if (length(value) == 0) character(0) else value,
        multiple = TRUE,
        options = list(
          create = TRUE,
          persist = TRUE,
          placeholder = if (!is.null(spec$placeholder)) as.character(spec$placeholder)[1] else "Add values and press Enter"
        )
      ),
      nullable_vector_text = selectizeInput(
        input_id,
        spec_label,
        choices = vector_choices,
        selected = if (length(value) == 0) character(0) else value,
        multiple = TRUE,
        options = list(
          create = TRUE,
          persist = TRUE,
          placeholder = if (!is.null(spec$placeholder)) as.character(spec$placeholder)[1] else "Optional list; leave empty for NULL"
        )
      ),
      selectize_text = selectizeInput(
        input_id,
        spec_label,
        choices = vector_choices,
        selected = {
          value_chr <- trimws(to_scalar(value, fallback = ""))
          if (!nzchar(value_chr)) as.character(spec$default) else value_chr
        },
        options = list(create = TRUE, persist = TRUE)
      ),
      sheet = textInput(input_id, spec_label, value = value),
      textInput(input_id, spec_label, value = value)
    )

    if (!is.null(spec$help)) {
      tagList(control, tags$p(class = "small-note", spec$help))
    } else {
      control
    }
  }

  safe_render_builder_control <- function(spec) {
    tryCatch(
      render_builder_control(spec),
      error = function(e) {
        fallback_key <- if (!is.null(spec$key) && length(spec$key) > 0) as.character(spec$key)[1] else "unknown_key"
        fallback_label <- if (!is.null(spec$label) && length(spec$label) > 0) as.character(spec$label)[1] else fallback_key
        fallback_default <- ""
        if (!is.null(spec$default) && length(spec$default) > 0) {
          fallback_default <- as.character(spec$default)[1]
        }
        fallback_value <- setting_display_value(initial_settings_text, fallback_key, default = fallback_default)

        textInput(
          setting_input_id(fallback_key),
          fallback_label,
          value = fallback_value
        )
      }
    )
  }

  # =========================================================================
  # Helper: Build section blocks for settings UI (builder or glossary)
  # =========================================================================
  # Consolidates the repeated logic of rendering form sections into structured
  # blocks with proper grid layout and titles.
  #
  # @param sections List of section objects with title and fields
  # @return Named list of tags$div blocks, keyed by section title
  build_section_blocks_ui <- function(sections) {
    blocks <- lapply(sections, function(section) {
      section_fields <- lapply(section$fields, safe_render_builder_control)
      tags$div(
        class = "settings-section-card",
        tags$h5(section$title),
        tags$div(
          class = "settings-fields-grid",
          do.call(tagList, section_fields)
        )
      )
    })
    names(blocks) <- vapply(sections, function(section) section$title, character(1))
    blocks
  }

  # =========================================================================
  # Helper: Extract all setting keys from sections
  # =========================================================================
  # Flattens the nested section/fields structure to extract all unique
  # setting keys for validation or indexing purposes.
  #
  # @param sections List of section objects with fields
  # @return Character vector of unique setting keys
  extract_settings_keys <- function(sections) {
    unique(unlist(lapply(sections, function(section) {
      vapply(section$fields, function(spec) spec$key, character(1))
    }), use.names = FALSE))
  }

  build_settings_builder_config <- function(current_text) {
    cfg <- current_text

    for (section in settings_form_sections) {
      for (spec in section$fields) {
        input_value <- input[[setting_input_id(spec$key)]]
        replacement <- switch(spec$type,
          checkbox = setting_value_logical(input_value),
          logical_select = setting_value_logical(identical(as.character(input_value)[1], "TRUE")),
          numeric = setting_value_numeric(input_value, default = spec$default),
          integer = setting_value_integer(input_value, default = spec$default),
          select = setting_value_text(input_value),
          multiselect = setting_value_vector_text(input_value),
          vector_numeric = setting_value_vector_numeric(input_value),
          vector_text = setting_value_vector_text(input_value),
          nullable_vector_text = setting_value_vector_text(input_value, allow_null = TRUE),
          selectize_text = setting_value_text(input_value),
          sheet = setting_value_sheet(input_value),
          setting_value_text(input_value)
        )

        cfg <- replace_or_append(cfg, spec$key, replacement)
      }
    }

    cfg
  }

  output$settings_builder_ui <- renderUI({
    section_blocks <- build_section_blocks_ui(settings_form_sections)
    make_section_card <- function(section_name) {
      if (!section_name %in% names(section_blocks)) {
        return(NULL)
      }

      section_blocks[[section_name]]
    }

    left_cards <- Filter(
      Negate(is.null),
      list(
        make_section_card("Statistics thresholds"),
        make_section_card("Normalization"),
        make_section_card("Feature filtering and naming")
      )
    )

    right_cards <- Filter(
      Negate(is.null),
      list(
        make_section_card("Heatmap"),
        make_section_card("Heatmap clustering"),
        tags$div(
          class = "settings-section-card",
          tags$h5("Save settings"),
          tags$p(
            class = "small-note",
            "Persist the current form values into config/settings.R."
          ),
          actionButton(
            "save_settings_form",
            "Save config/settings.R from form"
          )
        )
      )
    )

    bslib::page_fillable(
      # tags$p(
      #   class = "small-note",
      #   "Use this form to edit the main pipeline variables."
      # ),
      bslib::layout_columns(
        tags$div(class = "settings-column-stack", do.call(tagList, left_cards)),
        tags$div(class = "settings-column-stack", do.call(tagList, right_cards)),
        col_widths = c(6, 6)
      )
    )
  })

  output$settings_glossary_ui <- renderUI({
    glossary_text_for_key <- function(key) {
      if (key %in% names(settings_glossary_map)) {
        return(as.character(settings_glossary_map[key])[1])
      }
      "Controls this pipeline behavior."
    }

    section_blocks <- build_section_blocks_ui(settings_form_sections)
    settings_keys <- extract_settings_keys(settings_form_sections)
    glossary_keys <- settings_keys

    glossary_items <- lapply(glossary_keys, function(key) {
      text <- glossary_text_for_key(key)

      tags$li(
        tags$strong(key),
        ": ",
        text
      )
    })

    tagList(
      tags$div(
        class = "settings-builder-shell",
        tags$p(
          class = "small-note",
          "This guide explains what each configuration variable does. The editable form is in the other tab."
        ),
        tags$div(
          class = "settings-guide-card",
          tags$ul(
            class = "small-note settings-guide-list",
            do.call(tagList, glossary_items)
          )
        )
      )
    )
  })

  get_shiny_roots <- function() {
    normalize_root <- function(path) {
      normalizePath(path, winslash = "/", mustWork = FALSE)
    }

    roots <- c(Project = normalize_root(project_root))

    home_dir <- tryCatch(normalize_root("~"), error = function(e) "")
    if (nzchar(home_dir) && dir.exists(home_dir)) {
      roots <- c(roots, Home = home_dir)
    }

    onedrive_dir <- safe_trimws(Sys.getenv("OneDrive"))
    if (nzchar(onedrive_dir) && dir.exists(onedrive_dir)) {
      roots <- c(roots, OneDrive = normalize_root(onedrive_dir))
    }

    drive_paths <- vapply(LETTERS, function(letter) paste0(letter, ":/"), character(1))
    drive_paths <- drive_paths[dir.exists(drive_paths)]

    if (length(drive_paths) > 0) {
      drive_names <- paste0("Drive ", sub(":/$", ":", drive_paths))
      names(drive_paths) <- drive_names
      roots <- c(roots, drive_paths)
    }

    roots <- roots[!duplicated(roots)]
    roots <- roots[!is.na(roots) & nzchar(trimws(roots))]
    roots
  }

  if (requireNamespace("shinyFiles", quietly = TRUE)) {
    shinyFiles::shinyDirChoose(
      input,
      "browse_output_dir",
      roots = get_shiny_roots(),
      session = session
    )
  }

  find_missing_packages <- function() {
    installed <- rownames(installed.packages())
    setdiff(required_packages, installed)
  }

  replace_or_append <- function(text, key, value_expr) {
    pattern <- paste0("^\\s*", key, "\\s*<-")
    replacement <- paste0(key, " <- ", value_expr)
    lines <- strsplit(text, "\\n", fixed = FALSE)[[1]]

    idx <- grep(pattern, lines)
    if (length(idx) > 0) {
      lines[idx[1]] <- replacement
    } else {
      lines <- c(lines, replacement)
    }

    paste(lines, collapse = "\n")
  }


  parse_allowed_groups <- function(value) {
    raw <- safe_trimws(value)
    if (!nzchar(raw)) {
      return(character(0))
    }

    vals <- unlist(strsplit(raw, ",", fixed = TRUE), use.names = FALSE)
    vals <- trimws(vals)
    unique(vals[nzchar(vals)])
  }

  allowed_groups_hint_text <- function(value) {
    raw <- safe_trimws(value)
    if (!nzchar(raw)) {
      return("")
    }

    if (!grepl(",", raw, fixed = TRUE)) {
      return("Missing comma: separate the control and test group names with a comma, for example WT, TG.")
    }

    parsed <- parse_allowed_groups(raw)
    if (length(parsed) < 2) {
      return("Please provide two group names in the order control, test.")
    }

    ""
  }

  allowed_groups_missing_comma <- function(value) {
    raw <- safe_trimws(value)
    nzchar(raw) && !grepl(",", raw, fixed = TRUE)
  }

  make_output_subdir_from_data_file <- function(file_name) {
    stem <- tools::file_path_sans_ext(basename(file_name))
    stem <- gsub("[^A-Za-z0-9._-]+", "_", stem)
    stem <- gsub("^_+|_+$", "", stem)

    if (!nzchar(stem)) {
      stem <- format(Sys.time(), "run_%Y%m%d_%H%M%S")
    }

    file.path("output", stem)
  }

  metadata_column_mapping <- function() {
    list(
      sample = safe_trimws(input$metadata_col_sample),
      weight = safe_trimws(input$metadata_col_weight),
      group = safe_trimws(input$metadata_col_group),
      sex = safe_trimws(input$metadata_col_sex),
      model = safe_trimws(input$metadata_col_model)
    )
  }

  sanitize_input_id_fragment <- function(x) {
    x <- safe_trimws(x)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) {
      x <- "model"
    }
    x
  }

  metadata_model_groups_input_id <- function(model_name) {
    paste0("metadata_model_groups_", sanitize_input_id_fragment(model_name))
  }

  get_detected_metadata_models <- function() {
    md_path <- resolve_input_file("metadata")
    if (!is_valid_file_path(md_path)) {
      return(character(0))
    }

    mapping <- metadata_column_mapping()
    md <- tryCatch(
      read_metadata_with_mapping(md_path, mapping),
      error = function(e) NULL
    )

    if (is.null(md) || !("model" %in% names(md))) {
      return(character(0))
    }

    models <- unique(trimws(as.character(md$model)))
    models <- models[!is.na(models) & nzchar(models)]
    models <- models[!toupper(models) %in% c("NA", "N/A", "NULL")]
    sort(models)
  }

  get_detected_metadata_groups_by_model <- function() {
    md_path <- resolve_input_file("metadata")
    if (!is_valid_file_path(md_path)) {
      return(list())
    }

    mapping <- metadata_column_mapping()
    md <- tryCatch(
      read_metadata_with_mapping(md_path, mapping),
      error = function(e) NULL
    )

    if (is.null(md) || !("model" %in% names(md)) || !("group" %in% names(md))) {
      return(list())
    }

    md$model <- trimws(as.character(md$model))
    md$group <- trimws(as.character(md$group))

    keep <- !is.na(md$model) & nzchar(md$model) & !is_missing_like(md$model)
    md <- md[keep, , drop = FALSE]

    if (nrow(md) == 0) {
      return(list())
    }

    split_groups <- split(md$group, md$model)
    lapply(split_groups, function(values) {
      values <- unique(values[!is.na(values) & nzchar(values) & !is_missing_like(values)])
      sort(values)
    })
  }

  metadata_allowed_groups_by_model <- function() {
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

  format_named_character_vector <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return("NULL")
    }

    entries <- vapply(names(x), function(nm) {
      paste0(dQuote(nm), " = ", dQuote(as.character(x[[nm]])))
    }, character(1), USE.NAMES = FALSE)

    paste0("c(", paste(entries, collapse = ", "), ")")
  }

  has_metadata_mapping <- function(mapping) {
    any(vapply(mapping, nzchar, logical(1)))
  }

  metadata_mapped_rel_path <- function(uploaded_name) {
    stem <- tools::file_path_sans_ext(basename(uploaded_name))
    stem <- gsub("[^A-Za-z0-9._-]+", "_", stem)
    stem <- gsub("^_+|_+$", "", stem)
    if (!nzchar(stem)) {
      stem <- "metadata"
    }
    file.path("data", paste0(stem, "_mapped.csv"))
  }

  metadata_effective_rel_path <- function(source_path = NULL, uploaded_name = NULL, mapping = metadata_column_mapping()) {
    if (!has_metadata_mapping(mapping)) {
      return(NULL)
    }

    source_name <- uploaded_name
    if (!nzchar(safe_trimws(source_name))) {
      source_name <- basename(safe_trimws(source_path))
    }
    if (!nzchar(safe_trimws(source_name))) {
      source_name <- "metadata"
    }

    metadata_mapped_rel_path(source_name)
  }

  apply_metadata_mapping_to_df <- function(df, mapping) {
    if (!has_metadata_mapping(mapping)) {
      return(df)
    }

    actual <- names(df)
    actual_norm <- tolower(trimws(as.character(actual)))

    for (target in names(mapping)) {
      src <- tolower(trimws(as.character(mapping[[target]])))
      if (!nzchar(src)) {
        next
      }

      idx <- which(actual_norm == src)
      if (length(idx) == 0) {
        next
      }

      names(df)[idx[1]] <- target
      actual <- names(df)
      actual_norm <- tolower(trimws(as.character(actual)))
    }

    df
  }

  read_metadata_with_mapping <- function(path, mapping = metadata_column_mapping()) {
    md <- safe_read_table(path)
    apply_metadata_mapping_to_df(md, mapping)
  }

  persist_metadata_mapping <- function(source_path, uploaded_name = NULL, mapping = metadata_column_mapping()) {
    mapped_rel <- metadata_effective_rel_path(
      source_path = source_path,
      uploaded_name = uploaded_name,
      mapping = mapping
    )

    if (is.null(mapped_rel)) {
      return(NULL)
    }

    md <- read_metadata_with_mapping(source_path, mapping)
    utils::write.csv(md, file.path(project_root, mapped_rel), row.names = FALSE, na = "")
    mapped_rel
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
    cfg <- current_text
    allowed_groups <- parse_allowed_groups(input$allowed_metadata_groups)
    duplicate_strategy_effective <- safe_trimws(as.character(input$duplicate_name_strategy)[1])
    if (!nzchar(duplicate_strategy_effective)) {
      duplicate_strategy_effective <- "collapse_best_qc_rsd"
    }

    cfg <- replace_or_append(cfg, "output_dir", dQuote(input$output_dir))
    cfg <- replace_or_append(cfg, "use_weight_normalization", if (isTRUE(input$use_weight_normalization)) "TRUE" else "FALSE")
    normalization_mode_effective <- safe_trimws(input$normalization_mode)
    if (!normalization_mode_effective %in% c("none", "PQN", "QC_LOESS")) {
      normalization_mode_effective <- "PQN"
    }
    cfg <- replace_or_append(cfg, "normalization_mode", dQuote(normalization_mode_effective))
    cfg <- replace_or_append(cfg, "loess_min_qc_points", setting_value_integer(input$loess_min_qc_points, default = 5))
    cfg <- replace_or_append(cfg, "QC_LOESS_span", setting_value_numeric(input$QC_LOESS_span, default = 0.75))
    cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote(duplicate_strategy_effective))
    cfg <- replace_or_append(cfg, "p_value_cutoff", setting_value_numeric(input[[setting_input_id("p_value_cutoff")]], default = 0.05))
    cfg <- replace_or_append(cfg, "fdr_cutoff", setting_value_numeric(input[[setting_input_id("fdr_cutoff")]], default = 0.05))
    cfg <- replace_or_append(cfg, "run_metrics", dQuote(input$run_metrics))
    cfg <- replace_or_append(cfg, "heatmap_rank_metrics", dQuote(input$run_metrics))
    cfg <- replace_or_append(cfg, "statistical_test_type", dQuote(input$statistical_test_type))
    cfg <- replace_or_append(cfg, "test_is_paired", if (identical(as.character(input$test_is_paired)[1], "TRUE")) "TRUE" else "FALSE")
    # Derive p-value correction method from run_metrics selection so that
    # `run_metrics` becomes the primary choice for significance metric.
    # If the user selected any option including "FDR", default correction is "FDR" (BH),
    # otherwise use "raw" (no correction) to match p-value usage.
    derived_pvalue_correction <- if (grepl("FDR", as.character(input$run_metrics), ignore.case = TRUE)) {
      "FDR"
    } else {
      "raw"
    }
    cfg <- replace_or_append(cfg, "pvalue_correction_method", dQuote(derived_pvalue_correction))
    minimal_flag <- if (isTRUE(input$minimal_output)) "TRUE" else "FALSE"
    cfg <- replace_or_append(cfg, "minimal_output", minimal_flag)
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
    derived_alpha_sig <- if (grepl("FDR", as.character(input$run_metrics), ignore.case = TRUE)) {
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

  normalized_output_dir <- function() {
    out <- safe_trimws(input$output_dir)
    if (!nzchar(out)) {
      return("output")
    }
    out
  }

  build_expected_output_manifest <- function(minimal_mode = isTRUE(input$minimal_output)) {
    out_dir <- normalized_output_dir()
    cfg <- input$config_text
    strategy <- extract_config_value(cfg, "duplicate_name_strategy")
    strategy <- safe_trimws(if (is.null(strategy)) "" else strategy)

    volcano_enabled <- config_flag_value(cfg, "make_volcano_plots", default = TRUE)
    heatmap_all <- config_flag_value(cfg, "make_heatmap_by_model", default = TRUE)
    heatmap_sex <- config_flag_value(cfg, "make_heatmap_by_model_sex", default = TRUE)
    sig_heatmap_all <- config_flag_value(cfg, "make_sig_heatmap_by_model", default = TRUE)
    sig_heatmap_sex <- config_flag_value(cfg, "make_sig_heatmap_by_model_sex", default = TRUE)
    sig_heatmap_fvsm <- config_flag_value(cfg, "make_sig_heatmap_FvsM_within_group", default = TRUE)
    sig_metabolites_txt <- config_flag_value(cfg, "save_sig_metabolites_txt_per_model", default = TRUE)

    if (!nzchar(strategy)) {
      strategy <- "collapse_best_qc_rsd"
    }

    core_files <- c(
      file.path(out_dir, "PIPELINE_LOG.txt"),
      file.path(out_dir, "global", "audits_global", "filter_summary.csv"),
      file.path(out_dir, "global", "audits_global", "missing_exclusion_audit.csv"),
      file.path(out_dir, "global", "audits_global", "presence_filter_audit.csv"),
      file.path(out_dir, "global", "audits_global", "known_filter_audit.csv"),
      file.path(out_dir, "global", "audits_global", "qc_rsd_values_pre_variants.csv"),
      file.path(out_dir, "global", "audits_global", "low_variance_iqr_audit_ACTIVE.csv"),
      file.path(out_dir, "global", "audits_global", paste0("duplicate_name_audit_", strategy, ".csv")),
      file.path(out_dir, "global", "exports_global", "02_featureID_to_display_name_map.csv")
    )

    if (identical(strategy, "reference_or_best_qc_rsd")) {
      core_files <- c(
        core_files,
        file.path(
          out_dir,
          "global",
          "audits_global",
          "duplicate_name_reference_summary_reference_or_best_qc_rsd.csv"
        )
      )
    }

    if (!isTRUE(minimal_mode)) {
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", "06_pqn_factors_weight_then_PQN.csv"),
        file.path(out_dir, "global", "exports_global", "10_MATRIX_post_<ACTIVE_VARIANT>_postRSD_preLowVar_preDup_ALL.csv"),
        file.path(out_dir, "global", "exports_global", "10_TABLE_post_<ACTIVE_VARIANT>_postRSD_preLowVar_preDup_ALL_NAMED.csv")
      )
    }

    optional_files <- character(0)

    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "exports", "stats", "*.xlsx")
    )

    if (isTRUE(sig_metabolites_txt)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "exports", "stats", "significant_metabolites", "p_value", "*.txt"),
        file.path(out_dir, "<MODEL>", "exports", "stats", "significant_metabolites", "FDR", "*.txt")
      )
    }

    if (isTRUE(volcano_enabled)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "volcano", "*.png")
      )
    }

    if (isTRUE(heatmap_all)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "heatmap", "*.png")
      )
    }

    if (isTRUE(heatmap_sex)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "heatmap", "*.png")
      )
    }

    if (isTRUE(sig_heatmap_all)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
      )
    }

    if (isTRUE(sig_heatmap_sex)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
      )
    }

    if (isTRUE(sig_heatmap_fvsm)) {
      optional_files <- c(
        optional_files,
        file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
      )
    }

    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "pca", "*.png"),
      file.path(out_dir, "<MODEL>", "exports", "metaboanalyst", "*.csv")
    )

    list(
      out_dir = out_dir,
      core_files = unique(core_files),
      optional_files = unique(optional_files)
    )
  }

  render_manifest_list <- function(items) {
    tags$ul(lapply(items, function(item) tags$li(item)))
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
    pipeline_log_text("Starting pipeline...")

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
        pipeline_log_text(safe_read_file(output_log_file))
      } else if (file.exists(log_file)) {
        pipeline_log_text(safe_read_file(log_file))
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
    on.exit({
      clearing_inputs(FALSE)
    }, add = TRUE)
    reset_common_inputs()
    updateTextInput(session, "external_data_path", value = "")
    updateTextInput(session, "external_metadata_path", value = "")
    updateTextInput(session, "external_reference_path", value = "")
    cfg <- if (file.exists(active_config_path)) {
      safe_read_file(active_config_path)
    } else {
      read_initial_config()
    }
    cfg <- replace_or_append(cfg, "cd_file_path", dQuote(""))
    cfg <- replace_or_append(cfg, "metadata_path", dQuote(""))
    cfg <- replace_or_append(cfg, "reference_path", dQuote(""))
    updateTextAreaInput(session, "config_text", value = cfg)
    gallery_state$dir <- NULL
    gallery_state$prefix <- NULL
    process_state$log_file <- NULL
    process_state$pipeline_log_file <- NULL
    process_state$proc <- NULL
    process_state$running <- FALSE
    session_started_at <<- Sys.time()
    shinyjs::runjs("window.location.hash = '#data'; setTimeout(function() { window.location.hash = '#data'; }, 100);")
    status_message("All inputs cleared. Ready for a new analysis.")
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

  is_valid_file_path <- function(path) {
    !is.null(path) && nzchar(as.character(path)) && file.exists(as.character(path))
  }

  is_valid_image_path <- function(path) {
    !is.null(path) && nzchar(as.character(path)) && file.exists(as.character(path))
  }

  resolve_output_dir_abs <- function() {
    out <- trimws(input$output_dir)
    if (!nzchar(out)) {
      out <- "output"
    }

    if (is_absolute_path(out)) {
      return(normalizePath(out, winslash = "/", mustWork = FALSE))
    }

    normalizePath(file.path(project_root, out), winslash = "/", mustWork = FALSE)
  }

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

    info <- file.info(files)
    keep <- !is.na(info$mtime) & info$mtime >= session_started_at
    files <- files[keep]

    files[order(tolower(basename(files)))]
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

  rel_path_from_output <- function(path, out_dir) {
    path_norm <- gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
    out_norm <- gsub("\\\\", "/", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
    sub(paste0("^", out_norm, "/?"), "", path_norm)
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

  config_flag_value <- function(text, key, default = FALSE) {
    val <- extract_config_value(text, key)
    if (is.null(val) || !nzchar(trimws(as.character(val)))) {
      return(isTRUE(default))
    }

    val_norm <- toupper(trimws(as.character(val)))
    if (val_norm %in% c("TRUE", "T", "1", "YES")) {
      return(TRUE)
    }
    if (val_norm %in% c("FALSE", "F", "0", "NO")) {
      return(FALSE)
    }

    isTRUE(default)
  }

  current_config_text <- function(fallback_text = input$config_text) {
    if (file.exists(active_config_path)) {
      return(safe_read_file(active_config_path))
    }
    fallback_text
  }

  observeEvent(TRUE,
    {
      cfg <- input$config_text

      minimal_from_cfg <- config_flag_value(cfg, "minimal_output", default = FALSE)
      if (!identical(isTRUE(input$minimal_output), isTRUE(minimal_from_cfg))) {
        minimal_output_guard$updating <- TRUE
        updateCheckboxInput(session, "minimal_output", value = minimal_from_cfg)
      } else {
        minimal_output_guard$updating <- FALSE
      }
    },
    once = TRUE
  )

  observe({
    missing_pkgs <- find_missing_packages()
    missing_packages_state(missing_pkgs)

    if (length(missing_pkgs) == 0) {
      install_status_text("All required packages are installed.")
    } else {
      install_status_text(paste0(
        "Missing packages detected (",
        length(missing_pkgs),
        "): ",
        paste(missing_pkgs, collapse = ", ")
      ))
    }
  })

  observe({
    invalidateLater(2000, session)

    img_files <- get_result_image_files()

    if (length(img_files) > 0) {
      current <- selected_result_image()

      if (!is_valid_image_path(current)) {
        selected_result_image(img_files[1])
      }
    }
  })

  observeEvent(input$save_settings_form, {
    cfg <- build_settings_builder_config(current_config_text())
    clean_config <- normalize_config_text(cfg)
    writeLines(clean_config, active_config_path, useBytes = TRUE)
    updateTextAreaInput(session, "config_text", value = clean_config)
    status_message("Settings saved from the form to config/settings.R.")
  })

  observeEvent(input$data_file,
    {
      if (is.null(input$data_file) || is.null(input$data_file$name) || !nzchar(input$data_file$name)) {
        return()
      }

      suggested_output <- make_output_subdir_from_data_file(input$data_file$name)
      updateTextInput(session, "output_dir", value = suggested_output)
      status_message(paste("Output directory auto-updated based on data file:", suggested_output))
    },
    ignoreInit = TRUE
  )

  observeEvent(input$minimal_output,
    {
      if (isTRUE(minimal_output_guard$updating)) {
        minimal_output_guard$updating <- FALSE
        return()
      }

      if (isTRUE(input$minimal_output)) {
        manifest_min <- build_expected_output_manifest(minimal_mode = TRUE)
        manifest_full <- build_expected_output_manifest(minimal_mode = FALSE)
        skipped_files <- setdiff(manifest_full$core_files, manifest_min$core_files)

        showModal(modalDialog(
          title = "Confirm minimal output",
          size = "l",
          easyClose = TRUE,
          tags$p("Minimal output keeps your plots and statistics options. It only skips selected intermediate global exports."),
          tags$p(tags$strong("Output directory:"), paste0(" ", manifest_min$out_dir)),
          tags$h5("Files generated in minimal mode"),
          render_manifest_list(manifest_min$core_files),
          tags$h5("Files skipped when minimal mode is enabled"),
          if (length(skipped_files) > 0) {
            render_manifest_list(skipped_files)
          } else {
            tags$p("No additional core files are skipped with the current configuration.")
          },
          tags$p(
            class = "small-note",
            "Per-model plots and stats remain enabled according to your selected options."
          ),
          footer = tagList(
            actionButton("cancel_minimal_output", "Cancel"),
            actionButton("confirm_minimal_output", "Continue with minimal output", class = "btn-primary")
          )
        ))
      } else {
        set_minimal_output_status(FALSE)
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(input$confirm_minimal_output,
    {
      removeModal()
      set_minimal_output_status(TRUE)
    },
    ignoreInit = TRUE
  )

  observeEvent(input$cancel_minimal_output,
    {
      removeModal()
      minimal_output_guard$updating <- TRUE
      updateCheckboxInput(session, "minimal_output", value = FALSE)
      set_minimal_output_status(FALSE)
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
        updateSelectInput(session, "duplicate_name_strategy", selected = "reference_or_best_qc_rsd")
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
        updateSelectInput(session, "duplicate_name_strategy", selected = "collapse_best_qc_rsd")
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

  observeEvent(input$duplicate_name_strategy, {
    cfg <- input$config_text
    cfg <- replace_or_append(cfg, "duplicate_name_strategy", dQuote(safe_trimws(as.character(input$duplicate_name_strategy)[1])))
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
      if (!has_metadata_mapping(mapping)) {
        status_message("No metadata mapping values were provided.")
        return()
      }

      status_message("Metadata column mapping captured. It will be applied on save/run.")
    },
    ignoreInit = TRUE
  )

  output$metadata_model_alias_ui <- renderUI({
    models <- get_detected_metadata_models()
    model_groups <- get_detected_metadata_groups_by_model()

    if (length(models) == 0) {
      return(tags$p(
        class = "small-note",
        "Upload a metadata file with a model column to show one allowed-groups box per detected model."
      ))
    }

    alias_boxes <- lapply(models, function(model_name) {
      detected_groups <- model_groups[[model_name]]
      detected_groups_text <- if (length(detected_groups) == 0) {
        "None detected"
      } else {
        paste(detected_groups, collapse = ", ")
      }

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
          value = safe_trimws(input$allowed_metadata_groups)
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
      status_message("Output directory selected from browser.")
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
        "Metadata validation failed: please separate the control and test group names with a comma (for example, WT, TG)."
      )
      return()
    }

    if (length(missing_inputs) > 0) {
      status_message(paste0(
        "Missing required input(s) before running the pipeline: ",
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
      return()
    }

    cfg_to_run <- build_quick_config(current_config_text())

    if (isTRUE(action_confirm$ask_run)) {
      action_confirm$pending_run_cfg <- cfg_to_run

      showModal(modalDialog(
        title = "Confirm run",
        size = "m",
        easyClose = TRUE,
        tags$p("You are about to start the pipeline with the current settings."),
        tags$ul(
          tags$li(paste("Output directory:", safe_trimws(input$output_dir))),
          tags$li(paste("Minimal output:", if (isTRUE(input$minimal_output)) "Enabled" else "Disabled")),
          tags$li(paste("Use reference file:", if (isTRUE(input$use_reference_file)) "Enabled" else "Disabled"))
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
      showModal(modalDialog(
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
      showModal(modalDialog(
        title = "Clear all inputs?",
        size = "m",
        easyClose = TRUE,
        tags$p("This will clear all uploaded files, configuration fields, and the pipeline log."),
        tags$p("This action cannot be undone. Proceed?"),
        footer = tagList(
          actionButton("cancel_clear", "Cancel"),
          actionButton("confirm_clear", "Clear everything", class = "btn-warning")
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
    invalidateLater(700, session)

    live_log <- NULL

    if (!is.null(process_state$pipeline_log_file) && file.exists(process_state$pipeline_log_file)) {
      live_log <- process_state$pipeline_log_file
    } else if (!is.null(process_state$log_file) && file.exists(process_state$log_file)) {
      live_log <- process_state$log_file
    }

    if (!is.null(live_log)) {
      pipeline_log_text(safe_read_file(live_log))
    }

    if (isTRUE(process_state$running) &&
      !is.null(process_state$proc) &&
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
        pipeline_log_text(safe_read_file(process_state$pipeline_log_file))
      } else if (!is.null(process_state$log_file) && file.exists(process_state$log_file)) {
        pipeline_log_text(safe_read_file(process_state$log_file))
      }
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
        actionButton("clear_all", "Clear all", class = "btn-secondary", style = "padding:6px 12px; font-size:12px; background:#242424; border-color:#242424; color:#fff;")
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

    if (length(missing_pkgs) == 0) {
      return(NULL)
    }

    tags$div(
      tags$hr(),
      h4("Package Management"),
      tags$p(
        class = "small-note",
        "The pipeline relies on several R packages. Click the button to check if all required packages are installed and install any missing ones."
      ),
      tags$p(
        class = "small-note",
        paste("Missing packages:", paste(missing_pkgs, collapse = ", "))
      ),
      actionButton("install_missing_packages", "Install missing packages"),
      verbatimTextOutput("package_status"),
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
        if (length(groups) > 0 && length(allowed_groups_norm) > 0) {
          groups_norm <- toupper(trimws(groups))
          invalid_groups <- sort(unique(groups[!(groups_norm %in% allowed_groups_norm)]))
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

  output$results_gallery <- renderUI({
    invalidateLater(2000, session)

    out_dir <- resolve_output_dir_abs()

    if (!dir.exists(out_dir)) {
      return(tags$p("Output directory not found yet. Run the pipeline first."))
    }

    img_files <- get_result_image_files()

    if (length(img_files) == 0) {
      return(tags$p("No result images found yet. Run the pipeline to generate figures."))
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
      tags$div(class = "results-preview-title", "Selected Figure"),
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
