# Server-side helpers for the Settings Builder UI and config text generation.

render_builder_control <- function(spec, dynamic_choices = list()) {
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
    choices <- as.character(spec$choices)
    names(choices) <- names(spec$choices)
    choices
  } else {
    character(0)
  }

  if (spec_key %in% names(dynamic_choices)) {
    dynamic_spec_choices <- dynamic_choices[[spec_key]]
    dynamic_spec_choices <- as.character(dynamic_spec_choices)
    dynamic_spec_choices <- dynamic_spec_choices[!is.na(dynamic_spec_choices) & nzchar(dynamic_spec_choices)]
    spec_choices <- unique(c(spec_choices, dynamic_spec_choices))
  }

  input_id <- setting_input_id(spec_key)
  to_scalar <- function(x, fallback = "") {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(as.character(fallback))
    }
    as.character(x[1])
  }

  value <- switch(spec_type,
    checkbox = initial_setting_logical(spec_key, default = isTRUE(spec$default)),
    logical_select = if (initial_setting_logical(spec_key, default = isTRUE(spec$default))) "TRUE" else "FALSE",
    numeric = initial_setting_numeric(spec_key, default = spec$default),
    integer = initial_setting_numeric(spec_key, default = spec$default),
    multiselect = setting_default_vector(initial_setting_value(spec_key, default = spec$default)),
    vector_numeric = setting_default_numeric_vector(initial_setting_value(spec_key, default = spec$default)),
    vector_text = setting_default_vector(initial_setting_value(spec_key, default = spec$default)),
    detected_multiselect = setting_default_vector(initial_setting_value(spec_key, default = spec$default)),
    nullable_vector_text = setting_default_vector(initial_setting_value(spec_key, default = spec$default)),
    initial_setting_value(spec_key, default = spec$default)
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
          value_chr <- normalize_normalization_mode(value_chr, default = "none")
          if (identical(value_chr, "weight")) {
            value_chr <- "none"
          }
        }
        if (identical(spec_key, "active_variant") && identical(value_chr, "RSD_20")) {
          value_chr <- "RSD20"
        }
        if (identical(spec_key, "active_variant") && grepl("^QC_RSD[0-9.]+$", toupper(value_chr))) {
          value_chr <- "QC_RSD"
        }
        if (identical(spec_key, "active_variant") && grepl("^RSD[0-9.]+$", toupper(value_chr))) {
          value_chr <- "RSD"
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
    detected_multiselect = checkboxGroupInput(
      input_id,
      spec_label,
      choices = vector_choices,
      selected = if (length(value) == 0) character(0) else value,
      inline = TRUE
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

  control_with_help <- if (!is.null(spec$help)) {
    tagList(control, tags$p(class = "small-note", spec$help))
  } else {
    control
  }

  if (!is.null(spec$condition) && nzchar(as.character(spec$condition)[1])) {
    conditionalPanel(
      condition = as.character(spec$condition)[1],
      control_with_help
    )
  } else {
    control_with_help
  }
}

safe_render_builder_control <- function(spec, dynamic_choices = list()) {
  tryCatch(
    render_builder_control(spec, dynamic_choices = dynamic_choices),
    error = function(e) {
      fallback_key <- if (!is.null(spec$key) && length(spec$key) > 0) as.character(spec$key)[1] else "unknown_key"
      fallback_label <- if (!is.null(spec$label) && length(spec$label) > 0) as.character(spec$label)[1] else fallback_key
      fallback_default <- ""
      if (!is.null(spec$default) && length(spec$default) > 0) {
        fallback_default <- as.character(spec$default)[1]
      }
      fallback_value <- initial_setting_value(fallback_key, default = fallback_default)

      textInput(
        setting_input_id(fallback_key),
        fallback_label,
        value = fallback_value
      )
    }
  )
}

build_section_blocks_ui <- function(sections, dynamic_choices = list()) {
  section_field_columns <- function(section_title) {
    configured <- settings_builder_layout$field_columns
    columns <- 2L

    if (!is.null(configured) && length(configured) > 0) {
      if (section_title %in% names(configured)) {
        columns <- suppressWarnings(as.integer(configured[[section_title]]))
      } else if ("default" %in% names(configured)) {
        columns <- suppressWarnings(as.integer(configured[["default"]]))
      }
    }

    if (is.na(columns) || !columns %in% c(1L, 2L, 3L, 4L)) {
      columns <- 2L
    }

    columns
  }

  blocks <- lapply(sections, function(section) {
    section_fields <- lapply(section$fields, safe_render_builder_control, dynamic_choices = dynamic_choices)
    field_columns <- section_field_columns(section$title)
    block <- tags$div(
      class = "settings-section-card",
      tags$h5(section$title),
      section$note,
      tags$div(
        class = paste("settings-fields-grid", paste0("settings-fields-grid-", field_columns)),
        do.call(tagList, section_fields)
      )
    )

    if (!is.null(section$condition) && nzchar(as.character(section$condition)[1])) {
      conditionalPanel(condition = as.character(section$condition)[1], block)
    } else {
      block
    }
  })
  names(blocks) <- vapply(sections, function(section) section$title, character(1))
  blocks
}

extract_settings_keys <- function(sections) {
  unique(unlist(lapply(sections, function(section) {
    vapply(section$fields, function(spec) spec$key, character(1))
  }), use.names = FALSE))
}

apply_active_variant_config <- function(cfg, active_variant_value, rsd_threshold_value = NULL) {
  active_variant_value <- safe_trimws(active_variant_value)
  active_variant_upper <- toupper(active_variant_value)
  active_variant_upper <- sub("^QC_RSD[0-9.]+$", "QC_RSD", active_variant_upper)
  active_variant_upper <- sub("^RSD[0-9.]+$", "RSD", active_variant_upper)

  rsd_threshold <- suppressWarnings(as.numeric(as.character(rsd_threshold_value)[1]))
  if (length(rsd_threshold) == 0 || is.na(rsd_threshold) || rsd_threshold < 0) {
    rsd_threshold <- 20
  }
  rsd_threshold_expr <- paste0("c(", format(rsd_threshold, scientific = FALSE, trim = TRUE), ")")

  if (!nzchar(active_variant_value) || !active_variant_upper %in% c("NONE", "QC_RSD", "RSD")) {
    active_variant_value <- "none"
    active_variant_upper <- "NONE"
  }

  cfg <- replace_or_append(cfg, "rsd_thresholds", rsd_threshold_expr)

  if (identical(active_variant_upper, "RSD")) {
    cfg <- replace_or_append(cfg, "active_variant", dQuote("RSD"))
    cfg <- replace_or_append(cfg, "rsd_filter_type", dQuote("RSD"))
  } else if (identical(active_variant_upper, "QC_RSD")) {
    cfg <- replace_or_append(cfg, "active_variant", dQuote("QC_RSD"))
    cfg <- replace_or_append(cfg, "rsd_filter_type", dQuote("QC_RSD"))
  } else {
    cfg <- replace_or_append(cfg, "active_variant", dQuote("none"))
  }

  cfg
}

build_settings_builder_ui <- function(dynamic_choices = list()) {
  section_blocks <- build_section_blocks_ui(settings_form_sections, dynamic_choices = dynamic_choices)
  layout_tab_keys <- setdiff(names(settings_builder_layout), "field_columns")

  make_section_card <- function(section_name) {
    if (!section_name %in% names(section_blocks)) {
      return(NULL)
    }

    section_blocks[[section_name]]
  }

  make_cards <- function(...) {
    Filter(Negate(is.null), list(...))
  }

  layout_config <- function(tab_key) {
    config <- settings_builder_layout[[tab_key]]

    if (is.null(config)) {
      return(list(
        label = tab_key,
        sections = character(0),
        widths = NULL,
        narrow = FALSE
      ))
    }

    sections <- config$sections
    if (is.null(sections) || length(sections) == 0) {
      sections <- character(0)
    }

    widths <- config$widths
    if (is.null(widths) || length(widths) == 0) {
      widths <- NULL
    }

    label <- config$label
    if (is.null(label) || length(label) == 0 || !nzchar(as.character(label)[1])) {
      label <- tab_key
    }

    list(
      label = as.character(label)[1],
      sections = as.character(sections),
      widths = widths,
      narrow = isTRUE(config$narrow)
    )
  }

  cards_from_sections <- function(section_names) {
    do.call(make_cards, lapply(section_names, make_section_card))
  }

  make_columns <- function(cards, widths = NULL) {
    if (length(cards) == 0) {
      return(NULL)
    }

    if (!is.null(widths)) {
      widths <- suppressWarnings(as.integer(widths))
    }

    if (is.null(widths) || length(widths) != length(cards) || any(is.na(widths))) {
      widths <- rep(floor(12 / length(cards)), length(cards))
    }

    tags$div(
      class = "settings-layout-grid",
      do.call(
        tagList,
        Map(function(card, width) {
          width <- max(1L, min(12L, as.integer(width)))
          tags$div(
            class = "settings-column-stack",
            style = paste0("grid-column: span ", width, ";"),
            card
          )
        }, cards, widths)
      )
    )
  }

  tab_layouts <- lapply(layout_tab_keys, layout_config)

  save_bar <- tags$div(
    class = "settings-action-bar",
    tags$div(
      class = "settings-action-copy",
      tags$strong("Settings form"),
      tags$span("Current values are applied when running the pipeline.")
    ),
    actionButton(
      "save_settings_form", "Save config/settings.R from form", icon("save"), style = "background-color: #007bff; color: white;")
  )

  settings_tabs <- do.call(
    tabsetPanel,
    c(
      list(id = "settings_builder_group", type = "tabs"),
      lapply(tab_layouts, function(layout) {
        tabPanel(
          layout$label,
          tags$div(
            class = if (isTRUE(layout$narrow)) "settings-tab-pane settings-tab-pane-narrow" else "settings-tab-pane",
            make_columns(cards_from_sections(layout$sections), widths = layout$widths)
          )
        )
      })
    )
  )

  tags$div(
    class = "settings-builder-shell",
    save_bar,
    settings_tabs
  )
}

build_settings_glossary_ui <- function() {
  layout_tab_keys <- setdiff(names(settings_builder_layout), "field_columns")
  layout_tab_keys <- setdiff(layout_tab_keys, "exports")

  glossary_text_for_key <- function(key) {
    if (key %in% names(settings_glossary_map)) {
      return(as.character(settings_glossary_map[key])[1])
    }
    "Controls this pipeline behavior."
  }

  section_by_title <- setNames(settings_form_sections, vapply(settings_form_sections, function(section) section$title, character(1)))

  make_glossary_section <- function(section_name, width = 12L) {
    section <- section_by_title[[section_name]]
    if (is.null(section)) {
      return(NULL)
    }

    glossary_items <- lapply(section$fields, function(spec) {
      key <- spec$key
      text <- glossary_text_for_key(key)
      label <- if (!is.null(spec$label) && nzchar(as.character(spec$label)[1])) {
        as.character(spec$label)[1]
      } else {
        key
      }

      tags$li(
        tags$strong(label),
        tags$span(class = "settings-guide-key", paste0(" (", key, ")")),
        ": ",
        text
      )
    })

    guide_columns <- if (identical(section$title, "Output controls")) {
      1L
    } else if (isTRUE(width >= 12L) && length(glossary_items) >= 4L) {
      2L
    } else {
      1L
    }

    tags$div(
      class = "settings-guide-card",
      tags$h5(section$title),
      tags$ul(
        class = paste("small-note settings-guide-list", paste0("settings-guide-list-", guide_columns)),
        do.call(tagList, glossary_items)
      )
    )
  }

  make_guide_group <- function(tab_key) {
    layout <- settings_builder_layout[[tab_key]]
    title <- layout_label(tab_key)
    section_names <- layout$sections
    if (is.null(section_names) || length(section_names) == 0) {
      section_names <- character(0)
    }
    section_names <- as.character(section_names)
    if (length(section_names) == 0) {
      return(NULL)
    }

    widths <- suppressWarnings(as.integer(layout$widths))
    if (length(widths) != length(section_names) || any(is.na(widths))) {
      widths <- rep(floor(12 / length(section_names)), length(section_names))
    }
    widths <- pmax(1L, pmin(12L, as.integer(widths)))

    cards <- Map(make_glossary_section, section_names, widths)
    keep_cards <- !vapply(cards, is.null, logical(1))
    cards <- cards[keep_cards]
    widths <- widths[keep_cards]
    if (length(cards) == 0) {
      return(NULL)
    }

    tags$div(
      class = "settings-guide-group",
      tags$h4(title),
      tags$div(
        class = "settings-layout-grid settings-glossary-grid",
        do.call(
          tagList,
          Map(function(card, width) {
            tags$div(
              class = "settings-column-stack settings-glossary-column",
              style = paste0("grid-column: span ", width, ";"),
              card
            )
          }, cards, widths)
        )
      )
    )
  }

  layout_label <- function(tab_key) {
    label <- settings_builder_layout[[tab_key]]$label
    if (is.null(label) || length(label) == 0 || !nzchar(as.character(label)[1])) {
      return(tab_key)
    }
    as.character(label)[1]
  }

  glossary_groups <- Filter(
    Negate(is.null),
    lapply(layout_tab_keys, function(tab_key) {
      make_guide_group(tab_key)
    })
  )

  tagList(
    tags$div(
      class = "settings-builder-shell",
      tags$p(
        class = "small-note",
        "This guide follows the same tabs and cards as the editable settings form."
      ),
      do.call(tagList, glossary_groups)
    )
  )
}
