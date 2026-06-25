# Helpers for Data Overview and Pre-flight UI rendering.

preflight_status_item <- function(label, state, detail) {
  state <- match.arg(state, c("ok", "warning", "missing"))
  status_label <- switch(state,
    ok = "OK",
    warning = "Warning",
    missing = "Missing"
  )

  tags$div(
    class = paste("preflight-item", paste0("preflight-item-", state)),
    tags$span(class = "preflight-status", status_label),
    tags$div(
      class = "preflight-copy",
      tags$strong(label),
      tags$span(detail)
    )
  )
}

build_preflight_panel_ui <- function(data_info, md_info, comp_info, allowed_groups,
                                     use_reference_file = FALSE,
                                     requires_injection_order = FALSE,
                                     injection_order_ok = TRUE,
                                     output_dir_abs = "") {
  tags$div(
    class = "preflight-panel",
    tags$h5("Pre-flight check"),
    tags$div(
      class = "preflight-grid",
      preflight_status_item(
        "Data file",
        if (isTRUE(data_info$ok)) "ok" else "missing",
        if (isTRUE(data_info$ok)) data_info$path else data_info$msg
      ),
      preflight_status_item(
        "Metadata file",
        if (isTRUE(md_info$ok)) "ok" else "missing",
        if (isTRUE(md_info$ok)) md_info$path else md_info$msg
      ),
      preflight_status_item(
        "Allowed groups",
        if (length(allowed_groups) >= 2 && length(md_info$invalid_groups) == 0) {
          "ok"
        } else if (length(allowed_groups) >= 2) {
          "warning"
        } else {
          "missing"
        },
        if (length(allowed_groups) >= 2) {
          paste(allowed_groups, collapse = ", ")
        } else {
          "Select at least two groups: control first, then test."
        }
      ),
      preflight_status_item(
        "Reference file",
        if (!isTRUE(use_reference_file) || isTRUE(comp_info$ok)) "ok" else "missing",
        if (!isTRUE(use_reference_file)) {
          "Disabled for this run."
        } else if (isTRUE(comp_info$ok)) {
          comp_info$path
        } else {
          comp_info$msg
        }
      ),
      preflight_status_item(
        "Injection order",
        if (isTRUE(injection_order_ok)) "ok" else "missing",
        if (!isTRUE(requires_injection_order)) {
          "Not required for the selected normalization mode."
        } else if (isTRUE(injection_order_ok)) {
          "Available for QC-LOESS/QC-RSC normalization."
        } else {
          "Required for QC-LOESS or QC-RSC normalization."
        }
      ),
      preflight_status_item(
        "Output directory",
        if (nzchar(output_dir_abs)) "ok" else "missing",
        if (nzchar(output_dir_abs)) output_dir_abs else "No output directory resolved."
      )
    )
  )
}

build_data_matrix_overview_ui <- function(data_info) {
  tagList(
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
    }
  )
}

build_metadata_overview_ui <- function(md_info, allowed_groups) {
  tagList(
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
    }
  )
}

build_reference_overview_ui <- function(comp_info) {
  tagList(
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
          style = "color:#b91c1c;",
          paste("Reference file not selected. Please upload a reference file in the Inputs panel.")
        )
      } else {
        tags$p(style = "color:#b91c1c;", paste("Reference summary unavailable:", comp_info$msg))
      }
    }
  )
}

build_data_overview_ui <- function(data_info, md_info, comp_info, allowed_groups,
                                   use_reference_file = FALSE,
                                   requires_injection_order = FALSE,
                                   injection_order_ok = TRUE,
                                   output_dir_abs = "") {
  tags$div(
    style = "border: 1px solid #dbe4ef; border-radius: 10px; padding: 12px; background: #fff;",
    build_preflight_panel_ui(
      data_info = data_info,
      md_info = md_info,
      comp_info = comp_info,
      allowed_groups = allowed_groups,
      use_reference_file = use_reference_file,
      requires_injection_order = requires_injection_order,
      injection_order_ok = injection_order_ok,
      output_dir_abs = output_dir_abs
    ),
    build_data_matrix_overview_ui(data_info),
    build_metadata_overview_ui(md_info, allowed_groups),
    build_reference_overview_ui(comp_info)
  )
}
