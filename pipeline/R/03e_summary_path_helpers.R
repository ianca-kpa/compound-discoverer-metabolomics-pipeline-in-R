# =============================================================================
# 03e_summary_path_helpers.R
# Summary and path helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Filter summary helper
# -----------------------------------------------------------------------------
append_filter_summary <- function(summary_tbl, step, before, after, out_csv = NULL) {
  row <- tibble(
    step = step,
    n_features_before = as.integer(before),
    n_features_after = as.integer(after),
    n_removed = as.integer(before - after),
    pct_removed = round(100 * (before - after) / max(1, before), 2)
  )
  
  summary_tbl <- bind_rows(summary_tbl, row)
  
  if (!is.null(out_csv)) {
    readr::write_csv(summary_tbl, out_csv, na = "")
  }
  
  summary_tbl
}

format_qc_rsd_threshold_summary_lines <- function(qc_rsd_threshold_summary = NULL) {
  if (is.null(qc_rsd_threshold_summary) ||
      !is.data.frame(qc_rsd_threshold_summary) ||
      nrow(qc_rsd_threshold_summary) == 0) {
    return(character(0))
  }

  required_cols <- c("rsd_filter_type", "rsd_threshold", "kept")
  if (!all(required_cols %in% names(qc_rsd_threshold_summary))) {
    return(character(0))
  }

  qc_rsd_threshold_summary <- qc_rsd_threshold_summary %>%
    dplyr::arrange(rsd_filter_type, rsd_threshold)

  unlist(lapply(seq_len(nrow(qc_rsd_threshold_summary)), function(i) {
    row <- qc_rsd_threshold_summary[i, , drop = FALSE]
    effective_suffix <- ""
    if ("rsd_threshold_effective" %in% names(row) &&
        !is.na(row$rsd_threshold_effective) &&
        !identical(as.numeric(row$rsd_threshold), as.numeric(row$rsd_threshold_effective))) {
      effective_suffix <- paste0(
        " (effective cutoff: ",
        format(row$rsd_threshold_effective, scientific = FALSE, trim = TRUE),
        ")"
      )
    }
    paste0(
      row$rsd_filter_type,
      " <= ",
      format(row$rsd_threshold, scientific = FALSE, trim = TRUE),
      effective_suffix,
      ": ",
      row$kept
    )
  }), use.names = FALSE)
}

write_method_summary <- function(path,
                                 filter_summary = NULL,
                                 injection_order_source = NULL,
                                 qc_rsd_threshold_summary = NULL) {
  value_or_na <- function(expr) {
    value <- tryCatch(force(expr), error = function(e) NULL)
    if (is.null(value) || length(value) == 0 || all(is.na(value))) {
      return("not available")
    }
    if (is.logical(value)) {
      return(ifelse(isTRUE(value[1]), "TRUE", "FALSE"))
    }
    paste(as.character(value), collapse = ", ")
  }

  feature_count_lines <- "feature_counts: not available"
  if (!is.null(filter_summary) && is.data.frame(filter_summary) && nrow(filter_summary) > 0) {
    required_cols <- c("step", "n_features_before", "n_features_after", "n_removed", "pct_removed")
    if (all(required_cols %in% names(filter_summary))) {
      feature_count_lines <- unlist(lapply(seq_len(nrow(filter_summary)), function(i) {
        row <- filter_summary[i, , drop = FALSE]
        c(
          paste0("feature_counts.", row$step, ".before: ", value_or_na(row$n_features_before)),
          paste0("feature_counts.", row$step, ".after: ", value_or_na(row$n_features_after)),
          paste0("feature_counts.", row$step, ".removed: ", value_or_na(row$n_removed)),
          paste0("feature_counts.", row$step, ".pct_removed: ", value_or_na(row$pct_removed))
        )
      }), use.names = FALSE)
    }
  }

  rsd_threshold_lines <- format_qc_rsd_threshold_summary_lines(qc_rsd_threshold_summary)
  if (length(rsd_threshold_lines) > 0) {
    rsd_threshold_lines <- c("rsd_threshold_summary:", rsd_threshold_lines)
  }

  injection_order_source <- value_or_na(injection_order_source)
  real_injection_order_used <- if (grepl("^input_files_reference", injection_order_source)) {
    "TRUE"
  } else if (identical(injection_order_source, "not available")) {
    "not available"
  } else {
    "FALSE"
  }

  lines <- c(
    "method_summary",
    paste0("generated_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("normalization_mode: ", value_or_na(get0("normalization_mode", ifnotfound = NULL, inherits = TRUE))),
    paste0("use_weight_normalization: ", value_or_na(get0("use_weight_normalization", ifnotfound = NULL, inherits = TRUE))),
    paste0("rsd_filter_type: ", value_or_na(get0("rsd_filter_type", ifnotfound = NULL, inherits = TRUE))),
    paste0("rsd_thresholds: ", value_or_na(get0("rsd_thresholds", ifnotfound = NULL, inherits = TRUE))),
    paste0("active_variant: ", value_or_na(get0("active_variant", ifnotfound = NULL, inherits = TRUE))),
    paste0("low_variance_filter_method: ", value_or_na(get0("low_variance_filter_method", ifnotfound = NULL, inherits = TRUE))),
    paste0("low_variance_filter_fraction: ", value_or_na(get0("low_variance_filter_fraction", ifnotfound = NULL, inherits = TRUE))),
    paste0("duplicate_name_strategy: ", value_or_na(get0("duplicate_name_strategy", ifnotfound = NULL, inherits = TRUE))),
    paste0("log2_offset: ", value_or_na(get0("log2_offset", ifnotfound = NULL, inherits = TRUE))),
    paste0("pca_scaling: ", value_or_na(get0("pca_scaling", ifnotfound = NULL, inherits = TRUE))),
    paste0("real_injection_order_used: ", real_injection_order_used),
    paste0("injection_order_source: ", injection_order_source),
    rsd_threshold_lines,
    feature_count_lines
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

compact_path_component <- function(value, max_chars = 32L, fallback = "item") {
  value <- gsub("[^A-Za-z0-9_-]+", "_", as.character(value)[1])
  value <- gsub("_+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) value <- fallback

  max_chars <- max(12L, as.integer(max_chars)[1])
  if (nchar(value) <= max_chars) return(value)

  code_points <- utf8ToInt(value)
  checksum <- sum((as.numeric(code_points) * seq_along(code_points)) %% 10000019) %% 100000000
  suffix <- sprintf("_%08d", checksum)
  paste0(substr(value, 1, max_chars - nchar(suffix)), suffix)
}

make_compact_output_filename <- function(..., ext, max_chars = 96L, component_chars = 24L) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts)]
  parts <- vapply(parts, compact_path_component, character(1), max_chars = component_chars)
  stem <- paste(parts[nzchar(parts)], collapse = "_")
  stem <- compact_path_component(stem, max_chars = max_chars, fallback = "output")
  paste0(stem, ".", gsub("^\\.", "", ext))
}

# -----------------------------------------------------------------------------
# Directory helpers
# -----------------------------------------------------------------------------
normalize_model_names <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_squish(x)
  x[x %in% c("", "NA", "N/A")] <- NA_character_
  x
}

# Extract unique model names from metadata, ensuring they are cleaned and valid
get_models_from_metadata <- function(metadata_clean) {
  if (is.null(metadata_clean) || nrow(metadata_clean) == 0) {
    stop("metadata_clean is empty.")
  }
  
  if (!("model" %in% names(metadata_clean))) {
    stop("metadata_clean must contain a 'model' column.")
  }
  
  if (!("type" %in% names(metadata_clean))) {
    stop("metadata_clean must contain a 'type' column.")
  }
  
  models <- metadata_clean %>%
    dplyr::mutate(model = normalize_model_names(model)) %>%
    dplyr::filter(type == "Sample", !is.na(model), model != "") %>%
    dplyr::distinct(model) %>%
    dplyr::pull(model) %>%
    sort()
  
  if (length(models) == 0) {
    stop("No valid models found in metadata.")
  }
  
  models
}

# Create a structured list of paths for a given model
make_one_model_paths <- function(output_dir, model_name) {
  model_root <- file.path(output_dir, model_name)
  
  list(
    root = model_root,
    exports = list(
      root = model_root,
      stats = file.path(model_root, "stats")
    ),
    plots = list(
      root = file.path(model_root, "plots"),
      heatmap_global = file.path(model_root, "plots", "heatmap"),
      heatmap_all = file.path(model_root, "plots", "heatmap"),
      heatmap_by_sex = file.path(model_root, "plots", "heatmap"),
      heatmap_by_group = file.path(model_root, "plots", "heatmap"),
      heatmap_tg_vs_wt = file.path(model_root, "plots", "heatmap"),
      heatmap_tg_f_vs_tg_m = file.path(model_root, "plots", "heatmap"),
      heatmap_wt_f_vs_wt_m = file.path(model_root, "plots", "heatmap"),
      heatmap_significant_global = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_all = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_by_sex = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_by_group = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_tg_vs_wt = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_tg_f_vs_tg_m = file.path(model_root, "plots", "heatmap_significant"),
      heatmap_significant_wt_f_vs_wt_m = file.path(model_root, "plots", "heatmap_significant"),
      volcano_global = file.path(model_root, "plots", "volcano"),
      volcano_all = file.path(model_root, "plots", "volcano"),
      volcano_by_sex = file.path(model_root, "plots", "volcano"),
      volcano_by_group = file.path(model_root, "plots", "volcano"),
      volcano_tg_vs_wt = file.path(model_root, "plots", "volcano"),
      volcano_tg_f_vs_tg_m = file.path(model_root, "plots", "volcano"),
      volcano_wt_f_vs_wt_m = file.path(model_root, "plots", "volcano"),
      pca_global = file.path(model_root, "plots", "pca"),
      pca_all = file.path(model_root, "plots", "pca"),
      pca_by_sex = file.path(model_root, "plots", "pca"),
      pca_by_group = file.path(model_root, "plots", "pca"),
      pca_tg_vs_wt = file.path(model_root, "plots", "pca"),
      pca_tg_f_vs_tg_m = file.path(model_root, "plots", "pca"),
      pca_wt_f_vs_wt_m = file.path(model_root, "plots", "pca")
    )
  )
}

# Create model directories if they don't exist
create_model_dirs <- function(mp) {
  # Deliberately create only the model root here.
  # Subfolders are created on demand when a file is actually written.
  dir.create(mp$root, showWarnings = FALSE, recursive = TRUE)
}

setup_output_dirs <- function(output_dir, model_names = NULL) {
  if (!dir.exists(dirname(output_dir))) {
    stop("Parent dir does not exist: ", dirname(output_dir))
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  dirs <- list(
    root = output_dir,
    global = list(
      root = file.path(output_dir, "global"),
      exports = file.path(output_dir, "global", "exports_global"),
      audits = file.path(output_dir, "global", "audits_global")
    ),
    models = list()
  )
  
  if (!is.null(model_names)) {
    model_names <- normalize_model_names(model_names)
    model_names <- unique(model_names[!is.na(model_names) & model_names != ""])
    model_names <- sort(model_names)
    
    for (m in model_names) {
      mp <- make_one_model_paths(output_dir = output_dir, model_name = m)
      create_model_dirs(mp)
      dirs$models[[m]] <- mp
    }
  }
  
  dirs
}

setup_model_output_dirs <- function(paths, metadata_clean) {
  models <- get_models_from_metadata(metadata_clean)
  
  for (m in models) {
    mp <- make_one_model_paths(output_dir = paths$root, model_name = m)
    create_model_dirs(mp)
    paths$models[[m]] <- mp
  }
  
  paths
}

get_model_paths <- function(paths, model_name) {
  model_name <- normalize_model_names(model_name)
  
  if (length(model_name) != 1 || is.na(model_name) || !nzchar(model_name)) {
    stop("Invalid model_name supplied to get_model_paths().")
  }
  
  if (is.null(paths$models) || length(paths$models) == 0) {
    stop("No model path structure found inside 'paths$models'.")
  }
  
  if (!is.null(paths$models[[model_name]])) {
    return(paths$models[[model_name]])
  }
  
  available_models <- names(paths$models)
  
  if (is.null(available_models) || length(available_models) == 0) {
    stop("Model path structure exists, but no models were registered.")
  }
  
  stop(
    "Model path structure not found for model: ", model_name,
    ". Available models: ", paste(available_models, collapse = ", ")
  )
}

# Remove empty directories under a given root directory (non-recursive, only immediate subdirectories)
remove_empty_directories <- function(root_dir) {  
  if (!dir.exists(root_dir)) return(invisible(FALSE))
  
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  all_dirs <- all_dirs[order(nchar(all_dirs), decreasing = TRUE)]
  
  for (d in all_dirs) {
    if (!dir.exists(d)) next
    contents <- list.files(d, all.files = TRUE, no.. = TRUE)
    if (length(contents) == 0) {
      try(unlink(d, recursive = FALSE, force = TRUE), silent = TRUE)
    }
  }

  invisible(TRUE)
}
