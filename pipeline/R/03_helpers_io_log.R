# =============================================================================
# 03_helpers_io_log.R
# General helpers: logging / IO / text / directories
# =============================================================================

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
append_log_line <- function(log_path, line) {
  write(line, file = log_path, append = TRUE)
}

# Format dimensions of a data frame, matrix, or vector for logging
fmt_dims <- function(x) {
  if (is.null(x)) return("rows=NA cols=NA")
  if (is.data.frame(x)) return(paste0("rows=", nrow(x), " cols=", ncol(x)))
  if (is.matrix(x)) return(paste0("rows=", nrow(x), " cols=", ncol(x)))
  if (is.vector(x)) return(paste0("len=", length(x)))
  
  tryCatch(
    paste0("rows=", nrow(x), " cols=", ncol(x)),
    error = function(e) "rows=NA cols=NA"
  )
}

# Log a message about a written object, including its dimensions and file path
log_written_object <- function(log_path, file_path, object, note = NULL) {
  msg <- paste0("- ", basename(file_path), " -> ", fmt_dims(object), " | path: ", file_path)
  
  if (!is.null(note) && nzchar(note)) {
    msg <- paste0(msg, " | note: ", note)
  }
  
  append_log_line(log_path, msg)
}

# Run an expression while capturing all console output (messages, warnings, errors) to a log file
with_console_capture_to_file <- function(log_path, expr) {
  con <- file(log_path, open = "a", encoding = "UTF-8")
  
  on.exit({
    while (sink.number() > 0) {
      try(sink(), silent = TRUE)
    }
    try(close(con), silent = TRUE)
  }, add = TRUE)
  
  sink(con, append = TRUE, split = TRUE)
  
  tryCatch(
    withCallingHandlers(
      expr,
      message = function(m) {
        # Emit messages through stdout so sink(split=TRUE) mirrors to console + log once.
        cat(conditionMessage(m), sep = "")
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        cat("[WARNING] ", conditionMessage(w), sep = "")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat("[ERROR] ", conditionMessage(e), sep = "")
      stop(e)
    }
  )
}

# -----------------------------------------------------------------------------
# File read/write helpers
# -----------------------------------------------------------------------------
write_csv_safe <- function(df, path) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    stop("write_csv_safe() received an invalid path.")
  }
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path, na = "")
}

read_any_table <- function(path, sheet = 1) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx","xls","xlsm")) return(readxl::read_excel(path, sheet = sheet))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv","txt")) return(readr::read_tsv(path, show_col_types = FALSE))
  
  stop("Unsupported file extension: ", ext)
}

resolve_model_group_values <- function(model_name = NULL,
                                      control_label = NULL,
                                      treatment_label = NULL) {
  if (is.null(control_label)) {
    control_label <- get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE)
  }
  if (is.null(treatment_label)) {
    treatment_label <- get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)
  }

  labels <- list(
    control = as.character(control_label)[1],
    treatment = as.character(treatment_label)[1]
  )

  model_name <- trimws(as.character(model_name)[1])
  if (!nzchar(model_name)) {
    return(labels)
  }

  model_groups <- get0("model_allowed_groups_by_model", ifnotfound = NULL, inherits = TRUE)
  if (is.null(model_groups) || length(model_groups) == 0 || is.null(names(model_groups))) {
    return(labels)
  }

  idx <- match(model_name, trimws(names(model_groups)))
  if (is.na(idx)) {
    return(labels)
  }

  raw <- unlist(strsplit(as.character(model_groups[[idx]]), ",", fixed = TRUE), use.names = FALSE)
  raw <- trimws(raw)
  raw <- raw[nzchar(raw)]

  if (length(raw) >= 1) {
    labels$control <- raw[1]
  }
  if (length(raw) >= 2) {
    labels$treatment <- raw[2]
  } else if (length(raw) == 1) {
    labels$treatment <- raw[1]
  }

  labels
}

# Alias for compatibility (same as resolve_model_group_values)
get_comparison_group_labels_for_model <- function(model_name = NULL,
                                                 control_label = get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE),
                                                 treatment_label = get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)) {
  resolve_model_group_values(model_name, control_label, treatment_label)
}

# Map comparison group display values: converts global group names to display names if configured per-model
map_comparison_group_display_values <- function(values, model_name = NULL) {
  labels <- get_comparison_group_labels_for_model(model_name)
  control_global <- get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE)
  treatment_global <- get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)

  out <- as.character(values)
  out[out == control_global] <- labels$control
  out[out == treatment_global] <- labels$treatment
  out
}

# -----------------------------------------------------------------------------
# Text cleaning / sanitization helpers
# -----------------------------------------------------------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_trim(x)
  
  x_low <- tolower(x)
  
  x[x_low %in% c("not named","notnamed","unnamed","no name","noname")] <- NA_character_
  x[x %in% c("", "NA","N/A","n/a","-","Unknown","unknown","No results","no results")] <- NA_character_
  
  x
}

strip_v_suffix_end <- function(x) {
  x <- as.character(x)
  stringr::str_replace(x, "_v\\d+$", "")
}

sanitize_text_for_exports <- function(x, mode = c("greek_latin_ascii","ascii_translit")) {
  mode <- match.arg(mode)
  
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_squish(x)
  
  if (all(is.na(x))) return(x)
  
  if (mode == "greek_latin_ascii") {
    x <- stringi::stri_trans_general(x, "Greek-Latin; Latin-ASCII")
  } else {
    x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  }
  
  x <- stringr::str_replace_all(x, "[[:cntrl:]]+", "")
  x
}

# Normalize simple names (lowercase trim)
normalize_name <- function(x) {
  tolower(trimws(as.character(x)))
}

# Normalize model group pairs: map per-model alias values to control/treatment labels
normalize_model_group_pairs <- function(groups_vec, model_vec, pair_map, control_label, treatment_label) {
  if (is.null(pair_map) || length(pair_map) == 0) {
    return(groups_vec)
  }

  if (is.null(model_vec) || length(model_vec) == 0) {
    return(groups_vec)
  }

  out <- as.character(groups_vec)
  model_vec <- trimws(as.character(model_vec))
  model_keys <- names(pair_map)

  if (is.null(model_keys) || length(model_keys) == 0) {
    return(groups_vec)
  }

  for (model_name in model_keys) {
    pair_raw <- as.character(pair_map[[model_name]])
    pair_vals <- unlist(strsplit(pair_raw, ",", fixed = TRUE), use.names = FALSE)
    pair_vals <- trimws(pair_vals)
    pair_vals <- pair_vals[nzchar(pair_vals)]

    if (length(pair_vals) == 0) {
      next
    }

    model_idx <- !is.na(model_vec) & trimws(model_vec) == trimws(model_name)
    if (!any(model_idx)) {
      next
    }

    group_vals <- trimws(out[model_idx])
    group_norm <- toupper(group_vals)
    if (length(pair_vals) >= 1) {
      control_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[1])]
      if (length(control_idx) > 0) {
        out[control_idx] <- control_label
      }
    }

    if (length(pair_vals) >= 2) {
      treatment_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[2])]
      if (length(treatment_idx) > 0) {
        out[treatment_idx] <- treatment_label
      }
    } else {
      treatment_idx <- which(model_idx)[group_norm %in% toupper(pair_vals[1])]
      if (length(treatment_idx) > 0) {
        out[treatment_idx] <- treatment_label
      }
    }
  }

  out
}

normalize_config_text <- function(text) {
  text <- gsub("\u201C|\u201D", '"', text, perl = TRUE)
  text <- gsub("\u2018|\u2019", "'", text, perl = TRUE)
  text
}

# -----------------------------------------------------------------------------
# Helper: detect missing-like values
# -----------------------------------------------------------------------------
is_missing_like <- function(x) {
  xn <- toupper(trimws(as.character(x)))
  xn %in% c("", "NA", "N/A", "NULL")
}

# -----------------------------------------------------------------------------
# General helpers moved from app/global.R for sharing: trimming, parsing, encoding
# -----------------------------------------------------------------------------
safe_trimws <- function(value) {
  if (is.null(value) || length(value) == 0 || all(is.na(value))) {
    return("")
  }
  trimws(as.character(value)[1])
}

safe_read_table <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    return(utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel metadata files.")
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  stop("Unsupported metadata extension: ", ext)
}

# Parsing helpers
parse_as_text <- function(value) {
  value <- trimws(as.character(value)[1])
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  gsub("^\"|\"$|^'|'$", "", value)
}

parse_as_logical <- function(value) {
  tolower(trimws(as.character(value)[1])) %in% c("true", "t", "1", "yes")
}

parse_as_numeric <- function(value) {
  suppressWarnings(as.numeric(as.character(value)[1]))
}

parse_as_vector <- function(value) {
  if (length(value) == 1 && toupper(value) == "NULL") {
    return(character(0))
  }
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  parts <- strsplit(value, ",", fixed = TRUE)[[1]]
  parts <- gsub("^\"|\"$|^'|'$", "", trimws(parts))
  parts[nzchar(parts)]
}

extract_and_parse_setting <- function(config_text, key, parser_fn, default = NULL) {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(default)
  }
  parser_fn(value)
}

# Setting display helpers
setting_display_value <- function(config_text, key, default = "") {
  extract_and_parse_setting(config_text, key, parse_as_text, default)
}

setting_display_logical <- function(config_text, key, default = FALSE) {
  extract_and_parse_setting(config_text, key, parse_as_logical, default)
}

setting_display_numeric <- function(config_text, key, default = NA_real_) {
  numeric_value <- extract_and_parse_setting(config_text, key, parse_as_numeric, NA_real_)
  if (length(numeric_value) != 1 || is.na(numeric_value)) default else numeric_value
}

setting_display_integer <- function(config_text, key, default = NA_integer_) {
  value <- setting_display_numeric(config_text, key, default = NA_real_)
  if (length(value) != 1 || is.na(value)) default else as.integer(round(value))
}

setting_display_vector <- function(config_text, key) {
  extract_and_parse_setting(config_text, key, parse_as_vector, character(0))
}

setting_display_csv <- function(config_text, key, default = "") {
  vec <- setting_display_vector(config_text, key)
  if (length(vec) == 0) default else paste(vec, collapse = ", ")
}

setting_display_sheet <- function(config_text, key, default = "") {
  value <- setting_display_value(config_text, key, default = default)
  if (!setting_has_value(value)) {
    return(default)
  }
  suppressWarnings({
    numeric_value <- as.numeric(value)
  })
  if (!is.na(numeric_value) && abs(numeric_value - round(numeric_value)) < .Machine$double.eps^0.5) {
    as.character(as.integer(round(numeric_value)))
  } else {
    value
  }
}

# Encoding helpers
encode_text <- function(value, nullable = FALSE) {
  value <- safe_trimws(value)
  if (!nzchar(value) || (nullable && toupper(value) == "NULL")) {
    return("NULL")
  }
  dQuote(value)
}

encode_logical <- function(value) {
  val <- toupper(safe_trimws(value))
  if (val %in% c("TRUE", "T", "1", "YES")) "TRUE" else "FALSE"
}

encode_numeric <- function(value, default = 0) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) == 0 || is.na(numeric_value)) {
    numeric_value <- default
  }
  format(numeric_value, scientific = FALSE, trim = TRUE)
}

encode_integer <- function(value, default = 0L) {
  integer_value <- suppressWarnings(as.integer(round(as.numeric(value))))
  if (length(integer_value) == 0 || is.na(integer_value)) {
    integer_value <- default
  }
  as.character(integer_value)
}

encode_sheet <- function(value) {
  value <- safe_trimws(value)
  if (!nzchar(value)) {
    return('""')
  }
  if (grepl("^-?[0-9]+(\\.[0]+)?$", value)) {
    as.character(as.integer(round(as.numeric(value))))
  } else {
    dQuote(value)
  }
}

encode_vector_text <- function(value, allow_null = FALSE) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return(if (allow_null) "NULL" else "c()")
  }
  paste0("c(", paste(dQuote(items), collapse = ", "), ")")
}

encode_vector_numeric <- function(value) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]
  if (length(items) == 0) {
    return("c()")
  }
  numeric_items <- suppressWarnings(as.numeric(items))
  numeric_items <- numeric_items[!is.na(numeric_items)]
  if (length(numeric_items) == 0) {
    return("c()")
  }
  paste0("c(", paste(format(numeric_items, scientific = FALSE, trim = TRUE), collapse = ", "), ")")
}

# Setting value encoders
setting_value_text <- function(value) { encode_text(value, nullable = FALSE) }
setting_value_nullable_text <- function(value) { encode_text(value, nullable = TRUE) }
setting_value_logical <- function(value) { encode_logical(value) }
setting_value_numeric <- function(value, default = 0) { encode_numeric(value, default) }
setting_value_integer <- function(value, default = 0L) { encode_integer(value, default) }
setting_value_sheet <- function(value) { encode_sheet(value) }
setting_value_vector_text <- function(value, allow_null = FALSE) { encode_vector_text(value, allow_null) }
setting_value_vector_numeric <- function(value) { encode_vector_numeric(value) }

parse_num_robust <- function(x, decimal_mark = ".", grouping_mark = ",") {
  raw <- stringr::str_trim(as.character(x))
  
  suppressWarnings(
    readr::parse_number(
      raw,
      locale = readr::locale(decimal_mark = decimal_mark, grouping_mark = grouping_mark)
    )
  )
}

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
      root = file.path(model_root, "exports"),
      metaboanalyst = file.path(model_root, "exports", "metaboanalyst"),
      stats = file.path(model_root, "exports", "stats")
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
  
  dir.create(dirs$global$root, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$global$exports, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirs$global$audits, showWarnings = FALSE, recursive = TRUE)
  
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

# -----------------------------------------------------------------------------
# Runtime helper
# -----------------------------------------------------------------------------
run_step <- function(step_name, expr) {
  message("\n==================================================")
  message("START: ", step_name)
  message("TIME : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  t0 <- Sys.time()
  result <- force(expr)
  t1 <- Sys.time()
  
  elapsed <- round(as.numeric(difftime(t1, t0, units = "secs")), 2)
  
  message("DONE : ", step_name)
  message("ELAPSED: ", elapsed, " sec")
  message("==================================================")
  
  invisible(result)
}

# -----------------------------------------------------------------------------
# Console helpers
# -----------------------------------------------------------------------------
fmt_time_sec <- function(t0, t1 = Sys.time()) {
  round(as.numeric(difftime(t1, t0, units = "secs")), 2)
}

console_rule <- function(char = "=", n = 60) {
  paste(rep(char, n), collapse = "")
}

step_start <- function(step_no, step_total, title) {
  message(console_rule())
  message(sprintf("[STEP %02d/%02d] %s", step_no, step_total, title))
  message(console_rule())
  invisible(Sys.time())
}

step_info <- function(...) {
  message("[INFO] ", paste0(..., collapse = ""))
}

step_ok <- function(title, t0 = NULL) {
  if (is.null(t0)) {
    message("[OK] ", title)
  } else {
    message("[OK] ", title, " finished in ", fmt_time_sec(t0), " sec")
  }
}

step_warn <- function(...) {
  message("[WARN] ", paste0(..., collapse = ""))
}

step_fail <- function(...) {
  message("[ERROR] ", paste0(..., collapse = ""))
}
