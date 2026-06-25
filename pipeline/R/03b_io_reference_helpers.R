# =============================================================================
# 03b_io_reference_helpers.R
# IO and reference helpers
# =============================================================================

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

clean_input_file_sample_name <- function(x) {
  s <- as.character(x)
  s <- stringr::str_replace(s, "(?i)\\.raw.*$", "")
  s <- stringr::str_replace(s, "\\s*\\(.*\\)$", "")
  s <- stringr::str_trim(s)

  non_qc <- !stringr::str_detect(s, "^QC")
  s[non_qc] <- sub("_.*$", "", s[non_qc])

  s
}

parse_input_creation_datetime <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x))
  }
  if (inherits(x, "Date")) {
    return(as.POSIXct(x))
  }
  if (is.numeric(x)) {
    return(as.POSIXct(x * 86400, origin = "1899-12-30", tz = "UTC"))
  }

  x <- trimws(as.character(x))
  parsed <- suppressWarnings(as.POSIXct(
    x,
    format = "%m/%d/%Y %I:%M:%S %p",
    tz = "UTC"
  ))
  missing <- is.na(parsed) & nzchar(x)
  if (any(missing)) {
    parsed[missing] <- suppressWarnings(as.POSIXct(
      x[missing],
      format = "%m/%d/%Y %H:%M:%S",
      tz = "UTC"
    ))
  }
  missing <- is.na(parsed) & nzchar(x)
  if (any(missing)) {
    parsed[missing] <- suppressWarnings(as.POSIXct(x[missing], tz = "UTC"))
  }

  parsed
}

read_input_files_reference <- function(path = file.path("data", "Input Files.xlsx"), sheet = "InputFiles") {
  if (!file.exists(path)) {
    return(NULL)
  }

  ref <- tryCatch(
    read_any_table(path, sheet = sheet),
    error = function(e) NULL
  )

  if (is.null(ref) || nrow(ref) == 0) {
    return(NULL)
  }

  nms <- names(ref)
  nms_low <- tolower(trimws(nms))
  order_col <- nms[nms_low == "order"]
  sample_col <- nms[nms_low %in% c("samples", "sample", "file name", "filename")]
  file_col <- nms[nms_low %in% c("file name", "filename", "file path", "filepath", "path")]
  creation_col <- nms[nms_low %in% c("creation date", "creation datetime", "created", "created date", "date created")]
  sample_type_col <- nms[nms_low == "sample type"]

  if (length(order_col) > 0 && length(sample_col) > 0) {
    return(
      tibble::tibble(
        sample = clean_input_file_sample_name(ref[[sample_col[1]]]),
        input_order = suppressWarnings(as.numeric(ref[[order_col[1]]])),
        input_order_source = basename(path)
      ) %>%
        dplyr::filter(!is.na(sample), nzchar(sample), is.finite(input_order)) %>%
        dplyr::arrange(input_order) %>%
        dplyr::distinct(sample, .keep_all = TRUE)
    )
  }

  if (length(file_col) > 0 && length(creation_col) > 0) {
    out <- tibble::tibble(
      sample = clean_input_file_sample_name(basename(gsub("\\\\", "/", as.character(ref[[file_col[1]]])))),
      input_datetime = parse_input_creation_datetime(ref[[creation_col[1]]]),
      sample_type = if (length(sample_type_col) > 0) as.character(ref[[sample_type_col[1]]]) else NA_character_,
      input_order_source = basename(path)
    ) %>%
      dplyr::filter(!is.na(sample), nzchar(sample), !is.na(input_datetime))

    if (length(sample_type_col) > 0) {
      out <- out %>%
        dplyr::filter(tolower(trimws(sample_type)) %in% c("sample", "quality control", "qc"))
    }

    return(
      out %>%
        dplyr::arrange(input_datetime, sample) %>%
        dplyr::distinct(sample, .keep_all = TRUE) %>%
        dplyr::mutate(input_order = dplyr::row_number()) %>%
        dplyr::select(sample, input_order, input_order_source)
    )
  }

  NULL
}

# Allows each model to override the global control/treatment group labels.
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

order_pre_post_levels <- function(levels, preferred = NULL) {
  levels <- as.character(levels)
  levels <- unique(levels[!is.na(levels) & nzchar(levels)])

  preferred <- as.character(preferred)
  preferred <- unique(preferred[!is.na(preferred) & nzchar(preferred)])
  preferred <- preferred[preferred %in% levels]

  base_order <- c(preferred, setdiff(levels, preferred))
  norm <- tolower(trimws(base_order))
  pre_idx <- norm %in% "pre" | grepl("^pre([_ -]|$)", norm)
  post_idx <- norm %in% "post" | grepl("^post([_ -]|$)", norm)

  if (!any(pre_idx) && !any(post_idx)) {
    return(base_order)
  }

  c(base_order[pre_idx], base_order[post_idx], base_order[!(pre_idx | post_idx)])
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
