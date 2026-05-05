# Metabolomics Pipeline - global objects and helpers
library(shiny)
library(shinyjs)

# Ensure decimal separator is dot (.) not comma (,)
try(Sys.setlocale("LC_NUMERIC", "C"), silent = TRUE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pipeline_root <- file.path(project_root, "pipeline")
r_dir <- file.path(pipeline_root, "R")
config_dir <- file.path(pipeline_root, "config")
app_assets_dir <- file.path(project_root, "app", "assets")
example_config_path <- file.path(config_dir, "settings.example.R")
active_config_path <- file.path(config_dir, "settings.R")

if (dir.exists(app_assets_dir)) {
  addResourcePath("assets", app_assets_dir)
}

get_pipeline_required_packages <- function(packages_file = file.path(r_dir, "00_packages.R")) {
  fallback <- c(
    "tidyverse",
    "readr",
    "readxl",
    "openxlsx",
    "pheatmap",
    "ggrepel",
    "stringi",
    "RColorBrewer",
    "processx",
    "shinyFiles",
    "magick"
  )
  
  if (!file.exists(packages_file)) {
    return(fallback)
  }
  
  lines <- readLines(packages_file, warn = FALSE)
  start <- grep("^\\s*[A-Za-z0-9_.]+\\s*<-\\s*c\\(", lines)
  if (length(start) == 0) {
    return(fallback)
  }
  
  pkgs <- unique(unlist(lapply(start, function(idx) {
    end_rel <- grep("\\)", lines[(idx + 1):length(lines)])
    if (length(end_rel) == 0) {
      return(character(0))
    }

    end <- idx + end_rel[1]
    block <- paste(lines[idx:end], collapse = " ")
    matches <- unlist(regmatches(block, gregexpr("\"[^\"]+\"|'[^']+'", block, perl = TRUE)))
    gsub("^\"|\"$|^'|'$", "", matches)
  }), use.names = FALSE))
  pkgs <- unique(pkgs[nzchar(pkgs)])
  
  if (length(pkgs) == 0) {
    return(fallback)
  }
  
  pkgs
}

required_packages <- get_pipeline_required_packages()

script_paths <- c(file.path(pipeline_root, "run_pipeline.R"), sort(list.files(
  r_dir, pattern = "\\.R$", full.names = TRUE
)))

script_names <- vapply(script_paths, function(p) {
  gsub("\\\\", "/", sub(paste0(
    "^", gsub("\\\\", "/", project_root), "/?"
  ), "", gsub("\\\\", "/", p)))
}, character(1))

safe_read_file <- function(path) {
  if (!file.exists(path)) {
    return(paste("File not found:", path))
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

extract_config_value <- function(config_text, key) {
  pattern <- paste0("^\\s*", key, "\\s*<-\\s*(.+)$")
  lines <- strsplit(config_text, "\\n", fixed = FALSE)[[1]]
  hit <- grep(pattern, lines, perl = TRUE, value = TRUE)
  
  if (length(hit) == 0) {
    return(NULL)
  }
  
  raw <- sub(pattern, "\\1", hit[length(hit)], perl = TRUE)
  raw <- sub("\\s*#.*$", "", trimws(raw))
  
  if (grepl("^\".*\"$|^'.*'$", raw)) {
    return(gsub("^\"|\"$|^'|'$", "", raw))
  }
  
  raw
}

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\]|^/|^~", path)
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

validate_metadata_columns <- function(path, metadata_mapping = NULL, allowed_groups = c("TG", "WT")) {

  alias_map <- list(
    sample = c("sample", "sample_id", "sample_name", "id_sample", "id", "name"),
    weight = c("weight", "weight_mg", "mass", "mass_mg", "mg", "sample_weight", "sample_mass", "weight_g", "mass_g", "sample_weight_g", "sample_mass_g"),
    group = c("group", "treatment", "treat"),
    sex = c("sex", "gender"),
    model = c("model", "disease", "condition", "phenotype", "status")
  )

  md <- safe_read_table(path)
  actual <- tolower(trimws(as.character(names(md))))

  normalize_name <- function(x) {
    tolower(trimws(as.character(x)))
  }

  missing <- character(0)
  resolved_cols <- list()

  for (target in names(alias_map)) {
    mapped_col <- ""
    if (!is.null(metadata_mapping) && !is.null(metadata_mapping[[target]])) {
      mapped_col <- normalize_name(metadata_mapping[[target]])
    }

    if (nzchar(mapped_col)) {
      if (!(mapped_col %in% actual)) {
        missing <- c(missing, target)
      } else {
        resolved_cols[[target]] <- mapped_col
      }
    } else {
      aliases <- normalize_name(alias_map[[target]])
      found <- intersect(aliases, actual)
      if (length(found) == 0) {
        missing <- c(missing, target)
      } else {
        resolved_cols[[target]] <- found[1]
      }
    }
  }

  if (length(missing) > 0) {
    return(list(
      ok = FALSE,
      message = paste(
        "Metadata missing required columns or mappings:",
        paste(missing, collapse = ", ")
      )
    ))
  }

  allowed_groups <- unique(toupper(trimws(as.character(allowed_groups))))
  allowed_groups <- allowed_groups[nzchar(allowed_groups)]
  if (length(allowed_groups) < 2) {
    return(list(
      ok = FALSE,
      message = "Please provide at least two allowed group values in this order: control, test (e.g. WT, TG)."
    ))
  }

  group_col <- resolved_cols[["group"]]
  if (!is.null(group_col) && nzchar(group_col)) {
    is_missing_like <- function(x) {
      xn <- toupper(trimws(as.character(x)))
      xn %in% c("", "NA", "N/A", "NULL")
    }

    col_idx <- which(actual == group_col)
    if (length(col_idx) == 0) {
      return(list(ok = FALSE, message = paste0("Group column '", group_col, "' not found in metadata.")))
    }
    groups_raw <- as.character(md[[col_idx[1]]])
    groups_raw <- trimws(groups_raw)
    groups_raw <- groups_raw[!is.na(groups_raw)]
    groups_raw <- groups_raw[!is_missing_like(groups_raw)]

    groups_norm <- toupper(groups_raw)
    invalid_groups <- sort(unique(groups_raw[!(groups_norm %in% allowed_groups)]))

    if (length(invalid_groups) > 0) {
      return(list(
        ok = FALSE,
        message = paste0(
          "Invalid values in metadata group column ('",
          group_col,
          "'). Allowed values: ",
          paste(allowed_groups, collapse = ", "),
          ". Found: ",
          paste(invalid_groups, collapse = ", "),
          ". Please review your metadata file and fix the group column."
        )
      ))
    }
  }

  list(ok = TRUE, message = "Metadata validation passed.")
}

normalize_config_text <- function(text) {
  text <- gsub("\u201C|\u201D", '"', text, perl = TRUE)
  text <- gsub("\u2018|\u2019", "'", text, perl = TRUE)
  text
}

setting_input_id <- function(key) {
  paste0("settings_", key)
}

setting_has_value <- function(value) {
  !is.null(value) && length(value) > 0 && !all(is.na(value)) && nzchar(trimws(as.character(value)[1]))
}

setting_display_value <- function(config_text, key, default = "") {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(default)
  }

  value <- trimws(as.character(value)[1])
  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  value <- gsub("^\"|\"$|^'|'$", "", value)
  trimws(value)
}

setting_display_logical <- function(config_text, key, default = FALSE) {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(default)
  }

  tolower(trimws(as.character(value)[1])) %in% c("true", "t", "1", "yes")
}

setting_display_numeric <- function(config_text, key, default = NA_real_) {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(default)
  }

  numeric_value <- suppressWarnings(as.numeric(as.character(value)[1]))

  if (length(numeric_value) != 1 || is.na(numeric_value)) {
    default
  } else {
    numeric_value
  }
}

setting_display_integer <- function(config_text, key, default = NA_integer_) {
  value <- setting_display_numeric(config_text, key, default = default)
  if (length(value) != 1 || is.na(value)) {
    return(default)
  }
  as.integer(round(value))
}

setting_display_vector <- function(config_text, key) {
  value <- extract_config_value(config_text, key)
  if (!setting_has_value(value)) {
    return(character(0))
  }

  value <- trimws(value)
  if (length(value) == 1 && toupper(value) == "NULL") {
    return(character(0))
  }

  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  parts <- strsplit(value, ",", fixed = TRUE)[[1]]
  parts <- gsub("^\"|\"$|^'|'$", "", trimws(parts))
  parts[nzchar(parts)]
}

setting_display_csv <- function(config_text, key, default = "") {
  vec <- setting_display_vector(config_text, key)
  if (length(vec) == 0) {
    return(default)
  }
  paste(vec, collapse = ", ")
}

setting_default_vector <- function(default) {
  if (is.null(default) || length(default) == 0) {
    return(character(0))
  }

  value <- trimws(as.character(default)[1])
  if (!nzchar(value)) {
    return(character(0))
  }

  value <- sub("^c\\((.*)\\)$", "\\1", value, perl = TRUE)
  items <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  items <- gsub("^\"|\"$|^'|'$", "", items)
  items[nzchar(items)]
}

setting_default_numeric_vector <- function(default) {
  values <- setting_default_vector(default)
  numeric_values <- suppressWarnings(as.numeric(values))
  numeric_values <- numeric_values[!is.na(numeric_values)]
  if (length(numeric_values) == 0) {
    return(character(0))
  }
  format(numeric_values, scientific = FALSE, trim = TRUE)
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

setting_value_text <- function(value) {
  value <- safe_trimws(value)
  if (!nzchar(value)) {
    return('""')
  }
  dQuote(value)
}

setting_value_nullable_text <- function(value) {
  value <- safe_trimws(value)
  if (!nzchar(value) || toupper(value) == "NULL") {
    return("NULL")
  }
  dQuote(value)
}

setting_value_logical <- function(value) {
  val <- toupper(safe_trimws(value))
  if (val %in% c("TRUE", "T", "1", "YES")) "TRUE" else "FALSE"
}

setting_value_numeric <- function(value, default = 0) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) == 0 || is.na(numeric_value)) {
    numeric_value <- default
  }
  format(numeric_value, scientific = FALSE, trim = TRUE)
}

setting_value_integer <- function(value, default = 0L) {
  integer_value <- suppressWarnings(as.integer(round(as.numeric(value))))
  if (length(integer_value) == 0 || is.na(integer_value)) {
    integer_value <- default
  }
  as.character(integer_value)
}

setting_value_sheet <- function(value) {
  value <- safe_trimws(value)
  if (!nzchar(value)) {
    return('""')
  }

  if (grepl("^-?\\d+(\\.0+)?$", value)) {
    as.character(as.integer(round(as.numeric(value))))
  } else {
    dQuote(value)
  }
}

setting_value_vector_text <- function(value, allow_null = FALSE) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]

  if (length(items) == 0) {
    if (allow_null) {
      return("NULL")
    }
    return("c()")
  }

  paste0("c(", paste(dQuote(items), collapse = ", "), ")")
}

setting_value_vector_numeric <- function(value) {
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

settings_form_sections <- list(
  list(
    title = "Statistics thresholds",
    fields = list(
      list(key = "alpha_sig", label = "Alpha significance", type = "numeric", default = 0.05, step = 0.001, min = 0, max = 1),
      list(key = "fc_cutoff_log2", label = "FC cutoff (log2)", type = "numeric", default = 0, step = 0.1, min = 0, help = "Minimum absolute log2 fold-change used in significance filtering. Use 0 to disable this cutoff."),
      list(key = "heatmap_top_n", label = "Heatmap top N", type = "integer", default = 80, step = 1, min = 1, help = "Maximum number of top-ranked features shown in top heatmaps.")
    )
  ),
  list(
    title = "Heatmap clustering",
    fields = list(
      list(key = "heatmap_cluster_distance", label = "Heatmap cluster distance", type = "select", choices = c("euclidean", "manhattan"), default = "euclidean"),
      list(key = "heatmap_cluster_method", label = "Heatmap cluster method", type = "select", choices = c("ward.D2", "complete", "average"), default = "ward.D2")
    )
  ),
  list(
    title = "Plot generation",
    fields = list(
      list(key = "make_heatmap_by_model", label = "Heatmap by model", type = "logical_select", default = TRUE),
      list(key = "make_heatmap_by_model_sex", label = "Heatmap by model and sex", type = "logical_select", default = TRUE)
    )
  ),
  list(
    title = "PCA and heatmap style",
    fields = list(
      list(key = "pca_scaling", label = "PCA scaling", type = "select", choices = c("none", "pareto", "autoscale"), default = "pareto"),
      list(key = "heatmap_scale_method", label = "Heatmap scale method", type = "select", choices = c("none", "zscore", "pareto"), default = "zscore"),
      list(key = "sanitize_mode", label = "Sanitize mode", type = "select", choices = c("greek_latin_ascii", "ascii_translit"), default = "greek_latin_ascii")
    )
  )
)

settings_glossary_map <- c(
  use_reference_file = "Enables reference-table matching for duplicate handling.",
  output_dir = "Defines where all run outputs are written.",
  use_weight_normalization = "Applies sample-weight normalization before downstream analysis.",
  duplicate_name_strategy = "Sets how duplicate features are merged or kept.",
  run_metrics = "Selects the significance metric used in run-level decisions and rankings.",
  alpha_sig = "Significance threshold used for stats, volcano labels, and significant heatmaps.",
  fc_cutoff_log2 = "Minimum absolute log2 fold-change required where fold-change filtering is enabled; use 0 to disable.",
  pca_scaling = "Scaling mode applied before PCA: none, pareto, or autoscale.",
  heatmap_top_n = "Maximum number of ranked features shown in top heatmaps.",
  dup_mz_digits = "Rounding precision for m/z during duplicate detection.",
  dup_rt_digits = "Rounding precision for RT during duplicate detection.",
  sanitize_mode = "Chooses the text normalization strategy for feature names and labels.",
  make_heatmap_by_model = "TRUE generates top-ranked heatmaps per model.",
  make_heatmap_by_model_sex = "TRUE also generates top-ranked heatmaps split by sex.",
  heatmap_scale_method = "Scaling applied to heatmap matrices: none, zscore, or pareto.",
  heatmap_cluster_distance = "Recommended distance metric for heatmap clustering: euclidean or manhattan.",
  heatmap_cluster_method = "Recommended hierarchical clustering method: ward.D2, complete, or average."
)

read_initial_config <- function() {
  if (file.exists(active_config_path)) {
    return(safe_read_file(active_config_path))
  }
  if (file.exists(example_config_path)) {
    return(safe_read_file(example_config_path))
  }
  "# No settings file found in config/."
}

initial_settings_text <- read_initial_config()
