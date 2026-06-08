# Metabolomics Pipeline - global objects and helpers

# Compute project/pipeline paths early so we can ensure packages are
# installed before attempting to load libraries used by the app.
locate_project_root <- function(start_path = getwd()) {
  current_path <- normalizePath(start_path, winslash = "/", mustWork = FALSE)

  repeat {
    has_pipeline <- dir.exists(file.path(current_path, "pipeline", "R"))
    has_app <- dir.exists(file.path(current_path, "app"))
    has_app_entry <- file.exists(file.path(current_path, "app.R"))

    if (has_pipeline && has_app && has_app_entry) {
      return(normalizePath(current_path, winslash = "/", mustWork = TRUE))
    }

    parent_path <- normalizePath(file.path(current_path, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent_path, current_path)) {
      stop(
        "Could not locate the project root. Open the repository root or make sure app.R, app/, and pipeline/R/ exist."
      )
    }

    current_path <- parent_path
  }
}

project_root <- locate_project_root(getwd())
pipeline_root <- file.path(project_root, "pipeline")
r_dir <- file.path(pipeline_root, "R")
config_dir <- file.path(pipeline_root, "config")
app_assets_dir <- file.path(project_root, "app", "assets")

# Attempt to run the pipeline package installer script if present. This
# installs CRAN and Bioconductor packages (e.g. limma) so the Shiny app can
# load its dependencies later without error. Use try() so failures don't
# completely block the app UI from starting.
packages_script <- file.path(r_dir, "00_packages.R")
if (file.exists(packages_script)) {
  try(source(packages_script), silent = TRUE)
}

# Load app libraries after ensuring installation
library(shiny)
library(shinyjs)
library(bslib)

# Load shared pipeline helpers when available
helpers_path <- file.path(pipeline_root, "R", "03_helpers_io_log.R")
if (file.exists(helpers_path)) source(helpers_path)
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
    "magick",
    "limma",
    "ggplot2"
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
  r_dir,
  pattern = "\\.R$", full.names = TRUE
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

# Use shared implementations from pipeline/R/03_helpers_io_log.R when available
# (get_comparison_group_labels_for_model / map_comparison_group_display_values)

validate_metadata_columns <- function(path,
                                      metadata_mapping = NULL,
                                      allowed_groups = c("WT", "TG"),
                                      model_allowed_groups_by_model = NULL) {
  alias_map <- list(
    sample = c("sample", "sample_id", "sample_name", "id_sample", "id", "name"),
    weight = c("weight", "weight_mg", "mass", "mass_mg", "mg", "sample_weight", "sample_mass", "weight_g", "mass_g", "sample_weight_g", "sample_mass_g"),
    group = c("group", "treatment", "treat"),
    sex = c("sex", "gender"),
    model = c("model", "disease", "condition", "phenotype", "status")
  )

  md <- safe_read_table(path)
  actual <- tolower(trimws(as.character(names(md))))

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

  allowed_groups <- unique(trimws(as.character(allowed_groups)))
  allowed_groups <- allowed_groups[nzchar(allowed_groups)]
  model_group_values <- character(0)

  if (length(allowed_groups) < 2 && !is.null(model_allowed_groups_by_model) && length(model_allowed_groups_by_model) > 0) {
    inferred_groups <- unique(trimws(unlist(strsplit(as.character(model_allowed_groups_by_model), ",", fixed = TRUE), use.names = FALSE)))
    inferred_groups <- inferred_groups[nzchar(inferred_groups)]
    if (length(inferred_groups) >= 2) {
      allowed_groups <- inferred_groups
    }
  }

  if (!is.null(model_allowed_groups_by_model) && length(model_allowed_groups_by_model) > 0) {
    model_group_values <- unique(trimws(unlist(
      strsplit(as.character(model_allowed_groups_by_model), ",", fixed = TRUE),
      use.names = FALSE
    )))
    model_group_values <- model_group_values[nzchar(model_group_values)]
  }

  valid_group_values <- unique(c(allowed_groups, model_group_values))
  allowed_groups_norm <- toupper(valid_group_values)
  if (length(allowed_groups) < 2) {
    return(list(
      ok = FALSE,
      message = "Please provide at least two allowed group values in this order: control, test (e.g. WT, TG), or define them in model_allowed_groups_by_model."
    ))
  }

  group_col <- resolved_cols[["group"]]
  if (!is.null(group_col) && nzchar(group_col)) {
    col_idx <- which(actual == group_col)
    if (length(col_idx) == 0) {
      return(list(ok = FALSE, message = paste0("Group column '", group_col, "' not found in metadata.")))
    }
    groups_raw <- as.character(md[[col_idx[1]]])
    groups_raw <- trimws(groups_raw)

    if (!is.null(model_allowed_groups_by_model) && length(model_allowed_groups_by_model) > 0 && "model" %in% names(resolved_cols)) {
      model_idx <- which(actual == resolved_cols[["model"]])
      if (length(model_idx) > 0) {
        model_raw <- trimws(as.character(md[[model_idx[1]]]))
        valid_rows <- !is.na(groups_raw) & !is_missing_like(groups_raw) & !is.na(model_raw) & !is_missing_like(model_raw)

        if (any(valid_rows)) {
          groups_raw[valid_rows] <- normalize_model_group_pairs(
            groups_raw[valid_rows],
            model_raw[valid_rows],
            model_allowed_groups_by_model,
            allowed_groups[1],
            allowed_groups[2]
          )
        }
      }
    }

    groups_raw <- groups_raw[!is.na(groups_raw)]
    groups_raw <- groups_raw[!is_missing_like(groups_raw)]

    groups_norm <- toupper(groups_raw)
    invalid_groups <- sort(unique(groups_raw[!(groups_norm %in% allowed_groups_norm)]))

    if (length(invalid_groups) > 0) {
      return(list(
        ok = FALSE,
        message = paste0(
          "Invalid values in metadata group column ('",
          group_col,
          "'). Allowed values: ",
          paste(valid_group_values, collapse = ", "),
          ". Found: ",
          paste(invalid_groups, collapse = ", "),
          ". Please review your metadata file and fix the group column."
        )
      ))
    }
  }

  list(ok = TRUE, message = "Metadata validation passed.")
}


setting_input_id <- function(key) {
  paste0("settings_", key)
}

setting_has_value <- function(value) {
  !is.null(value) && length(value) > 0 && !all(is.na(value)) && nzchar(trimws(as.character(value)[1]))
}

settings_form_sections <- list(
  list(
    title = "Normalization",
    fields = list(
      list(key = "normalization_mode", label = "Main normalization", type = "select", choices = c("none", "PQN", "QC_LOESS"), default = "none"),
      list(key = "loess_min_qc_points", label = "LOESS minimum QC points", type = "integer", default = 5, step = 1, min = 5),
      list(key = "QC_LOESS_span", label = "QC-LOESS span", type = "numeric", default = 0.75, step = 0.05, min = 0.05, max = 1)
    )
  ),
  list(
    title = "RSD and IQR filters",
    fields = list(
      list(key = "rsd_filter_metric", label = "RSD filter metric", type = "select", choices = c("none", "qc_rsd", "rsd"), default = "none"),
      list(key = "rsd_thresholds", label = "RSD thresholds", type = "vector_numeric", default = c(20), placeholder = "Example: 10, 15, 20, 30"),
      list(key = "active_variant", label = "Active variant", type = "selectize_text", default = "BASE", choices = c("BASE", "QC_RSD10", "QC_RSD15", "QC_RSD20", "QC_RSD30", "RSD10", "RSD15", "RSD20", "RSD30")),
      list(key = "low_variance_filter_method", label = "Low-variance filter", type = "select", choices = c("none", "iqr"), default = "none"),
      list(key = "low_variance_filter_fraction", label = "IQR removal fraction", type = "numeric", default = 0.20, step = 0.05, min = 0, max = 1)
    )
  ),
  list(
    title = "Statistics thresholds",
    fields = list(
      list(key = "p_value_cutoff", label = "P-value cutoff", type = "numeric", default = 0.05, step = 0.001, min = 0, max = 1),
      list(key = "fdr_cutoff", label = "FDR cutoff", type = "numeric", default = 0.05, step = 0.001, min = 0, max = 1),
      list(key = "fc_cutoff_log2", label = "FC cutoff (log2)", type = "numeric", default = 0, step = 0.1, min = 0)
    )
  ),
  list(
    title = "PCA",
    fields = list(
      list(key = "pca_scaling", label = "PCA scaling", type = "select", choices = c("none", "pareto", "autoscale"), default = "pareto"),
      list(key = "ellipse_positive", label = "Enable group ellipses", type = "logical_select", default = TRUE)
    )
  ),
  list(
    title = "Heatmap",
    fields = list(
      list(key = "heatmap_cluster_distance", label = "Heatmap cluster distance", type = "select", choices = c("euclidean", "manhattan"), default = "euclidean"),
      list(key = "heatmap_cluster_method", label = "Heatmap cluster method", type = "select", choices = c("ward.D2", "complete", "average"), default = "ward.D2"),
      list(key = "heatmap_top_n", label = "Heatmap top N", type = "integer", default = 50, step = 1, min = 1)
    )
  ),
  list(
    title = "Feature filtering and naming",
    fields = list(
      list(key = "use_only_known", label = "Use only known features", type = "logical_select", default = FALSE),
      list(key = "sanitize_mode", label = "Sanitize mode", type = "select", choices = c("none", "greek_latin_ascii", "ascii_translit"), default = "none")
    )
  )
)

settings_builder_layout <- list(
  run_setup = list(
    label = "Run setup",
    sections = c(
      "Statistics thresholds",
      "Normalization",
      "Feature filters"
    ),
    widths = c(6, 6, 12)
  ),
  visual_outputs = list(
    label = "Visual outputs",
    sections = c(
      "Heatmap",
      "PCA"
    ),
    widths = c(6, 6)
  ),
  exports = list(
    label = "Exports",
    sections = c("Output controls"),
    widths = c(12),
    narrow = FALSE
  ),
  field_columns = c(
    default = 2,
    "Output controls" = 3,
    "PCA" = 2,
    "Normalization" = 2,
    "Feature filters" = 3,
    "RSD and IQR filters" = 3,
    "Feature filtering and naming" = 2,
    "Statistics thresholds" = 3,
    "Heatmap" = 3
  )
)

settings_glossary_map <- c(
  use_reference_file = "Enables reference-table matching for duplicate handling.",
  output_dir = "Defines where all run outputs are written.",
  use_weight_normalization = "Applies sample-weight normalization before downstream analysis.",
  normalization_mode = "Main normalization after optional weight normalization: none keeps the post-weight matrix, PQN applies QC-reference PQN, and QC_LOESS corrects signal drift using QC samples.",
  loess_min_qc_points = "Minimum number of valid QC values required for a feature to receive QC-LOESS drift correction.",
  QC_LOESS_span = "LOESS smoothing span for QC-LOESS; smaller values follow local drift more closely, larger values smooth more strongly.",
  rsd_filter_metric = "Controls RSD-based feature filtering: none keeps BASE, qc_rsd filters using QC sample RSD, and rsd filters using biological/sample RSD.",
  rsd_thresholds = "Numeric RSD cutoffs used to create variants such as QC_RSD20 or RSD20.",
  active_variant = "Selects which RSD-filtered variant is used downstream. BASE disables RSD filtering.",
  low_variance_filter_method = "Low-variance feature filtering method. none disables it; iqr removes the lowest-IQR fraction.",
  low_variance_filter_fraction = "Fraction of lowest-IQR features removed when low_variance_filter_method is iqr.",
  duplicate_name_strategy = "Sets how duplicate features are merged or kept.",
  run_metrics = "Selects the significance metric used in run-level decisions and rankings.",
  p_value_cutoff = "Raw p-value threshold used for significance decisions, volcano highlighting, and p-value ranked outputs.",
  fdr_cutoff = "Adjusted p-value threshold used for FDR significance decisions, volcano highlighting, and FDR ranked outputs.",
  fc_cutoff_log2 = "Minimum absolute log2 fold-change required when fold-change filtering is enabled; use 0 to disable the fold-change requirement.",
  use_only_known = "If TRUE, removes features without a known/canonical name before downstream statistics and plots.",
  pca_scaling = "Scaling applied before PCA: none leaves variables unchanged, pareto divides by sqrt(SD), and autoscale divides by SD.",
  heatmap_top_n = "Maximum number of top-ranked features shown in each top heatmap.",
  dup_mz_digits = "Rounding precision for m/z during duplicate detection.",
  dup_rt_digits = "Rounding precision for RT during duplicate detection.",
  sanitize_mode = "Text handling for exported feature names: none preserves characters, greek_latin_ascii transliterates Greek/Latin to ASCII, and ascii_translit uses general ASCII transliteration.",
  make_heatmap_by_model = "TRUE generates top-ranked heatmaps per model.",
  make_heatmap_by_model_sex = "TRUE also generates top-ranked heatmaps split by sex.",
  heatmap_scale_method = "Scaling applied to heatmap matrices: none, zscore, or pareto.",
  heatmap_cluster_distance = "Distance metric used for heatmap clustering: euclidean or manhattan.",
  heatmap_cluster_method = "Hierarchical clustering linkage used for heatmaps: ward.D2, complete, or average.",
  ellipse_positive = "Controls PCA ellipse rendering behavior; when TRUE, group ellipses are drawn only for eligible comparisons."
)

read_initial_config <- function() {
  if (file.exists(active_config_path)) {
    return(safe_read_file(active_config_path))
  }
  "# No settings file found in config/."
}

initial_settings_text <- read_initial_config()
