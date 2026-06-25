# Small server helpers shared across observers and renderers.

parse_allowed_groups <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(character(0))
  }

  raw <- trimws(as.character(value))
  raw <- raw[!is.na(raw) & nzchar(raw)]
  if (length(raw) == 0) {
    return(character(0))
  }

  vals <- unlist(strsplit(raw, ",", fixed = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  unique(vals[nzchar(vals)])
}

order_metadata_groups_for_control <- function(groups) {
  groups <- unique(trimws(as.character(groups)))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) < 2) {
    return(groups)
  }

  group_norm <- toupper(groups)
  control_idx <- which(group_norm %in% c("WT", "CONTROL", "CTRL", "UNTREATED"))
  if (length(control_idx) == 0) {
    return(groups)
  }

  c(groups[control_idx[1]], groups[-control_idx[1]])
}

allowed_groups_require_model_aliases <- function(value) {
  parsed <- toupper(parse_allowed_groups(value))
  !identical(parsed, c("WT", "TG"))
}

allowed_groups_hint_text <- function(value) {
  parsed <- parse_allowed_groups(value)
  if (length(parsed) == 0) {
    return("")
  }

  if (length(parsed) < 2) {
    return("Select at least two groups: control first, then test.")
  }

  ""
}

allowed_groups_missing_comma <- function(value) {
  parsed <- parse_allowed_groups(value)
  length(parsed) > 0 && length(parsed) < 2
}

make_output_subdir_from_data_file <- function(file_name) {
  stem <- tools::file_path_sans_ext(basename(as.character(file_name)[1]))
  stem <- gsub("[^A-Za-z0-9._-]+", "_", stem)
  stem <- gsub("^_+|_+$", "", stem)

  if (!nzchar(stem)) {
    stem <- format(Sys.time(), "run_%Y%m%d_%H%M%S")
  }

  file.path("output", stem)
}

strip_outer_quotes <- function(path) {
  path <- safe_trimws(path)
  while (grepl("^\".*\"$|^'.*'$|^“.*”$|^‘.*’$", path)) {
    path <- gsub("^\"|\"$|^'|'$|^“|”$|^‘|’$", "", path)
  }
  path
}

sanitize_output_dir_path <- function(path) {
  out <- strip_outer_quotes(path)
  if (!nzchar(out)) {
    return("output")
  }

  out <- gsub("\\\\", "/", out)
  prefix <- ""
  if (grepl("^[A-Za-z]:", out)) {
    prefix <- substr(out, 1, 2)
    out <- substring(out, 3)
  }

  parts <- strsplit(out, "/", fixed = TRUE)[[1]]
  parts <- vapply(parts, function(part) {
    if (!nzchar(part)) {
      return(part)
    }
    part <- gsub("[<>:\"|?*]", "_", part)
    part <- gsub("^\\s+|\\s+$", "", part)
    if (!nzchar(part)) "_" else part
  }, character(1))

  sanitized <- paste(parts, collapse = "/")
  paste0(prefix, sanitized)
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

render_manifest_list <- function(items) {
  tags$ul(lapply(items, function(item) tags$li(item)))
}

is_valid_file_path <- function(path) {
  !is.null(path) && nzchar(as.character(path)) && file.exists(as.character(path))
}

is_valid_image_path <- function(path) {
  !is.null(path) && nzchar(as.character(path)) && file.exists(as.character(path))
}

read_log_preview <- function(path, max_bytes = 512 * 1024) {
  if (!is_valid_file_path(path)) {
    return("")
  }

  info <- file.info(path)
  size <- suppressWarnings(as.numeric(info$size[1]))
  if (is.na(size) || size <= 0) {
    return("")
  }

  bytes_to_read <- min(size, max_bytes)
  start_at <- max(0, size - bytes_to_read)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  if (start_at > 0) {
    seek(con, where = start_at, origin = "start")
  }

  raw <- readBin(con, what = "raw", n = bytes_to_read)
  text <- rawToChar(raw)
  Encoding(text) <- "UTF-8"

  if (start_at > 0) {
    text <- sub("^[^\r\n]*(\r?\n)?", "", text)
    text <- paste0(
      "[Showing the last ",
      round(max_bytes / 1024),
      " KB of the log]\n\n",
      text
    )
  }

  text
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

get_shiny_roots <- function() {
  normalize_root <- function(path) {
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  add_root_if_exists <- function(roots, label, path) {
    if (nzchar(path) && dir.exists(path)) {
      roots <- c(roots, setNames(normalize_root(path), label))
    }
    roots
  }

  roots <- c(Project = normalize_root(project_root))

  home_dir <- tryCatch(normalize_root("~"), error = function(e) "")
  if (nzchar(home_dir) && dir.exists(home_dir)) {
    roots <- c(roots, Home = home_dir)
  }

  roots <- add_root_if_exists(roots, "Downloads", file.path(home_dir, "Downloads"))
  roots <- add_root_if_exists(roots, "Documents", file.path(home_dir, "Documents"))
  roots <- add_root_if_exists(roots, "Documentos", file.path(home_dir, "Documentos"))

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

find_missing_packages <- function() {
  installed <- rownames(installed.packages())
  setdiff(required_packages, installed)
}

metadata_qc_group_aliases <- function() {
  c("QC", "POOLED QC", "POOL QC", "POOLED_QC", "POOL_QC", "QUALITY CONTROL", "QUALITY_CONTROL")
}

metadata_missing_aliases <- function() {
  c("NA", "N/A", "NULL")
}

clean_metadata_values <- function(values, drop_qc = FALSE) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  values <- values[!is_missing_like(values)]

  if (isTRUE(drop_qc)) {
    values <- values[!(toupper(values) %in% metadata_qc_group_aliases())]
  }

  unique(values)
}

read_metadata_for_app <- function(path, mapping) {
  if (!is_valid_file_path(path)) {
    return(NULL)
  }

  tryCatch(
    read_metadata_with_mapping(path, mapping),
    error = function(e) NULL
  )
}

detect_metadata_models <- function(md) {
  if (is.null(md) || !("model" %in% names(md))) {
    return(character(0))
  }

  models <- clean_metadata_values(md$model)
  models <- models[!(toupper(models) %in% metadata_missing_aliases())]
  sort(unique(models))
}

detect_metadata_groups <- function(md, drop_qc = TRUE) {
  if (is.null(md) || !("group" %in% names(md))) {
    return(character(0))
  }

  order_metadata_groups_for_control(clean_metadata_values(md$group, drop_qc = drop_qc))
}

detect_metadata_groups_by_model <- function(md) {
  if (is.null(md) || !all(c("model", "group") %in% names(md))) {
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
  lapply(split_groups, function(groups) {
    order_metadata_groups_for_control(clean_metadata_values(groups, drop_qc = TRUE))
  })
}

make_multigroup_pair_choices <- function(groups) {
  groups <- unique(trimws(as.character(groups)))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) < 2) {
    return(character(0))
  }

  pairs <- utils::combn(groups, 2, simplify = FALSE)
  vapply(pairs, function(pair) paste(pair[1], "vs", pair[2]), character(1))
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

# Metadata mappings are persisted as generated CSVs so pipeline runs use canonical column names.
metadata_effective_rel_path <- function(source_path = NULL, uploaded_name = NULL, mapping) {
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

read_metadata_with_mapping <- function(path, mapping) {
  md <- safe_read_table(path)
  apply_metadata_mapping_to_df(md, mapping)
}

persist_metadata_mapping <- function(source_path, uploaded_name = NULL, mapping) {
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

rel_path_from_output <- function(path, out_dir) {
  path_norm <- gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
  out_norm <- gsub("\\\\", "/", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  sub(paste0("^", out_norm, "/?"), "", path_norm)
}

# Mirrors pipeline output decisions for the UI confirmation modal; keep in sync with pipeline writers.
build_expected_output_manifest_from_config <- function(cfg, out_dir, output_level = "standard") {
  strategy <- extract_config_value(cfg, "duplicate_name_strategy")
  strategy <- safe_trimws(if (is.null(strategy)) "" else strategy)

  output_level <- normalize_output_level(output_level)
  minimal_output_enabled <- identical(output_level, "minimal")
  standard_output_enabled <- identical(output_level, "standard")
  full_debug_output_enabled <- identical(output_level, "full_debug")

  volcano_enabled <- config_flag_value(cfg, "make_volcano_plots", default = TRUE)
  heatmap_all <- config_flag_value(cfg, "make_heatmap_by_model", default = TRUE)
  heatmap_sex <- config_flag_value(cfg, "make_heatmap_by_model_sex", default = TRUE)
  sig_heatmap_all <- config_flag_value(cfg, "make_sig_heatmap_by_model", default = FALSE)
  sig_heatmap_sex <- config_flag_value(cfg, "make_sig_heatmap_by_model_sex", default = FALSE)
  sig_heatmap_fvsm <- config_flag_value(cfg, "make_sig_heatmap_FvsM_within_group", default = FALSE)
  stats_excel_enabled <- config_flag_value(cfg, "save_stats_excel_per_model", default = TRUE)
  sig_metabolites_txt <- config_flag_value(cfg, "save_sig_metabolites_txt_per_model", default = TRUE)
  metaboanalyst_export_enabled <- config_flag_value(cfg, "export_metaboanalyst_ready", default = TRUE)
  full_output_enables <- function(value) isTRUE(value) || isTRUE(full_debug_output_enabled)
  volcano_enabled <- full_output_enables(volcano_enabled)
  heatmap_all <- full_output_enables(heatmap_all)
  heatmap_sex <- full_output_enables(heatmap_sex)
  sig_heatmap_all <- full_output_enables(sig_heatmap_all)
  sig_heatmap_sex <- full_output_enables(sig_heatmap_sex)
  sig_heatmap_fvsm <- full_output_enables(sig_heatmap_fvsm)
  stats_excel_enabled <- full_output_enables(stats_excel_enabled)
  sig_metabolites_txt <- full_output_enables(sig_metabolites_txt)
  metaboanalyst_export_enabled <- full_output_enables(metaboanalyst_export_enabled)
  metaboanalyst_duplicate_only_enabled <- full_output_enables(
    config_flag_value(cfg, "export_metaboanalyst_duplicate_only", default = TRUE)
  )
  normalization_mode <- safe_trimws(extract_config_value(cfg, "normalization_mode"))
  if (!nzchar(normalization_mode)) {
    normalization_mode <- "qcrsc"
  }
  normalization_mode_key <- normalize_normalization_mode(normalization_mode)
  qc_diagnostics <- config_flag_value(cfg, "make_qc_diagnostics", default = FALSE)
  comparison_mode <- tolower(safe_trimws(extract_config_value(cfg, "comparison_mode")))
  if (length(comparison_mode) == 0 || is.na(comparison_mode[1]) || !nzchar(comparison_mode[1])) comparison_mode <- "pairwise"

  if (!nzchar(strategy)) {
    strategy <- "collapse_best_qc_rsd"
  }

  core_files <- c(
    file.path(out_dir, "PIPELINE_LOG.txt"),
    file.path(out_dir, "README.txt")
  )

  if (standard_output_enabled || full_debug_output_enabled) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "audits_global", "qc_summary.csv"),
      file.path(out_dir, "global", "audits_global", "normalization_summary.csv")
    )
  }

  if (full_debug_output_enabled) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "audits_global", "filter_summary.csv"),
      file.path(out_dir, "global", "audits_global", "missing_exclusion_audit.csv"),
      file.path(out_dir, "global", "audits_global", "presence_filter_audit.csv"),
      file.path(out_dir, "global", "audits_global", "known_filter_audit.csv"),
      file.path(out_dir, "global", "audits_global", "qc_rsd_values_pre_variants.csv"),
      file.path(out_dir, "global", "audits_global", "low_variance_iqr_audit_ACTIVE.csv"),
      file.path(out_dir, "global", "audits_global", "matrix_exports_summary.csv"),
      file.path(out_dir, "global", "audits_global", paste0("duplicate_name_audit_", strategy, ".csv")),
      file.path(out_dir, "global", "exports_global", "02_featureID_to_display_name_map.csv")
    )
  }

  if (full_debug_output_enabled && isTRUE(metaboanalyst_export_enabled)) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "exports_global", "MA_ACTIVE_raw_GLOBAL_NO_QC.csv"),
      file.path(out_dir, "global", "exports_global", "MA_ACTIVE_raw_GLOBAL_WITH_QC.csv"),
      file.path(out_dir, "global", "exports_global", "MA_COMPATIBLE_duplicate_IQR_GLOBAL_NO_QC.csv")
    )
  }

  duplicate_only_export <- full_debug_output_enabled && isTRUE(metaboanalyst_export_enabled) && isTRUE(metaboanalyst_duplicate_only_enabled)
  if (isTRUE(duplicate_only_export)) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "exports_global", "MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv"),
      file.path(out_dir, "global", "exports_global", "MA_ACTIVE_duplicate_ONLY_GLOBAL_WITH_QC.csv")
    )
  }

  if (full_debug_output_enabled && identical(strategy, "reference_or_best_qc_rsd")) {
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

  if (full_debug_output_enabled) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "exports_global", "10_MATRIX_post_<ACTIVE_VARIANT>_postRSD_preLowVar_preDup_ALL.csv"),
      file.path(out_dir, "global", "exports_global", "10_TABLE_post_<ACTIVE_VARIANT>_postRSD_preLowVar_preDup_ALL_NAMED.csv")
    )

    if (normalization_mode_key %in% c("qc_loess", "qcrsc")) {
      correction_tag <- if (identical(normalization_mode_key, "qcrsc")) "QC_RSC" else "QC_LOESS"
      correction_file <- if (identical(normalization_mode_key, "qcrsc")) "qc_rsc" else "qc_loess"
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", paste0("11_MATRIX_ACTIVE_postFilters_post", correction_tag, "_WITH_QC.csv")),
        file.path(out_dir, "global", "exports_global", paste0("12_MATRIX_ACTIVE_post", correction_tag, "_postWeight_BIOLOGICAL_NO_QC.csv")),
        file.path(out_dir, "global", "exports_global", paste0("11_", correction_file, "_qc_correction_postFilters.csv")),
        file.path(out_dir, "global", "audits_global", paste0("qc_rsd_before_after_", correction_tag, "_summary.csv")),
        file.path(out_dir, "global", "audits_global", paste0("drift_spearman_before_after_", correction_tag, "_summary.csv")),
        file.path(out_dir, "global", "audits_global", "weight_normalization_audit_biological_samples.csv")
      )
    } else if (identical(normalization_mode_key, "weight")) {
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", "11_MATRIX_ACTIVE_postFilters_postWeight_ONLY_ALL.csv")
      )
    } else if (identical(normalization_mode_key, "pqn_qc")) {
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", "11_pqn_qc_normalization_factors_postFilters.csv"),
        file.path(out_dir, "global", "exports_global", "11_pqn_qc_reference_spectrum_postFilters.csv")
      )
    } else if (identical(normalization_mode_key, "pqn_sample")) {
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", "11_pqn_sample_normalization_factors_postFilters.csv"),
        file.path(out_dir, "global", "exports_global", "11_pqn_sample_reference_spectrum_postFilters.csv")
      )
    } else if (identical(normalization_mode_key, "cyclic_loess")) {
      core_files <- c(
        core_files,
        file.path(out_dir, "global", "exports_global", "11_MATRIX_ACTIVE_postFilters_postCyclicLOESS_WITH_QC.csv")
      )
    }
  }

  if ((standard_output_enabled || full_debug_output_enabled) && normalization_mode_key %in% c("qc_loess", "qcrsc")) {
    core_files <- c(
      core_files,
      file.path(out_dir, "global", "audits_global", "qc_pca_comparison_summary.csv")
    )
  }

  optional_files <- character(0)

  if (isTRUE(stats_excel_enabled)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "stats", "*.xlsx")
    )
  }

  if (isTRUE(stats_excel_enabled) && (standard_output_enabled || full_debug_output_enabled) && comparison_mode %in% c("multigroup", "both")) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "stats", "MULTIGROUP_README.txt"),
      file.path(out_dir, "<MODEL>", "plots", "pca", "*MULTIGROUP*.png")
    )
  }

  if (full_debug_output_enabled && isTRUE(sig_metabolites_txt)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "stats", "significant", "*.txt")
    )
  }

  if (isTRUE(volcano_enabled)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "volcano", "*.png")
    )
  }

  if (!minimal_output_enabled && isTRUE(heatmap_all)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "heatmap", "*.png")
    )
  }

  if (full_debug_output_enabled && isTRUE(heatmap_sex)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "heatmap", "*.png")
    )
  }

  if (full_debug_output_enabled && isTRUE(sig_heatmap_all)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
    )
  }

  if (full_debug_output_enabled && isTRUE(sig_heatmap_sex)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
    )
  }

  if (full_debug_output_enabled && isTRUE(sig_heatmap_fvsm)) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "<MODEL>", "plots", "heatmap_significant", "*.png")
    )
  }

  optional_files <- c(
    optional_files,
    file.path(out_dir, "<MODEL>", "plots", "pca", "*.png")
  )

  if (full_debug_output_enabled) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "global", "plots_global", "normalization", "*.png")
    )
  }

  if (full_debug_output_enabled && normalization_mode_key %in% c("qc_loess", "qcrsc")) {
    optional_files <- c(
      optional_files,
      file.path(out_dir, "global", "plots_global", "pca", "*.png"),
      file.path(out_dir, "global", "plots_global", "pca", "*.csv")
    )
  }

  list(
    out_dir = out_dir,
    core_files = unique(core_files),
    optional_files = unique(optional_files)
  )
}
