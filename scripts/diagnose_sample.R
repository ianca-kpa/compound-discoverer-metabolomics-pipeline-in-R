#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

default_input_csv <- "data/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv"
settings_path <- file.path("pipeline", "config", "settings.R")
if (file.exists(settings_path)) {
  settings_env <- new.env(parent = baseenv())
  try(source(settings_path, local = settings_env), silent = TRUE)
  configured_output_dir <- get0("output_dir", ifnotfound = "", envir = settings_env)
  configured_export <- file.path(configured_output_dir, "global", "exports_global", "MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv")
  if (nzchar(configured_output_dir) && file.exists(configured_export)) {
    default_input_csv <- configured_export
  }
}

sample_target <- ifelse(length(args) >= 1 && nzchar(args[1]), args[1], "PCA_OUTLIERS")
input_csv <- ifelse(
  length(args) >= 2 && nzchar(args[2]),
  args[2],
  default_input_csv
)
metadata_path <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], "")
out_dir <- ifelse(
  length(args) >= 4 && nzchar(args[4]),
  args[4],
  file.path(
    "output",
    if (toupper(sample_target) %in% c("PCA_OUTLIERS", "OUTLIERS", "AUTO")) {
      "diagnostics_pca_outliers"
    } else {
      paste0("diagnostics_", sample_target)
    }
  )
)
outlier_threshold <- ifelse(
  length(args) >= 5 && nzchar(args[5]),
  suppressWarnings(as.numeric(args[5])),
  3.5
)
if (!is.finite(outlier_threshold) || outlier_threshold <= 0) {
  outlier_threshold <- 3.5
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv)
}

helpers_path <- file.path("pipeline", "R", "03_helpers_io_log.R")
features_path <- file.path("pipeline", "R", "05_features_assay.R")
if (file.exists(helpers_path)) source(helpers_path)
if (file.exists(features_path)) source(features_path)

clean_sample_name <- function(x) {
  if (exists("clean_sample_from_area_col", mode = "function")) {
    return(clean_sample_from_area_col(x))
  }

  s <- sub("^Area\\s*:\\s*", "", x)
  s <- sub("(?i)\\.raw.*$", "", s, perl = TRUE)
  s <- sub("\\s*\\(.*\\)$", "", s)
  s <- str_trim(s)
  if (str_detect(s, "^QC")) return(s)
  if (str_detect(s, "_")) s <- str_split_fixed(s, "_", 2)[, 1]
  s
}

read_metadata <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("Metadata file not found: ", path)

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read metadata Excel files.")
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }

  as.data.frame(readr::read_csv(path, show_col_types = FALSE), stringsAsFactors = FALSE)
}

write_sample_metadata <- function(metadata, sample_name, out_dir) {
  if (is.null(metadata)) return(invisible(NULL))

  names(metadata) <- tolower(gsub("[^a-z0-9]+", "_", names(metadata)))
  possible_sample_cols <- c("sample", "sample_id", "id", "name", "sample_name")
  sample_col <- intersect(possible_sample_cols, names(metadata))[1]

  if (is.na(sample_col) || is.null(sample_col)) {
    cat("Could not find a sample/id column in metadata. Columns:", paste(names(metadata), collapse = ", "), "\n")
    return(invisible(NULL))
  }

  hit <- which(toupper(as.character(metadata[[sample_col]])) == toupper(sample_name))[1]
  if (is.na(hit)) {
    cat(sample_name, "not found in metadata column", sample_col, "\n")
    return(invisible(NULL))
  }

  row_sample <- metadata[hit, , drop = FALSE]
  out_path <- file.path(out_dir, paste0("metadata_", sample_name, "_row.csv"))
  write.csv(row_sample, out_path, row.names = FALSE)
  cat("Metadata row written to", out_path, "\n")
}

safe_file_part <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

prepare_sample_pca_matrix <- function(mat) {
  mat_pca <- mat

  finite_vals <- mat_pca[is.finite(mat_pca)]
  if (length(finite_vals) > 0 && min(finite_vals, na.rm = TRUE) >= 0 && stats::median(finite_vals, na.rm = TRUE) > 50) {
    mat_pca <- log2(mat_pca + 1)
  }

  keep_features <- apply(mat_pca, 1, function(v) {
    finite <- is.finite(v)
    sum(finite) >= 3 && stats::sd(v[finite], na.rm = TRUE) > 0
  })
  mat_pca <- mat_pca[keep_features, , drop = FALSE]

  if (nrow(mat_pca) < 2 || ncol(mat_pca) < 3) {
    return(NULL)
  }

  feature_medians <- apply(mat_pca, 1, stats::median, na.rm = TRUE)
  for (i in seq_len(nrow(mat_pca))) {
    bad <- !is.finite(mat_pca[i, ])
    if (any(bad)) {
      mat_pca[i, bad] <- feature_medians[i]
    }
  }

  sample_mat <- t(mat_pca)
  sample_mat <- scale(sample_mat, center = TRUE, scale = TRUE)
  sample_mat <- sample_mat[, apply(sample_mat, 2, function(v) all(is.finite(v)) && stats::sd(v) > 0), drop = FALSE]

  if (nrow(sample_mat) < 3 || ncol(sample_mat) < 2) {
    return(NULL)
  }

  sample_mat
}

detect_pca_outliers <- function(mat, samples, threshold = 3.5) {
  sample_mat <- prepare_sample_pca_matrix(mat)
  if (is.null(sample_mat)) {
    stop("Could not build PCA matrix with at least 3 samples and 2 variable features.")
  }

  pca <- stats::prcomp(sample_mat, center = FALSE, scale. = FALSE)
  if (is.null(pca$x) || ncol(pca$x) < 2) {
    stop("PCA did not return PC1 and PC2.")
  }

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$sample <- rownames(scores)
  scores <- scores[, c("sample", "PC1", "PC2")]

  center_pc1 <- stats::median(scores$PC1, na.rm = TRUE)
  center_pc2 <- stats::median(scores$PC2, na.rm = TRUE)
  scale_pc1 <- stats::mad(scores$PC1, center = center_pc1, constant = 1.4826, na.rm = TRUE)
  scale_pc2 <- stats::mad(scores$PC2, center = center_pc2, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale_pc1) || scale_pc1 == 0) scale_pc1 <- stats::sd(scores$PC1, na.rm = TRUE)
  if (!is.finite(scale_pc2) || scale_pc2 == 0) scale_pc2 <- stats::sd(scores$PC2, na.rm = TRUE)
  if (!is.finite(scale_pc1) || scale_pc1 == 0) scale_pc1 <- 1
  if (!is.finite(scale_pc2) || scale_pc2 == 0) scale_pc2 <- 1

  scores$robust_z_pc1 <- (scores$PC1 - center_pc1) / scale_pc1
  scores$robust_z_pc2 <- (scores$PC2 - center_pc2) / scale_pc2
  scores$pca_robust_distance <- sqrt(scores$robust_z_pc1^2 + scores$robust_z_pc2^2)
  scores$is_pca_outlier <- scores$pca_robust_distance >= threshold
  scores <- scores[order(scores$pca_robust_distance, decreasing = TRUE), ]

  scores
}

write_one_sample_diagnostic <- function(sample_idx,
                                        samples,
                                        mat,
                                        feature_names,
                                        input_csv,
                                        out_dir,
                                        metadata = NULL) {
  sample_name <- samples[sample_idx]
  vals <- as.numeric(mat[, sample_idx])

  cat("\nDiagnosing sample:", sample_name, "\n")
  cat("Input file:", input_csv, "\n")
  cat("Output directory:", out_dir, "\n\n")

  cat("Summary values:\n")
  print(summary(vals))
  cat("NA count:", sum(is.na(vals)), "/", length(vals), "\n")
  cat("Zeros or <=0 count:", sum(!is.finite(vals) | vals <= 0), "\n")

  sample_totals <- colSums(mat, na.rm = TRUE)
  missing_frac <- apply(mat, 2, function(v) mean(is.na(v) | !is.finite(v) | v <= 0))

  cat("\nTotal intensity:", sample_totals[sample_name], "\n")
  cat("Missing fraction:", missing_frac[sample_name], "\n")

  feature_medians <- apply(mat, 1, stats::median, na.rm = TRUE)
  feature_mean <- apply(mat, 1, mean, na.rm = TRUE)
  feature_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)
  threshold <- feature_mean + 3 * feature_sd
  n_extreme <- sum(vals > threshold, na.rm = TRUE)

  ord <- order(vals, decreasing = TRUE)
  topn <- head(ord, 50)
  top_tbl <- tibble::tibble(
    rank = seq_along(topn),
    feature_idx = topn,
    feature = feature_names[topn],
    sample_value = vals[topn],
    median_all = feature_medians[topn],
    fold_vs_median = ifelse(
      is.finite(feature_medians[topn]) & feature_medians[topn] > 0,
      vals[topn] / feature_medians[topn],
      NA_real_
    )
  )

  sample_file <- safe_file_part(sample_name)
  top_path <- file.path(out_dir, paste0("top_features_", sample_file, ".csv"))
  write.csv(top_tbl, top_path, row.names = FALSE)
  cat("\nTop features written to", top_path, "\n")
  print(head(top_tbl, 10))

  summary_path <- file.path(out_dir, paste0("diagnose_", sample_file, "_summary.csv"))
  summary_tbl <- data.frame(
    sample = sample_name,
    input_csv = input_csv,
    total_intensity = as.numeric(sample_totals[sample_name]),
    missing_fraction = as.numeric(missing_frac[sample_name]),
    na_count = sum(is.na(vals)),
    zero_or_nonpos = sum(!is.finite(vals) | vals <= 0),
    n_extreme = n_extreme,
    top_features_csv = top_path
  )
  write.csv(summary_tbl, summary_path, row.names = FALSE)

  write_sample_metadata(metadata, sample_name, out_dir)

  cat("Summary written to", summary_path, "\n")
  summary_tbl
}

input_df <- readr::read_csv(input_csv, show_col_types = FALSE)
area_cols <- grep("^Area", names(input_df), value = TRUE)
embedded_metadata <- NULL

if (length(area_cols) > 0) {
  sample_cols <- area_cols
  samples <- vapply(sample_cols, clean_sample_name, FUN.VALUE = character(1))
  value_df <- input_df %>%
    dplyr::select(dplyr::all_of(sample_cols)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ suppressWarnings(as.numeric(.x))))
  colnames(value_df) <- samples
  mat <- as.matrix(value_df)
  storage.mode(mat) <- "numeric"
  feature_names <- if ("Name" %in% names(input_df)) {
    as.character(input_df$Name)
  } else {
    paste0("row", seq_len(nrow(input_df)))
  }
} else {
  name_key <- tolower(names(input_df))
  sample_col_idx <- match(TRUE, name_key %in% c("sample", "sample_id", "id", "name", "sample_name"))

  is_numeric_like <- vapply(input_df, function(col) {
    if (is.numeric(col)) return(TRUE)
    vals <- suppressWarnings(as.numeric(as.character(col)))
    mean(is.na(vals)) < 0.5
  }, logical(1))

  if (!is.na(sample_col_idx)) {
    samples <- as.character(input_df[[sample_col_idx]])
    numeric_cols <- setdiff(which(is_numeric_like), sample_col_idx)
    value_df <- input_df[, numeric_cols, drop = FALSE]
    value_df[] <- lapply(value_df, function(col) suppressWarnings(as.numeric(as.character(col))))
    feature_names <- names(value_df)
    mat <- t(as.matrix(value_df))
    storage.mode(mat) <- "numeric"
    colnames(mat) <- samples
    rownames(mat) <- feature_names

    metadata_cols <- setdiff(seq_along(input_df), numeric_cols)
    embedded_metadata <- as.data.frame(input_df[, metadata_cols, drop = FALSE], stringsAsFactors = FALSE)
  } else {
    value_df <- input_df[, is_numeric_like, drop = FALSE]
    value_df[] <- lapply(value_df, function(col) suppressWarnings(as.numeric(as.character(col))))
    samples <- names(value_df)
    mat <- as.matrix(value_df)
    storage.mode(mat) <- "numeric"
    feature_names <- paste0("row", seq_len(nrow(input_df)))
  }
}

sample_totals <- colSums(mat, na.rm = TRUE)
missing_frac <- apply(mat, 2, function(v) mean(is.na(v) | !is.finite(v) | v <= 0))

sample_summary <- data.frame(
  sample = names(sample_totals),
  total = as.numeric(sample_totals),
  missing_frac = as.numeric(missing_frac)
)
write.csv(sample_summary, file.path(out_dir, "sample_totals_missing_frac.csv"), row.names = FALSE)

metadata <- read_metadata(metadata_path)
if (is.null(metadata) && !is.null(embedded_metadata)) {
  metadata <- embedded_metadata
}

outlier_mode <- toupper(sample_target) %in% c("PCA_OUTLIERS", "OUTLIERS", "AUTO")

if (isTRUE(outlier_mode)) {
  cat("Detecting PCA outliers using robust PC1/PC2 distance. Threshold:", outlier_threshold, "\n")
  pca_scores <- detect_pca_outliers(mat, samples, threshold = outlier_threshold)
  scores_path <- file.path(out_dir, "pca_outlier_candidates.csv")
  write.csv(pca_scores, scores_path, row.names = FALSE)
  cat("PCA outlier scores written to", scores_path, "\n")

  outlier_samples <- pca_scores$sample[pca_scores$is_pca_outlier]
  if (length(outlier_samples) == 0) {
    cat("No PCA outliers found at threshold", outlier_threshold, "\n")
    cat("Most distant samples:\n")
    print(head(pca_scores, 10))
    quit(status = 0)
  }

  cat("PCA outliers found:", paste(outlier_samples, collapse = ", "), "\n")
  sample_indices <- match(outlier_samples, samples)
} else {
  sample_matches <- which(toupper(samples) == toupper(sample_target))
  if (length(sample_matches) == 0) {
    sample_matches <- grep(sample_target, samples, ignore.case = TRUE)
  }

  if (length(sample_matches) == 0) {
    cat("Sample", sample_target, "not found. Samples available:\n")
    print(samples)
    quit(status = 0)
  }

  sample_indices <- sample_matches[1]
}

diagnostic_summaries <- lapply(
  sample_indices,
  write_one_sample_diagnostic,
  samples = samples,
  mat = mat,
  feature_names = feature_names,
  input_csv = input_csv,
  out_dir = out_dir,
  metadata = metadata
)

combined_summary <- dplyr::bind_rows(diagnostic_summaries)
combined_path <- file.path(out_dir, "diagnose_samples_summary.csv")
write.csv(combined_summary, combined_path, row.names = FALSE)

cat("\nCombined summary written to", combined_path, "\n")
cat("Done. Generated files:\n")
print(list.files(out_dir, full.names = TRUE))
