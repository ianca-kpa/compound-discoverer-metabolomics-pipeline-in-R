# =============================================================================
# 06a_normalization_core.R
# Normalization core helpers
# =============================================================================

# Normalization and filtering functions for the assay data, including:
# - Weight normalization (normalize_by_weight)
# - Median normalization for comparison plots (normalize_by_median)
# - QC-LOESS drift correction (normalize_qc_loess_ref)
# - Cyclic LOESS reference plot normalization (normalize_cyclic_loess_ref)
# - QC-RSC drift correction (normalize_qcrsc_qc_ref)
# - QC-based probabilistic quotient normalization (normalize_pqn_qc_ref)
# - Sample-reference probabilistic quotient normalization (normalize_pqn_sample_ref)
# - Missing value exclusion (filter_missing_exclusion)
# - Presence filtering and imputation (presence_filter_and_impute)
# - Known-only filtering (filter_known)
# - Low-variance filtering (filter_low_variance_deterministic
# - RSD calculation (calc_rsd, calc_qc_rsd)

# Preserve a user config value named apply_qcrsc_spectral_cleaning before defining
# the function with the same name.
if (exists("apply_qcrsc_spectral_cleaning", inherits = TRUE) &&
    !is.function(get("apply_qcrsc_spectral_cleaning", inherits = TRUE))) {
  assign(
    ".apply_qcrsc_spectral_cleaning_setting",
    get("apply_qcrsc_spectral_cleaning", inherits = TRUE),
    envir = .GlobalEnv
  )
}

normalize_by_weight <- function(assay_num_raw, metadata_aligned, sample_idx,
                                use_weight_normalization = TRUE,
                                stop_on_invalid_weight = TRUE,
                                invalid_weight_to_NA = TRUE) {
  if (!isTRUE(use_weight_normalization)) {
    return(assay_num_raw)
  }

  out <- assay_num_raw
  w <- metadata_aligned$weight

  bad <- sample_idx[is.na(w[sample_idx]) | w[sample_idx] <= 0]

  if (length(bad) > 0) {
    msg <- paste0(
      "Invalid weight for biological samples: ",
      paste(metadata_aligned$sample[bad], collapse = ", "),
      " | weights: ", paste(w[bad], collapse = ", ")
    )

    if (isTRUE(stop_on_invalid_weight)) stop(msg)
    warning(msg)

    if (isTRUE(invalid_weight_to_NA)) {
      out[bad, ] <- NA_real_
    }
  }

  good <- setdiff(sample_idx, bad)
  if (length(good) > 0) {
    out[good, ] <- sweep(out[good, , drop = FALSE], 1, w[good], "/")
  }

  out
}

normalize_by_median <- function(assay_num_raw) {
  out <- assay_num_raw

  row_medians <- apply(assay_num_raw, 1, function(v) {
    v <- v[is.finite(v) & v > 0]
    if (length(v) == 0) return(NA_real_)
    stats::median(v)
  })

  ref_median <- stats::median(row_medians, na.rm = TRUE)
  valid_rows <- is.finite(row_medians) & row_medians > 0 & is.finite(ref_median) & ref_median > 0

  if (any(valid_rows)) {
    out[valid_rows, ] <- sweep(out[valid_rows, , drop = FALSE], 1, row_medians[valid_rows], "/") * ref_median
  }

  out
}

normalize_cyclic_loess_ref <- function(assay_num_raw, log2_offset = 1,
                                       output_scale = c("log2", "linear")) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for cyclic LOESS normalization.")
  }

  output_scale <- match.arg(output_scale)
  mat_log <- log2_transform(t(assay_num_raw), log2_offset)
  mat_norm <- limma::normalizeCyclicLoess(mat_log, method = "fast")
  assay_num_cyclic_loess <- t(as.matrix(mat_norm))
  if (identical(output_scale, "linear")) {
    assay_num_cyclic_loess <- (2^assay_num_cyclic_loess) - log2_offset
  }
  colnames(assay_num_cyclic_loess) <- colnames(assay_num_raw)
  rownames(assay_num_cyclic_loess) <- rownames(assay_num_raw)

  assay_num_cyclic_loess
}

normalize_qc_loess_ref <- function(assay_num_weight, qc_idx,
                                   qc_loess_span = 0.75,
                                   min_qc_points = 4,
                                   injection_order = NULL) {
  assay_num_qc_loess <- assay_num_weight
  feature_ids <- colnames(assay_num_weight)
  if (is.null(injection_order)) {
    injection_order <- seq_len(nrow(assay_num_weight))
  } else {
    injection_order <- as.numeric(injection_order)
    if (length(injection_order) != nrow(assay_num_weight)) {
      stop("Injection_order must have one value per sample row.")
    }
  }

  audit <- vector("list", length(feature_ids))

  for (j in seq_along(feature_ids)) {
    values <- assay_num_weight[, j]
    qc_values <- values[qc_idx]
    qc_order <- injection_order[qc_idx]
    usable <- is.finite(qc_values) & qc_values > 0
    n_qc_used <- sum(usable)

    correction_factor <- rep(1, length(values))
    valid_qc_loess <- FALSE
    ref_intensity <- NA_real_

    if (n_qc_used >= max(3, min_qc_points) && length(unique(qc_order[usable])) >= 3) {
      fit_df <- data.frame(
        injection_order = qc_order[usable],
        log_intensity = log(qc_values[usable])
      )

      fit <- tryCatch(
        stats::loess(
          log_intensity ~ injection_order,
          data = fit_df,
          span = qc_loess_span,
          degree = 1,
          control = stats::loess.control(surface = "direct")
        ),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        predicted_log <- tryCatch(
          stats::predict(fit, newdata = data.frame(injection_order = injection_order)),
          error = function(e) rep(NA_real_, length(values))
        )

        if (sum(is.finite(predicted_log)) >= 3) {
          ref_log <- stats::median(predicted_log[qc_idx][is.finite(predicted_log[qc_idx])], na.rm = TRUE)
          if (is.finite(ref_log)) {
            correction_factor <- exp(predicted_log - ref_log)
            correction_factor[!is.finite(correction_factor) | correction_factor <= 0] <- 1
            positive_rows <- is.finite(values) & values > 0
            assay_num_qc_loess[positive_rows, j] <- values[positive_rows] / correction_factor[positive_rows]
            valid_qc_loess <- TRUE
            ref_intensity <- exp(ref_log)
          }
        }
      }
    }

    audit[[j]] <- tibble(
      featureID = feature_ids[j],
      n_qc_used = n_qc_used,
      valid_qc_loess = valid_qc_loess,
      qc_loess_ref_intensity = ref_intensity,
      median_qc_intensity = stats::median(qc_values[usable], na.rm = TRUE),
      correction_factor_min = min(correction_factor, na.rm = TRUE),
      correction_factor_median = stats::median(correction_factor, na.rm = TRUE),
      correction_factor_max = max(correction_factor, na.rm = TRUE)
    )
  }

  list(
    assay_num_qc_loess = assay_num_qc_loess,
    qc_loess_tbl = dplyr::bind_rows(audit)
  )
}

normalize_qcrsc_qc_ref <- function(assay_num_weight, qc_idx,
                                   min_qc_points = 4,
                                   injection_order = NULL,
                                   batch = NULL,
                                   max_iter = 5) {
  if (!requireNamespace("pmp", quietly = TRUE)) {
    stop("Package 'pmp' is required for QC-RSC normalization.")
  }

  if (is.null(injection_order)) {
    injection_order <- seq_len(nrow(assay_num_weight))
  } else {
    injection_order <- as.numeric(injection_order)
    if (length(injection_order) != nrow(assay_num_weight)) {
      stop("QC-RSC injection_order must have one value per sample row.")
    }
  }

  classes <- rep("Sample", nrow(assay_num_weight))
  classes[qc_idx] <- "QC"
  if (is.null(batch)) {
    batch <- rep(1L, nrow(assay_num_weight))
  } else {
    if (length(batch) != nrow(assay_num_weight)) {
      stop("QC-RSC batch must have one value per sample row.")
    }
    batch_labels <- trimws(as.character(batch))
    batch_missing <- is.na(batch_labels) | !nzchar(batch_labels)
    if (all(batch_missing)) {
      batch <- rep(1L, nrow(assay_num_weight))
    } else {
      if (any(batch_missing)) {
        stop("QC-RSC batch is missing for ", sum(batch_missing), " sample(s).")
      }
      batch <- as.integer(factor(batch_labels, levels = unique(batch_labels)))
    }
  }

  matrix_input <- t(as.matrix(assay_num_weight))
  mode(matrix_input) <- "numeric"

  corrected <- pmp::QCRSC(
    df = matrix_input,
    order = injection_order,
    batch = batch,
    classes = classes,
    spar = 0,
    log = TRUE,
    minQC = min_qc_points,
    qc_label = "QC"
  )

  corrected_mat <- t(as.matrix(corrected))
  colnames(corrected_mat) <- colnames(assay_num_weight)
  rownames(corrected_mat) <- rownames(assay_num_weight)

  qc_rsd_before <- calc_qc_rsd(assay_num_weight, qc_idx)

  qcrsc_tbl <- tibble::tibble(
    featureID = names(qc_rsd_before),
    qc_rsd_before = as.numeric(qc_rsd_before),
    qc_rsd_after = as.numeric(calc_qc_rsd(corrected_mat, qc_idx)),
    qc_rsc_applied = TRUE,
    qcrsc_source = "pmp::QCRSC",
    qcrsc_batch_count = length(unique(batch))
  )

  list(
    assay_num_qcrsc = corrected_mat,
    qcrsc_tbl = qcrsc_tbl
  )
}

# Removes features with strong QC/batch artifacts before QC-RSC correction.
apply_qcrsc_spectral_cleaning <- function(assay_num,
                                         metadata_aligned,
                                         qc_idx,
                                         batch = NULL,
                                         kruskal_p_cutoff = 1e-4,
                                         wilcoxon_p_cutoff = 1e-14,
                                         max_qc_rsd = 20) {
  if (is.null(batch) && "batch" %in% names(metadata_aligned)) {
    batch <- metadata_aligned$batch
  }

  if (is.null(batch)) {
    batch <- rep("1", nrow(metadata_aligned))
  }

  batch <- trimws(as.character(batch))
  batch[is.na(batch) | !nzchar(batch)] <- "1"

  is_qc <- seq_len(nrow(assay_num)) %in% qc_idx
  batch_levels <- unique(batch[!is.na(batch)])
  keep_feature <- rep(TRUE, ncol(assay_num))
  feature_ids <- colnames(assay_num)

  audit <- vector("list", length(feature_ids))

  for (j in seq_along(feature_ids)) {
    values <- assay_num[, j]
    qc_values <- values[is_qc]
    qc_rsd <- calc_rsd(qc_values)

    kruskal_p <- NA_real_
    kruskal_used <- FALSE
    kruskal_flag <- FALSE

    wilcoxon_p <- NA_real_
    wilcoxon_used <- FALSE
    wilcoxon_flag <- FALSE

    qc_rsd_flag <- is.finite(qc_rsd) && qc_rsd > max_qc_rsd

    if (length(batch_levels) >= 2) {
      bio_values <- values[!is_qc]
      bio_batches <- batch[!is_qc]
      bio_ok <- is.finite(bio_values) & !is.na(bio_batches)

      if (sum(bio_ok) > 0) {
        bio_values <- bio_values[bio_ok]
        bio_batches <- bio_batches[bio_ok]
        batch_counts <- table(bio_batches)
        if (sum(batch_counts > 0) >= 2) {
          kruskal_used <- TRUE
          kruskal_p <- tryCatch(
            stats::kruskal.test(bio_values, factor(bio_batches))$p.value,
            error = function(e) NA_real_
          )
          kruskal_flag <- is.finite(kruskal_p) && kruskal_p < kruskal_p_cutoff
        }

        if (!kruskal_flag) {
          batch_diffs <- vapply(batch_levels, function(b) {
            qc_b <- values[is_qc & batch == b]
            bio_b <- values[!is_qc & batch == b]
            qc_b <- qc_b[is.finite(qc_b)]
            bio_b <- bio_b[is.finite(bio_b)]
            if (length(qc_b) == 0 || length(bio_b) == 0) {
              return(NA_real_)
            }
            stats::median(qc_b, na.rm = TRUE) - stats::median(bio_b, na.rm = TRUE)
          }, numeric(1))

          batch_diffs <- batch_diffs[is.finite(batch_diffs)]
          if (length(batch_diffs) >= 2) {
            wilcoxon_used <- TRUE
            wilcoxon_p <- tryCatch(
              stats::wilcox.test(batch_diffs, mu = 0, exact = FALSE)$p.value,
              error = function(e) NA_real_
            )
            wilcoxon_flag <- is.finite(wilcoxon_p) && wilcoxon_p < wilcoxon_p_cutoff
          }
        }
      }
    }

    keep_feature[j] <- !(kruskal_flag || wilcoxon_flag || qc_rsd_flag)

    audit[[j]] <- tibble::tibble(
      featureID = feature_ids[j],
      batch_count = length(batch_levels),
      kruskal_used = kruskal_used,
      kruskal_p = kruskal_p,
      kruskal_flag = kruskal_flag,
      wilcoxon_used = wilcoxon_used,
      wilcoxon_p = wilcoxon_p,
      wilcoxon_flag = wilcoxon_flag,
      qc_rsd = qc_rsd,
      qc_rsd_flag = qc_rsd_flag,
      keep_feature = keep_feature[j],
      removal_reason = dplyr::case_when(
        kruskal_flag ~ "kruskal_wallis",
        wilcoxon_flag ~ "wilcoxon_signed_rank",
        qc_rsd_flag ~ "qc_rsd_gt_threshold",
        TRUE ~ "kept"
      )
    )
  }

  list(
    assay = assay_num[, keep_feature, drop = FALSE],
    keep_feature_ids = feature_ids[keep_feature],
    audit = dplyr::bind_rows(audit)
  )
}

normalize_pqn_qc_ref <- function(assay_num_weight, qc_idx, min_qc_points = 3) {
  assay_num_pqn <- assay_num_weight
  feature_ids <- colnames(assay_num_weight)
  sample_ids <- rownames(assay_num_weight)

  if (length(qc_idx) < min_qc_points) {
    stop(
      "PQN QC normalization requires at least ", min_qc_points,
      " QC samples, but found ", length(qc_idx), "."
    )
  }

  ref <- apply(assay_num_weight[qc_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  valid_ref <- is.finite(ref) & ref > 0

  if (sum(valid_ref) < 2) {
    stop("PQN QC normalization could not build a valid QC reference spectrum.")
  }

  factors <- rep(NA_real_, nrow(assay_num_weight))
  n_features_used <- integer(nrow(assay_num_weight))

  for (i in seq_len(nrow(assay_num_weight))) {
    values <- assay_num_weight[i, ]
    usable <- valid_ref & is.finite(values) & values > 0
    n_features_used[i] <- sum(usable)

    if (n_features_used[i] < 2) {
      next
    }

    quotient <- values[usable] / ref[usable]
    factor <- stats::median(quotient[is.finite(quotient) & quotient > 0], na.rm = TRUE)

    if (is.finite(factor) && factor > 0) {
      factors[i] <- factor
      assay_num_pqn[i, ] <- values / factor
    }
  }

  audit <- tibble(
    sample = sample_ids,
    pqn_factor = factors,
    n_features_used = n_features_used,
    normalized = is.finite(factors) & factors > 0
  )

  list(
    assay_num_pqn = assay_num_pqn,
    pqn_tbl = audit,
    pqn_reference_tbl = tibble(
      featureID = feature_ids,
      qc_reference_median = as.numeric(ref),
      valid_reference = valid_ref
    )
  )
}

normalize_pqn_sample_ref <- function(assay_num_weight, sample_idx, min_sample_points = 2) {
  assay_num_pqn <- assay_num_weight
  feature_ids <- colnames(assay_num_weight)
  sample_ids <- rownames(assay_num_weight)
  reference_idx <- intersect(sample_idx, seq_len(nrow(assay_num_weight)))

  if (length(reference_idx) < min_sample_points) {
    stop(
      "PQN sample normalization requires at least ", min_sample_points,
      " biological samples, but found ", length(reference_idx), "."
    )
  }

  ref <- apply(assay_num_weight[reference_idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  valid_ref <- is.finite(ref) & ref > 0

  if (sum(valid_ref) < 2) {
    stop("PQN sample normalization could not build a valid sample reference spectrum.")
  }

  factors <- rep(NA_real_, nrow(assay_num_weight))
  n_features_used <- integer(nrow(assay_num_weight))

  for (i in seq_len(nrow(assay_num_weight))) {
    values <- assay_num_weight[i, ]
    usable <- valid_ref & is.finite(values) & values > 0
    n_features_used[i] <- sum(usable)

    if (n_features_used[i] < 2) {
      next
    }

    quotient <- values[usable] / ref[usable]
    factor <- stats::median(quotient[is.finite(quotient) & quotient > 0], na.rm = TRUE)

    if (is.finite(factor) && factor > 0) {
      factors[i] <- factor
      assay_num_pqn[i, ] <- values / factor
    }
  }

  audit <- tibble(
    sample = sample_ids,
    pqn_factor = factors,
    n_features_used = n_features_used,
    normalized = is.finite(factors) & factors > 0
  )

  list(
    assay_num_pqn = assay_num_pqn,
    pqn_tbl = audit,
    pqn_reference_tbl = tibble(
      featureID = feature_ids,
      sample_reference_median = as.numeric(ref),
      valid_reference = valid_ref
    )
  )
}

# Lower scores indicate more stable sample distributions and QC behavior.
score_normalization_candidate <- function(assay_num,
                                          qc_idx = integer(0),
                                          injection_order = NULL,
                                          log2_offset = 1) {
  mat_log <- suppressWarnings(log2_transform(assay_num, log2_offset))
  finite_fraction <- mean(is.finite(mat_log))
  nonfinite_fraction <- 1 - finite_fraction

  sample_medians <- apply(mat_log, 1, stats::median, na.rm = TRUE)
  sample_iqrs <- apply(mat_log, 1, stats::IQR, na.rm = TRUE)
  sample_median_mad <- stats::median(abs(sample_medians - stats::median(sample_medians, na.rm = TRUE)), na.rm = TRUE)
  sample_iqr_mad <- stats::median(abs(sample_iqrs - stats::median(sample_iqrs, na.rm = TRUE)), na.rm = TRUE)

  qc_rsd_median <- NA_real_
  qc_drift_median_abs_cor <- NA_real_
  score <- sample_median_mad + (0.5 * sample_iqr_mad) + (100 * nonfinite_fraction)
  score_basis <- "sample_distribution"

  if (length(qc_idx) >= 3) {
    qc_rsd <- calc_qc_rsd(assay_num, qc_idx)
    qc_rsd_median <- stats::median(qc_rsd[is.finite(qc_rsd)], na.rm = TRUE)

    if (!is.null(injection_order)) {
      injection_order <- as.numeric(injection_order)
      if (length(injection_order) == nrow(assay_num)) {
        qc_order <- injection_order[qc_idx]
        drift <- apply(mat_log[qc_idx, , drop = FALSE], 2, function(v) {
          ok <- is.finite(v) & is.finite(qc_order)
          if (sum(ok) < 3 || length(unique(qc_order[ok])) < 3) return(NA_real_)
          suppressWarnings(stats::cor(v[ok], qc_order[ok], method = "spearman"))
        })
        qc_drift_median_abs_cor <- stats::median(abs(drift[is.finite(drift)]), na.rm = TRUE)
      }
    }

    if (!is.finite(qc_rsd_median)) qc_rsd_median <- 999
    if (!is.finite(qc_drift_median_abs_cor)) qc_drift_median_abs_cor <- 1

    score <- (qc_rsd_median / 100) +
      qc_drift_median_abs_cor +
      (0.25 * sample_median_mad) +
      (0.10 * sample_iqr_mad) +
      (100 * nonfinite_fraction)
    score_basis <- "qc_rsd_drift_distribution"
  }

  tibble::tibble(
    score = score,
    score_basis = score_basis,
    qc_rsd_median = qc_rsd_median,
    qc_drift_median_abs_cor = qc_drift_median_abs_cor,
    sample_median_mad = sample_median_mad,
    sample_iqr_mad = sample_iqr_mad,
    finite_fraction = finite_fraction
  )
}
