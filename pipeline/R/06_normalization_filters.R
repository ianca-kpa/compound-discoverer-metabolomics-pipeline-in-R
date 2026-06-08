# =============================================================================
# 06_normalization_filters.R
# Normalization + filters
# =============================================================================

# Normalization and filtering functions for the assay data, including:
# - Weight normalization (normalize_by_weight)
# - Probabilistic Quotient Normalization (normalize_pqn_qc_ref)
# - QC-LOESS signal drift correction (normalize_qc_loess)
# - Missing value exclusion (filter_missing_exclusion)
# - Presence filtering and imputation (presence_filter_and_impute)
# - Known-only filtering (filter_known)
# - Low-variance filtering (filter_low_variance_deterministic
# - RSD calculation (calc_rsd, calc_qc_rsd)

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

normalize_pqn_qc_ref <- function(assay_num_weight, qc_idx) {
  qc_ref <- apply(assay_num_weight[qc_idx, , drop = FALSE], 2, median, na.rm = TRUE)
  qc_ref[qc_ref == 0] <- NA_real_

  ratio <- sweep(assay_num_weight, 2, qc_ref, "/")
  pqn_factor <- apply(ratio, 1, median, na.rm = TRUE)

  valid <- is.finite(pqn_factor) & pqn_factor > 0
  pqn_used <- pqn_factor
  pqn_used[!valid] <- 1

  assay_num_pqn <- sweep(assay_num_weight, 1, pqn_used, "/")

  list(
    assay_num_pqn = assay_num_pqn,
    pqn_tbl = tibble(
      pqn_factor = pqn_factor,
      pqn_factor_used_for_norm = pqn_used,
      valid_pqn = valid
    )
  )
}

normalize_qc_loess <- function(assay_num_weight, metadata_aligned, qc_idx,
                               injection_order = NULL,
                               min_qc_points = 5,
                               loess_span = 0.75) {
  if (length(qc_idx) < min_qc_points) {
    stop(
      "QC-LOESS normalization requires at least ", min_qc_points,
      " QC samples. Found: ", length(qc_idx), "."
    )
  }

  if (is.null(injection_order)) {
    if ("injection_order" %in% names(metadata_aligned)) {
      injection_order <- metadata_aligned$injection_order
    } else {
      injection_order <- seq_len(nrow(assay_num_weight))
    }
  }

  injection_order <- suppressWarnings(as.numeric(injection_order))
  if (length(injection_order) != nrow(assay_num_weight)) {
    stop("QC-LOESS injection_order length must match the number of assay rows.")
  }

  missing_order <- !is.finite(injection_order)
  if (any(missing_order)) {
    injection_order[missing_order] <- seq_len(nrow(assay_num_weight))[missing_order]
  }

  assay_num_loess <- assay_num_weight
  audit_tbl <- tibble(
    featureID = colnames(assay_num_weight),
    qc_points_used = integer(ncol(assay_num_weight)),
    qc_median = numeric(ncol(assay_num_weight)),
    corrected = logical(ncol(assay_num_weight)),
    correction_reason = character(ncol(assay_num_weight))
  )

  for (j in seq_len(ncol(assay_num_weight))) {
    y_qc <- assay_num_weight[qc_idx, j]
    x_qc <- injection_order[qc_idx]
    good_qc <- is.finite(x_qc) & is.finite(y_qc) & y_qc > 0
    n_good <- sum(good_qc)
    audit_tbl$qc_points_used[j] <- n_good

    if (n_good < min_qc_points) {
      audit_tbl$corrected[j] <- FALSE
      audit_tbl$correction_reason[j] <- paste0("qc_points<", min_qc_points)
      audit_tbl$qc_median[j] <- NA_real_
      next
    }

    qc_median <- stats::median(y_qc[good_qc], na.rm = TRUE)
    audit_tbl$qc_median[j] <- qc_median

    if (!is.finite(qc_median) || qc_median <= 0) {
      audit_tbl$corrected[j] <- FALSE
      audit_tbl$correction_reason[j] <- "invalid_qc_median"
      next
    }

    fit_df <- data.frame(
      x = x_qc[good_qc],
      y = y_qc[good_qc]
    )

    fit <- tryCatch(
      stats::loess(
        y ~ x,
        data = fit_df,
        span = loess_span,
        degree = 1,
        family = "symmetric",
        na.action = stats::na.exclude,
        control = stats::loess.control(surface = "direct")
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      audit_tbl$corrected[j] <- FALSE
      audit_tbl$correction_reason[j] <- "loess_fit_failed"
      next
    }

    pred <- tryCatch(
      stats::predict(fit, newdata = data.frame(x = injection_order)),
      error = function(e) rep(NA_real_, length(injection_order))
    )

    valid_pred <- is.finite(pred) & pred > 0
    if (!any(valid_pred)) {
      audit_tbl$corrected[j] <- FALSE
      audit_tbl$correction_reason[j] <- "invalid_loess_prediction"
      next
    }

    factor <- qc_median / pred
    factor[!valid_pred] <- 1
    assay_num_loess[, j] <- assay_num_weight[, j] * factor
    audit_tbl$corrected[j] <- TRUE
    audit_tbl$correction_reason[j] <- "corrected"
  }

  list(
    assay_num_loess = assay_num_loess,
    loess_tbl = audit_tbl
  )
}

filter_missing_exclusion <- function(assay_num, feature_tbl, sample_idx, max_missing_fraction = 0.5, audit_path = NULL) {
  if (is.na(max_missing_fraction) || max_missing_fraction >= 1) {
    return(list(assay = assay_num, feature = feature_tbl))
  }

  miss <- apply(assay_num[sample_idx, , drop = FALSE], 2, function(v) mean(is.na(v) | v == 0))

  audit <- tibble(featureID = names(miss), missing_fraction = miss) %>%
    mutate(
      kept = missing_fraction <= max_missing_fraction,
      exclusion_reason = if_else(kept, "kept", paste0("missing_fraction>", max_missing_fraction))
    ) %>%
    arrange(desc(missing_fraction))

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  keep <- audit %>%
    filter(kept) %>%
    pull(featureID)

  list(
    assay = assay_num[, keep, drop = FALSE],
    feature = feature_tbl %>% filter(featureID %in% keep)
  )
}

presence_filter_and_impute <- function(assay_num, feature_tbl, sample_idx,
                                       min_fraction = 0,
                                       impute_half_min = TRUE,
                                       audit_path = NULL) {
  assay_work <- assay_num
  feature_work <- feature_tbl

  if (min_fraction > 0) {
    present <- apply(assay_work[sample_idx, , drop = FALSE], 2, function(v) mean(!is.na(v) & v != 0))

    audit <- tibble(featureID = names(present), present_rate = present) %>%
      mutate(
        kept = present_rate >= min_fraction,
        exclusion_reason = if_else(kept, "kept", paste0("present_rate<", min_fraction))
      ) %>%
      arrange(present_rate)

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    keep <- audit %>%
      filter(kept) %>%
      pull(featureID)
    assay_work <- assay_work[, keep, drop = FALSE]
    feature_work <- feature_work %>% filter(featureID %in% keep)
  }

  if (isTRUE(impute_half_min)) {
    for (j in seq_len(ncol(assay_work))) {
      v <- assay_work[sample_idx, j]
      nonmiss <- v[!is.na(v) & v > 0]
      if (length(nonmiss) == 0) next

      half_min <- 0.5 * min(nonmiss)
      miss_rows <- sample_idx[is.na(assay_work[sample_idx, j]) | assay_work[sample_idx, j] == 0]
      assay_work[miss_rows, j] <- half_min
    }
  }

  list(assay = assay_work, feature = feature_work)
}

filter_known <- function(assay, feature_tbl, use_only_known = TRUE, audit_path = NULL) {
  if (!isTRUE(use_only_known)) {
    return(list(assay = assay, feature = feature_tbl))
  }

  audit <- feature_tbl %>%
    mutate(
      kept = !is.na(Name_canon),
      exclusion_reason = if_else(kept, "kept", "Name missing/unknown")
    )

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  keep <- audit %>%
    filter(kept) %>%
    pull(featureID)

  list(
    assay = assay[, keep, drop = FALSE],
    feature = feature_tbl %>% filter(featureID %in% keep)
  )
}

calc_rsd <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)

  if (is.na(m) || m == 0) {
    return(NA_real_)
  }

  100 * s / m
}

calc_qc_rsd <- function(assay, qc_idx) {
  apply(assay[qc_idx, , drop = FALSE], 2, calc_rsd)
}

filter_low_variance_deterministic <- function(assay, feature_tbl, method = "none", frac = 0.2,
                                              sample_idx = NULL, audit_path = NULL) {
  if (method == "none" || is.na(frac) || frac <= 0) {
    return(list(assay = assay, feature = feature_tbl))
  }

  if (is.null(sample_idx)) sample_idx <- seq_len(nrow(assay))

  iqr_vals <- apply(assay[sample_idx, , drop = FALSE], 2, IQR, na.rm = TRUE)
  iqr_tbl <- tibble(featureID = names(iqr_vals), iqr = as.numeric(iqr_vals)) %>%
    mutate(iqr = if_else(is.na(iqr), -Inf, iqr))

  n <- nrow(iqr_tbl)
  k_remove <- floor(frac * n)

  if (k_remove < 1 || n < 2) {
    return(list(assay = assay, feature = feature_tbl))
  }

  ord <- order(iqr_tbl$iqr, decreasing = TRUE)
  keep_ids <- iqr_tbl$featureID[ord[seq_len(max(1, n - k_remove))]]

  audit <- iqr_tbl %>%
    mutate(
      rank_desc_iqr = rank(-iqr, ties.method = "first"),
      kept = featureID %in% keep_ids,
      exclusion_reason = if_else(kept, "kept", paste0("bottom_", round(100 * frac), "%_iqr"))
    ) %>%
    arrange(rank_desc_iqr)

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  list(
    assay = assay[, keep_ids, drop = FALSE],
    feature = feature_tbl %>% filter(featureID %in% keep_ids)
  )
}
