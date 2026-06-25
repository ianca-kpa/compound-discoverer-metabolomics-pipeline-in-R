# =============================================================================
# 06c_filter_helpers.R
# Filtering and RSD helpers
# =============================================================================

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

calc_rsd <- function(x, rsd_filter_type = "QC_RSD") {
  rsd_filter_type_local <- toupper(trimws(as.character(rsd_filter_type)[1]))
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)

  if (is.na(m) || m == 0) {
    return(NA_real_)
  }

  if (identical(rsd_filter_type_local, "RSD")) {
    return(s / m)
  }

  100 * s / m
}

calc_qc_rsd <- function(assay, qc_idx, rsd_filter_type = "QC_RSD") {
  apply(assay[qc_idx, , drop = FALSE], 2, calc_rsd, rsd_filter_type = rsd_filter_type)
}

filter_low_variance_deterministic <- function(assay, feature_tbl, method = "none", frac = 0.2,
                                              rounding = "floor", sample_idx = NULL,
                                              audit_path = NULL) {
  rounding <- tolower(trimws(as.character(rounding)[1]))
  if (!rounding %in% c("floor", "ceiling", "round")) {
    stop("Unsupported low-variance filter rounding: ", rounding,
         ". Use 'floor', 'ceiling', or 'round'.")
  }

  if (method == "none" || is.na(frac) || frac <= 0) {
    return(list(
      assay = assay,
      feature = feature_tbl,
      n_before = ncol(assay),
      n_removed = 0,
      n_after = ncol(assay),
      fraction = frac,
      rounding = rounding
    ))
  }

  if (is.null(sample_idx)) sample_idx <- seq_len(nrow(assay))

  iqr_vals <- apply(assay[sample_idx, , drop = FALSE], 2, IQR, na.rm = TRUE)
  iqr_tbl <- tibble(featureID = names(iqr_vals), iqr = as.numeric(iqr_vals)) %>%
    mutate(iqr = if_else(is.na(iqr), -Inf, iqr))

  n <- nrow(iqr_tbl)
  k_remove <- switch(
    rounding,
    floor = floor(frac * n),
    ceiling = ceiling(frac * n),
    round = round(frac * n)
  )
  k_remove <- min(k_remove, max(0, n - 1))

  if (k_remove < 1 || n < 2) {
    audit <- iqr_tbl %>%
      mutate(
        rank_desc_iqr = rank(-iqr, ties.method = "first"),
        kept = TRUE,
        exclusion_reason = "kept",
        low_variance_filter_fraction = frac,
        low_variance_filter_rounding = rounding,
        n_features_before_iqr = n,
        n_features_removed_iqr = 0,
        n_features_after_iqr = n
      ) %>%
      arrange(rank_desc_iqr)

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    return(list(
      assay = assay,
      feature = feature_tbl,
      n_before = n,
      n_removed = 0,
      n_after = n,
      fraction = frac,
      rounding = rounding
    ))
  }

  ord <- order(iqr_tbl$iqr, decreasing = TRUE)
  keep_ids <- iqr_tbl$featureID[ord[seq_len(max(1, n - k_remove))]]

  audit <- iqr_tbl %>%
    mutate(
      rank_desc_iqr = rank(-iqr, ties.method = "first"),
      kept = featureID %in% keep_ids,
      exclusion_reason = if_else(kept, "kept", paste0("bottom_", round(100 * frac), "%_iqr")),
      low_variance_filter_fraction = frac,
      low_variance_filter_rounding = rounding,
      n_features_before_iqr = n,
      n_features_removed_iqr = k_remove,
      n_features_after_iqr = n - k_remove
    ) %>%
    arrange(rank_desc_iqr)

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  list(
    assay = assay[, keep_ids, drop = FALSE],
    feature = feature_tbl %>% filter(featureID %in% keep_ids),
    n_before = n,
    n_removed = k_remove,
    n_after = n - k_remove,
    fraction = frac,
    rounding = rounding
  )
}
