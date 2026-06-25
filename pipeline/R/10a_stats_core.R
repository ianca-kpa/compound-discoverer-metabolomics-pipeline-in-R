# =============================================================================
# 10a_stats_core.R
# Statistics core helpers
# =============================================================================

# =============================================================================
# 10_stats_volcano.R
# Statistics + volcano + Excel export
# =============================================================================

# Summary:
  # This module contains functions for performing statistical tests (t-tests) between groups, 
  # generating volcano plots, and exporting results to Excel. The main function run_all_stats_
  # 5sets_per_model() computes statistics for predefined comparisons and creates volcano plots, 
  # while export_stats_excel_by_model() saves the results in an Excel workbook with conditional formatting.

resolve_volcano_styles <- function(volcano_style) {
  # Only the classic style is supported in this simplified workflow.
  # This function intentionally ignores other style values and
  # always returns a single-element character vector with "classic".
  "classic"
}

# normalize_metric_modes(): Normalizes the input for metrics to run, 
# expanding "FDR_and_p_value" into its components and validating the selection.
normalize_metric_modes <- function(metrics) {
  valid <- c("FDR", "p_value")

  if (is.null(metrics) || length(metrics) == 0) {
    return(valid)
  }

  expanded <- unlist(lapply(as.character(metrics), function(m) {
    if (identical(m, "FDR_and_p_value")) {
      return(c("FDR", "p_value"))
    }
    m
  }), use.names = FALSE)

  expanded <- unique(expanded)
  expanded <- expanded[expanded %in% valid]

  if (length(expanded) == 0) {
    stop("No valid metric selected. Use one of: FDR, p_value, FDR_and_p_value")
  }

  expanded
}

FLOOR_P <- 1e-300

resolve_statistical_test_type <- function(test_type) {
  test <- tolower(as.character(test_type)[1])

  if (is.na(test) || !nzchar(test)) {
    test <- "student"
  }

  if (!test %in% c("student", "welch", "wilcoxon", "limma")) {
    warning("Invalid statistical_test_type: ", test, ". Falling back to 'student'.")
    test <- "student"
  }

  test
}

resolve_pvalue_correction_method <- function(correction_method) {
  method <- tolower(as.character(correction_method)[1])

  if (is.na(method) || !nzchar(method)) {
    method <- "FDR"
  }

  valid_methods <- c("raw", "fdr", "bonferroni", "holm", "hochberg", "hommel", "by")
  if (!tolower(method) %in% valid_methods) {
    warning("Invalid pvalue_correction_method: ", method, ". Falling back to 'FDR'.")
    method <- "FDR"
  }

  # Map "raw" to appropriate method for p.adjust
  if (identical(tolower(method), "raw")) {
    return("none")
  }

  # Ensure proper capitalization for p.adjust
  switch(tolower(method),
    "fdr" = "BH",
    "bonferroni" = "bonferroni",
    "holm" = "holm",
    "hochberg" = "hochberg",
    "hommel" = "hommel",
    "by" = "BY",
    "BH"
  )
}

parse_multigroup_groups <- function(groups) {
  if (is.null(groups) || length(groups) == 0) {
    return(character(0))
  }

  values <- unlist(strsplit(as.character(groups), "[,;\\n]+", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  unique(values[!is.na(values) & nzchar(values)])
}

parse_multigroup_pair <- function(pair) {
  pair <- trimws(as.character(pair)[1])
  if (!nzchar(pair)) {
    return(character(0))
  }

  parts <- unlist(strsplit(pair, "\\s+(?i:vs)\\s+|\\s*/\\s*|\\s*,\\s*|\\s*;\\s*", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) < 2) {
    return(character(0))
  }
  parts[1:2]
}

make_pairwise_comparison_name <- function(den, num) {
  paste0("PAIR_", make.names(num), "_vs_", make.names(den))
}

build_multigroup_pairwise_configs <- function(groups,
                                             pairwise_mode = "selected",
                                             selected_pairs = NULL) {
  groups <- parse_multigroup_groups(groups)
  pairwise_mode <- tolower(trimws(as.character(pairwise_mode)[1]))
  if (!pairwise_mode %in% c("none", "all", "selected") || length(groups) < 2) {
    return(list())
  }

  pairs <- list()
  if (identical(pairwise_mode, "all")) {
    cmb <- utils::combn(groups, 2, simplify = FALSE)
    pairs <- lapply(cmb, function(x) c(den = x[1], num = x[2]))
  } else if (identical(pairwise_mode, "selected")) {
    pairs <- lapply(as.character(selected_pairs), parse_multigroup_pair)
    pairs <- Filter(function(x) length(x) == 2, pairs)
    pairs <- lapply(pairs, function(x) c(den = x[1], num = x[2]))
  }

  out <- list()
  for (pair in pairs) {
    den <- unname(pair[["den"]])
    num <- unname(pair[["num"]])
    if (!den %in% groups || !num %in% groups || identical(den, num)) {
      next
    }

    comp_name <- make_pairwise_comparison_name(den, num)
    out[[comp_name]] <- local({
      den_level <- den
      num_level <- num
      list(
        meta_filter = function(m, model_name = NULL) {
          dplyr::filter(m, group %in% c(den_level, num_level))
        },
        stats_compare_var = "group",
        stats_den = den_level,
        stats_num = num_level,
        use_model_groups = FALSE,
        pca_color_var = "group",
        pca_shape_var = "sex",
        pca_ellipse_color_var = "group",
        prefix = "MULTIGROUP_PAIRWISE",
        label = paste0(num_level, " vs ", den_level)
      )
    })
  }

  out
}

compute_limma_unpaired_pvalues <- function(mat_log2, groups, den_level, num_level,
                                           batch = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required when statistical_test_type = 'limma'.")
  }

  mat_log2 <- as.matrix(mat_log2)
  if (nrow(mat_log2) != length(groups)) {
    stop("Limma input must contain one matrix row per group value.")
  }

  keep_samples <- !is.na(groups) & groups %in% c(den_level, num_level)
  mat_work <- mat_log2[keep_samples, , drop = FALSE]
  group <- factor(groups[keep_samples], levels = c(den_level, num_level))

  if (sum(group == den_level) < 2 || sum(group == num_level) < 2) {
    return(stats::setNames(rep(NA_real_, ncol(mat_log2)), colnames(mat_log2)))
  }

  batch_adjusted <- FALSE
  batch_levels <- character(0)
  if (!is.null(batch)) {
    if (length(batch) != nrow(mat_log2)) {
      stop("Limma batch input must contain one value per matrix row.")
    }
    batch_work <- trimws(as.character(batch[keep_samples]))
    batch_missing <- is.na(batch_work) | !nzchar(batch_work)
    batch_levels <- unique(batch_work[!batch_missing])

    if (length(batch_levels) >= 2) {
      if (any(batch_missing)) {
        stop(
          "Batch adjustment requested implicitly because multiple batches are present, ",
          "but batch is missing for ", sum(batch_missing), " sample(s) in this comparison."
        )
      }
      batch_factor <- factor(batch_work)
      design <- stats::model.matrix(~ batch_factor + group)
      batch_adjusted <- TRUE
    } else {
      design <- stats::model.matrix(~ group)
    }
  } else {
    design <- stats::model.matrix(~ group)
  }

  if (qr(design)$rank < ncol(design)) {
    stop(
      "The limma design is rank-deficient: batch is confounded with the comparison ",
      den_level, " vs ", num_level, ". Ensure both biological groups are represented ",
      "across batches; a fully confounded design cannot be corrected statistically."
    )
  }

  # limma expects features in rows and samples in columns. Fitting the complete
  # feature matrix is also essential for empirical-Bayes variance moderation.
  fit <- limma::lmFit(t(mat_work), design)
  fit <- limma::eBayes(fit)

  group_coef <- grep("^group", colnames(design))
  if (length(group_coef) != 1) {
    stop("Could not identify the biological-group coefficient in the limma design.")
  }

  pvals <- as.numeric(fit$p.value[, group_coef])
  pvals <- stats::setNames(pvals, rownames(fit$p.value))
  attr(pvals, "batch_adjusted") <- batch_adjusted
  attr(pvals, "batch_levels") <- batch_levels
  pvals
}

perform_statistical_test <- function(x, y, test_type = "student", paired = FALSE) {
  test_type <- tolower(test_type)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(list(p_value = NA_real_, test_type_used = test_type, error = "Insufficient data"))
  }

  if (identical(test_type, "student")) {
    result <- tryCatch(
      {
        t <- t.test(y, x, var.equal = TRUE, paired = paired)
        list(p_value = as.numeric(t$p.value), test_type_used = "student", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "student", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "welch")) {
    result <- tryCatch(
      {
        t <- t.test(y, x, var.equal = FALSE, paired = paired)
        list(p_value = as.numeric(t$p.value), test_type_used = "welch", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "welch", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "wilcoxon")) {
    result <- tryCatch(
      {
        w <- wilcox.test(y, x, paired = paired, exact = FALSE)
        list(p_value = as.numeric(w$p.value), test_type_used = "wilcoxon", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "wilcoxon", error = as.character(e))
      }
    )
    return(result)
  }

  if (identical(test_type, "limma")) {
    if (isTRUE(paired)) {
      return(list(
        p_value = NA_real_,
        test_type_used = "limma",
        error = "Paired limma requires explicit pair identifiers and is not supported yet."
      ))
    }

    result <- tryCatch(
      {
        values <- matrix(c(x, y), ncol = 1, dimnames = list(NULL, "feature"))
        groups <- c(rep("den", length(x)), rep("num", length(y)))
        p_val <- compute_limma_unpaired_pvalues(values, groups, "den", "num")[[1]]
        list(p_value = as.numeric(p_val), test_type_used = "limma", error = NULL)
      },
      error = function(e) {
        list(p_value = NA_real_, test_type_used = "limma", error = as.character(e))
      }
    )
    return(result)
  }

  # Default fallback
  list(p_value = NA_real_, test_type_used = test_type, error = "Unknown test type")
}

# compute_ttest_stats_general(): Computes fold changes, log2 fold changes, p-values, FDR, 
# for a specified comparison between groups in the metadata. 
# It returns a data frame with these statistics along with feature information.
# Supports multiple statistical test types: Student's t-test, Welch's t-test, Wilcoxon rank-sum test, and limma.
compute_ttest_stats_general <- function(mat_log2, mat_raw, meta_sub, feature_tbl,
                                        compare_var = c("group", "sex"),
                                        num_level, den_level,
                                        statistical_test_type = "student",
                                        test_is_paired = FALSE,
                                        pvalue_correction_method = "FDR") {
  compare_var <- match.arg(compare_var)
  statistical_test_type <- resolve_statistical_test_type(statistical_test_type)
  pvalue_correction_method <- resolve_pvalue_correction_method(pvalue_correction_method)

  meta_sub <- meta_sub %>%
    filter(sample %in% rownames(mat_log2), sample %in% rownames(mat_raw)) %>%
    mutate(.ord = match(sample, rownames(mat_log2))) %>%
    arrange(.ord) %>%
    select(-.ord)

  if (nrow(meta_sub) < 4) {
    return(NULL)
  }

  v <- meta_sub[[compare_var]]
  if (sum(v == den_level, na.rm = TRUE) < 2 || sum(v == num_level, na.rm = TRUE) < 2) {
    return(NULL)
  }

  s <- meta_sub$sample
  sub_log2 <- mat_log2[s, , drop = FALSE]
  sub_raw <- mat_raw[s, , drop = FALSE]
  v2 <- meta_sub[[compare_var]]

  mean_den_pre <- colMeans(sub_raw[v2 == den_level, , drop = FALSE], na.rm = TRUE)
  mean_num_pre <- colMeans(sub_raw[v2 == num_level, , drop = FALSE], na.rm = TRUE)

  fc <- rep(NA_real_, length(mean_den_pre))
  ok <- is.finite(mean_den_pre) & is.finite(mean_num_pre) & mean_den_pre > 0 & mean_num_pre > 0
  fc[ok] <- mean_num_pre[ok] / mean_den_pre[ok]
  log2FC <- log2(fc)

  if (identical(statistical_test_type, "limma")) {
    if (isTRUE(test_is_paired)) {
      stop(
        "statistical_test_type = 'limma' cannot be combined with test_is_paired = TRUE yet. ",
        "Paired limma requires explicit pair identifiers in metadata."
      )
    }
    batch_values <- if ("batch" %in% names(meta_sub)) meta_sub$batch else NULL
    pvals <- compute_limma_unpaired_pvalues(
      sub_log2,
      v2,
      den_level,
      num_level,
      batch = batch_values
    )
    batch_adjusted <- isTRUE(attr(pvals, "batch_adjusted"))
    batch_levels_used <- attr(pvals, "batch_levels")
    pvals <- as.numeric(pvals[colnames(sub_log2)])
  } else {
    batch_values <- if ("batch" %in% names(meta_sub)) trimws(as.character(meta_sub$batch)) else character(0)
    batch_levels_used <- unique(batch_values[!is.na(batch_values) & nzchar(batch_values)])
    batch_adjusted <- FALSE
    if (length(batch_levels_used) >= 2) {
      stop(
        "Multiple batches are present in comparison ", num_level, " vs ", den_level,
        ", but statistical_test_type = '", statistical_test_type,
        "' cannot adjust for covariates. Use statistical_test_type = 'limma' for batch-adjusted inference."
      )
    }
    pvals <- rep(NA_real_, ncol(sub_log2))

    for (j in seq_len(ncol(sub_log2))) {
      x <- sub_log2[v2 == den_level, j]
      y <- sub_log2[v2 == num_level, j]

      if (all(is.na(x)) || all(is.na(y))) next
      if (length(unique(na.omit(c(x, y)))) < 2) next

      test_result <- perform_statistical_test(
        x = x,
        y = y,
        test_type = statistical_test_type,
        paired = test_is_paired
      )

      pvals[j] <- test_result$p_value
    }
  }

  # Apply p-value correction method
  adjusted_pvals <- if (identical(pvalue_correction_method, "none")) {
    pvals
  } else {
    p.adjust(pvals, method = pvalue_correction_method)
  }

  result_df <- tibble(
    featureID = colnames(sub_log2),
    FC_num_over_den = as.numeric(fc),
    log2FC_num_over_den = as.numeric(log2FC),
    p_value = pvals,
    FDR = adjusted_pvals,
    batch_adjusted = batch_adjusted,
    batch_levels = if (length(batch_levels_used) > 0) paste(batch_levels_used, collapse = " | ") else NA_character_
  )

  # If correction method is "none" (raw p-values), add FDR column for compatibility
  if (identical(pvalue_correction_method, "none")) {
    # Compute BH FDR for compatibility even when not used for significance
    result_df <- result_df %>%
      mutate(FDR_BH = p.adjust(p_value, method = "BH"))
  }

  result_df %>%
    left_join(
      feature_tbl %>%
        select(
          featureID,
          any_of(c(
            "display_name", "mz", "RT", "Name", "Name_canon",
            "Metabolika_pathways", "Formula"
          ))
        ),
      by = "featureID"
    )
}

perform_multigroup_test <- function(values, groups, test_type = "kruskal") {
  ok <- is.finite(values) & !is.na(groups) & nzchar(as.character(groups))
  values <- values[ok]
  groups <- droplevels(factor(groups[ok]))

  group_counts <- table(groups)
  groups_keep <- names(group_counts[group_counts >= 2])
  keep <- groups %in% groups_keep
  values <- values[keep]
  groups <- droplevels(groups[keep])

  if (length(levels(groups)) < 2 || length(values) < 4 || length(unique(values)) < 2) {
    return(list(p_value = NA_real_, test_type_used = test_type, error = "Insufficient data"))
  }

  test_type <- tolower(trimws(as.character(test_type)[1]))
  tryCatch(
    {
      if (identical(test_type, "anova")) {
        fit <- stats::aov(values ~ groups)
        p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
        return(list(p_value = as.numeric(p_val), test_type_used = "anova", error = NULL))
      }

      if (identical(test_type, "welch_anova")) {
        wt <- stats::oneway.test(values ~ groups, var.equal = FALSE)
        return(list(p_value = as.numeric(wt$p.value), test_type_used = "welch_anova", error = NULL))
      }

      kw <- stats::kruskal.test(values ~ groups)
      list(p_value = as.numeric(kw$p.value), test_type_used = "kruskal", error = NULL)
    },
    error = function(e) {
      list(p_value = NA_real_, test_type_used = test_type, error = conditionMessage(e))
    }
  )
}

compute_multigroup_stats_general <- function(mat_log2,
                                             mat_raw,
                                             meta_sub,
                                             feature_tbl,
                                             groups,
                                             multigroup_test = "kruskal",
                                             pvalue_correction_method = "FDR") {
  groups <- parse_multigroup_groups(groups)
  pvalue_correction_method <- resolve_pvalue_correction_method(pvalue_correction_method)

  meta_sub <- meta_sub %>%
    dplyr::filter(sample %in% rownames(mat_log2), sample %in% rownames(mat_raw)) %>%
    dplyr::mutate(.ord = match(sample, rownames(mat_log2))) %>%
    dplyr::arrange(.ord) %>%
    dplyr::select(-.ord)

  if (length(groups) == 0) {
    groups <- sort(unique(as.character(meta_sub$group[!is.na(meta_sub$group) & nzchar(as.character(meta_sub$group))])))
  }

  meta_sub <- meta_sub %>% dplyr::filter(group %in% groups)
  present_groups <- groups[groups %in% unique(as.character(meta_sub$group))]
  group_counts <- table(factor(meta_sub$group, levels = present_groups))
  valid_groups <- names(group_counts[group_counts >= 2])

  if (length(valid_groups) < 3 || nrow(meta_sub) < 6) {
    return(NULL)
  }

  meta_sub <- meta_sub %>% dplyr::filter(group %in% valid_groups)
  s <- meta_sub$sample
  sub_log2 <- mat_log2[s, , drop = FALSE]
  group_vec <- factor(meta_sub$group, levels = valid_groups)

  pvals <- rep(NA_real_, ncol(sub_log2))
  test_used <- rep(NA_character_, ncol(sub_log2))
  test_error <- rep(NA_character_, ncol(sub_log2))

  for (j in seq_len(ncol(sub_log2))) {
    test_result <- perform_multigroup_test(
      values = sub_log2[, j],
      groups = group_vec,
      test_type = multigroup_test
    )
    pvals[j] <- test_result$p_value
    test_used[j] <- test_result$test_type_used
    test_error[j] <- if (is.null(test_result$error)) NA_character_ else as.character(test_result$error)[1]
  }

  adjusted_pvals <- if (identical(pvalue_correction_method, "none")) {
    # MULTIGROUP_GLOBAL always exposes a genuine FDR column even when raw
    # p-values are the selected decision metric.
    stats::p.adjust(pvals, method = "BH")
  } else {
    p.adjust(pvals, method = pvalue_correction_method)
  }

  group_means <- lapply(valid_groups, function(g) {
    colMeans(mat_raw[meta_sub$sample[group_vec == g], , drop = FALSE], na.rm = TRUE)
  })
  names(group_means) <- valid_groups

  result_df <- tibble(
    featureID = colnames(sub_log2),
    comparison_type = "multigroup_global",
    groups_compared = paste(valid_groups, collapse = " | "),
    n_groups = length(valid_groups),
    n_samples = length(s),
    FC_num_over_den = NA_real_,
    log2FC_num_over_den = NA_real_,
    p_value = pvals,
    FDR = adjusted_pvals,
    test_type_used = test_used,
    test_error = test_error
  )

  for (g in valid_groups) {
    result_df[[paste0("mean_raw_", make.names(g))]] <- as.numeric(group_means[[g]])
  }

  result_df %>%
    dplyr::left_join(
      feature_tbl %>%
        dplyr::select(
          featureID,
          dplyr::any_of(c(
            "display_name", "mz", "RT", "Name", "Name_canon",
            "Metabolika_pathways", "Formula"
          ))
        ),
      by = "featureID"
    )
}

