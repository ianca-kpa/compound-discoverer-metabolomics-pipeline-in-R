# =============================================================================
# 11_heatmaps.R
# Heatmap functions
# Keeps the same palette + annotation colors
# =============================================================================

# -----------------------------------------------------------------------------
# Color palette
# -----------------------------------------------------------------------------
make_bwr_palette <- function(n = 101) {
  grDevices::colorRampPalette(c("#1134e6", "#F7F7F7", "#ea112a"))(n)
}

make_heatmap_breaks <- function(n = 101, symmetric = TRUE, limit = 3, mat = NULL) {
  if (isTRUE(symmetric)) {
    return(seq(-limit, limit, length.out = n))
  }

  if (is.null(mat)) {
    stop("mat must be provided when symmetric = FALSE")
  }

  rng <- range(mat, na.rm = TRUE)

  if (!all(is.finite(rng))) {
    rng <- c(-1, 1)
  }

  if (rng[1] == rng[2]) {
    rng <- c(rng[1] - 1, rng[2] + 1)
  }

  seq(rng[1], rng[2], length.out = n)
}

# -----------------------------------------------------------------------------
# Scaling helpers
# -----------------------------------------------------------------------------
scale_cols_zscore <- function(mat) {
  scale_matrix_columns_zscore(mat)
}

scale_cols_pareto <- function(mat) {
  scale_matrix_columns_pareto(mat)
}

apply_heatmap_scaling <- function(mat, method = c("none", "zscore", "pareto")) {
  method <- match.arg(method)

  if (method == "none") {
    return(replace_missing_matrix_values(mat))
  }

  if (method == "zscore") {
    return(scale_cols_zscore(mat))
  }

  scale_cols_pareto(mat)
}

# -----------------------------------------------------------------------------
# Annotation colors (same logic as before)
# -----------------------------------------------------------------------------
get_ma_annotation_colors <- function(model_name = NULL, group_values = NULL) {
  groups <- resolve_model_group_values(model_name)
  group_levels <- unique(trimws(as.character(group_values)))
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
  if (length(group_levels) == 0) {
    group_levels <- c(groups$control, groups$treatment)
  }
  group_levels <- order_pre_post_levels(
    group_levels,
    preferred = c(groups$control, groups$treatment)
  )

  group_colors <- setNames(
    rep_len(c("#4EEE94", "#FFA54F", "#14B8A6", "#A855F7", "#EF4444", "#22C55E"), length(group_levels)),
    group_levels
  )

  list(
    class = group_colors,
    sex   = c(F = "#CD0000", M = "#009ACD")
  )
}

# -----------------------------------------------------------------------------
# Small helper: safe feature labels
# -----------------------------------------------------------------------------
safe_feature_labels <- function(feature_ids, feature_tbl) {
  feat_map <- feature_tbl %>%
    dplyr::select(featureID, display_name)

  lbl <- feat_map$display_name[match(feature_ids, feat_map$featureID)]
  lbl[is.na(lbl) | trimws(lbl) == ""] <- feature_ids[is.na(lbl) | trimws(lbl) == ""]
  make.unique(lbl, sep = "__dup")
}

# -----------------------------------------------------------------------------
# Ranking helper for top heatmaps
# -----------------------------------------------------------------------------
compute_ttest_pvals_only <- function(mat_log2, meta_sub) {
  model_name <- if ("model" %in% names(meta_sub) && length(unique(meta_sub$model)) >= 1) unique(meta_sub$model)[1] else NULL
  groups <- resolve_model_group_values(model_name)

  meta_sub <- meta_sub %>%
    dplyr::filter(group %in% c(groups$control, groups$treatment), sample %in% rownames(mat_log2)) %>%
    dplyr::mutate(.ord = match(sample, rownames(mat_log2))) %>%
    dplyr::arrange(.ord) %>%
    dplyr::select(-.ord)

  if (nrow(meta_sub) < 4) {
    return(NULL)
  }

  if (sum(meta_sub$group == groups$control) < 2 || sum(meta_sub$group == groups$treatment) < 2) {
    return(NULL)
  }

  s <- meta_sub$sample
  g <- meta_sub$group
  sub <- mat_log2[s, , drop = FALSE]

  pvals <- rep(NA_real_, ncol(sub))

  for (j in seq_len(ncol(sub))) {
    x <- sub[g == groups$control, j]
    y <- sub[g == groups$treatment, j]

    if (all(is.na(x)) || all(is.na(y))) {
      next
    }

    if (length(unique(stats::na.omit(c(x, y)))) < 2) {
      next
    }

    pvals[j] <- tryCatch(
      stats::t.test(y, x)$p.value,
      error = function(e) NA_real_
    )
  }

  tibble::tibble(
    featureID = colnames(sub),
    p_value = pvals,
    FDR = stats::p.adjust(pvals, method = "BH")
  )
}

get_top_heatmap_rank_table <- function(stats_5sets_by_model,
                                       model_name,
                                       sex_level = NULL,
                                       comparison_mode = "pairwise") {
  if (is.null(stats_5sets_by_model) || is.null(model_name)) {
    return(NULL)
  }

  model_stats <- stats_5sets_by_model[[as.character(model_name)[1]]]
  if (is.null(model_stats)) {
    return(NULL)
  }

  use_multigroup <- tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")
  comp_name <- if (use_multigroup) "MULTIGROUP_GLOBAL" else "ALL_TGvsWT"
  if (!use_multigroup && !is.null(sex_level)) {
    sex_level <- toupper(trimws(as.character(sex_level)[1]))
    if (identical(sex_level, "F")) {
      comp_name <- "F_TGvsWT"
    } else if (identical(sex_level, "M")) {
      comp_name <- "M_TGvsWT"
    } else {
      return(NULL)
    }
  }

  st <- model_stats[[comp_name]]
  if (is.null(st) || nrow(st) == 0 || !"featureID" %in% names(st)) {
    return(NULL)
  }

  st %>%
    dplyr::select(dplyr::any_of(c("featureID", "p_value", "FDR")))
}

# -----------------------------------------------------------------------------
# Path resolver
# -----------------------------------------------------------------------------
resolve_heatmap_dir <- function(mp,
                                sig = FALSE,
                                by_sex = FALSE,
                                by_group = FALSE,
                                comp_type = NULL) {
  # comp_type: NULL (none), "tg_vs_wt", "tg_fvsm", "wt_fvsm"; legacy aliases are accepted.
  comp_key <- tolower(trimws(ifelse(is.null(comp_type), "", comp_type)))

  if (!sig && !by_sex && !by_group && is.null(comp_type)) return(mp$plots$heatmap_global)
  if (!sig &&  by_sex && !by_group && is.null(comp_type)) return(mp$plots$heatmap_by_sex)
  if (!sig && !by_sex &&  by_group && is.null(comp_type)) return(mp$plots$heatmap_by_group)
  if (!sig && !by_sex && !by_group && !is.null(comp_type)) {
    if (comp_key == "tg_vs_wt") return(mp$plots$heatmap_tg_vs_wt)
    if (comp_key %in% c("tg_fvsm", "tg_mvsf", "tg_f_vs_tg_m", "tg-f_vs_tg-m")) return(mp$plots$heatmap_tg_f_vs_tg_m)
    if (comp_key %in% c("wt_fvsm", "wt_mvsf", "wt_f_vs_wt_m", "wt-f_vs_wt-m")) return(mp$plots$heatmap_wt_f_vs_wt_m)
  }
  if ( sig && !by_sex && !by_group && is.null(comp_type)) return(mp$plots$heatmap_significant_global)
  if ( sig &&  by_sex && !by_group && is.null(comp_type)) return(mp$plots$heatmap_significant_by_sex)
  if ( sig && !by_sex &&  by_group && is.null(comp_type)) return(mp$plots$heatmap_significant_by_group)
  if ( sig && !by_sex && !by_group && !is.null(comp_type)) {
    if (comp_key == "tg_vs_wt") return(mp$plots$heatmap_significant_tg_vs_wt)
    if (comp_key %in% c("tg_fvsm", "tg_mvsf", "tg_f_vs_tg_m", "tg-f_vs_tg-m")) return(mp$plots$heatmap_significant_tg_f_vs_tg_m)
    if (comp_key %in% c("wt_fvsm", "wt_mvsf", "wt_f_vs_wt_m", "wt-f_vs_wt-m")) return(mp$plots$heatmap_significant_wt_f_vs_wt_m)
  }

  mp$plots$heatmap_global
}

# -----------------------------------------------------------------------------
# Internal drawer
# -----------------------------------------------------------------------------
draw_pheatmap_png <- function(mat_plot,
                              ann_col,
                              ann_colors,
                              out_png,
                              title_main,
                              cluster_samples = FALSE,
                              fontsize_row = 11,
                              fontsize_col = 8,
                              angle_col = 45) {
  hm_colors <- make_bwr_palette(heatmap_palette_n)
  hm_breaks <- make_heatmap_breaks(
    n = heatmap_palette_n,
    symmetric = heatmap_breaks_symmetric,
    limit = heatmap_breaks_limit,
    mat = mat_plot
  )

  grDevices::png(out_png, width = 2400, height = 1200, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  pheatmap::pheatmap(
    mat_plot,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = isTRUE(cluster_samples),
    clustering_distance_rows = heatmap_cluster_distance,
    clustering_method = heatmap_cluster_method,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = angle_col,
    color = hm_colors,
    breaks = hm_breaks,
    border_color = NA,
    show_colnames = TRUE,
    main = title_main
  )

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Top heatmaps by model
# paths = full pipeline paths object
# Requires get_model_paths(paths, model_name)
# -----------------------------------------------------------------------------
plot_heatmap_top_ttest_per_model <- function(mat_log2,
                                             metadata_aligned,
                                             feature_tbl,
                                             paths,
                                             stats_5sets_by_model = NULL,
                                             top_n = 50,
                                             rank_by = c("p_value", "FDR"),
                                             split_by_sex = FALSE,
                                             order_samples_by_group = TRUE,
                                             cluster_samples = FALSE,
                                             scale_method = c("none", "zscore", "pareto"),
                                             comparison_mode = "pairwise",
                                             multigroup_groups = character(0)) {
  rank_by <- match.arg(rank_by)
  scale_method <- match.arg(scale_method)
  use_multigroup <- tolower(as.character(comparison_mode)[1]) %in% c("multigroup", "both")
  configured_multigroup_groups <- parse_multigroup_groups(multigroup_groups)

  models <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))
  count <- 0

  sex_levels <- if (split_by_sex) {
    sort(unique(metadata_aligned$sex[metadata_aligned$type == "Sample"]))
  } else {
    "ALL"
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    groups <- resolve_model_group_values(m)
    out_dir <- resolve_heatmap_dir(
      mp = mp,
      sig = FALSE,
      by_sex = isTRUE(split_by_sex)
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    for (sx in sex_levels) {
      meta <- metadata_aligned %>%
        dplyr::filter(type == "Sample", model == m, sample %in% rownames(mat_log2))

      if (use_multigroup) {
        if (length(configured_multigroup_groups) > 0) {
          meta <- meta %>% dplyr::filter(group %in% configured_multigroup_groups)
        }
      } else {
        meta <- meta %>% dplyr::filter(group %in% c(groups$control, groups$treatment))
      }

      if (split_by_sex) {
        meta <- meta %>% dplyr::filter(sex == sx)
      }

      if (nrow(meta) < 4) {
        message("  - Skipping TOP heatmap for model=", m,
                if (split_by_sex) paste0(", sex=", sx) else "",
                " (fewer than 4 samples).")
        next
      }

      if (!use_multigroup && (sum(meta$group == groups$control) < 2 || sum(meta$group == groups$treatment) < 2)) {
        message("  - Skipping TOP heatmap for model=", m,
                if (split_by_sex) paste0(", sex=", sx) else "",
                " (need at least 2 ", groups$control, " and 2 ", groups$treatment, ").")
        next
      }

      if (use_multigroup && length(unique(stats::na.omit(meta$group))) < 2) {
        message("  - Skipping multigroup TOP heatmap for model=", m, " (fewer than 2 groups).")
        next
      }

      if (isTRUE(order_samples_by_group)) {
        group_levels <- order_pre_post_levels(meta$group, preferred = c(groups$control, groups$treatment))
        meta <- meta %>%
          dplyr::mutate(
            group = factor(group, levels = group_levels),
            sex   = factor(sex, levels = c("F", "M"))
          ) %>%
          dplyr::arrange(group, sex, sample)
      } else {
        meta <- meta %>%
          dplyr::arrange(sample)
      }

      st_rank <- get_top_heatmap_rank_table(
        stats_5sets_by_model = stats_5sets_by_model,
        model_name = m,
        sex_level = if (isTRUE(split_by_sex)) sx else NULL,
        comparison_mode = comparison_mode
      )

      if (is.null(st_rank) && !use_multigroup) {
        st_rank <- compute_ttest_pvals_only(mat_log2, meta)
      }

      if (is.null(st_rank)) {
        message("  - Skipping TOP heatmap for model=", m,
                if (split_by_sex) paste0(", sex=", sx) else "",
                " (ranking table is NULL).")
        next
      }

      st_rank <- st_rank %>%
        dplyr::filter(!is.na(.data[[rank_by]])) %>%
        dplyr::arrange(.data[[rank_by]])

      if (nrow(st_rank) < 2) {
        message("  - Skipping TOP heatmap for model=", m,
                if (split_by_sex) paste0(", sex=", sx) else "",
                " (fewer than 2 ranked features).")
        next
      }

      feat_ids <- unique(st_rank$featureID)
      feat_ids <- feat_ids[feat_ids %in% colnames(mat_log2)]

      n_use <- min(top_n, length(feat_ids))
      if (n_use < 2) {
        message("  - Skipping TOP heatmap for model=", m,
                if (split_by_sex) paste0(", sex=", sx) else "",
                " (fewer than 2 usable features).")
        next
      }

      feat_ids <- feat_ids[seq_len(n_use)]
      mat_use <- mat_log2[meta$sample, feat_ids, drop = FALSE]
      colnames(mat_use) <- safe_feature_labels(colnames(mat_use), feature_tbl)

      ann_col <- meta %>%
        dplyr::select(class = group, sex) %>%
        as.data.frame()
      rownames(ann_col) <- meta$sample
      ann_colors <- get_ma_annotation_colors(m, group_values = meta$group)

      mat_scaled <- apply_heatmap_scaling(mat_use, method = scale_method)
      mat_plot <- t(mat_scaled)

      sex_tag <- if (split_by_sex) paste0("sex_", sx) else "sex_ALL"
      analysis_tag <- if (use_multigroup) "multigroup" else "pairwise"

      out_png <- file.path(
        out_dir,
        make_compact_output_filename(
          "HM_TOP", paste0("m", m), analysis_tag, sex_tag,
          paste0("top", n_use), paste0("rank", rank_by), paste0("s", scale_method),
          ext = "png"
        )
      )

      title_main <- paste0(
        "Heatmap - model=", m,
        " | ", sex_tag,
        " | analysis=", analysis_tag,
        " | Top ", n_use,
        " ranked by ", rank_by,
        " | scale=", scale_method
      )

      draw_pheatmap_png(
        mat_plot = mat_plot,
        ann_col = ann_col,
        ann_colors = ann_colors,
        out_png = out_png,
        title_main = title_main,
        cluster_samples = cluster_samples,
        fontsize_row = 11,
        fontsize_col = 8,
        angle_col = 45
      )

      message("  ✓ TOP heatmap saved: ", out_png)
      count <- count + 1
    }
  }

  message("  ✓ TOP heatmaps created: ", count)
}

# -----------------------------------------------------------------------------
# Significant heatmaps from stats tables
# -----------------------------------------------------------------------------
plot_sig_heatmap_from_stats <- function(mat_log2,
                                        meta,
                                        feature_tbl,
                                        stats_df,
                                        sig_metric = c("FDR", "p_value"),
                                        alpha_sig = 0.05,
                                        fc_cutoff_log2 = 1,
                                        require_fc_cutoff = TRUE,
                                        sig_max = 50,
                                        scale_method = c("none", "zscore", "pareto"),
                                        order_samples_by_group = TRUE,
                                        cluster_samples = FALSE,
                                        out_png,
                                        title_main) {
  sig_metric <- match.arg(sig_metric)
  scale_method <- match.arg(scale_method)
  model_name <- if ("model" %in% names(meta) && length(unique(meta$model)) >= 1) unique(meta$model)[1] else NULL
  groups <- resolve_model_group_values(model_name)

  if (is.null(stats_df) || nrow(stats_df) < 2) {
    message("  - Significant heatmap skipped: stats table is NULL or too small.")
    return(FALSE)
  }

  sig_tbl <- stats_df %>%
    dplyr::filter(!is.na(.data[[sig_metric]]), .data[[sig_metric]] < alpha_sig)

  if (isTRUE(require_fc_cutoff)) {
    sig_tbl <- sig_tbl %>%
      dplyr::filter(
        !is.na(log2FC_num_over_den),
        abs(log2FC_num_over_den) >= fc_cutoff_log2
      )
  }

  sig_feats <- sig_tbl %>%
    dplyr::arrange(.data[[sig_metric]]) %>%
    dplyr::pull(featureID) %>%
    unique()

  if (length(sig_feats) == 0) {
    message("  - Significant heatmap skipped: no significant features found.")
    return(FALSE)
  }

  if (length(sig_feats) == 1) {
    message("  - Significant heatmap skipped: only 1 significant feature found.")
    return(FALSE)
  }

  if (length(sig_feats) > sig_max) {
    sig_feats <- head(sig_feats, sig_max)
  }

  meta <- meta %>% dplyr::filter(sample %in% rownames(mat_log2))
  if (nrow(meta) < 4) {
    message("  - Significant heatmap skipped: fewer than 4 samples.")
    return(FALSE)
  }

  if (isTRUE(order_samples_by_group) && all(c("group", "sex") %in% names(meta))) {
    group_levels <- order_pre_post_levels(meta$group, preferred = c(groups$control, groups$treatment))
    meta <- meta %>%
      dplyr::mutate(
        group = factor(group, levels = group_levels),
        sex   = factor(sex, levels = c("F", "M"))
      ) %>%
      dplyr::arrange(group, sex, sample)
  } else {
    meta <- meta %>% dplyr::arrange(sample)
  }

  feat_keep <- intersect(colnames(mat_log2), sig_feats)
  mat_m <- mat_log2[meta$sample, feat_keep, drop = FALSE]

  if (ncol(mat_m) < 2) {
    message("  - Significant heatmap skipped: fewer than 2 plottable features.")
    return(FALSE)
  }

  colnames(mat_m) <- safe_feature_labels(colnames(mat_m), feature_tbl)

  ann_colors <- get_ma_annotation_colors(model_name, group_values = meta$group)

  ann_col <- meta %>%
    dplyr::select(dplyr::any_of(c("group", "sex"))) %>%
    as.data.frame()

  if (ncol(ann_col) > 0) {
    if ("group" %in% names(ann_col)) {
      names(ann_col)[names(ann_col) == "group"] <- "class"
    }
    rownames(ann_col) <- meta$sample
  } else {
    ann_col <- NULL
  }

  mat_scaled <- apply_heatmap_scaling(mat_m, method = scale_method)
  mat_plot <- t(mat_scaled)

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  draw_pheatmap_png(
    mat_plot = mat_plot,
    ann_col = ann_col,
    ann_colors = ann_colors,
    out_png = out_png,
    title_main = title_main,
    cluster_samples = cluster_samples,
    fontsize_row = 11,
    fontsize_col = 8,
    angle_col = 45
  )

  message("  ✓ Significant heatmap saved: ", out_png)
  TRUE
}
