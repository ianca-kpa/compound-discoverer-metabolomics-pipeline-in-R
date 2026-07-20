# =============================================================================
# 09_pca.R
# PCA functions
# =============================================================================

# Summary:
  # This module contains functions for performing PCA on the log2-transformed assay matrix, with different
  # scaling options (none, pareto, autoscale) as specified in the settings. It includes a main driver function 
  # that generates PCA plots for each model and subsets of samples, and saves them to the appropriate output directories.

# scale_uv(): Autoscaling (unit variance scaling) - centers to mean and scales to unit variance.
scale_uv <- function(mat) {
  scale_matrix_columns_zscore(mat)
}

# scale_pareto(): Pareto scaling - centers to mean and scales to square root of standard deviation.
scale_pareto <- function(mat) {
  scale_matrix_columns_pareto(mat)
}

# apply_pca_scaling(): Apply the specified scaling method to the matrix before PCA, 
# and return both the scaled matrix and a label for the method used.
apply_pca_scaling <- function(mat, method = c("pareto", "autoscale", "none")) {
  method <- match.arg(method)

  if (method == "none") {
    out <- replace_missing_matrix_values(mat)
    return(list(mat = out, label = "none"))
  }

  if (method == "pareto") {
    return(list(mat = scale_pareto(mat), label = "pareto"))
  }

  list(mat = scale_uv(mat), label = "autoscale")
}

# -----------------------------------------------------------------------------
# Remove invalid / constant columns before PCA
# -----------------------------------------------------------------------------
prepare_matrix_for_pca <- function(mat) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) {
    return(NULL)
  }

  keep_cols <- apply(mat, 2, function(v) {
    vv <- stats::na.omit(v)
    length(vv) > 0 && length(unique(vv)) > 1
  })

  mat2 <- mat[, keep_cols, drop = FALSE]

  if (ncol(mat2) < 2) {
    return(NULL)
  }

  mat2 <- as.matrix(mat2)
  mat2[is.na(mat2)] <- 0
  mat2
}

pca_legend_title <- function(var_name) {
  switch(
    var_name,
    group = "Group",
    sex = "Sex",
    tools::toTitleCase(gsub("_", " ", as.character(var_name)))
  )
}

pca_color_values <- function(var_name, levels, model_name = NULL) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]

  if (length(levels) == 0) {
    return(c(ALL = "#6B7280"))
  }

  vals <- setNames(rep("#6B7280", length(levels)), levels)

  if (identical(var_name, "group")) {
    groups <- resolve_model_group_values(model_name)
    group_levels <- order_pre_post_levels(c(groups$control, groups$treatment))

    # Blue for control and orange for treatment, keeping the order stable across models.
    ctrl_idx <- which(tolower(levels) == tolower(group_levels[1]))[1]
    trt_idx <- which(tolower(levels) == tolower(group_levels[2]))[1]

    if (!is.na(ctrl_idx)) vals[ctrl_idx] <- "#4EEE94"
    if (!is.na(trt_idx)) vals[trt_idx] <- "#FFA54F"

    remaining <- which(vals == "#6B7280")
    fallback <- c("#14B8A6", "#A855F7", "#EF4444", "#22C55E")
    if (length(remaining) > 0) {
      vals[remaining] <- rep_len(fallback, length(remaining))
    }

    return(vals)
  }

  if (identical(var_name, "sex")) {
    female_idx <- which(tolower(levels) == "f")[1]
    male_idx <- which(tolower(levels) == "m")[1]

    if (!is.na(female_idx)) vals[female_idx] <- "#CD0000"
    if (!is.na(male_idx)) vals[male_idx] <- "#009ACD"

    remaining <- which(vals == "#6B7280")
    fallback <- c("#6366F1", "#EC4899", "#16A34A")
    if (length(remaining) > 0) {
      vals[remaining] <- rep_len(fallback, length(remaining))
    }

    return(vals)
  }

  fallback <- c("#2563EB", "#F97316", "#14B8A6", "#A855F7")
  vals[] <- rep_len(fallback, length(levels))
  vals
}

pca_shape_values <- function(var_name, levels) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & nzchar(levels)]

  if (length(levels) == 0) {
    return(c(ALL = 16))
  }

  vals <- setNames(rep(16, length(levels)), levels)

  if (identical(var_name, "sex")) {
    female_idx <- which(tolower(levels) == "f")[1]
    male_idx <- which(tolower(levels) == "m")[1]
    if (!is.na(female_idx)) vals[female_idx] <- 16
    if (!is.na(male_idx)) vals[male_idx] <- 17
  } else {
    shape_seq <- c(16, 17, 15, 18, 8, 3, 7)
    vals <- setNames(rep_len(shape_seq, length(levels)), levels)
  }

  vals
}

pca_apply_legend_theme <- function(p) {
  p + ggplot2::theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9)
  ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 2)
    )
}

# -----------------------------------------------------------------------------
# Draw and save one PCA plot
# -----------------------------------------------------------------------------
plot_one_pca_subset <- function(mat_log2 = mat_log2,
                                meta = meta,
                                out_png = out_png,
                                title_main = title_main,
                                pca_scaling = pca_scaling,
                                color_var = "group",
                                shape_var = "sex",
                                ellipse_color_var = NULL,
                                draw_ellipse = TRUE,
                                ellipse_positive = TRUE,
                                pca_label_samples = TRUE,
                                log_path = NULL) {
  meta <- meta %>%
    dplyr::filter(sample %in% rownames(mat_log2)) %>%
    dplyr::distinct(sample, .keep_all = TRUE)

  if (nrow(meta) < 3) {
    message("  - PCA skipped: fewer than 3 samples for ", basename(out_png))
    return(FALSE)
  }

  mat_sub <- mat_log2[meta$sample, , drop = FALSE]
  mat_sub <- prepare_matrix_for_pca(mat_sub)

  if (is.null(mat_sub)) {
    message("  - PCA skipped: matrix too small or without variable features for ", basename(out_png))
    return(FALSE)
  }

  sc <- apply_pca_scaling(mat_sub, method = pca_scaling)

  pca <- tryCatch(
    stats::prcomp(sc$mat, center = FALSE, scale. = FALSE),
    error = function(e) NULL
  )

  if (is.null(pca)) {
    message("  - PCA failed for ", basename(out_png))
    return(FALSE)
  }

  if (is.null(pca$x) || ncol(pca$x) < 2) {
    message("  - PCA skipped: PC1/PC2 not available for ", basename(out_png))
    return(FALSE)
  }

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1 <- round(100 * var_exp[1], 1)
  pc2 <- round(100 * var_exp[2], 1)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta, by = "sample")

  if (!(color_var %in% names(scores))) {
    scores[[color_var]] <- "ALL"
  }
  if (!is.null(shape_var) && length(shape_var) > 0 && !(shape_var %in% names(scores))) {
    scores[[shape_var]] <- "ALL"
  }

  model_name <- if ("model" %in% names(meta) && length(unique(meta$model)) >= 1) unique(meta$model)[1] else NULL

  scores[[color_var]] <- as.factor(scores[[color_var]])
  if (!is.null(shape_var) && length(shape_var) > 0) {
    scores[[shape_var]] <- as.factor(scores[[shape_var]])
  }

  if (identical(color_var, "group")) {
    groups <- resolve_model_group_values(model_name)
    lvl <- levels(scores[[color_var]])
    ordered <- order_pre_post_levels(lvl, preferred = c(groups$control, groups$treatment))
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[color_var]] <- factor(scores[[color_var]], levels = ordered)
  } else if (identical(color_var, "sex")) {
    lvl <- levels(scores[[color_var]])
    ordered <- c("F", "M")
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[color_var]] <- factor(scores[[color_var]], levels = ordered)
  }

  color_breaks <- levels(scores[[color_var]])
  color_values <- pca_color_values(color_var, color_breaks, model_name = model_name)

  # determine which variable should provide ellipse colors (allow separate mapping)
  ellipse_var <- if (!is.null(ellipse_color_var) && nzchar(as.character(ellipse_color_var))) ellipse_color_var else color_var
  if (!(ellipse_var %in% names(scores))) {
    scores[[ellipse_var]] <- scores[[color_var]]
  }
  scores[[ellipse_var]] <- as.factor(scores[[ellipse_var]])

  if (identical(ellipse_var, "group")) {
    groups <- resolve_model_group_values(model_name)
    lvl <- levels(scores[[ellipse_var]])
    ordered <- order_pre_post_levels(lvl, preferred = c(groups$control, groups$treatment))
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[ellipse_var]] <- factor(scores[[ellipse_var]], levels = ordered)
  } else if (identical(ellipse_var, "sex")) {
    lvl <- levels(scores[[ellipse_var]])
    ordered <- c("F", "M")
    ordered <- ordered[ordered %in% lvl]
    ordered <- c(ordered, setdiff(lvl, ordered))
    scores[[ellipse_var]] <- factor(scores[[ellipse_var]], levels = ordered)
  }

  ellipse_breaks <- levels(scores[[ellipse_var]])
  ellipse_values <- pca_color_values(ellipse_var, ellipse_breaks, model_name = model_name)

  # If ellipses should only be drawn for 'positive' case and the flag is FALSE,
  # when coloring by sex we prefer black points and rely on shape to distinguish sex.
  if (identical(color_var, "sex") && isFALSE(ellipse_positive)) {
    color_values[] <- "#000000"
  }

  if (!is.null(shape_var) && length(shape_var) > 0) {
    if (identical(shape_var, "group")) {
      groups <- resolve_model_group_values(model_name)
      lvl <- levels(scores[[shape_var]])
      ordered <- order_pre_post_levels(lvl, preferred = c(groups$control, groups$treatment))
      ordered <- ordered[ordered %in% lvl]
      ordered <- c(ordered, setdiff(lvl, ordered))
      scores[[shape_var]] <- factor(scores[[shape_var]], levels = ordered)
    } else if (identical(shape_var, "sex")) {
      lvl <- levels(scores[[shape_var]])
      ordered <- c("F", "M")
      ordered <- ordered[ordered %in% lvl]
      ordered <- c(ordered, setdiff(lvl, ordered))
      scores[[shape_var]] <- factor(scores[[shape_var]], levels = ordered)
    }
  }

  if (is.null(shape_var) || length(shape_var) == 0) {
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]]
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::labs(
        title = title_main,
        x = paste0("PC1 (", pc1, "%)"),
        y = paste0("PC2 (", pc2, "%)"),
        color = pca_legend_title(color_var)
      ) +
      ggplot2::theme_minimal()

    p <- p + ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = TRUE
    )
    if (isTRUE(pca_label_samples)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),
        show.legend = FALSE,
        max.overlaps = 20
      )
    }
    # add filled group ellipses if requested and group sizes allow
    if (isTRUE(draw_ellipse) && isTRUE(ellipse_positive)) {
      grp_counts <- table(scores[[color_var]])
      valid_grps <- names(grp_counts[grp_counts >= 3])
      if (length(valid_grps) > 0) {
        ell_data <- scores[scores[[color_var]] %in% valid_grps, , drop = FALSE]
        p <- p + ggplot2::stat_ellipse(
          data = ell_data,
          mapping = ggplot2::aes(x = PC1, y = PC2, fill = .data[[ellipse_var]], group = .data[[ellipse_var]]),
          inherit.aes = FALSE,
          geom = "polygon",
          alpha = 0.15,
          level = 0.95,
          show.legend = TRUE,
          type = "norm"
        )
        p <- p + ggplot2::scale_fill_manual(values = ellipse_values, breaks = ellipse_breaks, name = paste0(pca_legend_title(ellipse_var), " (ellipse)"))
      }
    }
  } else {
    p <- ggplot2::ggplot(
      scores,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]],
        shape = .data[[shape_var]]
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::labs(
        title = title_main,
        x = paste0("PC1 (", pc1, "%)"),
        y = paste0("PC2 (", pc2, "%)"),
        color = pca_legend_title(color_var),
        shape = pca_legend_title(shape_var)
      ) +
      ggplot2::theme_minimal()

    shape_breaks <- levels(scores[[shape_var]])
    shape_values <- pca_shape_values(shape_var, shape_breaks)

    p <- p + ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = TRUE
    ) +
      ggplot2::scale_shape_manual(
        values = shape_values,
        breaks = shape_breaks,
        drop = TRUE
      )
    if (isTRUE(pca_label_samples)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),
        show.legend = FALSE,
        max.overlaps = 20
      )
    }
    # add filled group ellipses if requested and group sizes allow
    if (isTRUE(draw_ellipse) && isTRUE(ellipse_positive)) {
      grp_counts <- table(scores[[color_var]])
      valid_grps <- names(grp_counts[grp_counts >= 3])
      if (length(valid_grps) > 0) {
        ell_data <- scores[scores[[color_var]] %in% valid_grps, , drop = FALSE]
        p <- p + ggplot2::stat_ellipse(
          data = ell_data,
          mapping = ggplot2::aes(x = PC1, y = PC2, fill = .data[[ellipse_var]], group = .data[[ellipse_var]]),
          inherit.aes = FALSE,
          geom = "polygon",
          alpha = 0.12,
          level = 0.95,
          show.legend = TRUE,
          type = "norm"
        )
        p <- p + ggplot2::scale_fill_manual(values = ellipse_values, breaks = ellipse_breaks, name = paste0(pca_legend_title(ellipse_var), " (ellipse)"))
      }
    }
  }

  p <- pca_apply_legend_theme(p)

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )

  if (!is.null(log_path)) {
    log_written_object(
      log_path = log_path,
      file_path = out_png,
      object = scores,
      note = paste0(
        "PCA plot | title=", title_main,
        " | scaling=", sc$label,
        " | color=", color_var,
        " | shape=", shape_var
      )
    )
  }

  message("  ✓ PCA saved: ", out_png)
  TRUE
}

# -----------------------------------------------------------------------------
# PCA driver
# - one PCA per model for each comparison defined in COMPARISON_CONFIGS
# -----------------------------------------------------------------------------
plot_pca_per_model <- function(mat_log2 = mat_log2,
                               metadata_aligned = metadata_aligned,
                               paths = paths,
                               pca_scaling = pca_scaling,
                               ellipse_positive = TRUE,
                               pca_label_samples = TRUE,
                               log_path = NULL,
                               comparison_names = NULL,
                               include_secondary_pca = TRUE) {
  
  models <- metadata_aligned %>%
    dplyr::filter(type == "Sample") %>%
    dplyr::pull(model) %>%
    unique() %>%
    sort()
  n_done <- 0

  get_pca_dir_for_comparison <- function(prefix, mp) {
    # accept both legacy lowercase-style prefixes and comparison config prefixes
    if (prefix %in% c("tg_vs_wt", "tg-vs-wt", "ALL_TGvsWT", "ALL_TGvsWT")) {
      return(mp$plots$pca_tg_vs_wt)
    }
    if (prefix %in% c("tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m", "F_TGvsWT", "M_TGvsWT")) {
      return(mp$plots$pca_by_sex)
    }
    if (prefix %in% c("tg-f_vs_tg-m", "tg_f_vs_tg_m", "TG_FvsM", "TG_MvsF")) {
      return(mp$plots$pca_tg_f_vs_tg_m)
    }
    if (prefix %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "WT_FvsM", "WT_MvsF")) {
      return(mp$plots$pca_wt_f_vs_wt_m)
    }

    # try lowercase-normalized fallback
    np <- tolower(prefix)
    if (np %in% c("tg_vs_wt", "tg-vs-wt")) return(mp$plots$pca_tg_vs_wt)
    if (np %in% c("tg-f_vs_wt-f", "tg_f_vs_wt_f", "tg-m_vs_wt-m", "tg_m_vs_wt_m")) return(mp$plots$pca_by_sex)
    if (np %in% c("tg-f_vs_tg-m", "tg_f_vs_tg_m", "tg_fvsm", "tg_mvsf")) return(mp$plots$pca_tg_f_vs_tg_m)
    if (np %in% c("wt-f_vs_wt-m", "wt_f_vs_wt_m", "wt_fvsm", "wt_mvsf")) return(mp$plots$pca_wt_f_vs_wt_m)

    stop("Unknown comparison prefix for PCA: ", prefix)
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    model_groups <- resolve_model_group_values(m)
    message("  - PCA model groups: model=", m, " | control=", model_groups$control, " | treatment=", model_groups$treatment)

    # ensure only the PCA directories required by the comparison set exist
    comparison_configs_active <- COMPARISON_CONFIGS
    if (!is.null(comparison_names)) {
      comparison_configs_active <- comparison_configs_active[intersect(names(comparison_configs_active), comparison_names)]
    }

    pca_dirs <- vapply(
      comparison_configs_active,
      function(cfg) get_pca_dir_for_comparison(cfg$prefix, mp),
      character(1)
    )

    for (d in unique(pca_dirs)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    meta_model <- metadata_aligned %>%
      dplyr::filter(type == "Sample", model == m)

    for (comp_name in names(comparison_configs_active)) {
      cfg <- comparison_configs_active[[comp_name]]
      model_groups <- resolve_model_group_values(m)
      meta_sub <- cfg$meta_filter(meta_model, model_name = m)
      out_dir <- get_pca_dir_for_comparison(cfg$prefix, mp)

      comp_label <- switch(
        comp_name,
        tg_vs_wt = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=ALL"),
        "tg-f_vs_wt-f" = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=F"),
        "tg-m_vs_wt-m" = paste0(model_groups$treatment, " vs ", model_groups$control, " | sex=M"),
        "tg-f_vs_tg-m" = paste0("FvsM within ", model_groups$treatment),
        "wt-f_vs_wt-m" = paste0("FvsM within ", model_groups$control),
        cfg$label
      )

      # For ALL_TGvsWT produce two PCA plots: ellipse colored by group, and ellipse colored by sex
      if (identical(comp_name, "ALL_TGvsWT") || identical(cfg$prefix, "ALL_TGvsWT")) {
        out_base_group <- file.path(
          out_dir,
          make_compact_output_filename("PCA", paste0("m", m), comp_name, paste0("s", pca_scaling), "egroup", ext = "png")
        )
        ok1 <- plot_one_pca_subset(
          mat_log2 = mat_log2,
          meta = meta_sub,
          out_png = out_base_group,
          title_main = paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ") | ellipse=group"),
          pca_scaling = pca_scaling,
          color_var = cfg$pca_color_var,
          shape_var = cfg$pca_shape_var,
          ellipse_color_var = "group",
          draw_ellipse = TRUE,
          ellipse_positive = ellipse_positive,
          pca_label_samples = pca_label_samples,
          log_path = log_path
        )
        if (isTRUE(ok1)) n_done <- n_done + 1

        if (isTRUE(include_secondary_pca)) {
          out_base_sex <- file.path(
            out_dir,
            make_compact_output_filename("PCA", paste0("m", m), comp_name, paste0("s", pca_scaling), "esex", ext = "png")
          )
          ok2 <- plot_one_pca_subset(
            mat_log2 = mat_log2,
            meta = meta_sub,
            out_png = out_base_sex,
            title_main = paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ") | ellipse=sex"),
            pca_scaling = pca_scaling,
            color_var = cfg$pca_color_var,
            shape_var = cfg$pca_shape_var,
            ellipse_color_var = "sex",
            draw_ellipse = TRUE,
            ellipse_positive = ellipse_positive,
            pca_label_samples = pca_label_samples,
            log_path = log_path
          )
          if (isTRUE(ok2)) n_done <- n_done + 1
        }
      } else {
        title_base <- paste0("PCA - model=", m, " | ", comp_label, " (", pca_scaling, ")")

        ok <- plot_one_pca_subset(
          mat_log2 = mat_log2,
          meta = meta_sub,
          out_png = file.path(
            out_dir,
            make_compact_output_filename("PCA", paste0("m", m), comp_name, paste0("s", pca_scaling), ext = "png")
          ),
          title_main = title_base,
          pca_scaling = pca_scaling,
          color_var = cfg$pca_color_var,
          shape_var = cfg$pca_shape_var,
          ellipse_color_var = cfg$pca_ellipse_color_var,
          draw_ellipse = TRUE,
          ellipse_positive = ellipse_positive,
          pca_label_samples = pca_label_samples,
          log_path = log_path
        )
        if (isTRUE(ok)) n_done <- n_done + 1
      }
    }
  }

  message("  ✓ PCA plots created: ", n_done)
}

# Exploratory PCA across all configured multi-group levels. This complements
# the five primary comparison-specific PCA outputs and has no directional
# numerator/denominator interpretation.
plot_pca_multigroup_per_model <- function(mat_log2,
                                          metadata_aligned,
                                          paths,
                                          multigroup_groups = character(0),
                                          pca_scaling = "pareto",
                                          ellipse_positive = TRUE,
                                          pca_label_samples = TRUE,
                                          log_path = NULL) {
  configured_groups <- parse_multigroup_groups(multigroup_groups)
  models <- sort(unique(as.character(metadata_aligned$model[metadata_aligned$type == "Sample"])))
  n_done <- 0L

  for (m in models) {
    mp <- get_model_paths(paths, m)
    meta_model <- metadata_aligned %>%
      dplyr::filter(type == "Sample", model == m)

    if (length(configured_groups) > 0) {
      meta_model <- meta_model %>% dplyr::filter(group %in% configured_groups)
    }

    present_groups <- unique(stats::na.omit(as.character(meta_model$group)))
    if (length(present_groups) < 3) {
      message("  - Multi-group PCA skipped for model=", m, ": fewer than 3 groups are available.")
      next
    }

    out_png <- file.path(
      mp$plots$pca_global,
      make_compact_output_filename("PCA", paste0("m", m), "MULTIGROUP", paste0("s", pca_scaling), ext = "png")
    )
    ok <- plot_one_pca_subset(
      mat_log2 = mat_log2,
      meta = meta_model,
      out_png = out_png,
      title_main = paste0("Exploratory multi-group PCA | model=", m, " | groups=", paste(present_groups, collapse = ", ")),
      pca_scaling = pca_scaling,
      color_var = "group",
      shape_var = "sex",
      ellipse_color_var = "group",
      draw_ellipse = TRUE,
      ellipse_positive = ellipse_positive,
      pca_label_samples = pca_label_samples,
      log_path = log_path
    )
    if (isTRUE(ok)) n_done <- n_done + 1L
  }

  message("  ✓ Exploratory multi-group PCA plots created: ", n_done)
  invisible(n_done)
}

plot_global_pca_exports <- function(mat_log2_pre,
                                    mat_log2_post,
                                    metadata_aligned,
                                    paths,
                                    pca_scaling = "pareto",
                                    ellipse_positive = TRUE,
                                    log_path = NULL,
                                    file_tag = NULL,
                                    title_label = NULL,
                                    stages = NULL,
                                    include_per_model = TRUE) {
  out_dir <- file.path(paths$global$root, "plots_global", "pca")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  file_prefix <- if (is.null(file_tag) || !nzchar(file_tag)) {
    "PCA_GLOBAL"
  } else {
    paste0("PCA_GLOBAL_", file_tag)
  }

  title_prefix <- if (is.null(title_label) || !nzchar(title_label)) {
    "Global PCA"
  } else {
    paste0("Global PCA ", title_label)
  }

  stage_specs <- list(
    pre_normalization = list(mat = mat_log2_pre, label = "before normalization"),
    post_normalization = list(mat = mat_log2_post, label = "after normalization")
  )
  if (!is.null(stages)) {
    stage_specs <- stage_specs[intersect(names(stage_specs), stages)]
  }

  model_names <- metadata_aligned %>%
    dplyr::filter(type == "Sample") %>%
    dplyr::pull(model) %>%
    unique() %>%
    sort()

  n_done <- 0

  for (stage_name in names(stage_specs)) {
    stage_spec <- stage_specs[[stage_name]]

    # One all-sample overview: QC and Sample rows together, colored by model and shaped by type.
    ok_overview <- plot_one_pca_subset(
      mat_log2 = stage_spec$mat,
      meta = metadata_aligned,
      out_png = file.path(out_dir, make_compact_output_filename(file_prefix, stage_name, "model_type", ext = "png")),
      title_main = paste0(title_prefix, " ", stage_spec$label, " | color=model | shape=type | ", pca_scaling),
      pca_scaling = pca_scaling,
      color_var = "model",
      shape_var = "type",
      ellipse_color_var = "model",
      draw_ellipse = TRUE,
      ellipse_positive = ellipse_positive,
      pca_label_samples = FALSE,
      log_path = log_path
    )
    if (isTRUE(ok_overview)) n_done <- n_done + 1

    # One group-oriented PCA per model, with QC rows kept in the view and labels suppressed.
    if (isTRUE(include_per_model)) for (m in model_names) {
      meta_model <- metadata_aligned %>%
        dplyr::filter(type %in% c("QC", "Sample"), model == m)

      ok_model <- plot_one_pca_subset(
        mat_log2 = stage_spec$mat,
        meta = meta_model,
        out_png = file.path(out_dir, make_compact_output_filename(file_prefix, stage_name, paste0("m", m), "group_type", ext = "png")),
        title_main = paste0(title_prefix, " ", stage_spec$label, " | model=", m, " | color=group | shape=type | ", pca_scaling),
        pca_scaling = pca_scaling,
        color_var = "group",
        shape_var = "type",
        ellipse_color_var = "group",
        draw_ellipse = TRUE,
        ellipse_positive = ellipse_positive,
        pca_label_samples = FALSE,
        log_path = log_path
      )
      if (isTRUE(ok_model)) n_done <- n_done + 1
    }
  }

  message("  ✓ Global PCA plots created: ", n_done)
}

plot_pca_pre_post_per_model <- function(mat_log2_pre,
                                        mat_log2_post,
                                        metadata_biological_final,
                                        paths,
                                        pca_scaling = "pareto",
                                        pca_label_samples = TRUE,
                                        correction_label = "NORMALIZATION",
                                        log_path = NULL) {
  correction_label <- toupper(correction_label)
  correction_display <- gsub("_", "-", correction_label)

  common_features <- intersect(colnames(mat_log2_pre), colnames(mat_log2_post))
  if (length(common_features) < 2) {
    message("  - Biological pre/post PCA skipped: fewer than 2 shared features.")
    return(FALSE)
  }

  common_samples <- intersect(rownames(mat_log2_pre), rownames(mat_log2_post))
  common_samples <- intersect(common_samples, metadata_biological_final$sample)
  if (length(common_samples) < 2) {
    message("  - Biological pre/post PCA skipped: fewer than 2 shared biological samples.")
    return(FALSE)
  }

  meta <- metadata_biological_final %>%
    dplyr::filter(sample %in% common_samples) %>%
    dplyr::distinct(sample, .keep_all = TRUE)

  model_names <- sort(unique(meta$model))
  model_names <- model_names[!is.na(model_names) & nzchar(model_names)]
  if (length(model_names) == 0) {
    message("  - Biological pre/post PCA skipped: no models found.")
    return(FALSE)
  }

  n_done <- 0
  stage_levels <- c(paste0("pre_", correction_label), paste0("post_", correction_label))
  stage_colors <- setNames(
    c("#B91C1C", "#2563EB"),
    stage_levels
  )

  for (m in model_names) {
    meta_model <- meta %>%
      dplyr::filter(model == m)

    model_samples <- meta_model$sample
    if (length(model_samples) < 2) {
      message("  - Biological pre/post PCA skipped for model=", m, ": fewer than 2 samples.")
      next
    }

    pre <- mat_log2_pre[model_samples, common_features, drop = FALSE]
    post <- mat_log2_post[model_samples, common_features, drop = FALSE]
    rownames(pre) <- paste0(model_samples, "__pre")
    rownames(post) <- paste0(model_samples, "__post")

    mat <- rbind(pre, post)
    mat <- prepare_matrix_for_pca(mat)
    if (is.null(mat)) {
      message("  - Biological pre/post PCA skipped for model=", m, ": matrix too small or without variable features.")
      next
    }

    sc <- apply_pca_scaling(mat, method = pca_scaling)
    pca <- tryCatch(
      stats::prcomp(sc$mat, center = FALSE, scale. = FALSE),
      error = function(e) NULL
    )
    if (is.null(pca) || is.null(pca$x) || ncol(pca$x) < 2) {
      message("  - Biological pre/post PCA failed for model=", m, ".")
      next
    }

    var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
      tibble::rownames_to_column("audit_sample") %>%
      dplyr::mutate(
        stage = dplyr::if_else(
          grepl("__pre$", audit_sample),
          paste0("pre_", correction_label),
          paste0("post_", correction_label)
        ),
        stage = factor(stage, levels = stage_levels),
        sample = sub("__(pre|post)$", "", audit_sample)
      ) %>%
      dplyr::left_join(
        meta_model %>% dplyr::select(sample, group, sex, model),
        by = "sample"
      )

    stage_counts <- table(scores$stage)
    ellipse_stages <- names(stage_counts[stage_counts >= 3])
    ellipse_stages <- stage_levels[stage_levels %in% ellipse_stages]

    p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = stage, shape = stage))

    if (length(ellipse_stages) > 0) {
      ellipse_data <- scores[scores$stage %in% ellipse_stages, , drop = FALSE]
      p <- p +
        ggplot2::stat_ellipse(
          data = ellipse_data,
          mapping = ggplot2::aes(x = PC1, y = PC2, fill = stage, group = stage),
          inherit.aes = FALSE,
          geom = "polygon",
          alpha = 0.12,
          level = 0.95,
          type = "norm",
          show.legend = TRUE
        ) +
        ggplot2::scale_fill_manual(values = stage_colors, breaks = stage_levels, name = "Stage")
    }

    p <- p +
      ggplot2::geom_point(size = 3, alpha = 0.9) +
      ggplot2::scale_color_manual(values = stage_colors, breaks = stage_levels) +
      ggplot2::scale_shape_discrete(breaks = stage_levels) +
      ggplot2::labs(
        title = paste0("Biological PCA before vs after ", correction_display, " | model=", m, " (", pca_scaling, ")"),
        x = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(100 * var_exp[2], 1), "%)"),
        color = "Stage",
        shape = "Stage"
      ) +
      ggplot2::theme_minimal()

    if (isTRUE(pca_label_samples)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), show.legend = FALSE, max.overlaps = 40)
    }

    mp <- get_model_paths(paths, m)
    out_dir <- mp$plots$pca_global
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_png <- file.path(out_dir, make_compact_output_filename("PCA", paste0("m", m), "prepost", correction_label, paste0("s", pca_scaling), ext = "png"))
    out_scores <- file.path(out_dir, paste0("PCA_ACTIVE_model_", m, "_pre_vs_post_", correction_label, "_scores.csv"))
    out_variance <- file.path(out_dir, paste0("PCA_ACTIVE_model_", m, "_pre_vs_post_", correction_label, "_variance.csv"))

    ggplot2::ggsave(out_png, p, width = 7, height = 5, dpi = 300)
    write_csv_safe(scores, out_scores)
    write_csv_safe(
      tibble::tibble(
        PC = paste0("PC", seq_along(var_exp)),
        variance_fraction = as.numeric(var_exp),
        variance_percent = round(100 * as.numeric(var_exp), 3)
      ),
      out_variance
    )

    if (!is.null(log_path)) {
      log_written_object(log_path, out_png, scores, note = paste("Biological PCA before vs after", correction_display, "for model", m))
      log_written_object(log_path, out_scores, scores, note = "Biological pre/post PCA scores")
      log_written_object(log_path, out_variance, var_exp, note = "Biological pre/post PCA variance")
    }

    n_done <- n_done + 1
  }

  message("  ✓ Biological pre/post PCA plots created: ", n_done)
  n_done > 0
}

plot_qc_loess_audit_pca <- function(assay_num_pre,
                                    assay_num_post,
                                    metadata_aligned,
                                    qc_idx,
                                    out_dir,
                                    pca_scaling = "pareto",
                                    log2_offset = 1,
                                    pca_label_samples = TRUE,
                                    injection_order = NULL,
                                    correction_label = "QC-LOESS",
                                    log_path = NULL) {
  correction_label <- toupper(correction_label)
  correction_display <- gsub("_", "-", correction_label)
  if (length(qc_idx) < 3) {
    message("  - QC ", correction_display, " PCA audit skipped: fewer than 3 QC samples.")
    return(FALSE)
  }

  if (is.null(injection_order)) {
    injection_order <- seq_len(nrow(assay_num_pre))
  } else {
    injection_order <- as.numeric(injection_order)
    if (length(injection_order) != nrow(assay_num_pre)) {
      stop("QC ", correction_display, " PCA audit injection_order must have one value per sample row.")
    }
  }

  common_features <- intersect(colnames(assay_num_pre), colnames(assay_num_post))
  if (length(common_features) < 2) {
    message("  - QC ", correction_display, " PCA audit skipped: fewer than 2 shared features.")
    return(FALSE)
  }

  pre <- log2_transform(assay_num_pre[qc_idx, common_features, drop = FALSE], log2_offset)
  post <- log2_transform(assay_num_post[qc_idx, common_features, drop = FALSE], log2_offset)
  rownames(pre) <- paste0(rownames(assay_num_pre)[qc_idx], "__pre")
  rownames(post) <- paste0(rownames(assay_num_post)[qc_idx], "__post")

  mat <- rbind(pre, post)
  mat <- prepare_matrix_for_pca(mat)
  if (is.null(mat)) {
    message("  - QC ", correction_display, " PCA audit skipped: matrix too small or without variable features.")
    return(FALSE)
  }

  sc <- apply_pca_scaling(mat, method = pca_scaling)
  pca <- tryCatch(
    stats::prcomp(sc$mat, center = FALSE, scale. = FALSE),
    error = function(e) NULL
  )
  if (is.null(pca) || is.null(pca$x) || ncol(pca$x) < 2) {
    message("  - QC ", correction_display, " PCA audit failed.")
    return(FALSE)
  }

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  stage_levels <- c(paste0("pre_", correction_label), paste0("post_", correction_label))
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("audit_sample") %>%
    dplyr::mutate(
      stage = if_else(grepl("__pre$", audit_sample), paste0("pre_", correction_label), paste0("post_", correction_label)),
      stage = factor(stage, levels = stage_levels),
      sample = sub("__(pre|post)$", "", audit_sample)
    )

  variance <- tibble(
    PC = paste0("PC", seq_along(var_exp)),
    variance_fraction = as.numeric(var_exp),
    variance_percent = round(100 * as.numeric(var_exp), 3)
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_png <- file.path(out_dir, make_compact_output_filename("PCA_QC", "prepost", correction_label, paste0("s", pca_scaling), ext = "png"))
  out_scores <- file.path(out_dir, paste0("PCA_QC_pre_vs_post_", correction_label, "_scores.csv"))
  out_variance <- file.path(out_dir, paste0("PCA_QC_pre_vs_post_", correction_label, "_variance.csv"))

  stage_colors <- setNames(c("#B91C1C", "#2563EB"), stage_levels)
  stage_counts <- table(scores$stage)
  ellipse_stages <- names(stage_counts[stage_counts >= 3])
  ellipse_stages <- stage_levels[stage_levels %in% ellipse_stages]

  p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = stage, shape = stage))

  if (length(ellipse_stages) > 0) {
    ellipse_data <- scores[scores$stage %in% ellipse_stages, , drop = FALSE]
    p <- p +
      ggplot2::stat_ellipse(
        data = ellipse_data,
        mapping = ggplot2::aes(x = PC1, y = PC2, fill = stage, group = stage),
        inherit.aes = FALSE,
        geom = "polygon",
        alpha = 0.12,
        level = 0.95,
        type = "norm",
        show.legend = TRUE
      ) +
      ggplot2::scale_fill_manual(values = stage_colors, breaks = stage_levels, name = "Stage")
  }

  p <- p +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = stage_colors, breaks = stage_levels) +
    ggplot2::scale_shape_discrete(breaks = stage_levels) +
    ggplot2::labs(
      title = paste0("QC PCA before vs after ", correction_display, " (", pca_scaling, ")"),
      x = paste0("PC1 (", round(100 * var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * var_exp[2], 1), "%)"),
      color = "Stage",
      shape = "Stage"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(pca_label_samples)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), show.legend = FALSE, max.overlaps = 40)
  }

  ggplot2::ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  write_csv_safe(scores, out_scores)
  write_csv_safe(variance, out_variance)

  if (!is.null(log_path)) {
    log_written_object(log_path, out_png, scores, note = paste("QC PCA audit before vs after", correction_display))
    log_written_object(log_path, out_scores, scores, note = "QC PCA audit scores")
    log_written_object(log_path, out_variance, variance, note = "QC PCA audit variance")
  }

  TRUE
}

plot_qc_loess_audit_metrics <- function(assay_num_pre,
                                        assay_num_post,
                                        metadata_aligned,
                                        qc_idx,
                                        out_dir,
                                        audit_dir = NULL,
                                        log2_offset = 1,
                                        injection_order = NULL,
                                        correction_label = "QC-LOESS",
                                        log_path = NULL,
                                        export_plots = TRUE,
                                        export_detailed_table = TRUE) {
  correction_label <- toupper(correction_label)
  correction_display <- gsub("_", "-", correction_label)
  if (length(qc_idx) < 3) {
    message("  - QC ", correction_display, " metric audit skipped: fewer than 3 QC samples.")
    return(FALSE)
  }

  if (is.null(injection_order)) {
    injection_order <- seq_len(nrow(assay_num_pre))
  } else {
    injection_order <- as.numeric(injection_order)
    if (length(injection_order) != nrow(assay_num_pre)) {
      stop("QC ", correction_display, " metric audit injection_order must have one value per sample row.")
    }
  }

  common_features <- intersect(colnames(assay_num_pre), colnames(assay_num_post))
  if (length(common_features) < 1) {
    message("  - QC ", correction_display, " metric audit skipped: no shared features.")
    return(FALSE)
  }

  pre_qc <- assay_num_pre[qc_idx, common_features, drop = FALSE]
  post_qc <- assay_num_post[qc_idx, common_features, drop = FALSE]
  injection_order <- injection_order[qc_idx]

  safe_cor <- function(v) {
    ok <- is.finite(v) & is.finite(injection_order)
    if (sum(ok) < 3 || length(unique(v[ok])) < 2) return(NA_real_)
    stats::cor(injection_order[ok], v[ok], method = "spearman")
  }

  rsd_pre <- apply(pre_qc, 2, calc_rsd)
  rsd_post <- apply(post_qc, 2, calc_rsd)
  log_pre <- log2_transform(pre_qc, log2_offset)
  log_post <- log2_transform(post_qc, log2_offset)
  drift_pre <- apply(log_pre, 2, safe_cor)
  drift_post <- apply(log_post, 2, safe_cor)

  metrics <- tibble(
    featureID = common_features,
    qc_rsd_pre = as.numeric(rsd_pre),
    qc_rsd_post = as.numeric(rsd_post),
    qc_rsd_delta = qc_rsd_post - qc_rsd_pre,
    drift_cor_pre = as.numeric(drift_pre),
    drift_cor_post = as.numeric(drift_post),
    abs_drift_cor_pre = abs(drift_cor_pre),
    abs_drift_cor_post = abs(drift_cor_post),
    abs_drift_cor_delta = abs_drift_cor_post - abs_drift_cor_pre,
    rsd_improved = is.finite(qc_rsd_pre) & is.finite(qc_rsd_post) & qc_rsd_post < qc_rsd_pre,
    drift_improved = is.finite(abs_drift_cor_pre) & is.finite(abs_drift_cor_post) & abs_drift_cor_post < abs_drift_cor_pre
  )

  summary_tbl <- tibble(
    metric = c(
      "QC samples",
      "Features evaluated",
      "Median QC RSD",
      "Median QC RSD",
      "QC RSD improved features",
      "Median absolute drift correlation",
      "Median absolute drift correlation",
      "Drift correlation improved features"
    ),
    value = c(
      length(qc_idx),
      length(common_features),
      round(stats::median(metrics$qc_rsd_pre, na.rm = TRUE), 3),
      round(stats::median(metrics$qc_rsd_post, na.rm = TRUE), 3),
      round(100 * mean(metrics$rsd_improved, na.rm = TRUE), 3),
      round(stats::median(metrics$abs_drift_cor_pre, na.rm = TRUE), 3),
      round(stats::median(metrics$abs_drift_cor_post, na.rm = TRUE), 3),
      round(100 * mean(metrics$drift_improved, na.rm = TRUE), 3)
    ),
    unit = c("n", "n", "%", "%", "%", "rho", "rho", "%"),
    stage = c("QC", "feature", paste0("pre_", correction_label), paste0("post_", correction_label), "post_vs_pre", paste0("pre_", correction_label), paste0("post_", correction_label), "post_vs_pre"),
    note = c(
      paste("QC rows used for", correction_display, "audit."),
      "Shared features evaluated in pre/post matrices.",
      "Lower values indicate tighter QC reproducibility.",
      "Lower values indicate tighter QC reproducibility.",
      paste("Percent of features with lower QC RSD after", correction_display, "."),
      "Absolute Spearman correlation with injection order.",
      "Absolute Spearman correlation with injection order.",
      paste("Percent of features with lower absolute drift correlation after", correction_display, ".")
    )
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_metrics <- file.path(out_dir, paste0("QC_", correction_label, "_audit_feature_metrics.csv"))
  if (isTRUE(export_detailed_table)) {
    write_csv_safe(metrics, out_metrics)
  }

  rsd_long <- metrics %>%
    dplyr::select(featureID, pre = qc_rsd_pre, post = qc_rsd_post) %>%
    tidyr::pivot_longer(-featureID, names_to = "stage", values_to = "qc_rsd") %>%
    dplyr::mutate(stage = factor(stage, levels = order_pre_post_levels(stage)))

  drift_long <- metrics %>%
    dplyr::select(featureID, pre = abs_drift_cor_pre, post = abs_drift_cor_post) %>%
    tidyr::pivot_longer(-featureID, names_to = "stage", values_to = "abs_drift_cor") %>%
    dplyr::mutate(stage = factor(stage, levels = order_pre_post_levels(stage)))

  box_long <- rbind(
    data.frame(stage = "pre", value = as.vector(log_pre)),
    data.frame(stage = "post", value = as.vector(log_post))
  ) %>%
    dplyr::mutate(stage = factor(stage, levels = order_pre_post_levels(stage)))

  out_rsd <- file.path(out_dir, make_compact_output_filename("QC", correction_label, "RSD_prepost", ext = "png"))
  out_drift <- file.path(out_dir, make_compact_output_filename("QC", correction_label, "drift_prepost", ext = "png"))
  out_box <- file.path(out_dir, make_compact_output_filename("QC", correction_label, "boxplot_prepost", ext = "png"))

  if (isTRUE(export_plots)) ggplot2::ggsave(
    out_rsd,
    ggplot2::ggplot(rsd_long, ggplot2::aes(x = stage, y = qc_rsd, fill = stage)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.25) +
      ggplot2::labs(title = paste("QC RSD before vs after", correction_display), x = NULL, y = "QC RSD (%)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none"),
    width = 6,
    height = 4,
    dpi = 300
  )

  if (isTRUE(export_plots)) ggplot2::ggsave(
    out_drift,
    ggplot2::ggplot(drift_long, ggplot2::aes(x = stage, y = abs_drift_cor, fill = stage)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.25) +
      ggplot2::labs(title = paste("QC drift correlation before vs after", correction_display), x = NULL, y = "Absolute Spearman rho") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none"),
    width = 6,
    height = 4,
    dpi = 300
  )

  if (isTRUE(export_plots)) ggplot2::ggsave(
    out_box,
    ggplot2::ggplot(box_long, ggplot2::aes(x = stage, y = value, fill = stage)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.1) +
      ggplot2::labs(title = paste("QC log2 intensity before vs after", correction_display), x = NULL, y = "log2 intensity") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none"),
    width = 6,
    height = 4,
    dpi = 300
  )

  if (!is.null(audit_dir)) {
    qc_rsd_summary <- metrics %>%
      dplyr::transmute(
        featureID,
        qc_rsd_pre = qc_rsd_pre,
        qc_rsd_post = qc_rsd_post,
        qc_rsd_delta = qc_rsd_delta,
        qc_rsd_improved = rsd_improved
      )

    drift_summary <- metrics %>%
      dplyr::transmute(
        featureID,
        spearman_pre = drift_cor_pre,
        spearman_post = drift_cor_post,
        abs_spearman_pre = abs_drift_cor_pre,
        abs_spearman_post = abs_drift_cor_post,
        abs_spearman_delta = abs_drift_cor_delta,
        drift_improved = drift_improved
      )

    write_csv_safe(summary_tbl, file.path(audit_dir, "qc_pca_comparison_summary.csv"))
    write_csv_safe(qc_rsd_summary, file.path(audit_dir, paste0("qc_rsd_before_after_", correction_label, "_summary.csv")))
    write_csv_safe(drift_summary, file.path(audit_dir, paste0("drift_spearman_before_after_", correction_label, "_summary.csv")))
  }

  if (!is.null(log_path)) {
    if (isTRUE(export_detailed_table)) log_written_object(log_path, out_metrics, metrics, note = paste("QC", correction_display, "feature audit metrics"))
    if (isTRUE(export_plots)) {
      log_written_object(log_path, out_rsd, rsd_long, note = "QC RSD audit plot")
      log_written_object(log_path, out_drift, drift_long, note = "QC drift correlation audit plot")
      log_written_object(log_path, out_box, box_long, note = "QC log2 intensity audit boxplot")
    }
  }

  TRUE
}
