# =============================================================================
# 02_comparisons.R
# Comparison definitions
# =============================================================================

comparison_group_control <- get0("comparison_group_control", ifnotfound = "WT", inherits = TRUE)
comparison_group_treatment <- get0("comparison_group_treatment", ifnotfound = "TG", inherits = TRUE)

comparison_group_levels <- c(comparison_group_control, comparison_group_treatment)

COMPARISON_CONFIGS <- list(
  tg_vs_wt = list(
    meta_filter = function(m) dplyr::filter(m, group %in% comparison_group_levels),
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = "sex",
    prefix = "tg_vs_wt",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=ALL")
  ),

  f_vs_m = list(
    meta_filter = function(m) dplyr::filter(m, sex %in% c("F", "M")),
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = "group",
    prefix = "f_vs_m",
    label = "F vs M | group=ALL"
  ),
  
  "tg-f_vs_wt-f" = list(
    meta_filter = function(m) dplyr::filter(m, sex == "F", group %in% comparison_group_levels),
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "tg-f_vs_wt-f",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=F")
  ),
  
  "tg-m_vs_wt-m" = list(
    meta_filter = function(m) dplyr::filter(m, sex == "M", group %in% comparison_group_levels),
    stats_compare_var = "group",
    stats_den = comparison_group_control,
    stats_num = comparison_group_treatment,
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "tg-m_vs_wt-m",
    label = paste0(comparison_group_treatment, " vs ", comparison_group_control, " | sex=M")
  ),
  
  "tg-f_vs_tg-m" = list(
    meta_filter = function(m) dplyr::filter(m, group == comparison_group_treatment, sex %in% c("F", "M")),
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "tg-f_vs_tg-m",
    label = paste0("F vs M within ", comparison_group_treatment)
  ),
  
  "wt-f_vs_wt-m" = list(
    meta_filter = function(m) dplyr::filter(m, group == comparison_group_control, sex %in% c("F", "M")),
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "wt-f_vs_wt-m",
    label = paste0("F vs M within ", comparison_group_control)
  )
)

COMPARISON_NAMES <- names(COMPARISON_CONFIGS)
