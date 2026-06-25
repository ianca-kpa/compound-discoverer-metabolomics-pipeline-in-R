# Metabolomics Pipeline - global objects and helpers

# Compute project/pipeline paths early so we can ensure packages are
# installed before attempting to load libraries used by the app.
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
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

app_helpers_path <- file.path(project_root, "app", "helpers.R")
if (file.exists(app_helpers_path)) source(app_helpers_path)

server_util_helpers_path <- file.path(project_root, "app", "server_util_helpers.R")
if (file.exists(server_util_helpers_path)) source(server_util_helpers_path)

server_settings_helpers_path <- file.path(project_root, "app", "server_settings_helpers.R")
if (file.exists(server_settings_helpers_path)) source(server_settings_helpers_path)

metadata_helpers_path <- file.path(project_root, "app", "metadata_helpers.R")
if (file.exists(metadata_helpers_path)) source(metadata_helpers_path)

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


settings_form_sections <- list(
  list(
    title = "Normalization",
    fields = list(
      list(
        key = "normalization_mode",
        label = "Preprocessing scenario",
        type = "select",
        choices = c("QC-RSC" = "qcrsc", "QC-LOESS" = "qc_loess", "Cyclic LOESS" = "cyclic_loess", "PQN-QC" = "pqn_qc", "PQN Sample" = "pqn_sample", "None" = "none"),
        default = "none"
      ),
      list(key = "make_qc_diagnostics", label = "Run QC/normalization audit", type = "logical_select", default = FALSE)
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
    title = "Statistical analysis",
    note = tags$div(
      class = "small-note settings-section-note",
      tags$strong("Multi-group analysis: "),
      "exploratory and complementary. The five primary pairwise comparisons remain enabled. The global test has no directional FC or volcano plot; select only biologically interpretable follow-up pairs."
    ),
    fields = list(
      list(
        key = "comparison_mode",
        label = "How many groups will be compared?",
        type = "select",
        choices = c(
          "Primary pairwise comparisons" = "pairwise",
          "Primary + exploratory multigroup" = "multigroup",
          "Primary + multigroup (legacy both)" = "both"
        ),
        default = "pairwise"
      ),
      list(
        key = "statistical_test_type",
        label = "Pairwise test",
        type = "select",
        choices = c("student", "welch", "wilcoxon", "limma"),
        default = "student"
      ),
      list(
        key = "test_is_paired",
        label = "Pairwise design",
        type = "select",
        choices = c("Unpaired" = "FALSE", "Paired" = "TRUE"),
        default = "FALSE"
      ),
      list(key = "run_metrics", label = "Significance metric", type = "select", choices = c("FDR", "p_value", "FDR_and_p_value"), default = "FDR_and_p_value")
    )
  ),
  list(
    title = "Multi-group statistics",
    condition = "input.settings_comparison_mode == 'multigroup' || input.settings_comparison_mode == 'both'",
    fields = list(
      list(
        key = "multigroup_groups",
        label = "Groups to compare",
        type = "detected_multiselect",
        default = character(0)
      ),
      list(key = "multigroup_test", label = "Global test (3+ groups)", type = "select", choices = c("kruskal", "anova", "welch_anova"), default = "anova"),
      list(
        key = "multigroup_pairwise_mode",
        label = "Pairwise follow-up",
        type = "select",
        choices = c("Selected biologically interpretable pairs" = "selected", "No additional follow-up" = "none"),
        default = "selected"
      ),
      list(
        key = "multigroup_pairwise_pairs",
        label = "Selected pairs",
        type = "nullable_vector_text",
        default = NULL,
        condition = "input.settings_multigroup_pairwise_mode == 'selected'"
      )
    )
  ),
  list(
    title = "PCA",
    fields = list(
      list(key = "pca_scaling", label = "PCA scaling", type = "select", choices = c("none", "pareto", "autoscale"), default = "pareto"),
      list(key = "pca_label_samples", label = "Label PCA samples", type = "logical_select", default = TRUE)
    )
  ),
  list(
    title = "Heatmap",
    fields = list(
      list(key = "heatmap_top_n", label = "Heatmap top N", type = "integer", default = 50, step = 1, min = 1),
      list(key = "make_heatmap_by_model", label = "Heatmap by model", type = "logical_select", default = TRUE),
      list(key = "make_heatmap_by_model_sex", label = "Heatmap by model and sex", type = "logical_select", default = TRUE),
      list(key = "heatmap_scale_method", label = "Heatmap scale method", type = "select", choices = c("none", "zscore", "pareto"), default = "zscore")
    )
  ),
  list(
    title = "Output controls",
    note = tags$div(
      class = "small-note settings-section-note",
      tags$strong("Output level: "),
      tags$br(),
      "Standard is recommended.",
      tags$br(),
      "Minimal keeps final statistics, principal/combined PCA, primary volcano plots, README, and log.",
      tags$br(),
      "Full / Debug adds every intermediate and technical artifact."
    ),
    fields = list(
      list(
        key = "output_level",
        label = "Output level",
        type = "select",
        choices = c("Minimal" = "minimal", "Standard" = "standard", "Full / Debug" = "full_debug"),
        default = "standard"
      )
    )
  ),
  list(
    title = "Volcano",
    fields = list(
      list(key = "make_volcano_plots", label = "Pairwise volcano plots", type = "logical_select", default = TRUE),
      list(key = "volcano_add_labels", label = "Add volcano labels", type = "logical_select", default = TRUE),
      list(key = "volcano_add_cutoff_lines", label = "Add cutoff lines", type = "logical_select", default = TRUE)
    )
  ),
  list(
    title = "Feature filters",
    fields = list(
      list(key = "active_variant", label = "RSD filter", type = "select", choices = c("none", "QC_RSD", "RSD"), default = "none"),
      list(key = "rsd_thresholds", label = "RSD threshold", type = "numeric", default = 20, step = 1, min = 0),
      list(key = "low_variance_filter_method", label = "Low-variance filter", type = "select", choices = c("iqr", "none"), default = "none"),
      list(key = "low_variance_filter_fraction", label = "IQR fraction", type = "numeric", default = 0.20, step = 0.01, min = 0, max = 1),
      list(key = "use_only_known", label = "Use only known features", type = "logical_select", default = TRUE),
      list(
        key = "duplicate_name_strategy",
        label = "Duplicate handling",
        type = "select",
        choices = c("reference_or_best_qc_rsd", "keep_separate", "collapse_mean", "collapse_sum", "collapse_best_qc_rsd"),
        default = "collapse_best_qc_rsd"
      )
    )
  )
)

settings_builder_layout <- list(
  run_setup = list(
    label = "Run setup",
    sections = c(
      "Statistical analysis",
      "Multi-group statistics",
      "Statistics thresholds",
      "Normalization",
      "Feature filters"
    ),
    widths = c(12, 12, 6, 6, 12)
  ),
  visual_outputs = list(
    label = "Visual outputs",
    sections = c(
      "Heatmap",
      "Volcano",
      "PCA"
    ),
    widths = c(6, 6, 12)
  ),
  exports = list(
    label = "Exports",
    sections = c("Output controls"),
    widths = c(12),
    narrow = FALSE
  ),
  field_columns = c(
    default = 2,
    "Output controls" = 2,
    "PCA" = 2,
    "Normalization" = 2,
    "Feature filters" = 3,
    "Statistics thresholds" = 3,
    "Statistical analysis" = 2,
    "Multi-group statistics" = 2
    # "Heatmap" = 4,
    # "Volcano" = 3
  )
)

settings_glossary_map <- c(
  normalization_mode = "Preprocessing scenario applied before downstream filtering and statistics. Use PQN Sample when QC samples are absent or too few.",
  make_qc_diagnostics = "When TRUE, runs the optional QC/normalization audit plots and summary tables. Leave FALSE for faster routine runs.",
  active_variant = "Chooses the RSD filter family: none disables RSD filtering, QC_RSD uses QC samples, and RSD uses raw SD/mean.",
  rsd_thresholds = "Numeric RSD cutoff used with QC_RSD or RSD. Example: QC_RSD with 30 is saved as rsd_thresholds <- c(30).",
  low_variance_filter_method = "Controls the low-variance filter. Use iqr to remove the lowest-variance features, or none to keep them.",
  low_variance_filter_fraction = "Fraction of features removed by the IQR low-variance filter. Example: 0.20 removes the lowest 20%.",
  use_only_known = "When TRUE, keeps only features with known identities.",
  p_value_cutoff = "P-value threshold used for significance decisions and volcano/heatmap filtering.",
  fdr_cutoff = "Adjusted p-value/FDR threshold used when the selected metric includes FDR.",
  fc_cutoff_log2 = "Minimum absolute log2 fold-change required where fold-change filtering is enabled. Use 0 to disable this cutoff.",
  statistical_test_type = "Test used for two-group comparisons and optional pairwise follow-ups. It is separate from the global multi-group test.",
  test_is_paired = "Controls whether two-group comparisons use paired or unpaired observations.",
  run_metrics = "Controls whether significance outputs use raw p-values, FDR-adjusted values, or both.",
  duplicate_name_strategy = "Controls how duplicated metabolite names are retained or collapsed.",
  output_level = "Controls exported files without changing calculations: Minimal keeps final stats/PCA/main volcano/README/log; Standard adds main heatmaps and summarized QC/normalization outputs; Full / Debug exports all intermediate and technical artifacts.",
  comparison_mode = "The five primary pairwise comparisons always run. Multi-group mode adds an exploratory global test and optional selected follow-up pairs.",
  multigroup_groups = "Editable group names used for multi-group analysis. Leave empty to use all detected biological groups per model.",
  multigroup_test = "Exploratory global test for three or more groups. It detects whether groups differ, but does not estimate a directional fold change.",
  multigroup_pairwise_mode = "Controls additional follow-up beyond the five primary comparisons. Prefer selected biologically interpretable pairs; the global test itself is never used for volcano plots.",
  multigroup_pairwise_pairs = "Selected follow-up pairs. Use entries like pre vs post, pre/post, or pre,post.",
  pca_scaling = "Scaling applied before PCA: none, pareto, or autoscale.",
  pca_label_samples = "When TRUE, sample labels are printed on PCA plots.",
  heatmap_top_n = "Maximum number of ranked features shown in standard heatmaps.",
  make_heatmap_by_model = "When TRUE, generates heatmaps for each model.",
  make_heatmap_by_model_sex = "When TRUE, generates heatmaps split by model and sex.",
  heatmap_scale_method = "Scaling used inside heatmap matrices: none, zscore, or pareto.",
  make_volcano_plots = "When TRUE, generates volcano plots only for pairwise comparisons. MULTIGROUP_GLOBAL is always excluded because it has no directional effect.",
  volcano_add_labels = "When TRUE, labels significant points in volcano plots.",
  volcano_add_cutoff_lines = "When TRUE, draws cutoff guide lines on volcano plots.",
  export_metaboanalyst_ready = "When TRUE, writes the main MetaboAnalyst-compatible export.",
  save_stats_excel_per_model = "When TRUE, writes one Excel stats workbook per model.",
  save_sig_metabolites_txt_per_model = "When TRUE, writes plain-text lists of significant metabolites per model/comparison."
)

initial_settings_text <- read_initial_config()

initial_setting_value <- function(key, default = "") {
  setting_display_value(initial_settings_text, key, default = default)
}

initial_setting_logical <- function(key, default = FALSE) {
  setting_display_logical(initial_settings_text, key, default = default)
}

initial_setting_numeric <- function(key, default = NA_real_) {
  setting_display_numeric(initial_settings_text, key, default = default)
}
