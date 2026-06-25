#!/usr/bin/env Rscript

# Manual runner for QC-LOESS/QC-RSC weight comparison plots.
# The main pipeline generates the same plots when make_qc_diagnostics is TRUE.

args <- commandArgs(trailingOnly = TRUE)

settings_path <- file.path("pipeline", "config", "settings.R")
if (file.exists(settings_path)) {
  source(settings_path, local = .GlobalEnv)
}

default_run_dir <- get0("output_dir", ifnotfound = "output", inherits = TRUE)
exports_dir <- ifelse(
  length(args) >= 1 && nzchar(args[1]),
  args[1],
  file.path(default_run_dir, "global", "exports_global")
)
plots_dir <- ifelse(
  length(args) >= 2 && nzchar(args[2]),
  args[2],
  file.path(dirname(exports_dir), "plots_global", "normalization")
)

if (length(args) == 0) {
  message("No arguments supplied; using exports directory from output_dir/default: ", exports_dir)
}

source("pipeline/R/00_packages.R")
source("pipeline/R/06_normalization_filters.R")

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Read matrices
assay_raw_path <- file.path(exports_dir, "MA_ACTIVE_raw_GLOBAL_WITH_QC.csv")
metadata_path <- file.path(exports_dir, "04_sampleData_aligned.csv")
inj_order_path <- file.path(exports_dir, "00_drift_injection_order_used.csv")

if (!file.exists(assay_raw_path) || !file.exists(metadata_path)) {
  stop(
    "Required exported matrices not found in: ", exports_dir,
    "\nExpected files: ",
    basename(assay_raw_path), ", ", basename(metadata_path)
  )
}

assay_raw_df <- read.csv(assay_raw_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# ensure numeric matrix
assay_raw <- as.matrix(as.data.frame(lapply(assay_raw_df, function(col) as.numeric(as.character(col))), check.names = FALSE))
rownames(assay_raw) <- rownames(assay_raw_df)
metadata_aligned <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

# identify sample and qc indices using metadata
sample_idx <- which(metadata_aligned$type == "Sample" & metadata_aligned$sample %in% rownames(assay_raw))
qc_idx <- which(metadata_aligned$type == "QC" & metadata_aligned$sample %in% rownames(assay_raw))

inj_order <- NULL
if (file.exists(inj_order_path)) {
  inj_df <- read.csv(inj_order_path, stringsAsFactors = FALSE)
  # inj_df expected to have columns sample, injection_order
  if ("injection_order" %in% names(inj_df)) {
    ord <- inj_df$injection_order
    # reindex to rows of assay_raw
    ord_rows <- match(rownames(assay_raw), inj_df$sample)
    if (all(!is.na(ord_rows))) {
      inj_order <- ord[ord_rows]
    }
  }
}

out_png <- file.path(plots_dir, "qc_loess_weight_comparison.png")

plot_qc_loess_weight_comparison(
  assay_num_base = assay_raw,
  metadata_aligned = metadata_aligned,
  sample_idx = sample_idx,
  qc_idx = qc_idx,
  out_png = out_png,
  qc_loess_span = get0("qc_loess_span", ifnotfound = 0.75),
  qc_loess_min_qc_points = get0("qc_loess_min_qc_points", ifnotfound = 4),
  injection_order = inj_order,
  log2_offset = get0("log2_offset", ifnotfound = 1)
)

cat("Generated:", out_png, "\n")

# QC-RSC comparison
out_png_qcrsc <- file.path(plots_dir, "qc_qcrsc_weight_comparison.png")
tryCatch(
  plot_qc_qcrsc_weight_comparison(
    assay_num_base = assay_raw,
    metadata_aligned = metadata_aligned,
    sample_idx = sample_idx,
    qc_idx = qc_idx,
    out_png = out_png_qcrsc,
    min_qc_points = get0("qc_loess_min_qc_points", ifnotfound = 4),
    injection_order = inj_order,
    log2_offset = get0("log2_offset", ifnotfound = 1)
  ),
  error = function(e) {
    cat("QC-RSC comparison skipped: ", conditionMessage(e), "\n")
  }
)

cat("Generated:", out_png_qcrsc, "\n")
