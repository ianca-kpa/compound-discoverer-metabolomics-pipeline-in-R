#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || !nzchar(args[1])) {
  stop("Usage: Rscript scripts/generate_sample_report.R <diagnostics_dir> [sample] [output_file]")
}

diagnostics_dir <- args[1]
sample <- ifelse(length(args) >= 2 && nzchar(args[2]), args[2], "")
output_file <- ifelse(
  length(args) >= 3 && nzchar(args[3]),
  args[3],
  file.path(
    diagnostics_dir,
    if (nzchar(sample)) {
      paste0("sample_diagnostic_report_", gsub("[^A-Za-z0-9_.-]+", "_", sample), ".html")
    } else {
      "sample_diagnostic_report.html"
    }
  )
)

if (!dir.exists(diagnostics_dir)) {
  stop("Diagnostics directory not found: ", diagnostics_dir)
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

table_to_html <- function(tbl, max_rows = 25) {
  if (is.null(tbl) || nrow(tbl) == 0) return("<p>No rows available.</p>")
  tbl <- head(tbl, max_rows)
  header <- paste0("<th>", html_escape(names(tbl)), "</th>", collapse = "")
  rows <- apply(tbl, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

render_fallback_html <- function() {
  summary_tbl <- read_optional_csv(file.path(diagnostics_dir, "diagnose_samples_summary.csv"))
  if (!is.null(summary_tbl) && nzchar(sample) && "sample" %in% names(summary_tbl)) {
    summary_tbl <- summary_tbl[toupper(summary_tbl$sample) == toupper(sample), , drop = FALSE]
  }
  pca_tbl <- read_optional_csv(file.path(diagnostics_dir, "pca_outlier_candidates.csv"))
  totals_tbl <- read_optional_csv(file.path(diagnostics_dir, "sample_totals_missing_frac.csv"))

  samples_to_show <- character(0)
  if (!is.null(summary_tbl) && "sample" %in% names(summary_tbl)) samples_to_show <- summary_tbl$sample
  if (length(samples_to_show) == 0 && nzchar(sample)) samples_to_show <- sample

  top_sections <- character(0)
  for (sample_name in samples_to_show) {
    top_path <- file.path(diagnostics_dir, paste0("top_features_", gsub("[^A-Za-z0-9_.-]+", "_", sample_name), ".csv"))
    top_sections <- c(
      top_sections,
      paste0("<h3>", html_escape(sample_name), "</h3>", table_to_html(read_optional_csv(top_path), max_rows = 25))
    )
  }
  if (length(top_sections) == 0) {
    top_sections <- "<p>No selected sample was available for top-feature tables.</p>"
  }

  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Sample diagnostic report</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1180px;margin:32px auto;padding:0 18px;color:#17202a}",
    "table{border-collapse:collapse;width:100%;margin:14px 0 28px}th,td{border:1px solid #d5d8dc;padding:6px 8px;text-align:left;font-size:13px}",
    "th{background:#eef2f5}code{background:#f4f6f7;padding:2px 4px}</style></head><body>",
    "<h1>Sample diagnostic report</h1>",
    "<p><strong>Diagnostics directory:</strong> <code>", html_escape(diagnostics_dir), "</code></p>",
    "<h2>Summary</h2>", table_to_html(summary_tbl),
    "<h2>PCA outlier scores</h2>", table_to_html(pca_tbl),
    "<h2>Sample totals</h2>", table_to_html(totals_tbl),
    "<h2>Top features</h2>", paste(top_sections, collapse = "\n"),
    "</body></html>"
  )
  writeLines(html, output_file, useBytes = TRUE)
}

template <- file.path("scripts", "sample_report.Rmd")
can_render_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE) &&
  isTRUE(tryCatch(rmarkdown::pandoc_available("1.12.3"), error = function(e) FALSE)) &&
  file.exists(template)

if (isTRUE(can_render_rmarkdown)) {
  rmarkdown::render(
    input = template,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    params = list(diagnostics_dir = diagnostics_dir, sample = sample),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
} else {
  render_fallback_html()
}

cat("Report written to", output_file, "\n")
