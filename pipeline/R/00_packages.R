# =============================================================================
# 00_packages.R
# Package installation/loading
# =============================================================================

message("Step 1: Checking / installing dependencies...")

cran_packages <- c(
  "tidyverse", "readr", "readxl", "openxlsx",
  "pheatmap", "ggrepel", "stringi", "RColorBrewer"
)

app_packages <- c(
  "processx",
  "shinyFiles",
  "magick"
)

required_packages <- unique(c(cran_packages, app_packages))

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(pheatmap)
  library(ggrepel)
  library(stringi)
  library(RColorBrewer)
})

message("R version: ", R.version.string)