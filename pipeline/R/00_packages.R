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

bioc_packages <- c(
  "pmp",
  "limma"
)

required_packages <- unique(c(cran_packages, app_packages, bioc_packages))

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
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
  library(pmp)
  library(limma)
})

message("R version: ", R.version.string)
