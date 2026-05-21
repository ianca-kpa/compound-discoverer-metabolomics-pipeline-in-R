# Compound Discoverer Refinement and Metabolomics Pipeline

A modular R pipeline for untargeted metabolomics data processing, focused on:

- Compound Discoverer table refinement,
- normalization and filtering,
- statistical analysis and visualization,
- MetaboAnalyst-ready exports,
- execution via script or Shiny app.

---

## Overview

This project processes Compound Discoverer feature tables together with experimental metadata and produces:

- cleaned assay matrices,
- filter audit reports,
- PCA outputs,
- statistical test results,
- volcano plots,
- heatmaps,
- downstream export tables.

The workflow is orchestrated by `pipeline/run_pipeline.R` and modularized in `pipeline/R/00` through `pipeline/R/12`.

---

## Repository Structure

```text
.
|-- app.R
|-- app/
|   |-- global.R
|   |-- server.R
|   |-- ui.R
|   `-- assets/
|-- pipeline/
|   |-- run_pipeline.R
|   |-- config/
|   |   |-- settings.example.R
|   |   `-- settings.R (local, not versioned)
|   `-- R/
|       |-- 00_packages.R
|       |-- 01_validation.R
|       |-- 02_comparisons.R
|       |-- 03_helpers_io_log.R
|       |-- 04_metadata.R
|       |-- 05_features_assay.R
|       |-- 06_normalization_filters.R
|       |-- 07_duplicates.R
|       |-- 08_exports.R
|       |-- 09_pca.R
|       |-- 10_stats_volcano.R
|       |-- 11_heatmaps.R
|       `-- 12_main_pipeline.R
|-- images/
`-- output/ (generated at runtime)
```

---

## Requirements

- R ≥ 4.5.3 (recommended)
- RStudio Desktop (recommended)
- Required CRAN packages are installed automatically by `pipeline/R/00_packages.R`

---

# Installation

## 1. Install R

Download and install R from:

https://cran.r-project.org/

---

## 2. Install RStudio Desktop

Download and install RStudio Desktop from:

https://posit.co/download/rstudio-desktop/

Recommended order:

1. Install R
2. Install RStudio Desktop
3. Open the project `.Rproj` file in RStudio

---

## 3. Download or Clone the Repository

### Option A — Download ZIP

Repository:

https://github.com/ianca-kpa/compound-discoverer-metabolomics-pipeline-in-R

1. Click **Code**
2. Click **Download ZIP**
3. Extract the archive
4. Open the extracted folder

---

### Option B — Clone with Git

```powershell
git clone https://github.com/ianca-kpa/compound-discoverer-metabolomics-pipeline-in-R.git
cd compound-discoverer-metabolomics-pipeline-in-R
```

---

# Initial Setup

## 1. Open the Project in RStudio

Open the `.Rproj` file located in the repository root directory.

Using the `.Rproj` file ensures that relative paths work correctly.

---

## 2. Create the Local Settings File

In the RStudio Console:

```r
file.copy(
  "pipeline/config/settings.example.R",
  "pipeline/config/settings.R",
  overwrite = FALSE
)
```

Or in PowerShell:

```powershell
Copy-Item pipeline/config/settings.example.R pipeline/config/settings.R
```

The `settings.R` file is local and should contain your input paths and analysis parameters.

---

# Configuring the Analysis

Main fields commonly edited in `pipeline/config/settings.R`:

- Input files:
  - `cd_file_path`
  - `metadata_path`
  - `comparison_path`

- Comparison groups:
  - `comparison_group_control`
  - `comparison_group_treatment`

- Output:
  - `output_dir`

- Filters and thresholds:
  - `missing_exclusion_max_fraction`
  - `presence_filter_min_fraction`
  - `alpha_sig`
  - `fc_cutoff_log2`

- Export options:
  - `export_metaboanalyst_ready`
  - `save_stats_excel_per_model`

---

# Running the Pipeline

## Option A — Run Through the Shiny App

Launch the app:

```r
shiny::runApp(".")
```

The app allows you to:

- load Compound Discoverer tables and metadata,
- configure groups and filters,
- generate `settings.R`,
- run the pipeline,
- monitor logs,
- review generated plots and exports.

To start execution:

1. Upload the required files
2. Configure parameters if needed
3. Choose the output directory
4. Click **Run pipeline**

---

## Option B — Run Directly from RStudio

```r
source("pipeline/run_pipeline.R")
```

---

## Option C — Run from the Command Line

```powershell
Rscript pipeline/run_pipeline.R
```

---

# Expected Outputs

Results are written to `output/` (or the folder defined by `output_dir`).

Main outputs include:

- `PIPELINE_LOG.txt`
- filter audit tables,
- processed matrices,
- statistical outputs,
- PCA plots,
- volcano plots,
- heatmaps,
- MetaboAnalyst-ready exports.

Example structure:

```text
output/
|-- PIPELINE_LOG.txt
|-- global/
|   |-- audits_global/
|   `-- exports_global/
`-- <MODEL>/
    |-- exports/
    |   |-- metaboanalyst/
    |   `-- stats/
    `-- plots/
        |-- pca/
        |-- volcano/
        |-- heatmap/
        `-- heatmap_significant/
```

---

# Troubleshooting

### `pipeline/config/settings.R not found`

Create the local settings file from `settings.example.R`.

---

### Missing package error

Re-run the pipeline to trigger automatic package installation.

---

### Spreadsheet read error

Check:

- file extension,
- sheet name/index,
- file paths in `settings.R`.

---

### `Rscript` not recognized

Run the pipeline from RStudio:

```r
source("pipeline/run_pipeline.R")
```

Or add R to the system PATH.

---

# Contributing

Contributions are welcome through issues or pull requests with a clear description of:

- the proposed change,
- the experimental context,
- the expected impact on the pipeline.
