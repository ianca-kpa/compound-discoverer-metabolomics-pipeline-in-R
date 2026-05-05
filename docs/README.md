# Compound Discoverer Refinement and Metabolomics Pipeline

A modular R pipeline for untargeted metabolomics data processing, focused on:
- Compound Discoverer table refinement,
- robust normalization and filtering,
- statistical analysis and visualization,
- MetaboAnalyst-ready exports,
- execution via script and Shiny app.

## Overview

This project processes Compound Discoverer feature tables together with experimental metadata and produces:
- cleaned assay matrices,
- filter audit reports,
- PCA outputs,
- statistical test results,
- volcano plots,
- heatmaps,
- downstream export tables.

The main workflow is orchestrated by `pipeline/run_pipeline.R` and modularized in `pipeline/R/00` through `pipeline/R/12`.

## Repository structure

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

## Requirements

- R 4.5.3 (recommended)
- RStudio Desktop (recommended)
- CRAN packages installed automatically by `pipeline/R/00_packages.R`

## Guided Tutorial For First-Time Users

This section explains how to install the required programs, obtain the scripts, open the project, and launch the app.

### 1. Install R

R is the programming language used by the pipeline and the Shiny app.

1. Go to https://cran.r-project.org/
2. Click the download link for your operating system.
   - Windows: click **Download R for Windows**, then **base**, and download the installer.
   - macOS: click **Download R for macOS** and choose the recommended installer.
3. Run the installer and keep the default options.

After installation, you do not need to open R directly. The recommended interface for using the project is RStudio.

### 2. Install RStudio Desktop

RStudio is the program used to open, edit, and run the project more easily.

1. Go to https://posit.co/download/rstudio-desktop/
2. Download **RStudio Desktop** for your operating system.
3. Run the installer and keep the default options.

Recommended order:
1. Install R first.
2. Install RStudio Desktop second.
3. Open the project file `project/metabolomics-pipeline_git.Rproj` in RStudio.

### 3. Obtain The Project Scripts

You can obtain the project in two ways.

#### Option A: download as ZIP

Use this option if you do not use Git.

1. Open the GitHub repository page:
   https://github.com/ianca-kpa/metabolomics-pipeline_git
2. Click **Code**.
3. Click **Download ZIP**.
4. Extract the ZIP file to a folder on your computer.
5. Open the extracted folder.

Important: run the project from the extracted folder, not from inside the compressed ZIP file.

#### Option B: clone with Git

Use this option if Git is installed on your computer.

```powershell
git clone https://github.com/ianca-kpa/metabolomics-pipeline_git.git
cd metabolomics-pipeline_git
```

This command creates a local copy of all scripts and folders.

### 4. Open The Project In RStudio

1. Open RStudio Desktop.
2. Go to **File > Open Project...**
3. Select:

```text
project/metabolomics-pipeline_git.Rproj
```

Opening the `.Rproj` file helps RStudio use the correct project folder. This makes relative paths such as `pipeline/run_pipeline.R`, `app.R`, and `output/` work correctly.

### 5. Start The App

In RStudio, open `app.R` and click **Run App**.

You can also run this command in the RStudio Console:

```r
shiny::runApp(".")
```

The first run may take longer because R checks and installs missing packages.

### 6. Create The Local Settings File

Before running the pipeline, create your personal settings file from the example file.

In the RStudio Console:

```r
file.copy("pipeline/config/settings.example.R", "pipeline/config/settings.R", overwrite = FALSE)
```

Or in PowerShell:

```powershell
Copy-Item pipeline/config/settings.example.R pipeline/config/settings.R
```

The file `settings.R` is local and should contain your own input paths and analysis parameters.

### 7. Choose Input Files And Parameters

You can configure the analysis in two ways:

- Through the Shiny app interface.
- By editing `pipeline/config/settings.R` directly in RStudio.

Key fields to review:
- Input files: `cd_file_path`, `metadata_path`, `comparison_path`
- Comparison groups: `comparison_group_control`, `comparison_group_treatment`
- Output directory: `output_dir`
- Filters and thresholds: `missing_exclusion_max_fraction`, `presence_filter_min_fraction`, `alpha_sig`, `fc_cutoff_log2`
- Export options: `export_metaboanalyst_ready`, `save_stats_excel_per_model`

### 8. Run The Pipeline

Recommended option for new users: run it through the Shiny app interface.

In the left side panel:
- upload the main file under **Data file**,
- upload the metadata spreadsheet under **Metadata file**,
- enable **Manually map metadata columns and configure groups** if you need to provide column names manually,
- enable **Weight normalization** if weight-based normalization is part of the analysis,
- enable **Use reference file for duplicate matching** if you will use a reference file for duplicate handling,
- choose the folder under **Output directory**,
- adjust **Duplicate handling strategy**, **Run metrics**, and **Use only known metabolites**.

In the main panel:
- use **Data - Overview** to review the loaded data,
- use **Settings Builder** to adjust advanced parameters and save `settings.R`,
- click **Run pipeline** to start execution,
- monitor execution in **Pipeline Log**,
- open **Gallery Results** to view generated figures and access the results folder,
- use **Stop pipeline** only if you need to interrupt a running execution.

You can also run the pipeline directly.

In the RStudio Console:

```r
source("pipeline/run_pipeline.R")
```

Or in PowerShell:

```powershell
Rscript pipeline/run_pipeline.R
```

### 9. Find The Results

By default, results are saved in `output/`, unless another output folder is configured in `settings.R` or in the app.

The main execution log is:

```text
output/PIPELINE_LOG.txt
```

Use this log to check which steps were completed and to diagnose errors.

## Quick Start From The Command Line

### 1. Clone the repository

```powershell
git clone https://github.com/ianca-kpa/metabolomics-pipeline_git.git
cd metabolomics-pipeline_git
```

### 2. Create your local settings file

```powershell
Copy-Item pipeline/config/settings.example.R pipeline/config/settings.R
```

### 3. Configure input paths and parameters

Edit `pipeline/config/settings.R` with your file paths and experimental rules.

Key fields to review:
- Input files: `cd_file_path`, `metadata_path`, `comparison_path`
- Comparison groups: `comparison_group_control`, `comparison_group_treatment`
- Output directory: `output_dir`
- Filters and thresholds: `missing_exclusion_max_fraction`, `presence_filter_min_fraction`, `alpha_sig`, `fc_cutoff_log2`
- Export options: `export_metaboanalyst_ready`, `save_stats_excel_per_model`

### 4. Run the pipeline

Option A (terminal):

```powershell
Rscript pipeline/run_pipeline.R
```

Option B (R or RStudio):

```r
source("pipeline/run_pipeline.R")
```

### 5. Run through the Shiny control panel (optional)

```r
shiny::runApp(".")
```

From the app you can:
- review scripts,
- edit settings,
- validate inputs,
- launch the pipeline,
- monitor logs in real time.

The app separates quick run controls from advanced settings:
- the left control panel defines input files, output folder, comparison groups, duplicate handling, run metrics, and whether to keep only known metabolites,
- the Settings tab keeps advanced numerical and plotting parameters such as significance thresholds, PCA/heatmap scaling, heatmap generation, and heatmap clustering.

Comparison groups are configured in the left panel as a comma-separated pair. The first value is treated as control and the second as treatment, for example `WT, TG` or `Control, Treated`.

Recommended heatmap clustering options are intentionally limited in the app:
- distance: `euclidean` or `manhattan`,
- method: `ward.D2`, `complete`, or `average`.

## Expected outputs

Outputs are written to `output/` (or the folder configured by `output_dir`), including:
- `PIPELINE_LOG.txt` for execution traceability,
- filter audit tables,
- processed matrices by step,
- statistical outputs by model/comparison,
- PCA, volcano plots, and heatmaps,
- MetaboAnalyst-compatible exports (when enabled).

The output tree is intentionally compact to keep results easy to browse:

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
    |       |-- STATS_ACTIVE_model_<MODEL>.xlsx
    |       `-- significant_metabolites/
    |           |-- p_value/
    |           `-- FDR/
    `-- plots/
        |-- pca/
        |-- volcano/
        |-- heatmap/
        `-- heatmap_significant/
```

Plot folders are grouped by plot type, not by comparison. The comparison is encoded in the filename, for example `tg_vs_wt`, `f_vs_m`, `tg-f_vs_wt-f`, `tg-m_vs_wt-m`, `tg-f_vs_tg-m`, or `wt-f_vs_wt-m`.

Significant-metabolite TXT exports are metric-specific:
- files under `significant_metabolites/p_value/` use `p_value < alpha_sig`,
- files under `significant_metabolites/FDR/` use `FDR < alpha_sig`.

Significant heatmaps can also require the fold-change cutoff when `sig_heatmap_require_fc_cutoff = TRUE`, using `abs(log2FC) >= fc_cutoff_log2`.

## Troubleshooting

- Error: `pipeline/config/settings.R not found`
  - Copy `pipeline/config/settings.example.R` to `pipeline/config/settings.R`.
- Missing package error
  - Re-run the pipeline to trigger automatic package installation.
- Spreadsheet read error
  - Check file extension, sheet name/index, and file paths in `settings.R`.
- `Rscript` not recognized in terminal
  - Run with `source("pipeline/run_pipeline.R")` in RStudio, or add R to your Windows PATH.

## Contributing

Contributions are welcome through issues or pull requests with a clear description of the change, experimental context, and expected impact on the pipeline.
