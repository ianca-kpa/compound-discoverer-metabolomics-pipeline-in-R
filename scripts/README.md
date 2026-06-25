# Scripts

This folder contains small runners and diagnostics around the main pipeline.

## Active scripts

- `diagnose_sample.R`: canonical sample diagnostic runner. Use `PCA_OUTLIERS` to diagnose every sample row that is very distant from the rest in PC1/PC2 space, or pass a specific sample name such as `SAMPLE_ID`.
- `generate_sample_report.R`: renders an HTML report from a diagnostics directory produced by `diagnose_sample.R`. It uses `rmarkdown` when Pandoc is available and falls back to a simple self-contained HTML report otherwise.
- `generate_qc_loess_weight_plot.R`: manual runner for the QC-LOESS/QC-RSC weight/no-weight comparison plots. The main pipeline can also generate these plots when QC diagnostics are enabled.

Usage:

```powershell
Rscript scripts/diagnose_sample.R PCA_OUTLIERS data/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv
Rscript scripts/diagnose_sample.R PCA_OUTLIERS data/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv "" output/diagnostics_pca_outliers 3.5
Rscript scripts/diagnose_sample.R SAMPLE_ID data/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv
Rscript scripts/diagnose_sample.R SAMPLE_ID output/<run>/global/exports_global/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv path/to/metadata.xlsx output/diagnostics_SAMPLE_ID_metadata
Rscript scripts/generate_sample_report.R output/diagnostics_pca_outliers SAMPLE_ID
Rscript scripts/generate_qc_loess_weight_plot.R output/<run>/global/exports_global output/<run>/global/plots_global/normalization
```

The fifth argument is the robust PCA distance threshold. The default is `3.5`.
The script supports both Compound Discoverer-style matrices with `Area: sample` columns and MetaboAnalyst-style matrices with a `sample` column plus feature columns.
When no input CSV is provided, `diagnose_sample.R` first tries the configured `output_dir` export and then falls back to `data/MA_ACTIVE_duplicate_ONLY_GLOBAL_NO_QC.csv`.

## Legacy scripts

Files under `scripts/legacy/` are retained for reference while older diagnostics are consolidated. Prefer the active scripts above for new runs.
