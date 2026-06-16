ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;700&family=Space+Grotesk:wght@500;700&display=swap"
    ),
    tags$link(
      rel = "stylesheet",
      href = "assets/styles.css"
    ),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('scrollPipelineLog', function(message) {
        setTimeout(function() {
          var logBox = document.getElementById('pipeline_log_box');

          if (logBox) {
            logBox.scrollTop = logBox.scrollHeight;
          }
        }, 100);
      });
    "))
  ),
  div(
    class = "app-shell",
    div(
      style = "text-align: center;",
      class = "hero",
      tags$h2("Compound Discoverer Refinement and Metabolomics Pipeline"),
      tags$p(
        "Control panel for data upload, normalization, filtering, statistics, MetaboAnalyst-ready exports, and real-time pipeline logs."
      )
    ),
    div(
      class = "layout-shell",
      sidebarLayout(
        sidebarPanel(
          width = 5,
          class = "control-card",
          uiOutput("package_management_ui"),
          div(
            style = "text-align: center;",
            class = "install-status-box",
            textOutput("install_status_text")
          ),
          tags$hr(),
          h4("Inputs"),
          tags$p(class = "small-note", "Upload the Compound Discoverer export, metadata, and optional reference table. Accepted formats are CSV, TSV, TXT, XLSX, and XLS. Uploaded files are copied to the project 'data' directory and written into settings.R."),
          div(
            class = "big-file-input",
            fileInput(
              "data_file",
              "Data file",
              accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
            )
          ),
          div(
            class = "big-file-input",
            fileInput(
              "metadata_file",
              "Metadata file",
              accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
            )
          ),
          conditionalPanel(
            condition = "input.settings_normalization_mode == 'QC_LOESS'",
            selectInput(
              "existing_injection_order_path",
              "Existing input_order file",
              choices = c("None" = "", available_injection_order_files()),
              selected = ""
            ),
            div(
              class = "big-file-input",
              fileInput(
                "injection_order_file",
                "Injection order file (QC-LOESS)",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
              )
            ),
            tags$p(
              class = "small-note",
              "Supports sample + order columns, or Compound Discoverer File Name + Creation Date exports. Upload overrides the existing-file choice."
            )
          ),
          checkboxInput(
            "manual_metadata_cols",
            "Manual metadata setup",
            value = FALSE
          ),
          conditionalPanel(
            condition = "input.manual_metadata_cols == true",
            tags$div(
              style = "margin-top: 10px;",
              tags$h5("Allowed metadata groups by model"),
              tags$p(
                class = "small-note",
                "Use one box per model. The first value is treated as control and the second as test."
              ),
              uiOutput("metadata_model_alias_ui")
            ),
            tags$div(
              class = "metadata-mapping-panel",
              checkboxInput(
                "show_metadata_column_fields",
                "Show metadata column fields",
                value = FALSE
              ),
              conditionalPanel(
                condition = "input.show_metadata_column_fields == true",
                div(
                  class = "grid-3x2",
                  textInput("metadata_col_sample", "Sample column", value = ""),
                  textInput("metadata_col_weight", "Weight column", value = ""),
                  textInput("metadata_col_group", "Group column", value = ""),
                  textInput("metadata_col_sex", "Sex column", value = ""),
                  textInput("metadata_col_model", "Model column", value = "")
                ),
                actionButton("apply_metadata_cols", "Apply metadata column mapping")
              )
            )
          ),
          tags$hr(),
          checkboxInput(
            "use_weight_normalization",
            tagList(
              "Weight normalization",
              tags$span(
                "   -> divides each sample by its metadata weight before the main normalization step.",
                class = "small-note"
              )
            ),
            value = setting_display_logical(initial_settings_text, "use_weight_normalization", default = FALSE)
          ),
          checkboxInput(
            "use_reference_file",
            "Use reference file for duplicate matching",
            value = setting_display_logical(initial_settings_text, "use_reference_file", default = FALSE)
          ),
          conditionalPanel(
            condition = "input.use_reference_file == true",
            div(
              class = "big-file-input",
              fileInput(
                "reference_file",
                "Reference file",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
              )
            ),
            checkboxInput(
              "manual_reference_cols",
              "Manually set reference column names",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.manual_reference_cols == true",
              tags$p(
                class = "small-note",
                "If reference column names vary, set them manually below. Leave blank to auto-detect metabolite/name, reference ion, m/z, and RT columns."
              ),
              div(
                class = "grid-2x2",
                textInput("reference_col_metabolite", "Reference metabolite column", value = ""),
                textInput("reference_col_ref_ion", "Reference ion column", value = ""),
                textInput("reference_col_mz", "Reference m/z column", value = ""),
                textInput("reference_col_rt", "Reference RT column", value = "")
              ),
              actionButton("apply_reference_cols", "Apply reference column names")
            )
          ),
          conditionalPanel(
            condition = "input.use_reference_file == false",
            tags$p(
              class = "small-note",
              "Reference file disabled. Duplicate handling will use the selected non-reference strategy, usually best QC RSD."
            )
          ),
          tags$hr(),
          h4("Configuration"),
          div(
            class = "output-dir-block",
            tags$label(`for` = "output_dir", "Output directory"),
            div(
              class = "output-dir-row",
              div(class = "output-dir-button-wrap", uiOutput("output_dir_browser_ui")),
              div(
                class = "output-dir-input-wrap",
                textInput(
                  "output_dir",
                  label = NULL,
                  value = setting_display_value(initial_settings_text, "output_dir", default = "output"),
                  width = "100%"
                )
              )
            ),
            textOutput("output_dir_status")
          ),
          tags$p(
            class = "small-note",
            "Analysis options are configured in the Settings Builder cards."
          )
        ),
        mainPanel(
          width = 7,
          class = "content-card",
          uiOutput("top_status_banner"),
          fluidRow(
            column(6, actionButton("stop_pipeline", "Stop pipeline", icon = icon("stop"))),
            column(6, actionButton("run_pipeline", "Run pipeline", icon = icon("play")))
          ),
          tags$hr(),
          tabsetPanel(
            id = "main_tabs",
            tabPanel(
              "Data - Overview",
              h4("Data Overview"),
              uiOutput("data_overview")
            ),
            tabPanel(
              "Settings Builder",
              tabsetPanel(
                id = "settings_subtabs",
                tabPanel(
                  "Settings form",
                  tags$p("Use the form below to configure the most common run settings. The variable guide explains what each field controls."),
                  uiOutput("settings_builder_ui"),
                  tags$div(
                    style = "display:none;",
                    class = "save_settings_form",
                    textAreaInput("config_text", "settings.R content", value = initial_settings_text, rows = 35, width = "auto")
                  )
                ),
                tabPanel(
                  "Variable guide",
                  tags$p("Use this guide to check normalization, filtering, duplicate-handling, statistics, PCA, and heatmap settings before editing the form."),
                  uiOutput("settings_glossary_ui")
                )
              )
            ),
            tabPanel(
              "Pipeline Log",
              tags$p("Run the pipeline from the app and inspect the latest execution log as it updates."),
              tags$pre(
                id = "pipeline_log_box",
                style = "max-height: 850px; overflow-y: auto; background: #111; color: #f2f2f2; padding: 12px;",
                textOutput("pipeline_log")
              )
            ),
            tabPanel(
              "Results Gallery",
              tags$div(
                style = "display:flex; align-items:center; justify-content:space-between; gap:10px; margin-bottom:10px;",
                h4(style = "margin:0;", "Results Gallery"),
                tags$div(
                  style = "display:flex; align-items:center; gap:8px;",
                  actionButton(
                    "refresh_results_gallery",
                    "Refresh",
                    icon = icon("refresh"),
                    style = "background-color: #007bff; color: white;"
                  ),
                  actionButton("open_output_dir_gallery", "Open output folder", icon = icon("folder-open"))
                )
              ),
              h4("QC/PCA comparison"),
              tags$p(
                class = "small-note",
                "Retained features, QC/sample RSD, IQR filtering, QC-LOESS drift audit, and PCA figure counts from the latest run."
              ),
              tableOutput("qc_pca_comparison_summary"),
              tags$hr(),
              h4("QC-LOESS adjustment summary"),
              tags$p(
                class = "small-note",
                "Pre/post QC-RSD and drift-correlation metrics showing how much the QC-LOESS correction changed technical variability."
              ),
              tableOutput("loess_adjustment_summary"),
              tags$hr(),
              h4("Figures"),
              uiOutput("results_gallery")
            ),
          )
        )
      )
    )
  )
)
