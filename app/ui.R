ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;700&family=Space+Grotesk:wght@500;700&display=swap"
    ),
    tags$link(
      rel = "stylesheet",
      href = paste0(
        "assets/styles.css?v=",
        as.integer(file.info(file.path(app_assets_dir, "styles.css"))$mtime)
      )
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
      tags$h2("CD-MetaboRefine"),
      tags$p(
        "An R/Shiny pipeline for quality control, normalization, statistical analysis, and visualization of Compound Discoverer metabolomics data."
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
          textOutput("install_status_text", container = span)
        ),
        tags$hr(),

        h4("Inputs"),
        tags$p(class = "small-note", "Upload your data, metadata, and reference files here. Accepted formats are CSV, TSV, TXT, XLSX, and XLS. Uploaded files will be copied to the 'data' directory in the project root."),
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
          condition = "input.settings_normalization_mode == 'qc_loess' || input.settings_normalization_mode == 'qcrsc'",
          div(
            class = "big-file-input",
            fileInput(
              "injection_order_file",
              "Injection order file",
              accept = c(".xlsx")
            )
          )
        ),
        uiOutput("allowed_metadata_groups_ui"),
        checkboxInput(
          "manual_metadata_cols",
          "Manually map metadata columns and configure groups",
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.manual_metadata_cols != true",
          tags$div(
            class = "metadata-mapping-alert",
            tags$strong("Model-specific groups need this option."),
            tags$span("Turn it on to detect/configure groups within each model before running the pipeline.")
          )
        ),
        uiOutput("allowed_metadata_groups_hint"),
        tags$hr(),
        conditionalPanel(
          condition = "input.manual_metadata_cols == true",
          tags$div(
            class = "metadata-mapping-panel",
            uiOutput("metadata_model_alias_ui"),
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
              )
            ),
            actionButton("apply_metadata_cols", "Apply metadata/group settings")
          )
        ),
        tags$hr(),
        checkboxInput(
          "use_weight_normalization",
          "Weight normalization",
          value = initial_setting_logical("use_weight_normalization", default = FALSE)
        ),
        checkboxInput(
          "use_reference_file",
          "Use reference file for duplicate matching",
          value = initial_setting_logical("use_reference_file", default = FALSE)
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
              "If column names vary, set them manually below. Leave blank to auto-detect."
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
            "Reference file disabled. Reference input and manual reference-column settings are hidden."
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
                value = initial_setting_value("output_dir", default = "output"),
                width = "100%"
              )
            )
          )
        )
      ),
      
      mainPanel(
        width = 7,
        class = "content-card",
        uiOutput("top_status_banner"),
        fluidRow(
          column(6, actionButton("stop_pipeline", "Stop pipeline", icon = icon("stop"), style = "background-color: #dc3545; color: white;")),
          column(6, actionButton("run_pipeline", "Run pipeline", icon = icon("play")))
        ),
        uiOutput("run_readiness_hint"),
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
                tags$p("Use the form below to configure key variables with standardized input fields.  The variable guide is in the other tab."),
                uiOutput("settings_builder_ui"),
                tags$div(
                  style = "display:none;",
                  class = "save_settings_form",
                  textAreaInput(
                    "config_text",
                    "settings.R content",
                    value = initial_settings_text,
                    rows = 35,
                    width = "auto"
                  )
                )
              ),
              tabPanel(
                "Variable guide",
                tags$p("Use this guide to check what each variable controls before editing the form."),
                uiOutput("settings_glossary_ui")
              )
            )
          ),
          tabPanel(
            "Pipeline Log",
            tags$p("Run pipeline from UI and inspect the latest execution log."),
            uiOutput("pipeline_log_summary"),
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
                style = "display:flex; align-items:center; gap:8px; flex-wrap:wrap;",
                tags$div(
                  class = "results-filter-control",
                  selectInput(
                    "result_gallery_filter",
                    label = NULL,
                    choices = c(
                      "All figures" = "all",
                      "PCA" = "pca",
                      "Volcano" = "volcano",
                      "Heatmap" = "heatmap",
                      "QC / Normalization" = "qc_norm",
                      "Other" = "other"
                    ),
                    selected = "all",
                    width = "180px"
                  )
                ),
                actionButton(style = "background-color: #007bff; color: white; ", "refresh_results_gallery", "Refresh", icon = icon("refresh")),
                actionButton("open_output_dir_gallery", "Open output folder", icon = icon("folder-open"))
              )
            ),
            uiOutput("results_gallery_summary"),
            tags$hr(),
            h4("QC/PCA comparison"),
            tags$p(
              class = "small-note",
              "RSD, drift residual, technical/biological PCA, and retained feature counts from the latest run."
            ),
            tableOutput("qc_pca_comparison_summary"),
            tags$hr(),
            h4("Figures"),
            uiOutput("results_gallery")
          )
          ,
          
        )
      )
      )
    )
  )
)
