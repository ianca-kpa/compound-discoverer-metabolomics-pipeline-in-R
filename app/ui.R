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
        "Control panel for inputs, settings, MetaboAnalyst-ready exports, and real-time pipeline execution logs."
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
        checkboxInput(
          "manual_metadata_cols",
          "Manually map metadata columns and configure groups",
          value = FALSE
        ),
        # textInput(
        #   "allowed_metadata_groups",
        #   "Allowed metadata groups (control first, test second)",
        #   value = "WT, TG"
        # ),
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
              ),
              actionButton("apply_metadata_cols", "Apply metadata column mapping")
            )
          )
        ),
        checkboxInput(
          "use_weight_normalization",
          tagList(
            "Weight normalization",
            tags$span(
              "   → Weight normalization controls weight-based scaling.",
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
                value = setting_display_value(initial_settings_text, "output_dir", default = "output"),
                width = "100%"
              )
            )
          ),
          textOutput("output_dir_status")
        ),
        checkboxInput(
          "minimal_output",
          tagList(
            "Minimal output",
            tags$span(
              "   → keeps plots/statistics based on selected options; skips selected intermediate global exports.",
              class = "small-note"
            )
          ),
          value = setting_display_logical(initial_settings_text, "minimal_output", default = FALSE)
        ),
        conditionalPanel(
          condition = "input.minimal_output == true",
          tags$div(
            class = "small-note",
            style = "margin-top:6px; color:#9a3412;",
            "Minimal output enabled: plots and statistics are kept according to your selected options. Only selected intermediate global exports are skipped."
          )
        ),
        div(
          class = "grid-3x1",
          selectInput(
            "duplicate_name_strategy",
            "Duplicate handling strategy",
            choices = c(
              "reference_or_best_qc_rsd",
              "keep_separate",
              "collapse_mean",
              "collapse_sum",
              "collapse_best_qc_rsd"
            ),
            selected = setting_display_value(initial_settings_text, "duplicate_name_strategy", default = "collapse_best_qc_rsd")
          ),
          selectInput(
            "run_metrics",
            "Run metrics",
            choices = c("FDR", "p_value", "FDR_and_p_value"),
            selected = setting_display_value(initial_settings_text, "run_metrics", default = "FDR_and_p_value")
          ),
          selectInput(
            "use_only_known",
            "Use only known metabolites",
            choices = c("TRUE", "FALSE"),
            selected = if (setting_display_logical(initial_settings_text, "use_only_known", default = TRUE)) "TRUE" else "FALSE"
          )
        ),
        tags$hr(),
        h4("Significance Testing"),
        tags$p(class = "small-note", "Configure statistical test type, pairing, and p-value correction method."),
        div(
          class = "grid-3x1",
          selectInput(
            "statistical_test_type",
            "Statistical test",
            choices = c("student", "welch", "wilcoxon", "limma"),
            selected = setting_display_value(initial_settings_text, "statistical_test_type", default = "student")
          ),
          selectInput(
            "test_is_paired",
            "Test type",
            choices = c("Unpaired" = "FALSE", "Paired" = "TRUE"),
            selected = if (setting_display_logical(initial_settings_text, "test_is_paired", default = FALSE)) "TRUE" else "FALSE"
          )
        )

      ),
      
      mainPanel(
        width = 7,
        class = "content-card",
        uiOutput("top_status_banner"),
        fluidRow(
          column(6, actionButton("stop_pipeline", "Stop pipeline")),
          column(6, actionButton("run_pipeline", "Run pipeline"))
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
                tags$p("Use the form below to configure key variables with standardized input fields.  The variable guide is in the other tab."),
                uiOutput("settings_builder_ui"),
                tags$div(
                  style = "display:none;",
                  class = "save_settings_form",
                  textAreaInput(
                    "config_text",
                    "settings.R content",
                    value = read_initial_config(),
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
            tags$pre(
              id = "pipeline_log_box",
              style = "max-height: 850px; overflow-y: auto; background: #111; color: #f2f2f2; padding: 12px;",
              textOutput("pipeline_log")
            )
          ),
          tabPanel(
            "Gallery Results",
            tags$div(
              style = "display:flex; align-items:center; justify-content:space-between; gap:10px; margin-bottom:10px;",
              h4(style = "margin:0;", "Results Gallery"),
              actionButton("open_output_dir_gallery", "Open output folder")
            ),
            uiOutput("results_gallery")
          )
          ,
          
        )
      )
      )
    )
  )
)
