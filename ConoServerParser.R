# ====================================================================================================
# ConoServerParser: 
# A Flexible Shiny Application for Filtering, Exploring, and Exporting Conotoxin Sequences from ConoServer
# ====================================================================================================
# Author: Dany Domínguez-Pérez

# -------------------------
# LOAD REQUIRED LIBRARIES
# -------------------------
library(shiny)
library(seqinr)
library(DT)
library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(plotly)
library(readr)
library(writexl)
library(zip)
library(htmlwidgets)
library(shinycssloaders)
library(stringr)

# --- Example File Paths (Adjust as needed) ---
EXAMPLE_CONOSERVER_FASTA_PATH <- "data/conoserver_protein.fa"
EXAMPLE_CONOPREC_CSV_PATH <- "data/Y_superfamily.csv"

# -------------------------
# PARSE FASTA FUNCTION 
# -------------------------
parse_fasta <- function(filepath) {
  lines <- readLines(filepath)
  if (length(lines) == 0) { stop("FASTA file is empty or could not be read.") }
  headers_idx <- grep("^>", lines)
  if (length(headers_idx) == 0) { stop("No FASTA headers found. Check file format.") }
  sequences <- sapply(seq_along(headers_idx), function(i) {
    start <- headers_idx[i] + 1
    end <- if (i < length(headers_idx)) headers_idx[i + 1] - 1 else length(lines)
    if (start > end) return("")
    paste(lines[start:end], collapse = "")
  })
  headers <- lines[headers_idx]
  metadata <- strsplit(sub("^>", "", headers), "\\|")
  if (length(metadata) == 0) { 
    max_fields <- 0
    df_meta <- data.frame(matrix(ncol = 0, nrow = length(headers)))
  } else {
    max_fields <- tryCatch(max(sapply(metadata, length)), error = function(e) 0)
    if (!is.numeric(max_fields) || max_fields < 0) max_fields <- 0
    if (max_fields > 0) {
      df_meta_list <- lapply(metadata, function(m) { 
        length(m) <- max_fields
        m[is.na(m)] <- "Unclassified"
        m[m == ""] <- "Unclassified"
        m 
      })
      df_meta <- do.call(rbind, df_meta_list)
      df_meta <- as.data.frame(df_meta, stringsAsFactors = FALSE)
    } else { 
      df_meta <- data.frame(matrix(ncol = 0, nrow = length(headers))) 
    }
  }
  base_col_names <- c("ConoID", "Name", "Organism", "ProteinType", "ConopeptideClass", 
                      "Superfamily", "Pharmacological", "Framework", "SequenceEvidence")
  col_names <- character(max_fields)
  if (max_fields > 0) {
    common_length <- min(length(base_col_names), max_fields)
    col_names[1:common_length] <- base_col_names[1:common_length]
    if (max_fields > length(base_col_names)) { 
      col_names[(length(base_col_names) + 1):max_fields] <- paste0("ExtraField_", seq_len(max_fields - length(base_col_names))) 
    }
    if (ncol(df_meta) > 0) { 
      if (length(col_names) == ncol(df_meta)) { 
        colnames(df_meta) <- col_names 
      } else { 
        warning("Mismatch: col names vs metadata.") 
      } 
    }
  }
  if (nrow(df_meta) != length(headers)) {
    warning("Mismatch: headers vs metadata rows.")
    if (ncol(df_meta) == 0) { 
      df_meta <- data.frame(matrix(ncol = 0, nrow = length(headers)))
    } else { 
      if (nrow(df_meta) > length(headers)) { 
        df_meta <- df_meta[1:length(headers), , drop = FALSE] 
      } 
    }
  }
  df <- data.frame(Header = headers, Sequence = sequences, stringsAsFactors = FALSE)
  if (ncol(df_meta) > 0 && nrow(df_meta) == nrow(df)) { 
    df <- cbind(df, df_meta) 
  }
  standard_cols <- c("Header", "Sequence", "ConoID", "Name", "Organism", "ProteinType", 
                     "ConopeptideClass", "Superfamily", "Pharmacological", "Framework", "SequenceEvidence")
  for (col in standard_cols) { 
    if (!col %in% names(df)) { 
      df[[col]] <- "Unclassified" 
    } 
  }
  if ("Superfamily" %in% names(df)) { 
    df$Superfamily <- trimws(df$Superfamily)
    df$Superfamily <- ifelse(is.na(df$Superfamily) | df$Superfamily == "" | df$Superfamily == "NA", 
                             "Unclassified", df$Superfamily) 
  }
  if ("Pharmacological" %in% names(df)) {
    df$Pharmacological <- trimws(df$Pharmacological)
    df$Pharmacological <- ifelse(is.na(df$Pharmacological) | df$Pharmacological == "" | 
                                   df$Pharmacological == "NA", "Unclassified", df$Pharmacological)
    known_pharmas <- c("alpha conotoxin", "chi conotoxin", "delta conotoxin", "epsilon conotoxin", 
                       "gamma conotoxin", "iota conotoxin", "kappa conotoxin", "mu conotoxin", 
                       "omega conotoxin", "rho conotoxin", "sigma conotoxin", "tau conotoxin")
    df$Pharmacological <- sapply(df$Pharmacological, function(ph) { 
      ph_trimmed <- trimws(ph)
      if (is.na(ph_trimmed) || ph_trimmed == "" || ph_trimmed == "Unclassified" || ph_trimmed == "NA") { 
        return("Unclassified") 
      }
      ph_lower <- tolower(ph_trimmed)
      match_idx <- match(ph_lower, tolower(known_pharmas))
      if (!is.na(match_idx)) { 
        known_pharmas[match_idx] 
      } else { 
        ph_trimmed 
      } 
    }, USE.NAMES = FALSE)
  }
  if ("ProteinType" %in% names(df)) { 
    df$ProteinType <- tolower(trimws(df$ProteinType))
    df$ProteinType <- ifelse(is.na(df$ProteinType) | df$ProteinType == "", "Unclassified", df$ProteinType) 
  }
  if ("Level" %in% names(df)) { 
    names(df)[names(df) == "Level"] <- "SequenceEvidence" 
    df$SequenceEvidence <- trimws(df$SequenceEvidence)
    df$SequenceEvidence <- ifelse(is.na(df$SequenceEvidence) | df$SequenceEvidence == "" | df$SequenceEvidence == "NA", "Unclassified", df$SequenceEvidence) 
  }
  extra_cols <- setdiff(names(df), standard_cols)
  present_standard_cols <- intersect(standard_cols, names(df))
  final_col_order <- c(present_standard_cols, extra_cols)
  df <- df[, final_col_order, drop = FALSE]
  return(df)
}

# -------------------------
# UI DEFINITION 
# -------------------------
ui <- fluidPage(
  titlePanel("ConoServerParser: Conotoxin FASTA Analysis & Exploration Suite"),
  tags$head(tags$script(HTML("
    Shiny.addCustomMessageHandler('selectFilteredRows', function(data) {
      var table = $('#seq_table table').DataTable();
      if (table) {
        setTimeout(function() {
          if (data.deselect) { table.rows({ search: 'applied' }).deselect(); } 
          else { table.rows({ search: 'applied' }).select(); }
        }, 100);
      } else { console.error('DataTable instance not found for #seq_table'); }
    });
  "))),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.main_tabs == 'Metadata Table' || input.main_tabs == 'Summary Statistics' || 
                     input.main_tabs == 'FASTA Preview' || input.main_tabs == 'Plots'",
        h4("ConoServer Data Upload & Processing"),
        fileInput("fasta_file", "Upload ConoServer FASTA File:", accept = c(".fa", ".fasta")),
        actionButton("load_example_precursor_fasta", "Load Example ConoServer FASTA", 
                     class = "btn-block btn-sm btn-info"),
        br(),
        h4("Column Selection for Table"),
        uiOutput("column_selector"),
        h4("Filters for Data"),
        uiOutput("column_filters"),
        splitLayout(
          actionButton("clear_filters", "Clear All Filter Values", class = "btn btn-secondary btn-sm"),
          actionButton("reset_filters", "Reset All Filters", class = "btn btn-secondary btn-sm")
        ),
        h4("FASTA Export Options"),
        uiOutput("fasta_header_selector"),
        selectInput("header_sep", "FASTA Header Separator:", 
                    choices = c("|", " ", ",", ";", "_", "~", ":"), selected = "|"),
        downloadButton("download_filtered_fasta", "Download Selected Sequences (.fasta)"),
        br(), br(),
        downloadButton("download_zip_by_superfamily", "Download FASTA by Gene Superfamily (.zip)"),
        br(),
        downloadButton("download_zip_by_pharmacological", "Download FASTA by Pharmacological Family (.zip)"),
        br(),
        tags$div(
          style = "max-height: 250px; overflow-y: auto; border: 1px solid #ddd; padding: 5px; 
                   margin-top: 5px; background-color: #f9f9f9;",
          verbatimTextOutput("fasta_stats_message")
        )
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Mature Peptide FASTA Generator'",
        h4("Mature Peptide Data (from ConoPrec)"),
        fileInput("conoprec_csv_files", "Upload ConoPrec CSV File(s):", multiple = TRUE, 
                  accept = c("text/csv", ".csv")),
        actionButton("load_example_conoprec_csv", "Load Example ConoPrec CSV", 
                     class = "btn-block btn-sm btn-info"),
        br(),
        uiOutput("mature_peptide_summary_ui"),
        downloadButton("download_mature_fasta_zip", "Download Mature Peptides FASTA (.zip)")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Metadata Table", 
                 DT::DTOutput("seq_table"), 
                 br(),
                 actionButton("select_filtered", "Select Filtered Rows", 
                              class = "btn btn-primary btn-sm"),
                 actionButton("deselect_filtered", "Deselect Filtered Rows", 
                              class = "btn btn-secondary btn-sm"), 
                 br(), br(),
                 splitLayout(
                   downloadButton("download_filtered_tsv", "Filtered (.tsv)"),
                   downloadButton("download_filtered_csv", "Filtered (.csv)"),
                   downloadButton("download_filtered_xlsx", "Filtered (.xlsx)")
                 ), 
                 br(),
                 splitLayout(
                   downloadButton("download_metadata_tsv", "All (.tsv)"),
                   downloadButton("download_metadata_csv", "All (.csv)"),
                   downloadButton("download_metadata_xlsx", "All (.xlsx)")
                 )
        ),
        tabPanel("Summary Statistics", 
                 verbatimTextOutput("summary_input_statistics"),
                 splitLayout(
                   downloadButton("download_stats_tsv", ".tsv"),
                   downloadButton("download_stats_csv", ".csv)"),
                   downloadButton("download_stats_xlsx", ".xlsx")
                 )
        ),
        tabPanel("FASTA Preview", 
                 tags$div(
                   style = "max-height: 500px; overflow-y: auto; border: 1px solid #ddd; 
                            padding: 10px; background-color: #fdfdfd; font-family: monospace; 
                            white-space: pre;",
                   verbatimTextOutput("fasta_preview")
                 )
        ),
        tabPanel("Plots",
                 fluidRow(
                   column(6, uiOutput("plot_axis_selectors_x")),
                   column(6, uiOutput("plot_axis_selectors_stack"))
                 ), 
                 hr(),
                 plotlyOutput("sequence_plots", height = "600px") %>% withSpinner(color = "#0dc5c1"), 
                 br(),
                 actionButton("open_plot_window", "Open Plot in New Window", 
                              icon = icon("external-link-alt"))
        ),
        tabPanel("Mature Peptide FASTA Generator",
                 h3("Generate Mature Peptide FASTA from ConoPrec Output"),
                 p("This section allows you to upload CSV file(s) generated by tools like ConoPrec, 
                   which predict mature conotoxin peptide regions from precursor sequences."),
                 p("The application will parse these CSVs, extract the mature peptide sequences, 
                   and modify the original FASTA headers (replacing 'precursor' with 'mature' in 
                   the description) to generate new FASTA files containing only the mature peptides."),
                 uiOutput("conoprec_upload_status"),
                 DT::DTOutput("mature_peptide_table_preview")
        )
      )
    )
  )
)

# -------------------------
# SERVER DEFINITION
# -------------------------
server <- function(input, output, session) {
  fasta_data <- reactiveVal(NULL)
  processed_conoprec_data <- reactiveVal(list())
  
  process_precursor_fasta_file <- function(filepath, filename_for_notification = NULL) {
    processed_filename <- if (!is.null(filename_for_notification)) filename_for_notification else basename(filepath)
    id <- showNotification(paste("Parsing FASTA file:", processed_filename), duration = NULL, 
                           closeButton = FALSE, type = "message")
    on.exit(removeNotification(id), add = TRUE)
    tryCatch({
      df <- parse_fasta(filepath)
      fasta_data(df)
      removeNotification(id)
      showNotification(paste("Successfully parsed", nrow(df), "sequences from", processed_filename), 
                       duration = 5, type = "message")
      session$sendCustomMessage("selectFilteredRows", list(deselect = TRUE))
    }, error = function(e) { 
      removeNotification(id)
      showNotification(paste("Error parsing FASTA", processed_filename, ":", e$message), 
                       type = "error", duration = 15)
      fasta_data(NULL) 
    })
  }
  
  observeEvent(input$fasta_file, { 
    req(input$fasta_file)
    process_precursor_fasta_file(input$fasta_file$datapath[1], input$fasta_file$name[1]) 
  })
  
  observeEvent(input$load_example_precursor_fasta, {
    if (file.exists(EXAMPLE_CONOSERVER_FASTA_PATH)) { 
      process_precursor_fasta_file(EXAMPLE_CONOSERVER_FASTA_PATH, basename(EXAMPLE_CONOSERVER_FASTA_PATH))
    } else { 
      showNotification(paste("Example FASTA not found:", EXAMPLE_CONOSERVER_FASTA_PATH, 
                             "(Ensure it's in 'data/' or update path. CWD:", getwd(), ")"), 
                       type = "error", duration = 10) 
    }
  })
  
  process_conoprec_csv_files_logic <- function(files_df) {
    showNotification("Processing ConoPrec CSV file(s)...", duration = NULL, closeButton = FALSE, 
                     id = "csv_proc_id", type = "message")
    on.exit(removeNotification("csv_proc_id"), add = TRUE)
    all_dfs <- list()
    error_occurred <- FALSE
    processed_count <- 0
    total_mature_peptides_found <- 0
    for (i in 1:nrow(files_df)) {
      file_info_name <- files_df[i, "name"]
      file_info_datapath <- files_df[i, "datapath"]
      tryCatch({
        df <- read_csv(file_info_datapath, show_col_types = FALSE)
        req_cols <- c("Name", "Mature")
        if (!all(req_cols %in% names(df))) { 
          stop(paste("CSV '", file_info_name, "' missing required columns (Name, Mature).")) 
        }
        df_mature_current_file <- df %>% filter(!is.na(Mature) & Mature != "" & str_trim(Mature) != "")
        if (nrow(df_mature_current_file) == 0) { 
          showNotification(paste("No valid (non-empty) mature peptides in", file_info_name), 
                           type = "warning", duration = 7)
          next 
        }
        df_mature_processed <- df_mature_current_file %>% mutate(
          NewName = sapply(Name, function(name_str) { 
            parts <- str_split(name_str, "\\|", n = 3)
            if (length(parts[[1]]) >= 2) { 
              parts[[1]][2] <- str_replace_all(parts[[1]][2], regex("precursor", ignore_case = TRUE), "mature")
              return(paste(parts[[1]], collapse = "|")) 
            } else { 
              return(name_str) 
            } 
          }, USE.NAMES = FALSE),
          OriginalFile = file_info_name 
        ) %>% select(OriginalFile, NewName, MatureSequence = Mature)
        all_dfs[[file_info_name]] <- df_mature_processed
        processed_count <- processed_count + 1
        total_mature_peptides_found <- total_mature_peptides_found + nrow(df_mature_processed)
      }, error = function(e) { 
        showNotification(paste("Error processing CSV", file_info_name, ":", e$message), 
                         type = "error", duration = 10)
        error_occurred <<- TRUE 
      })
    }
    removeNotification("csv_proc_id")
    if (!error_occurred && processed_count > 0) { 
      showNotification(paste("Processed", processed_count, "CSV file(s), found", 
                             total_mature_peptides_found, "mature peptides."), 
                       type = "message", duration = 7)
    } else if (processed_count == 0 && !error_occurred) { 
      showNotification("No CSVs processed or no mature peptides found.", type = "warning", duration = 7) 
    }
    processed_conoprec_data(all_dfs)
  }
  
  observeEvent(input$conoprec_csv_files, { 
    req(input$conoprec_csv_files)
    process_conoprec_csv_files_logic(input$conoprec_csv_files) 
  })
  
  observeEvent(input$load_example_conoprec_csv, {
    if (file.exists(EXAMPLE_CONOPREC_CSV_PATH)) {
      example_file_df <- data.frame(name = basename(EXAMPLE_CONOPREC_CSV_PATH), 
                                    size = file.size(EXAMPLE_CONOPREC_CSV_PATH), 
                                    type = "text/csv", 
                                    datapath = EXAMPLE_CONOPREC_CSV_PATH, 
                                    stringsAsFactors = FALSE)
      process_conoprec_csv_files_logic(example_file_df)
    } else { 
      showNotification(paste("Example ConoPrec CSV not found at:", EXAMPLE_CONOPREC_CSV_PATH, 
                             "(Ensure it's in 'data/' or update path. CWD:", getwd(), ")"), 
                       type = "error", duration = 10) 
    }
  })
  
  output$fasta_header_selector <- renderUI({ 
    req(fasta_data())
    df <- fasta_data()
    cols_for_header <- setdiff(names(df), c("Header", "Sequence"))
    if (length(cols_for_header) == 0) return(p("No metadata for header."))
    default_selection <- intersect(c("ConoID", "Name", "Organism", "Superfamily", 
                                     "Pharmacological", "SequenceEvidence"), cols_for_header)
    if (length(default_selection) == 0) default_selection <- head(cols_for_header, min(length(cols_for_header), 5))
    checkboxGroupInput("fasta_header_fields", "Select Fields for FASTA Header:", 
                       choices = cols_for_header, selected = default_selection) 
  })
  
  plot_suitable_fields <- reactive({ 
    req(fasta_data())
    df <- fasta_data()
    potential_fields <- setdiff(names(df), c("Header", "Sequence", "ConoID", "Name"))
    suitable_fields <- Filter(function(col) { 
      if (!col %in% names(df) || all(is.na(df[[col]]))) return(FALSE)
      count <- length(unique(na.omit(df[[col]])))
      count > 1 && count < 100 
    }, potential_fields)
    if (length(suitable_fields) == 0) suitable_fields <- Filter(function(col) { 
      if (!col %in% names(df) || all(is.na(df[[col]]))) return(FALSE)
      length(unique(na.omit(df[[col]]))) > 0 
    }, potential_fields)
    if (length(suitable_fields) == 0) return(NULL)
    suitable_fields 
  })
  
  output$plot_axis_selectors_x <- renderUI({ 
    suitable_fields <- plot_suitable_fields()
    if (is.null(suitable_fields)) return(p("No columns for X-axis."))
    default_x <- if ("Superfamily" %in% suitable_fields) "Superfamily" else if ("SequenceEvidence" %in% suitable_fields) "SequenceEvidence" else suitable_fields[1]
    selectInput("x_axis_field", "Select X-Axis Field:", choices = suitable_fields, selected = default_x) 
  })
  
  output$plot_axis_selectors_stack <- renderUI({ 
    suitable_fields <- plot_suitable_fields()
    if (is.null(suitable_fields)) return(p("No columns for Stacking."))
    current_x <- input$x_axis_field
    default_x_fallback <- if ("Superfamily" %in% suitable_fields) "Superfamily" else if ("SequenceEvidence" %in% suitable_fields) "SequenceEvidence" else if (length(suitable_fields) > 0) suitable_fields[1] else NULL
    x_to_avoid <- if (!is.null(current_x) && current_x %in% suitable_fields) current_x else default_x_fallback
    default_stack <- if ("Pharmacological" %in% suitable_fields && "Pharmacological" != x_to_avoid) { 
      "Pharmacological" 
    } else if ("SequenceEvidence" %in% suitable_fields && "SequenceEvidence" != x_to_avoid) {
      "SequenceEvidence"
    } else { 
      next_choice <- setdiff(suitable_fields, x_to_avoid)
      if (length(next_choice) > 0) next_choice[1] else suitable_fields[1] 
    }
    if (is.null(default_stack) && length(suitable_fields) > 0) default_stack <- suitable_fields[1]
    selectInput("stack_field", "Select Stacking Field:", choices = suitable_fields, selected = default_stack) 
  })
  
  output$column_selector <- renderUI({ 
    req(fasta_data())
    cols <- colnames(fasta_data())
    standard_visible <- c("ConoID", "Name", "Organism", "Superfamily", "Pharmacological", "ProteinType", "SequenceEvidence")
    default_visible <- intersect(standard_visible, cols)
    if ("Sequence" %in% cols) default_visible <- c(default_visible, "Sequence")
    extra_cols <- setdiff(cols, c(standard_visible, "Header", "Sequence"))
    default_visible <- unique(c(default_visible, head(extra_cols, 2)))
    if (length(default_visible) == 0) default_visible <- head(cols, min(length(cols), 7))
    pickerInput("visible_columns", "Select Columns to Show:", choices = cols, selected = default_visible, 
                multiple = TRUE, options = list(`actions-box` = TRUE, `selected-text-format` = "count > 3", 
                                                `count-selected-text` = "{0}/{1}")) 
  })
  
  filter_choices <- reactive({ 
    req(fasta_data())
    df <- fasta_data()
    cols_to_filter <- intersect(c("Organism", "ProteinType", "ConopeptideClass", "Superfamily", 
                                  "Pharmacological", "Framework", "SequenceEvidence"), names(df))
    choices_list <- list()
    for (col in cols_to_filter) { 
      unique_vals <- na.omit(unique(df[[col]]))
      unique_vals <- as.character(sort(unique_vals))
      if (length(unique_vals) > 1 || (length(unique_vals) == 1 && unique_vals[1] != "Unclassified")) { 
        if (length(unique_vals) <= 500) { 
          choices_list[[col]] <- unique_vals 
        } 
      } 
    }
    choices_list 
  })
  
  output$column_filters <- renderUI({ 
    req(fasta_data())
    choices_list <- filter_choices()
    if (length(choices_list) == 0) return(p("No applicable filters."))
    filter_elements <- lapply(names(choices_list), function(col) { 
      unique_vals <- choices_list[[col]]
      col_sanitized <- make.names(col)
      tagList(
        pickerInput(inputId = paste0("filter_", col_sanitized), label = paste("Filter", col), 
                    choices = unique_vals, selected = unique_vals, multiple = TRUE, 
                    options = list(`actions-box` = TRUE, `selected-text-format` = "count > 3", 
                                   `count-selected-text` = "{0} selected", 
                                   `none-selected-text` = "Select values...", 
                                   `live-search` = length(unique_vals) > 10, 
                                   `virtual-scroll` = if (length(unique_vals) > 50) 20 else FALSE)),
        checkboxInput(inputId = paste0("exclude_", col_sanitized), label = paste("Exclude selected", col), 
                      value = FALSE)
      ) 
    })
    do.call(tagList, filter_elements) 
  })
  
  filtered_data_raw <- reactive({ 
    req(fasta_data())
    df <- fasta_data()
    choices_list <- filter_choices()
    cols_with_filters <- names(choices_list)
    for (col in cols_with_filters) { 
      col_sanitized <- make.names(col)
      filter_input_id <- paste0("filter_", col_sanitized)
      exclude_input_id <- paste0("exclude_", col_sanitized)
      if (filter_input_id %in% names(input)) { 
        filter_values <- input[[filter_input_id]]
        exclude_logic <- isTRUE(input[[exclude_input_id]])
        if (!is.null(filter_values)) { 
          current_unique_vals <- choices_list[[col]]
          is_filter_active <- !( (!exclude_logic && length(filter_values) == length(current_unique_vals)) || 
                                   (exclude_logic && length(filter_values) == 0) )
          if (is_filter_active) { 
            if (exclude_logic) { 
              df <- df[!(df[[col]] %in% filter_values), , drop = FALSE]
              if ("NA" %in% filter_values) { 
                df <- df[!is.na(df[[col]]), , drop = FALSE] 
              } 
            } else { 
              if ("NA" %in% filter_values) { 
                filter_values_no_na <- setdiff(filter_values, "NA")
                df <- df[(df[[col]] %in% filter_values_no_na) | is.na(df[[col]]), , drop = FALSE] 
              } else { 
                df <- df[df[[col]] %in% filter_values, , drop = FALSE] 
              }
              if (length(filter_values) == 0) { 
                df <- df[0, , drop = FALSE] 
              } 
            } 
          } 
        } 
      }
      if (nrow(df) == 0) break 
    }
    df 
  })
  
  filtered_data <- debounce(filtered_data_raw, millis = 400)
  
  dt_initComplete_js <- "function(settings,json){
    var table=this.api();
    var shinyInputDebounce=null;
    var updateShinySelection=function(){
      var selectedRows=table.rows({selected:true}).indexes().toArray();
      Shiny.setInputValue('seq_table_rows_selected',selectedRows,{priority:'event'});
    };
    table.on('draw.dt select.dt deselect.dt',function(){
      if(shinyInputDebounce!==null){clearTimeout(shinyInputDebounce);}
      shinyInputDebounce=setTimeout(updateShinySelection,100);
    });
    setTimeout(updateShinySelection,200);
  }"
  
  output$seq_table <- DT::renderDT({ 
    req(filtered_data(), input$visible_columns)
    df_display_full <- filtered_data()
    visible_cols_exist <- intersect(input$visible_columns, names(df_display_full))
    if (length(visible_cols_exist) == 0) { shiny::validate("No valid columns selected.") }
    df_display <- df_display_full[, visible_cols_exist, drop = FALSE]
    col_display_names <- setNames(names(df_display), names(df_display))
    datatable(df_display, 
              colnames = col_display_names, 
              options = list(pageLength = 10, 
                             lengthMenu = list(c(10, 100, 1000, -1), c('10', '100', '1000', 'All')), 
                             searching = TRUE, 
                             info = TRUE, 
                             lengthChange = TRUE, 
                             autoWidth = TRUE, 
                             scrollX = TRUE, 
                             select = list(style = 'multi', selector = 'td'), 
                             initComplete = DT::JS(dt_initComplete_js)), 
              filter = "none", 
              selection = 'none', 
              rownames = FALSE, 
              extensions = c('Select')) 
  }, server = FALSE)
  
  observeEvent(input$clear_filters, { 
    choices_list <- filter_choices()
    for (col in names(choices_list)) { 
      col_sanitized <- make.names(col)
      updatePickerInput(session, paste0("filter_", col_sanitized), selected = character(0))
      updateCheckboxInput(session, paste0("exclude_", col_sanitized), value = FALSE) 
    }
    session$sendCustomMessage("selectFilteredRows", list(deselect = TRUE)) 
  })
  
  observeEvent(input$reset_filters, { 
    choices_list <- filter_choices()
    for (col in names(choices_list)) { 
      col_sanitized <- make.names(col)
      updatePickerInput(session, paste0("filter_", col_sanitized), selected = choices_list[[col]])
      updateCheckboxInput(session, paste0("exclude_", col_sanitized), value = FALSE) 
    }
    session$sendCustomMessage("selectFilteredRows", list(deselect = TRUE)) 
  })
  
  observeEvent(input$select_filtered, { 
    session$sendCustomMessage("selectFilteredRows", list(deselect = FALSE)) 
  })
  
  observeEvent(input$deselect_filtered, { 
    session$sendCustomMessage("selectFilteredRows", list(deselect = TRUE)) 
  })
  
  summary_stats_reactive <- reactive({ 
    req(fasta_data())
    df <- fasta_data()
    output_lines <- c(paste("Total input sequences:", nrow(df), "\n"))
    stats_list <- list()
    exclude_cols <- c("Header", "Sequence", "ConoID", "Name")
    cols_to_summarize <- setdiff(names(df), c(exclude_cols, names(df)[startsWith(names(df), "ExtraField_")]))
    for (col in cols_to_summarize) {
      if (col %in% names(df)) {
        counts <- table(df[[col]], useNA = "ifany")
        if (length(counts) == 0) next
        sorted_counts <- sort(counts, decreasing = TRUE)
        col_header <- paste0(col, ": ", length(sorted_counts), " unique")
        output_lines <- c(output_lines, col_header)
        count_details <- sapply(seq_along(sorted_counts), function(i) {
          entry_name <- names(sorted_counts)[i]
          if (is.na(entry_name)) entry_name <- "NA/Missing"
          paste(" -", entry_name, ":", sorted_counts[i])
        })
        output_lines <- c(output_lines, count_details, "")
        stats_df <- data.frame(Category = names(sorted_counts), Count = as.integer(sorted_counts), 
                               stringsAsFactors = FALSE)
        stats_df$Category[is.na(stats_df$Category)] <- "NA/Missing"
        stats_list[[col]] <- stats_df
      }
    }
    list(text = paste(output_lines, collapse = "\n"), data = stats_list) 
  })
  
  output$summary_input_statistics <- renderPrint({ 
    result <- summary_stats_reactive()
    cat(result$text) 
  })
  
  output$fasta_stats_message <- renderPrint({ 
    req(fasta_data(), filtered_data())
    df_filtered <- filtered_data()
    df_all <- fasta_data()
    selected_indices <- input$seq_table_rows_selected
    total_filtered_count <- nrow(df_filtered)
    total_input_count <- nrow(df_all)
    cat(paste("Filtered sequences:", total_filtered_count, "/", total_input_count, "\n"))
    selected_df <- NULL
    context_msg <- ""
    is_selection_valid <- FALSE
    if (!is.null(selected_indices) && length(selected_indices) > 0) {
      max_index <- max(selected_indices)
      if (max_index < total_filtered_count) {
        selected_df <- df_filtered[selected_indices + 1, , drop = FALSE]
        cat(paste("Selected sequences:", nrow(selected_df), "/", total_filtered_count, "\n\n"))
        context_msg <- "Stats for selected (vs filtered):"
        is_selection_valid <- TRUE
      } else {
        cat("Selection invalid. Stats for filtered.\n")
        selected_df <- df_filtered
        context_msg <- "Stats for filtered data:"
      }
    } else {
      cat("No sequences selected.\n\n")
      selected_df <- df_filtered
      context_msg <- "Stats for filtered data:"
    }
    if (is.null(selected_df)) { 
      cat("Error determining data for stats.\n")
      return() 
    }
    if (nrow(selected_df) == 0 && is_selection_valid) { 
      cat("Selected data is empty.\n")
      return() 
    }
    if (nrow(selected_df) == 0 && !is_selection_valid && total_filtered_count == 0) { 
      cat("Filtered data is empty.\n")
      return() 
    }
    cat(context_msg, "\n")
    generate_stats_string <- function(data_subset, data_total, column_name) {
      if (!column_name %in% names(data_subset) || !column_name %in% names(data_total)) return(NULL)
      counts_subset <- table(data_subset[[column_name]], useNA = "ifany")
      counts_total <- table(data_total[[column_name]], useNA = "ifany")
      if (length(counts_subset) == 0) return(paste0(column_name, ": 0 unique\n"))
      sorted_counts_subset <- sort(counts_subset, decreasing = TRUE)
      lines <- c(paste0(column_name, " (", length(sorted_counts_subset), " unique):"))
      subset_names <- names(sorted_counts_subset)
      details <- sapply(seq_along(sorted_counts_subset), function(i) {
        name <- subset_names[i]
        display_name <- if (is.na(name)) "NA/Missing" else name
        subset_count <- sorted_counts_subset[[i]]
        total_count_val <- counts_total[name]
        if (is.na(name)) { 
          total_count_val <- counts_total[is.na(names(counts_total))]
          if (length(total_count_val) == 0) total_count_val <- 0
        } else if (is.na(total_count_val)) { 
          total_count_val <- 0 
        }
        if (length(total_count_val) != 1) total_count_val <- 0
        paste0(" - ", display_name, ": ", subset_count, " / ", total_count_val)
      })
      c(lines, details, "")
    }
    for (col_stat in c("Superfamily", "Pharmacological", "ProteinType", "SequenceEvidence", "Organism")) {
      if (col_stat %in% names(selected_df)) {
        stats_output <- generate_stats_string(selected_df, df_filtered, col_stat)
        if (!is.null(stats_output)) cat(paste(stats_output, collapse = "\n"))
      }
    }
  })
  
  output$fasta_preview <- renderPrint({ 
    req(filtered_data(), input$fasta_header_fields)
    df_filtered <- filtered_data()
    selected_indices <- input$seq_table_rows_selected
    selected_fields <- input$fasta_header_fields
    sep <- input$header_sep
    if (!"Sequence" %in% names(df_filtered)) { 
      cat("Error: 'Sequence' column missing.")
      return() 
    }
    valid_fields <- intersect(selected_fields, names(df_filtered))
    if (length(valid_fields) == 0) { 
      cat("Error: No valid FASTA header fields.")
      return() 
    }
    if (nrow(df_filtered) == 0) { 
      cat("No sequences match filters.")
      return() 
    }
    max_preview <- 50
    preview_df <- NULL
    context_msg <- ""
    if (!is.null(selected_indices) && length(selected_indices) > 0) {
      max_index <- max(selected_indices)
      if (max_index >= nrow(df_filtered)) {
        cat("Selection invalid. Previewing filtered.\n")
        preview_df <- head(df_filtered, max_preview)
        context_msg <- paste("Showing first", nrow(preview_df), "of", nrow(df_filtered), "filtered:\n")
      } else {
        selected_df <- df_filtered[selected_indices + 1, , drop = FALSE]
        preview_df <- head(selected_df, max_preview)
        context_msg <- paste("Showing", nrow(preview_df), "of", nrow(selected_df), "selected (out of", 
                             nrow(df_filtered), "filtered):\n")
      }
    } else {
      preview_df <- head(df_filtered, max_preview)
      context_msg <- paste("No selection. Previewing first", nrow(preview_df), "of", 
                           nrow(df_filtered), "filtered:\n")
    }
    if (!is.null(preview_df) && nrow(preview_df) > 0) {
      fasta_str <- create_fasta_content(preview_df, valid_fields, sep)
      cat(context_msg)
      cat(fasta_str)
    } else { 
      cat("No data for preview.") 
    }
  })
  
  # --- Plotting Logic ---
  plot_object <- reactive({
    req(filtered_data(), input$x_axis_field, input$stack_field)
    df <- filtered_data()
    x_field <- input$x_axis_field
    stack_field <- input$stack_field
    if (!(x_field %in% names(df)) || !(stack_field %in% names(df))) { 
      return(plotly::plot_ly() %>% layout(title = "Error: Plot field(s) not found.")) 
    }
    if (x_field == stack_field && length(unique(df[[x_field]])) > 1) { 
      return(plotly::plot_ly() %>% layout(title = "Error: X & Stack fields must be different.")) 
    }
    if (nrow(df) == 0) { 
      return(plotly::plot_ly() %>% layout(title = "No data for plotting.")) 
    }
    df_plot <- df %>% 
      mutate(across(all_of(c(x_field, stack_field)), 
                    ~ if_else(is.na(.) | trimws(as.character(.)) == "", "Unclassified", as.character(.)))) %>% 
      count(!!sym(x_field), !!sym(stack_field), name = "Count") %>% 
      ungroup()
    if (nrow(df_plot) == 0) { 
      return(plotly::plot_ly() %>% layout(title = "No category combinations.")) 
    }
    x_order <- df_plot %>% 
      group_by(!!sym(x_field)) %>% 
      summarise(Total = sum(Count)) %>% 
      arrange(desc(Total)) %>% 
      pull(!!sym(x_field))
    df_plot[[x_field]] <- factor(df_plot[[x_field]], levels = x_order)
    stack_order <- sort(unique(df_plot[[stack_field]]))
    df_plot[[stack_field]] <- factor(df_plot[[stack_field]], levels = stack_order)
    max_count <- if (nrow(df_plot) > 0) max(df_plot$Count, na.rm = TRUE) else 0
    min_count <- if (nrow(df_plot) > 0) min(df_plot$Count[df_plot$Count > 0], na.rm = TRUE) else 1
    use_log_scale <- max_count > 0 && (max_count / max(1, min_count)) > 50
    
    p <- plot_ly(
      data = df_plot, 
      x = ~get(x_field), 
      y = ~Count, 
      type = 'bar',
      color = ~get(stack_field), 
      colors = "viridis",
      text = ~paste("<b>", x_field, ": ", get(x_field), "</b><br>", 
                    stack_field, ": ", get(stack_field), "<br>", "Count: ", Count),
      hoverinfo = "text"
    ) %>%
      layout(
        title = list(text = paste("Sequence Counts by", x_field, "and", stack_field), x = 0.05),
        xaxis = list(title = x_field, type = 'category', tickangle = -45, automargin = TRUE),
        yaxis = list(title = paste("Number of Sequences", if (use_log_scale) "(Log Scale)" else ""), 
                     type = if (use_log_scale) "log" else "linear", autorange = TRUE, tickformat = ',d'),
        barmode = "stack",
        legend = list(title = list(text = stack_field), orientation = "v", traceorder = "normal", 
                      itemsizing = 'constant', yanchor = "top", y = 1, xanchor = "left", x = 1.02),
        margin = list(l = 60, r = 180, b = 130, t = 60, pad = 5),
        hovermode = "closest"
      ) %>%
      config(
        toImageButtonOptions = list(
          format = "png",
          filename = paste0("ConoToxin_Plot_", Sys.Date()),
          width = 1600,
          height = 900,
          scale = 2.5
        ),
        displaylogo = FALSE,
        modeBarButtons = list(list("toImage", "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d", 
                                   "pan2d", "zoom2d"))
      )
    return(p)
  })
  
  output$sequence_plots <- renderPlotly({ plot_object() })
  
  observeEvent(input$open_plot_window, { 
    showModal(modalDialog(
      title = "Plot (Expanded View). Use plot toolbar (camera icon hover) to download as higher-resolution PNG.",
      plotlyOutput("modalPlotOutput", height = "600px") %>% withSpinner(color = "#0dc5c1"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Dismiss")
    )) 
  })
  
  output$modalPlotOutput <- renderPlotly({ plot_object() })
  
  # --- Download Handlers ---
  get_data_for_download <- function(filtered_df, selected_idx) {
    if (!is.null(selected_idx) && length(selected_idx) > 0) {
      max_index <- if (length(selected_idx) > 0) max(selected_idx) else -1
      if (max_index >= nrow(filtered_df)) {
        showNotification("Selection invalid. Downloading filtered.", type = "warning", duration = 7)
        return(filtered_df)
      }
      return(filtered_df[selected_idx + 1, , drop = FALSE])
    } else {
      return(filtered_df)
    }
  }
  
  create_fasta_content <- function(df_subset, header_fields, separator) {
    if (!"Sequence" %in% names(df_subset) || nrow(df_subset) == 0 || is.null(header_fields) || 
        length(header_fields) == 0) return("")
    valid_fields <- intersect(header_fields, names(df_subset))
    if (length(valid_fields) == 0) return("")
    header_data <- df_subset[, valid_fields, drop = FALSE]
    sequences <- df_subset$Sequence
    fasta_headers <- apply(header_data, 1, function(row) {
      row_values <- ifelse(is.na(row), "NA", as.character(row))
      paste0(">", paste(row_values, collapse = separator))
    })
    paste(paste(fasta_headers, sequences, sep = "\n"), collapse = "\n")
  }
  
  output$download_filtered_fasta <- downloadHandler(
    filename = function() paste0("ConoServerParser_precursors_", Sys.Date(), ".fasta"),
    content = function(file) {
      req(filtered_data(), input$fasta_header_fields, input$header_sep)
      df_to_download <- get_data_for_download(filtered_data(), input$seq_table_rows_selected)
      if (nrow(df_to_download) == 0 || !"Sequence" %in% names(df_to_download) || 
          length(intersect(input$fasta_header_fields, names(df_to_download))) == 0) {
        showNotification("Cannot generate FASTA: No data/Sequence/valid headers.", 
                         type = "error", duration = 8)
        writeLines("", file)
        return()
      }
      fasta_str <- create_fasta_content(df_to_download, input$fasta_header_fields, input$header_sep)
      if (nzchar(fasta_str)) { 
        writeLines(fasta_str, file) 
      } else { 
        writeLines("", file) 
      }
    }
  )
  
  # MODIFICATION: Changed default prefix to "" (empty string) to remove "precursors_" from zipped FASTA files.
  # Also updated file_rel_path logic to handle empty prefix correctly.
  create_zip_archive <- function(file, df_source, group_by_col, header_fields, separator, prefix = "") {
    if (!group_by_col %in% names(df_source)) stop(paste("Grouping column '", group_by_col, "' not found.", sep = ""))
    if (nrow(df_source) == 0) { 
      showNotification("No data for zip.", type = "warning")
      writeLines("", file) # Create an empty file for the download handler
      return() 
    }
    
    # Ensure the grouping column is character and handles NA/empty strings
    df_source[[group_by_col]] <- ifelse(is.na(df_source[[group_by_col]]), "Unclassified", 
                                        as.character(df_source[[group_by_col]]))
    df_source[[group_by_col]] <- ifelse(trimws(df_source[[group_by_col]]) == "", "Unclassified", 
                                        df_source[[group_by_col]])
    
    split_data <- split(df_source, df_source[[group_by_col]])
    
    temp_dir <- tempdir()
    fasta_files_full <- character() # To store full paths of created FASTA files for later cleanup
    # files_to_zip_relative <- character() # Not strictly needed if using basename for zipping
    
    # Ensure temp_dir is cleaned up even if errors occur
    on.exit({ unlink(fasta_files_full, recursive = TRUE, force = TRUE) }, add = TRUE) 
    
    for (group_name in names(split_data)) {
      df_subset <- split_data[[group_name]]
      fasta_content <- create_fasta_content(df_subset, header_fields, separator)
      
      if (nzchar(fasta_content)) {
        # Sanitize group_name for use in filename
        safe_group_name <- gsub("[^A-Za-z0-9_.-]+", "_", group_name)
        safe_group_name <- ifelse(safe_group_name == "", "Unclassified", safe_group_name) # Handle empty after sanitization
        
        # MODIFICATION: Construct filename based on whether prefix is empty or not
        if (nzchar(prefix)) {
          file_rel_path <- paste0(prefix, "_", safe_group_name, ".fasta")
        } else {
          file_rel_path <- paste0(safe_group_name, ".fasta")
        }
        
        file_full_path <- file.path(temp_dir, file_rel_path)
        
        tryCatch({
          writeLines(fasta_content, file_full_path)
          fasta_files_full <- c(fasta_files_full, file_full_path)
          # files_to_zip_relative <- c(files_to_zip_relative, file_rel_path) # Storing relative path if zipping from a different CWD
        }, error = function(e) { 
          warning(paste("Failed to write FASTA file for group:", group_name, "-", e$message)) 
        })
      }
    }
    
    if (length(fasta_files_full) == 0) { 
      showNotification("No FASTA files were generated to include in the zip.", type = "warning")
      writeLines("", file) # Create an empty file for the download handler
      return() 
    }
    
    # Zipping files
    # It's often safer to change directory to the temp_dir to avoid issues with full paths in zip archives
    # especially if some zip utilities/platforms don't handle them well or create unwanted directory structures.
    wd_orig <- getwd()
    setwd(temp_dir)
    # Ensure working directory is restored and files are cleaned up
    on.exit({ setwd(wd_orig); unlink(fasta_files_full, recursive = TRUE, force = TRUE) }, add = TRUE, after = FALSE) 
    
    tryCatch({
      # Using basename() ensures that only the file names are used, not their full paths from temp_dir
      zip::zip(zipfile = file, files = basename(fasta_files_full), mode = "cherry-pick") 
    }, error = function(e) { 
      # setwd(wd_orig) # Already handled by on.exit
      stop(paste("Failed to create zip archive:", e$message)) 
    })
  }
  
  output$download_zip_by_superfamily <- downloadHandler(
    filename = function() paste0("split_by_Gene_Superfamily_", Sys.Date(), ".zip"),
    content = function(file) {
      req(filtered_data(), input$fasta_header_fields, input$header_sep)
      df_process <- get_data_for_download(filtered_data(), input$seq_table_rows_selected)
      # MODIFICATION: Removed prefix argument, will use the function's new default (prefix = "")
      create_zip_archive(file, df_process, "Superfamily", input$fasta_header_fields, input$header_sep)
    },
    contentType = "application/zip"
  )
  
  output$download_zip_by_pharmacological <- downloadHandler(
    filename = function() paste0("split_by_Pharmacological_Family_", Sys.Date(), ".zip"), # MODIFIED: Zip filename for consistency
    content = function(file) {
      req(filtered_data(), input$fasta_header_fields, input$header_sep)
      df_process <- get_data_for_download(filtered_data(), input$seq_table_rows_selected)
      # MODIFICATION: Removed prefix argument, will use the function's new default (prefix = "")
      create_zip_archive(file, df_process, "Pharmacological", input$fasta_header_fields, input$header_sep)
    },
    contentType = "application/zip"
  )
  
  output$download_filtered_tsv <- downloadHandler(
    filename = function() paste0("filtered_metadata_precursors_", Sys.Date(), ".tsv"),
    content = function(file) { 
      req(filtered_data())
      write_tsv(filtered_data(), file, na = "NA") 
    }
  )
  
  output$download_filtered_csv <- downloadHandler(
    filename = function() paste0("filtered_metadata_precursors_", Sys.Date(), ".csv"),
    content = function(file) { 
      req(filtered_data())
      write_csv(filtered_data(), file, na = "NA") 
    }
  )
  
  output$download_filtered_xlsx <- downloadHandler(
    filename = function() paste0("filtered_metadata_precursors_", Sys.Date(), ".xlsx"),
    content = function(file) { 
      req(filtered_data())
      df_write <- filtered_data()
      tryCatch(
        writexl::write_xlsx(list(FilteredData = df_write), file),
        error = function(e) { 
          showNotification(paste("XLSX Error:", e$message), type = "error")
          stop(e) 
        }
      ) 
    }
  )
  
  output$download_metadata_tsv <- downloadHandler(
    filename = function() paste0("all_metadata_precursors_", Sys.Date(), ".tsv"),
    content = function(file) { 
      req(fasta_data())
      write_tsv(fasta_data(), file, na = "NA") 
    }
  )
  
  output$download_metadata_csv <- downloadHandler(
    filename = function() paste0("all_metadata_precursors_", Sys.Date(), ".csv"),
    content = function(file) { 
      req(fasta_data())
      write_csv(fasta_data(), file, na = "NA") 
    }
  )
  
  output$download_metadata_xlsx <- downloadHandler(
    filename = function() paste0("all_metadata_precursors_", Sys.Date(), ".xlsx"),
    content = function(file) { 
      req(fasta_data())
      df_write <- fasta_data()
      tryCatch(
        writexl::write_xlsx(list(AllMetadata = df_write), file),
        error = function(e) { 
          showNotification(paste("XLSX Error:", e$message), type = "error")
          stop(e) 
        }
      ) 
    }
  )
  
  output$download_stats_tsv <- downloadHandler(
    filename = function() paste0("summary_stats_precursors_", Sys.Date(), ".tsv"),
    content = function(file) { 
      stats_list <- summary_stats_reactive()$data
      req(length(stats_list) > 0)
      bind_rows(stats_list, .id = "SummaryColumn") %>% write_tsv(file, na = "NA") 
    }
  )
  
  output$download_stats_csv <- downloadHandler(
    filename = function() paste0("summary_stats_precursors_", Sys.Date(), ".csv"),
    content = function(file) { 
      stats_list <- summary_stats_reactive()$data
      req(length(stats_list) > 0)
      bind_rows(stats_list, .id = "SummaryColumn") %>% write_csv(file, na = "NA") 
    }
  )
  
  output$download_stats_xlsx <- downloadHandler(
    filename = function() paste0("summary_stats_precursors_", Sys.Date(), ".xlsx"),
    content = function(file) { 
      stats_list <- summary_stats_reactive()$data
      req(length(stats_list) > 0)
      tryCatch(
        writexl::write_xlsx(stats_list, file),
        error = function(e) { 
          showNotification(paste("XLSX Error:", e$message), type = "error")
          stop(e) 
        }
      ) 
    }
  )
  
  output$mature_peptide_summary_ui <- renderUI({ 
    data_list <- processed_conoprec_data()
    if (length(data_list) == 0) { 
      return(p("No ConoPrec CSVs processed.")) 
    }
    summary_tags <- tagList()
    for (filename in names(data_list)) {
      num_mature <- nrow(data_list[[filename]])
      summary_tags <- tagAppendChild(summary_tags, p(strong(filename), ": ", num_mature, " mature peptides found."))
    }
    summary_tags 
  })
  
  output$conoprec_upload_status <- renderUI({ 
    req(input$conoprec_csv_files)
    data_list <- processed_conoprec_data()
    if (length(data_list) > 0) {
      total_peptides <- sum(sapply(data_list, nrow))
      p(strong(paste("Total mature peptides processed from", length(data_list), "file(s):", total_peptides)))
    } else if (!is.null(input$conoprec_csv_files) && length(data_list) == 0) {
      p("Processing or no mature peptides found...")
    } else {
      p("Upload ConoPrec CSVs.")
    }
  })
  
  output$mature_peptide_table_preview <- DT::renderDT({ 
    data_list <- processed_conoprec_data()
    if (length(data_list) == 0) {
      return(datatable(data.frame(Message = "No mature peptide data."), rownames = FALSE, 
                       options = list(dom = 't')))
    }
    combined_df <- bind_rows(data_list, .id = "SourceFile")
    datatable(combined_df, 
              caption = "Preview of Processed Mature Peptides", 
              options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE), 
              rownames = FALSE) 
  })
  
  output$download_mature_fasta_zip <- downloadHandler(
    filename = function() { paste0("ConoPrec_mature_peptides_", Sys.Date(), ".zip") },
    content = function(file) {
      data_list <- processed_conoprec_data()
      if (length(data_list) == 0) {
        showNotification("No processed mature peptide data.", type = "warning", duration = 5)
        zz <- zip_list(character(0)) # Create an empty zip archive
        writeBin(zz$zip_archive, file)
        return()
      }
      
      temp_dir <- tempdir()
      fasta_files_to_zip_paths <- character()
      # Ensure created files are cleaned up
      on.exit(unlink(list.files(temp_dir, full.names = TRUE, pattern = "\\.fasta$"), recursive = TRUE, force = TRUE), add = TRUE) 
      
      for (original_csv_name in names(data_list)) {
        df_mature <- data_list[[original_csv_name]]
        if (nrow(df_mature) > 0) {
          fasta_content_mature <- character(nrow(df_mature) * 2)
          for (j in 1:nrow(df_mature)) {
            fasta_content_mature[(j * 2) - 1] <- paste0(">", df_mature$NewName[j])
            fasta_content_mature[j * 2] <- df_mature$MatureSequence[j]
          }
          fasta_str_mature <- paste(fasta_content_mature, collapse = "\n")
          
          base_name <- tools::file_path_sans_ext(original_csv_name)
          mature_fasta_filename <- paste0(base_name, "_mature.fasta")
          mature_fasta_filepath <- file.path(temp_dir, mature_fasta_filename)
          
          tryCatch({
            writeLines(fasta_str_mature, mature_fasta_filepath)
            fasta_files_to_zip_paths <- c(fasta_files_to_zip_paths, mature_fasta_filepath)
          }, error = function(e) { 
            warning(paste("Could not write FASTA for mature peptides from:", original_csv_name, "-", e$message)) 
          })
        }
      }
      
      if (length(fasta_files_to_zip_paths) == 0) {
        showNotification("No FASTA files generated for mature peptides.", type = "warning", duration = 5)
        zz <- zip_list(character(0)) # Create an empty zip archive
        writeBin(zz$zip_archive, file)
        return()
      }
      
      wd_orig <- getwd()
      setwd(temp_dir)
      # Ensure working directory is restored and files are cleaned up
      on.exit({ setwd(wd_orig); unlink(fasta_files_to_zip_paths, recursive = TRUE, force = TRUE) }, add = TRUE, after = FALSE) 
      
      tryCatch({
        zip::zip(zipfile = file, files = basename(fasta_files_to_zip_paths), mode = "cherry-pick")
      }, error = function(e) { 
        # setwd(wd_orig) # Already handled by on.exit
        showNotification(paste("Error creating zip for mature peptides:", e$message), type = "error", duration = 10)
        stop(e) # Propagate error to stop download
      })
    },
    contentType = "application/zip"
  )
}

# -------------------------
# LAUNCH APPLICATION
# -------------------------
shinyApp(ui = ui, server = server)