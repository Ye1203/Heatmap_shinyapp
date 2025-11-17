source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/heatmap_09d.R")

tab_09d_ui <- fluidPage(
  useShinyjs(),  
  tags$style(HTML("
    #condition_ranklist .rank-list-item {
      background-color: #e6f2ff !important;
      border: 1px solid #cce0ff;
      padding: 6px 10px;
      margin-bottom: 5px;
      border-radius: 6px;
    }
  ")),
  
  # ---------- step1 UI ----------
  div(
    id = "step1_ui",
    fluidRow(
      column(4,
             h4("Sample Selection"),
             textInput(
               "input_folder_09d",
               "Input Folder",
               value = "",
               placeholder = "eg. /../G216_G221/Scripts"
             ),
             actionButton("search_09d", "Search for 09d files"),
             br(),br(),
             downloadButton("download_gene_list", "Download Full Gene list"),
             br(),br(),
             fileInput("edited_gene_list", "Upload Custom Gene List", accept = c(".xlsx", ".xls")),
             conditionalPanel(
               condition = "output.gene_list_uploaded",
               checkboxInput("segex_transfer","transfer segex to 09d")
             ),
             conditionalPanel(
               condition = "output.gene_list_uploaded",
               actionButton("clear_gene_list","Clear Gene List")
             ),
             numericInput("set_all_row", "Set Row cluster number for All", value = 1, min = 1),
             numericInput("set_all_col", "Set Column cluster number for All", value = 1, min = 1),
             actionButton("apply_all_row_col", "Apply to All"),
             actionButton("next_09d", "Next Step")
      ),
      column(8,
             h4("Drag folders from left to right to exclude samples from execution"),
             uiOutput("bucket_ui_09d")
      )
    )
  ),
  
  # ---------- step2 UI ----------
  div(
    id = "step2_ui",
    style = "display:none;",  
    fluidRow(
      column(6,
             h4("Additional Settings"),
             textInput("color_input_09d", "Custom Colors (comma separated, #RGB code):", 
                       "red,yellow,#3fdc04"),
             selectInput("genic_region_09d", "Genic Region", 
                         choices = c("ExonCollapsed", "FullGeneBody"), selected = "ExonCollapsed"),
             selectInput("using_package_f_09d", "Using Package", 
                         choices = c("EdgeR", "DESeq2"), selected = "EdgeR"),
             selectInput("scale_function_09d", "Scale Function", 
                         choices = c("maxscale", "zscale", "both"), selected = "both"),
             textInput("abs_log2FC_09d", "Absolute log2 Fold Change (comma-separated)", value = "1,2"),
             textInput("FDR_09d", "FDR Thresholds (comma-separated)", value = "0.05"),
             textInput("tpm_cutoff_09d", "TPM Cutoff (comma-separated)", value = "0,1"),
             br(),
             actionButton("back_09d", "Go Back"),
             actionButton("run_heatmap_09d", "Generate Heatmaps", class = "btn-primary"),
             br(), br(),
             downloadButton("download_heatmap_09d", "Download Results")
      ),
      column(6,
             h4("Order of Conditions"),
             uiOutput("condition_ranklist_ui")  
      )
    )
  )
)

tab_09d_server <- function(input, output, session) {
  
  step_09d <- reactiveVal("step1")
  folder_list <- reactiveVal(character(0))
  selected_folders <- reactiveVal(character(0))
  input_folder_stored <- reactiveVal(NULL)
  folder_settings <- reactiveValues()
  condition_order <- reactiveVal(character(0))
  output_dir_09d <- reactiveVal(tempdir())
  
  special_settings <- reactiveValues(
    row_all_individual = 1,
    col_all_individual = 1,
    row_all_group = 1,
    col_all_group = 1
  )
  
  gene_list_uploaded <- reactiveVal(FALSE)
  
  observeEvent(input$search_09d, {
    req(input$input_folder_09d)
    input_folder_stored(input$input_folder_09d)
    all_dirs <- list.dirs(path = input$input_folder_09d, full.names = FALSE, recursive = FALSE)
    matched <- sort(all_dirs[grepl("^09d", all_dirs)])
    folder_list(matched)
    selected_folders(matched)
  })
  
  observe({
    gene_list_uploaded(!is.null(input$edited_gene_list))
  })
  
  output$gene_list_uploaded <- reactive({
    gene_list_uploaded()
  })
  outputOptions(output, "gene_list_uploaded", suspendWhenHidden = FALSE)
  
  output$condition_ranklist_ui <- renderUI({
    rank_list(
      text = NULL,
      labels = condition_order(),
      input_id = "condition_order_09d",
      options = sortable_options(animation = 150),
      css_id = "condition_ranklist"
    )
  })
  
  observeEvent(input$condition_order_09d, {
    condition_order(input$condition_order_09d)
  })
  
  output$download_gene_list <- downloadHandler(
    filename = function() {
      paste0("Original_gene_list_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(input$input_folder_09d)
      input_folder <- input$input_folder_09d
      
      all_dirs <- list.dirs(path = input_folder, full.names = TRUE, recursive = FALSE)
      dirs_09d <- all_dirs[grepl("/09d", all_dirs)]
      if(length(dirs_09d) == 0) stop("No folder starting with '09d' found in input folder.")
      target_09d <- dirs_09d[1]
      sub_dirs <- list.dirs(target_09d, full.names = TRUE, recursive = FALSE)
      target_sub <- sub_dirs[grepl("FullGeneBody$", sub_dirs)]
      if(length(target_sub) == 0) {
        target_sub <- sub_dirs[grepl("ExonCollapsed$", sub_dirs)]
        if(length(target_sub) == 0) stop("Neither FullGeneBody nor ExonCollapsed folder found.")
      }
      target_sub <- target_sub[1]
      diffexp_files <- list.files(target_sub, pattern = "^DiffExp", full.names = TRUE)
      if(length(diffexp_files) == 0) stop("No DiffExp file found.")
      diffexp_file <- diffexp_files[1]
      original_data <- read_data(diffexp_file)
      if(!("id" %in% colnames(original_data))) stop("The read data does not contain an 'id' column.")
      gene_list <- data.frame(id = original_data$id)
      openxlsx::write.xlsx(gene_list, file = file, colNames = FALSE)
    }
  )
  
  output$download_heatmap_09d <- downloadHandler(
    filename = function() {
      paste0("Heatmap_Result_09d_", Sys.Date(), ".zip")
    },
    content = function(file) {
      zip::zipr(zipfile = file, files = list.files(output_dir_09d(), full.names = TRUE))
    }
  )
  
  output$bucket_ui_09d <- renderUI({
    folders <- folder_list()
    selected <- input$selected_09d %||% selected_folders()
    selected_folders(selected)
    if (length(folders) == 0) return(h6("No folders found."))
    
    bucket_ui <- bucket_list(
      header = NULL,
      group_name = "bucket_09d_group",
      orientation = "horizontal",
      add_rank_list("Included samples", selected, input_id = "selected_09d"),
      add_rank_list("Excluded samples", setdiff(folders, selected), input_id = "available_09d")
    )
    
    param_ui <- lapply(selected, function(folder) {
      folder_path <- file.path(input_folder_stored(), folder)
      condition_1 <- "Not Found"
      condition_2 <- "Not Found"
      
      try({
        cond1_file <- file.path(folder_path, "Condition_1.txt")
        cond2_file <- file.path(folder_path, "Condition_2.txt")
        
        if (file.exists(cond1_file)) {
          dat1 <- read.table(cond1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
          if ("Description" %in% colnames(dat1)) condition_1 <- dat1$Description[1]
        }
        
        if (file.exists(cond2_file)) {
          dat2 <- read.table(cond2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
          if ("Description" %in% colnames(dat2)) condition_2 <- dat2$Description[1]
        }
      }, silent = TRUE)
      
      row_id <- paste0("row_", folder)
      col_id <- paste0("col_", folder)
      if (is.null(folder_settings[[row_id]])) folder_settings[[row_id]] <- 1
      if (is.null(folder_settings[[col_id]])) folder_settings[[col_id]] <- 1
      
      observeEvent(input[[row_id]], {
        folder_settings[[row_id]] <- input[[row_id]]
      }, ignoreInit = TRUE)
      
      observeEvent(input[[col_id]], {
        folder_settings[[col_id]] <- input[[col_id]]
      }, ignoreInit = TRUE)
      
      fluidRow(
        column(3, p(strong("Condition 1:"), condition_1)),
        column(3, p(strong("Condition 2:"), condition_2)),
        column(2, textInput(row_id, "Row", value = folder_settings[[row_id]])),
        column(2, textInput(col_id, "Column", value = folder_settings[[col_id]]))
      )
    })
    
    param_ui <- c(param_ui, list(
      tags$hr(),
      fluidRow(
        column(3, strong("All Individual")),
        column(3, ""),
        column(2, textInput("row_all_individual", "Row", value = special_settings$row_all_individual)),
        column(2, textInput("col_all_individual", "Column", value = special_settings$col_all_individual))
      ),
      fluidRow(
        column(3, strong("All Group")),
        column(3, ""),
        column(2, textInput("row_all_group", "Row", value = special_settings$row_all_group)),
        column(2, textInput("col_all_group", "Column", value = special_settings$col_all_group))
      )
    ))
    
    tagList(bucket_ui, br(), param_ui)
  })
  
  observeEvent(input$apply_all_row_col, {
    selected <- selected_folders()
    for (folder in selected) {
      folder_settings[[paste0("row_", folder)]] <- input$set_all_row
      folder_settings[[paste0("col_", folder)]] <- input$set_all_col
    }
    special_settings$row_all_individual <- input$set_all_row
    special_settings$col_all_individual <- input$set_all_col
    special_settings$row_all_group <- input$set_all_row
    special_settings$col_all_group <- input$set_all_col
  })
  
  observe({
    if (!is.null(input$row_all_individual)) special_settings$row_all_individual <- input$row_all_individual
    if (!is.null(input$col_all_individual)) special_settings$col_all_individual <- input$col_all_individual
    if (!is.null(input$row_all_group)) special_settings$row_all_group <- input$row_all_group
    if (!is.null(input$col_all_group)) special_settings$col_all_group <- input$col_all_group
  })
  
  observeEvent(input$clear_gene_list, {
    shinyjs::reset("edited_gene_list")
    gene_list_uploaded(FALSE)
    
    showModal(modalDialog(
      title = "Success",
      "Gene list has been cleared.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  })
  
  # ------------------ Step switching ------------------
  observeEvent(input$next_09d, {
    input_folder <- input_folder_stored()
    
      selected <- selected_folders()
      conds <- character()
      for (folder in selected) {
        folder_path <- file.path(input_folder_stored(), folder)
        cond1_file <- file.path(folder_path, "Condition_1.txt")
        cond2_file <- file.path(folder_path, "Condition_2.txt")
        cond1 <- NA
        cond2 <- NA
        try({
          if (file.exists(cond1_file)) {
            dat1 <- read.table(cond1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            if ("Description" %in% colnames(dat1)) cond1 <- dat1$Description[1]
          }
          if (file.exists(cond2_file)) {
            dat2 <- read.table(cond2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
            if ("Description" %in% colnames(dat2)) cond2 <- dat2$Description[1]
          }
        }, silent = TRUE)
        conds <- c(conds, cond1, cond2)
      }
      condition_order_vals <- unique(na.omit(conds))
    
    condition_order(condition_order_vals)
    
    shinyjs::hide("step1_ui")
    shinyjs::show("step2_ui")
  })
  
  
  observeEvent(input$back_09d, {
    shinyjs::show("step1_ui")
    shinyjs::hide("step2_ui")
  })
  
  observeEvent(input$run_heatmap_09d, {
    showModal(modalDialog(
      title = "Heatmap Generation",
      "Generating heatmaps... Please wait.",
      footer = NULL
    ))
    
    input_folder_list <- selected_folders()
    
    cluster_row_vec <- unname(c(
      sapply(input_folder_list, function(folder) {
        as.numeric(folder_settings[[paste0("row_", folder)]]) %||% 1
      }),
      as.numeric(special_settings$row_all_individual),
      as.numeric(special_settings$row_all_group)
    ))
    
    cluster_col_vec <- unname(c(
      sapply(input_folder_list, function(folder) {
        as.numeric(folder_settings[[paste0("col_", folder)]]) %||% 1
      }),
      as.numeric(special_settings$col_all_individual),
      as.numeric(special_settings$col_all_group)
    ))
    
    out_dir <- file.path(tempdir(), paste0("heatmap_09d_", as.integer(Sys.time())))
    dir.create(out_dir)
    output_dir_09d(out_dir)
    
    edited_gene_path <- if (gene_list_uploaded()) input$edited_gene_list$datapath else NULL
    tryCatch({
      do_heatmap_09d(
        input_folder = input_folder_stored(),
        output_folder = out_dir,
        input_folder_list = input_folder_list,
        segex_transfer = input$segex_transfer,
        edited_gene_list = edited_gene_path,
        cluster_number_row = cluster_row_vec,
        cluster_number_col = cluster_col_vec,
        groups_order = input$condition_order_09d,
        custom_color = strsplit(input$color_input_09d, ",")[[1]],
        genic_region = input$genic_region_09d,
        using_package_f = input$using_package_f_09d,
        scale_function = input$scale_function_09d,
        abs_log2FC = as.numeric(strsplit(input$abs_log2FC_09d, ",")[[1]]),
        FDR = as.numeric(strsplit(input$FDR_09d, ",")[[1]]),
        tpm_cutoff = as.numeric(strsplit(input$tpm_cutoff_09d, ",")[[1]])
      )
      
      removeModal()
      showModal(modalDialog(
        title = "Success",
        "✅ Heatmap generation complete!",
        "Click 'Download Results' to download all result.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("❌ Heatmap generation failed:", e$message),
        easyClose = TRUE
      ))
    })
  })
}
