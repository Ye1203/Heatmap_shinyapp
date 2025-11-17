source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/heatmap_ChipSeq.R")

tab_chipseq_ui <- function() {
  fluidPage(
    useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        h4("CHIP-Seq"),
        fileInput("datafile_chipseq", "Upload CSV or Excel file", accept = c(".csv", ".xlsx", ".xls")),
        numericInput("cluster_col_chipseq", "Cluster Number (Columns):", value = 1, min = 1),
        numericInput("cluster_row_chipseq", "Cluster Number (Rows):", value = 1, min = 1),
        textInput("color_input_chipseq", "Custom Colors (comma separated, #RGB code):", "red,yellow,#3fdc04"),
        selectInput("original_datatype", "Original Data Type", 
                    choices = c("log2", "linear"), selected = "log2"),
        selectInput("datatype_want", "The Type of Output Data You Want", 
                    choices = c("both", "log2", "linear"), selected = "both"),
        actionButton("run_btn_chipseq", "Generate Heatmaps"),
        br(),
        br(),
        downloadButton("download_zip_chipseq", "Download Results")
      ),
      mainPanel(
        verbatimTextOutput("status_chipseq")
      )
    )
  )
}

tab_chipseq_server <- function(input, output, session) {
  output_dir_chipseq <- reactiveVal(tempdir())
  
  disable("download_zip_chipseq")
  
  observeEvent(input$run_btn_chipseq, {
    req(input$datafile_chipseq)
    
    shinyjs::disable("run_btn_chipseq")
    shinyjs::html("run_btn_chipseq", "Processing...")
    showModal(modalDialog(
      title = NULL,
      "Generating heatmaps... Please wait.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    disable("download_zip_chipseq")
    
    out_dir <- file.path(tempdir(), paste0("results_chipseq_", as.integer(Sys.time())))
    dir.create(out_dir)
    output_dir_chipseq(out_dir)
    
    cluster_col <- input$cluster_col_chipseq
    cluster_row <- input$cluster_row_chipseq
    custom_colors <- trimws(unlist(strsplit(input$color_input_chipseq, ",")))
    
    tryCatch({
      result <- do_heatmap_chipseq(
        original_data = input$datafile_chipseq$datapath,
        original_name = input$datafile_chipseq$name,   
        output_folder = out_dir,
        cluster_number_col = cluster_col,
        cluster_number_row = cluster_row,
        custom_color = custom_colors,
        original_datatype = input$original_datatype,
        datatype_want = input$datatype_want
      )
      
      enable("download_zip_chipseq")
      
      output$status_chipseq <- renderText({
        if (result > 6500) {
          paste0(
            "✅ Heatmap generation completed.\n",
            "Files saved to temporary folder.\n",
            "Click \"Download Results\" to download all files.\n",
            "❗ Input row is larger than 6500.\n",
            "Adobe Acrobat may not open the file correctly.\n",
            "Please use other apps like 'Preview' or split the heatmap into multiple pages."
          )
        } else {
          paste0(
            "✅ Heatmap generation completed.\n",
            "Files saved to temporary folder.\n",
            "Click \"Download Results\" to download all files."
          )
        }
      })
    }, error = function(e) {
      output$status_chipseq <- renderText({
        paste("❌ Error:", e$message)
      })
    }, finally = {
      removeModal()
      shinyjs::enable("run_btn_chipseq")
      shinyjs::html("run_btn_chipseq", "Generate Heatmaps")
    })
  })
  
  output$download_zip_chipseq <- downloadHandler(
    filename = function() {
      paste0(file_path_sans_ext(basename(input$datafile_chipseq$name)), "_Heatmap_result.zip")
    },
    content = function(file) {
      zip::zipr(zipfile = file, files = list.files(output_dir_chipseq(), full.names = TRUE))
    }
  )
}
