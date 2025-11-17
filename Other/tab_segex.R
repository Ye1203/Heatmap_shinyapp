
source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/heatmap_segex.R")
tab_segex_ui <- function() {
  fluidPage(
    useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        h4("segex"),
        fileInput("datafile", "Upload CSV or Excel file", accept = c(".csv", ".xlsx", ".xls")),
        numericInput("cluster_col", "Cluster Number (Columns):", value = 1, min = 1),
        numericInput("cluster_row", "Cluster Number (Rows):", value = 1, min = 1),
        textInput("color_input", "Custom Colors (comma separated, #RGB code):", "red,yellow,#3fdc04"),
        actionButton("run_btn", "Generate Heatmaps"),
        br(),
        br(),
        downloadButton("download_zip", "Download Results")
      ),
      mainPanel(
        verbatimTextOutput("status")
      )
    )
  )
}

tab_segex_server <- function(input, output, session) {
  output_dir <- reactiveVal(tempdir())
  
  disable("download_zip")
  
  observeEvent(input$run_btn, {
    req(input$datafile)
    
    shinyjs::disable("run_btn")
    shinyjs::html("run_btn", "Processing...")
    showModal(modalDialog(
      title = NULL,
      "Generating heatmaps... Please wait.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    disable("download_zip")
    
    out_dir <- file.path(tempdir(), paste0("results_", as.integer(Sys.time())))
    dir.create(out_dir)
    output_dir(out_dir)
    
    cluster_col <- input$cluster_col
    cluster_row <- input$cluster_row
    custom_colors <- trimws(unlist(strsplit(input$color_input, ",")))
    
    tryCatch({
      result <- do_heatmap_segex(
        original_data = input$datafile$datapath,
        original_name = input$datafile$name,   
        output_folder = out_dir,
        cluster_number_col = cluster_col,
        cluster_number_row = cluster_row,
        custom_color = custom_colors
      )
      
      enable("download_zip")
      
      output$status <- renderText({
        if (result > 6500) {
          paste0(
            "✅ Heatmap generation completed.\n",
            "Files saved to temporary folder.\n",
            "Click \"Download Results\" to download all files.\n",
            "❗ Input row is larger than 6500.\n",
            "Adobe Acrobat may not open the file correctly.\n",
            "Please use other apps like 'Preview' or upload data in batches."
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
      output$status <- renderText({
        paste("❌ Error:", e$message)
      })
    }, finally = {
      removeModal()
      shinyjs::enable("run_btn")
      shinyjs::html("run_btn", "Generate Heatmaps")
    })
  })
  
  output$download_zip <- downloadHandler(
    filename = function() {
      paste0(file_path_sans_ext(basename(input$datafile$name)), "_Heatmap_result.zip")
    },
    content = function(file) {
      zip::zipr(zipfile = file, files = list.files(output_dir(), full.names = TRUE))
    }
  )
}
