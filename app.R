library(shiny)
library(shinyjs)
library(openxlsx)
library(shinyjqui)
library(sortable)
library(zip)
library(tidyverse) #2.0.0
library(RColorBrewer) #1.1.3
library(pheatmap) #1.0.13
library(tools) #4.4.3
library(gridExtra) #2.3
library(openxlsx) #4.2.8
library(grid)
library(readxl)
options(shiny.maxRequestSize = 100 * 1024^2)

source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/tab_09d.R")
source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/tab_segex.R")
source("/projectnb/wax-es/00_shinyapp/Heatmap/Other/tab_ChipSeq.R")
options(warn = -1)

ui <- fluidPage(
  titlePanel("Heatmap Generator"),
  tabsetPanel(id = "tabs",
              tabPanel("09d", tab_09d_ui),
              tabPanel("segex", tab_segex_ui()),
              tabPanel("CHIP-seq", tab_chipseq_ui())
  )
)

server <- function(input, output, session) {
  tab_09d_server(input, output, session)
  tab_segex_server(input, output, session)
  tab_chipseq_server(input, output, session)
  
}

shinyApp(ui, server)
