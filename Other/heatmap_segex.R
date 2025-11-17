
heatmap_group_segex <- function(original_name,
                          data,
                          output_folder,
                          scale_function,
                          cluster_number_col,
                          cluster_number_row,
                          custom_color){
  
  # Get data start with tpm
  tpm_cols <- grep("^tpm_", colnames(data), value = TRUE)
  df <- data[, c("Probe ID", tpm_cols)]
  # Get unique column
  base_names <- sub("\\.\\d+$", "", colnames(df))
  unique_indices <- !duplicated(base_names)
  df <- df[, unique_indices]
  
  if (scale_function == "maxscale") {
    min_scale <- apply(select(df, starts_with("tpm_")), 1, min, na.rm = TRUE)
    max_scale <- apply(select(df, starts_with("tpm_")), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    df <- df %>%
      mutate(across(starts_with("tpm_"), 
                    ~ case_when(
                      .x != 0  ~ as.numeric(.x) / max_norm,
                      .x == 0 ~ 0
                    ), .names = "maxscale_{.col}")) %>%
      as.data.frame()}
  else if (scale_function == "zscale") {
    mean_x <- apply(select(df, starts_with("tpm_")), 1, mean, na.rm = TRUE)
    sd_x <- apply(select(df, starts_with("tpm_")), 1, sd, na.rm = TRUE)
    df <- df %>%
      mutate(across(starts_with("tpm_"), 
                    ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                    .names = "zscale_{.col}")) %>%
      as.data.frame()  
    min_scale <- apply(select(df, starts_with("zscale_")), 1, min, na.rm = TRUE)
    max_scale <- apply(select(df, starts_with("zscale_")), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    df <- df %>%
      mutate(across(starts_with("zscale_"), 
                    ~ case_when(
                      .x != 0  ~ (.x / max_norm),
                      .x == 0 ~ 0
                    )))
  }
  
  df <- df[, !grepl("^tpm_", colnames(df))]
  
  data_save <- df[, colnames(df) != "Probe ID"]
  
  # remove scale_tpm_ on df
  colnames(df) <- ifelse(
    colnames(df) == "Probe ID",
    "Probe ID",
    sub(".*tpm_", "", colnames(df))
  )
  
  rownames(df) <- df$`Probe ID`
  
  df$`Probe ID` <- NULL
  
  # HM clustered long
  HM_origin <- pheatmap(df, 
                        cellwidth = 20, 
                        cellheight = 2 , 
                        fontsize_col = 10,
                        fontsize_row = 2,
                        border_color = NA, 
                        cutree_cols = cluster_number_col, 
                        cutree_rows = cluster_number_row, 
                        cluster_rows = TRUE, 
                        cluster_cols = TRUE,
                        treeheight_row = 200,
                        treeheight_col = 80,
                        angle_col = 90, 
                        color = custom_color,
                        silent = TRUE)
  
  
  row_order_origin <- HM_origin$tree_row$order
  col_order_origin <- HM_origin$tree_col$order
  row_clusters <- cutree(HM_origin$tree_row, k = cluster_number_row)
  col_clusters <- cutree(HM_origin$tree_col, k = cluster_number_col)
  
  # change cluster number
  ordered_row_clusters <- row_clusters[row_order_origin]
  unique_ordered_clusters <- unique(ordered_row_clusters)
  new_cluster_mapping <- setNames(1:length(unique_ordered_clusters), unique_ordered_clusters)
  row_clusters_ordered <- new_cluster_mapping[as.character(row_clusters)]
  names(row_clusters_ordered) <- names(row_clusters)
  row_clusters <- row_clusters_ordered
  
  #save col order
  col_order_table <- data.frame(
    COL_NUMBER = seq_along(col_order_origin),  
    CLUSTER_NUMBER = col_clusters[col_order_origin],
    COL_NAME = colnames(df)[col_order_origin],
    stringsAsFactors = FALSE
  )
  
  heatmap_width_ratio <- 0.75  
  heatmap_height_ratio <- 0.9 
  
  cellwidth <- (8.27 * heatmap_width_ratio * 72-275 - 5 *cluster_number_col) / ncol(df)    
  cellheight <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72-100 - 5 * cluster_number_row) / nrow(df)}else{(11.69 * heatmap_height_ratio * 72- 110 - 5 * cluster_number_row) / nrow(df)}
  cellwidth_noclustering <- (8.27 * heatmap_width_ratio * 72-260) / ncol(df)    
  cellheight_noclustering <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72 -15 - 5 * cluster_number_row) / nrow(df)}else{(11.69 * heatmap_height_ratio * 72- 25 - 5 * cluster_number_row) / nrow(df)}
  
  # HM noclustering long
  HM_noclustering_origin <- pheatmap(df, 
                                     cellwidth = 20, 
                                     cellheight = 2 , 
                                     fontsize_col = 10,
                                     fontsize_row = 2,
                                     border_color = NA, 
                                     cutree_cols = cluster_number_col, 
                                     cutree_rows = cluster_number_row, 
                                     cluster_rows = TRUE, 
                                     cluster_cols = FALSE,
                                     treeheight_row = 200,
                                     angle_col = 90, 
                                     color = custom_color,
                                     silent = TRUE)
  
  col_order_table_noclustering <- data.frame(
    COL_NUMBER = seq_along(colnames(df)),  
    COL_NAME = colnames(df),
    stringsAsFactors = FALSE
  )
  
  colnames(df) <- seq_len(ncol(df))
  
  HM_noclustering_scale <- pheatmap(df, 
                                    cellwidth = cellwidth_noclustering, 
                                    cellheight = cellheight_noclustering,
                                    fontsize_col = 100/ncol(df),
                                    fontsize_row = 2,
                                    border_color = NA, 
                                    cutree_cols = cluster_number_col, 
                                    cutree_rows = cluster_number_row, 
                                    cluster_rows = TRUE, 
                                    cluster_cols = FALSE,
                                    treeheight_row = 130,
                                    angle_col = 90, 
                                    show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                                    show_rownames = FALSE,
                                    color = custom_color,
                                    silent = TRUE)
  
  df <- df[, col_order_origin] 
  colnames(df) <- seq_len(ncol(df))
  HM_scale <- pheatmap(df, 
                       cellwidth = cellwidth, 
                       cellheight = cellheight,
                       fontsize_col = 100/ncol(df),
                       fontsize_row = 2,
                       border_color = NA, 
                       cutree_cols = cluster_number_col, 
                       cutree_rows = cluster_number_row, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE,
                       treeheight_row = 130,
                       treeheight_col = 80,
                       angle_col = 90, 
                       show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                       show_rownames = FALSE,
                       color = custom_color,
                       silent = TRUE)
  
  
  origin_pdf_width <- (ncol(df) * 20) / 72 + 6 + cluster_number_row
  origin_pdf_height <- (nrow(df) * 2) / 72 + 10 +  cluster_number_col
  origin_pdf_width_noclustering <- (ncol(df) * 20) / 72 + 6
  origin_pdf_height_noclustering <- (nrow(df) * 2) / 72 + 8
  
  cluster_long_path <- paste0("Heatmap_SamplesClustered_Long_",scale_function,".pdf")
  nocluster_long_path <- paste0("Heatmap_NoSamplesClustering_Long_",scale_function,".pdf")
  cluster_short_path <- paste0("Heatmap_SamplesClustered_Short_",scale_function,".pdf")
  nocluster_short_path <- paste0("Heatmap_NoSamplesClustering_Short_",scale_function,".pdf")
  
  # SAVE integrated long cluster
  pdf(file.path(output_folder, cluster_long_path), width = origin_pdf_width, height = origin_pdf_height)
  grid.newpage()
  title <- paste0(file_path_sans_ext(basename(original_name)),"_SampleClustered",
                  "\t", scale_function,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_origin$gtable) 
  popViewport()
  dev.off()
  
  
  
  # SAVE integrated long noclustered
  pdf(file.path(output_folder, nocluster_long_path), width = origin_pdf_width_noclustering, height = origin_pdf_height_noclustering)
  grid.newpage()
  title <- paste0(file_path_sans_ext(basename(original_name)),"_NoSampleClustering",
                  "\t", scale_function,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_noclustering_origin$gtable) 
  popViewport()
  dev.off()
  rm(HM_noclustering_origin)
  
  # SAVE integrated short cluster
  pdf(file.path(output_folder, cluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(file_path_sans_ext(basename(original_name)),"_SampleClustered\t", scale_function)
  if (nchar(title) > 80) {
    title <- paste0(substr(title, 1, 80), "-\n", substr(title, 81, nchar(title)))
  }
  info <- paste0("GENE NUMBER: ", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.35, y = 0.95, width = heatmap_width_ratio, height = heatmap_height_ratio, just = c("center","top")))
  print(HM_scale, newpage = FALSE)
  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                padding = unit(c(1.5, 1.5), "mm")), 
    colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm")), 
    rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm"))
  )
  grid.text("Gene clusters are numbered from\ntop to bottom", x = 1, y = 1, gp = gpar(fontsize = 8, col = "red"), just = c("left", "top"))
  grid.text("First column: Index (->)\nSecond column: Group cluster number\nThird column: Sample Group", x = 1, y = 0.01, gp = gpar(fontsize = 8, col = "red"), just = c("left", "top"))
  grid.text(info, x = 1, y = 0.98, gp = gpar(fontsize = 8), just = c("left","top"))
  pushViewport(viewport(x = 1, y = 0.8, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  # SAVE integrated short noclustering
  pdf(file.path(output_folder, nocluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(file_path_sans_ext(basename(original_name)),"_NoSampleClustering\t", scale_function)
  if (nchar(title) > 80) {
    title <- paste0(substr(title, 1, 80), "-\n", substr(title, 81, nchar(title)))
  }
  info <- paste0("GENE NUMBER: ", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.35, y = 0.95, width = heatmap_width_ratio, height = heatmap_height_ratio, just = c("center","top")))
  print(HM_noclustering_scale, newpage = FALSE)
  table_theme <- ttheme_default(
    core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                padding = unit(c(1.5, 1.5), "mm")), 
    colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm")), 
    rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                   padding = unit(c(1.5, 1.5), "mm"))
  )
  grid.text("Gene clusters are numbered from\ntop to bottom", x = 1, y = 1, gp = gpar(fontsize = 8, col = "red"), just = c("left", "top"))
  grid.text("First column: Index (->)\nSecond column: Sample Group", x = 1, y = 0.01, gp = gpar(fontsize = 8, col = "red"), just = c("left", "top"))
  
  grid.text(info, x = 1, y = 0.98, gp = gpar(fontsize = 8), just = c("left","top"))
  pushViewport(viewport(x = 1, y = 0.8, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table_noclustering, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  return(list(
    col_order_table = col_order_table,
    col_order_table_noclustering = col_order_table_noclustering,
    row_order_origin = row_order_origin,
    row_clusters = row_clusters,
    data_save = data_save
  ))
}

do_heatmap_segex <- function(original_data,
                       original_name,
                       output_folder,
                       cluster_number_col = 1,
                       cluster_number_row = 1,
                       custom_color = NULL){
  
  if (is.null(custom_color) || length(custom_color) == 0 || length(custom_color) == 1) {
    custom_color <- c("red", "yellow", "#3fdc04")
  }
  custom_color <- colorRampPalette(custom_color)(200)
  if (tolower(file_ext(original_data)) == "csv") {
    lines <- readLines(original_data)
    down_line <- which(grepl("^.*,\\s*Down:", lines))
    info <- lines[down_line+1]
    data <- read.table(original_data, sep = ",", skip = down_line + 1, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    matches <- regmatches(info, gregexpr("[0-9]+:[^,]+/[^,]+", info))[[1]]
    probe_col <- which(colnames(data) == "Probe ID")
    gene_col <- which(colnames(data) == "Gene Symbol")
    data <- data[, c(probe_col, (gene_col + 1):ncol(data))]
    numeric_cols <- setdiff(colnames(data), "Probe ID")
    data[numeric_cols] <- lapply(data[numeric_cols], function(x) as.numeric(as.character(x)))
  }
  else if (tolower(file_ext(original_data)) %in% c("xlsx", "xls")) {
    all_data <- read.xlsx(original_data, colNames = FALSE, check.names = FALSE)
    down_line <- which(all_data[[2]] == "Down:")
    info_line <- all_data[down_line + 1, ]
    info_string <- paste(as.character(info_line), collapse = ",")
    info <- unlist(strsplit(info_string, split = ","))
    data <- read.xlsx(original_data, startRow = down_line + 3, colNames = FALSE, check.names = FALSE)
    colnames(data) <- as.character(unlist(data[1, ]))
    data <- data[-1, ]
    valid_info <- info[!info %in% c("NA", "", NA)]
    matches <- regmatches(valid_info, gregexpr("[0-9]+:[^,]+/[^,]+", valid_info))
    matches <- unlist(matches)
    probe_col <- which(colnames(data) == "Probe ID")
    gene_col <- which(colnames(data) == "Gene Symbol")
    data <- data[, c(probe_col, (gene_col + 1):ncol(data))]
    numeric_cols <- setdiff(colnames(data), "Probe ID")
    data[numeric_cols] <- lapply(data[numeric_cols], function(x) as.numeric(as.character(x)))
  }
  
  new_colnames <- colnames(data)
  
  for (match in matches) {
    parts <- strsplit(match, ":|/")[[1]] 
    id <- parts[1]
    groupA <- parts[2]
    groupB <- parts[3]
    col_A <- paste0(id, ": Intensity-2")
    col_B <- paste0(id, ": Intensity-1")
    new_colnames[new_colnames == col_A] <- paste0("tpm_",groupA)
    new_colnames[new_colnames == col_B] <- paste0("tpm_",groupB)
  }
  
  colnames(data) <- new_colnames
  
  # MAXSCALE
  maxscale_result <- heatmap_group_segex(original_name,
                                   data = data,
                                   output_folder = output_folder,
                                   scale_function = "maxscale",
                                   cluster_number_col = cluster_number_col,
                                   cluster_number_row = cluster_number_row,
                                   custom_color = custom_color)  
  # ZSCALE
  zscale_result <- heatmap_group_segex(original_name,
                                 data = data,
                                 output_folder = output_folder,
                                 scale_function = "zscale",
                                 cluster_number_col = cluster_number_col,
                                 cluster_number_row = cluster_number_row,
                                 custom_color = custom_color)
  
  excel_sheets <- list(
    Maxscale_NoCluster_ColOrder = maxscale_result[["col_order_table_noclustering"]],  
    Maxscale_Cluster_ColOrder = maxscale_result[["col_order_table"]],
    Zscale_NoCluster_ColOrder = zscale_result[["col_order_table_noclustering"]],  
    Zscale_Cluster_ColOrder = zscale_result[["col_order_table"]]
  )
  
  write.xlsx(excel_sheets, file = file.path(output_folder, "Summary_SampleOrder.xlsx"), rowNames = FALSE)
  
  first_sheet <- cbind(
    data[, 1, drop = FALSE],                   
    max_order = match(seq_len(nrow(data)), maxscale_result[["row_order_origin"]]),
    max_cluster = maxscale_result[["row_clusters"]],
    zscale_order = match(seq_len(nrow(data)), zscale_result[["row_order_origin"]]),
    zscale_cluster = zscale_result[["row_clusters"]],
    data[, -1, drop = FALSE]
  )
  
  second_sheet <- cbind(
    data[, 1, drop = FALSE],                   
    max_order = match(seq_len(nrow(data)), maxscale_result[["row_order_origin"]]),
    max_cluster = maxscale_result[["row_clusters"]],
    zscale_order = match(seq_len(nrow(data)), zscale_result[["row_order_origin"]]),
    zscale_cluster = zscale_result[["row_clusters"]],
    maxscale_result[["data_save"]],
    zscale_result[["data_save"]]
  )
  
  
  wb <- createWorkbook()
  addWorksheet(wb, "TFS_OriginalData")
  writeData(wb, sheet = "TFS_OriginalData", first_sheet)
  addWorksheet(wb, "TFS_MaxScale_Zscale")
  writeData(wb, sheet = "TFS_MaxScale_Zscale", second_sheet)
  saveWorkbook(wb, file = file.path(output_folder, "HeatmapData.xlsx"), overwrite = TRUE)
  return(nrow(data))
  
}