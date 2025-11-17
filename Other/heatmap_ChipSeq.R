determine_group<- function(n) {
  if (n <= 26) {
    return(LETTERS[1:n])
  }
  s <- character(n)
  for (i in seq_len(n)) {
    num <- i
    name <- ""
    while (num > 0) {
      num <- num - 1
      name <- paste0(LETTERS[(num %% 26) + 1], name)
      num <- num %/% 26
    }
    s[i] <- name
  }
  s
}
heatmap_group_chipseq <- function(original_name,
                                  data,
                                  output_folder,
                                  scale_function,
                                  cluster_number_col,
                                  cluster_number_row,
                                  custom_color,
                                  col_annotate,
                                  col_color,
                                  from,
                                  original_datatype,
                                  datatype_want){
  
  df <- data
  col_dataprocessing <- colnames(data)[-1]
  if(original_datatype == "log2" & datatype_want == "linear"){
    df[, col_dataprocessing] <- lapply(df[, col_dataprocessing], function(x) 2^as.numeric(x))
  }else if(original_datatype == "linear" & datatype_want == "log2"){
    df[, col_dataprocessing] <- lapply(df[, col_dataprocessing], function(x) log2(as.numeric(x)+1))
  }
  if (scale_function == "maxscale") {
    min_scale <- apply(select(df, col_dataprocessing), 1, min, na.rm = TRUE)
    max_scale <- apply(select(df, col_dataprocessing), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    
    df <- df %>%
      mutate(across(all_of(col_dataprocessing), 
                    ~ case_when(
                      .x != 0 ~ as.numeric(.x) / max_norm,
                      TRUE ~ 0
                    ),
                    .names = paste0(datatype_want, "_maxscale_{.col}")
      ))
    
  } else if (scale_function == "zscale") {
    mean_x <- apply(select(df, col_dataprocessing), 1, mean, na.rm = TRUE)
    sd_x <- apply(select(df, col_dataprocessing), 1, sd, na.rm = TRUE)
    df <- df %>%
      mutate(across(col_dataprocessing, 
                    ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                    .names = paste0(datatype_want, "_zscale_{.col}"))) %>%
      as.data.frame()  
    min_scale <- apply(select(df, starts_with(paste0(datatype_want, "_zscale"))), 1, min, na.rm = TRUE)
    max_scale <- apply(select(df, starts_with(paste0(datatype_want, "_zscale"))), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    df <- df %>%
      mutate(across(starts_with(paste0(datatype_want, "_zscale_")), 
                    ~ case_when(
                      .x != 0  ~ (.x / max_norm),
                      .x == 0 ~ 0
                    )))
  }
  first_column <- colnames(df)[1]
  if (scale_function == "maxscale") {
    df <- df[, c(first_column, grep("*maxscale_", colnames(df), value = TRUE)), drop = FALSE]
  } else if (scale_function == "zscale") {
    df <- df[, c(first_column, grep("*zscale_", colnames(df), value = TRUE)), drop = FALSE]
  }
  data_save <- df[, colnames(df) != first_column]
  
  # remove scale_tpm_ on df
  colnames(df) <- ifelse(
    colnames(df) == first_column,
    first_column,
    sub(".*scale_", "", colnames(df))
  )
  row_names_to_set <- df[[first_column]]
  df[[first_column]] <- NULL
  rownames(df) <- row_names_to_set
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
                        annotation_col = col_annotate,
                        annotation_colors = col_color,
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
  col_order_table$Group <- col_annotate[col_order_table$COL_NAME, "Group"]
  col_order_table <- col_order_table %>%
    select(COL_NUMBER,CLUSTER_NUMBER,Group,COL_NAME)
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
                                     annotation_col = col_annotate,
                                     annotation_colors = col_color,
                                     treeheight_row = 200,
                                     angle_col = 90, 
                                     color = custom_color,
                                     silent = TRUE)
  col_order_table_noclustering <- data.frame(
    COL_NUMBER = seq_along(colnames(df)),  
    COL_NAME = colnames(df),
    stringsAsFactors = FALSE
  )
  col_order_table_noclustering$Group <- col_annotate[col_order_table_noclustering$COL_NAME, "Group"]
  col_order_table_noclustering <- col_order_table_noclustering %>%
    select(COL_NUMBER,Group,COL_NAME)
  colnames(df) <- seq_len(ncol(df))
  rownames(col_annotate) <- as.character(seq_len(ncol(df)))
  
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
                                    annotation_col = col_annotate,
                                    annotation_colors = col_color,
                                    treeheight_row = 130,
                                    angle_col = 90, 
                                    show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                                    show_rownames = FALSE,
                                    color = custom_color,
                                    silent = TRUE)
  
  df <- df[, col_order_origin] 
  col_annotate <- col_annotate[col_order_origin, , drop = FALSE]
  #FIX HERE
  colnames(df) <- seq_len(ncol(df))
  rownames(col_annotate) <- as.character(seq_len(ncol(df)))
  
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
                       annotation_col = col_annotate,
                       annotation_colors = col_color,
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
  
  cluster_long_path <- paste0("Heatmap_SamplesClustered_Long_",scale_function,"_",datatype_want,".pdf")
  nocluster_long_path <- paste0("Heatmap_NoSamplesClustering_Long_",scale_function,"_",datatype_want,".pdf")
  cluster_short_path <- paste0("Heatmap_SamplesClustered_Short_",scale_function,"_",datatype_want,".pdf")
  nocluster_short_path <- paste0("Heatmap_NoSamplesClustering_Short_",scale_function,"_",datatype_want,".pdf")
  # SAVE integrated long cluster
  pdf(file.path(output_folder, cluster_long_path), width = origin_pdf_width, height = origin_pdf_height)
  grid.newpage()
  title <- paste0(file_path_sans_ext(basename(original_name)),"_SampleClustered",
                  "\nSCALE FUNCTION: ", scale_function,",\t\tNORMALIZED: ", from,
                  "\n", original_datatype, "->", datatype_want,
                  "\nGENOMIC REGION (Row): ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_origin$gtable) 
  popViewport()
  dev.off()
  
  
  
  # SAVE integrated long noclustered
  pdf(file.path(output_folder, nocluster_long_path), width = origin_pdf_width_noclustering, height = origin_pdf_height_noclustering)
  grid.newpage()
  title <- paste0(file_path_sans_ext(basename(original_name)),"_NoSampleClustering",
                  "\nSCALE FUNCTION: ", scale_function,",\t\tNORMALIZED: ", from,
                  "\n", original_datatype, "->", datatype_want,
                  "\nGENOMIC REGION (Row): ", nrow(df))
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
  info <- paste0("\nSCALE FUNCTION: ", scale_function,"\nNORMALIZED: ", from, "\n", original_datatype, "->", datatype_want, "\n\nGENOMIC REGION (Row):", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
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
  pushViewport(viewport(x = 1, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  # SAVE integrated short noclustering
  pdf(file.path(output_folder, nocluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(file_path_sans_ext(basename(original_name)),"_NoSampleClustering\t", scale_function)
  if (nchar(title) > 80) {
    title <- paste0(substr(title, 1, 80), "-\n", substr(title, 81, nchar(title)))
  }
  info <- paste0("\nSCALE FUNCTION: ", scale_function,"\nNORMALIZED: ", from, "\n", original_datatype, "->", datatype_want, "\n\nGENOMIC REGION (Row):", nrow(df))
  grid.newpage()
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
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
  pushViewport(viewport(x = 1, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
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

do_heatmap_chipseq <- function(original_data,
                               original_name,
                               output_folder,
                               cluster_number_col = 1,
                               cluster_number_row = 1,
                               custom_color = NULL,
                               original_datatype = "log2",
                               datatype_want = "both"){
  
  if (is.null(custom_color) || length(custom_color) == 0 || length(custom_color) == 1) {
    custom_color <- c("red", "yellow", "#3fdc04")
  }
  custom_color <- colorRampPalette(custom_color)(200)
  for(s in excel_sheets(original_data)){
    new_output_folder <- file.path(output_folder, s)
    dir.create(new_output_folder, recursive = TRUE)
    data <- read_excel(original_data, sheet = s)
    all_cols <- colnames(data)
    sorted_cols <- c(all_cols[1], sort(setdiff(all_cols, all_cols[1])))
    data <- data[, sorted_cols]
    
    col_annotate <- colnames(data)[-1]
    group_names <- sub("([(:]).*", "", col_annotate)
    group_map <- setNames(determine_group(length(unique(group_names))), unique(group_names))
    col_annotate <- data.frame(
      original_colname = col_annotate,
      Group = group_map[group_names],
      stringsAsFactors = FALSE
    )
    
    #from heamap_09d.R
    group_list <- unique(col_annotate$Group)
    group_colors <- color_selected(length(group_list))
    names(group_colors) <- group_list
    col_color <- list(Group = group_colors)
    rownames(col_annotate) <- col_annotate$original_colname
    col_annotate$original_colname <- NULL
    if(datatype_want == "both"){
      maxscale_result_log2 <- heatmap_group_chipseq(original_name,
                                                    data = data,
                                                    output_folder = new_output_folder,
                                                    scale_function = "maxscale",
                                                    cluster_number_col = cluster_number_col,
                                                    cluster_number_row = cluster_number_row,
                                                    custom_color = custom_color,
                                                    col_annotate = col_annotate,
                                                    col_color = col_color,
                                                    from = s,
                                                    original_datatype = original_datatype,
                                                    datatype_want = "log2")  
      # ZSCALE
      zscale_result_log2 <- heatmap_group_chipseq(original_name,
                                                  data = data,
                                                  output_folder = new_output_folder,
                                                  scale_function = "zscale",
                                                  cluster_number_col = cluster_number_col,
                                                  cluster_number_row = cluster_number_row,
                                                  custom_color = custom_color,
                                                  col_annotate = col_annotate,
                                                  col_color = col_color,
                                                  from = s,
                                                  original_datatype = original_datatype,
                                                  datatype_want = "log2")
      # MAXSCALE linear
      maxscale_result_linear <- heatmap_group_chipseq(original_name,
                                                      data = data,
                                                      output_folder = new_output_folder,
                                                      scale_function = "maxscale",
                                                      cluster_number_col = cluster_number_col,
                                                      cluster_number_row = cluster_number_row,
                                                      custom_color = custom_color,
                                                      col_annotate = col_annotate,
                                                      col_color = col_color,
                                                      from = s,
                                                      original_datatype = original_datatype,
                                                      datatype_want = "linear")  
      # ZSCALE linear
      zscale_result_linear <- heatmap_group_chipseq(original_name,
                                                    data = data,
                                                    output_folder = new_output_folder,
                                                    scale_function = "zscale",
                                                    cluster_number_col = cluster_number_col,
                                                    cluster_number_row = cluster_number_row,
                                                    custom_color = custom_color,
                                                    col_annotate = col_annotate,
                                                    col_color = col_color,
                                                    from = s,
                                                    original_datatype = original_datatype,
                                                    datatype_want = "linear")
      
      excel_sheets <- list(
        Log2_Maxscale_NoCluster = maxscale_result_log2[["col_order_table_noclustering"]],  
        Log2_Maxscale_Cluster = maxscale_result_log2[["col_order_table"]],
        Linear_Maxscale_NoCluster = maxscale_result_linear[["col_order_table_noclustering"]],  
        Linear_Maxscale_Cluster = maxscale_result_linear[["col_order_table"]],
        Log2_Zscale_NoCluster = zscale_result_log2[["col_order_table_noclustering"]],  
        Log2_Zscale_Cluster = zscale_result_log2[["col_order_table"]],
        Linear_Zscale_NoCluster = zscale_result_linear[["col_order_table_noclustering"]],  
        Linear_Zscale_Cluster = zscale_result_linear[["col_order_table"]]
      )
      
      write.xlsx(excel_sheets, file = file.path(new_output_folder, "Summary_SampleOrder.xlsx"), rowNames = FALSE)
      
      first_sheet <- cbind(
        data[, 1, drop = FALSE],                   
        log2_max_order = match(seq_len(nrow(data)), maxscale_result_log2[["row_order_origin"]]),
        log2_max_cluster = maxscale_result_log2[["row_clusters"]],
        log2_zscale_order = match(seq_len(nrow(data)), zscale_result_log2[["row_order_origin"]]),
        log2_zscale_cluster = zscale_result_log2[["row_clusters"]],
        linear_max_order = match(seq_len(nrow(data)), maxscale_result_linear[["row_order_origin"]]),
        linear_max_cluster = maxscale_result_linear[["row_clusters"]],
        linear_zscale_order = match(seq_len(nrow(data)), zscale_result_linear[["row_order_origin"]]),
        linear_zscale_cluster = zscale_result_linear[["row_clusters"]],
        data[, -1, drop = FALSE]
      )
      
      second_sheet <- cbind(
        data[, 1, drop = FALSE],                   
        log2_max_order = match(seq_len(nrow(data)), maxscale_result_log2[["row_order_origin"]]),
        log2_max_cluster = maxscale_result_log2[["row_clusters"]],
        log2_zscale_order = match(seq_len(nrow(data)), zscale_result_log2[["row_order_origin"]]),
        log2_zscale_cluster = zscale_result_log2[["row_clusters"]],
        linear_max_order = match(seq_len(nrow(data)), maxscale_result_linear[["row_order_origin"]]),
        linear_max_cluster = maxscale_result_linear[["row_clusters"]],
        linear_zscale_order = match(seq_len(nrow(data)), zscale_result_linear[["row_order_origin"]]),
        linear_zscale_cluster = zscale_result_linear[["row_clusters"]],
        maxscale_result_log2[["data_save"]],
        zscale_result_log2[["data_save"]],
        maxscale_result_linear[["data_save"]],
        zscale_result_linear[["data_save"]]
      )
      
      wb <- createWorkbook()
      addWorksheet(wb, paste0("OriginalData_",toTitleCase(original_datatype)))
      writeData(wb, sheet = paste0("OriginalData_",toTitleCase(original_datatype)), first_sheet)
      addWorksheet(wb, "MaxScale_Zscale_Values")
      writeData(wb, sheet = "MaxScale_Zscale_Values", second_sheet)
      saveWorkbook(wb, file = file.path(new_output_folder,  "HeatmapData.xlsx"), overwrite = TRUE)
      }
    else{maxscale_result <- heatmap_group_chipseq(original_name,
                                                       data = data,
                                                       output_folder = new_output_folder,
                                                       scale_function = "maxscale",
                                                       cluster_number_col = cluster_number_col,
                                                       cluster_number_row = cluster_number_row,
                                                       custom_color = custom_color,
                                                       col_annotate = col_annotate,
                                                       col_color = col_color,
                                                       from = s,
                                                       original_datatype = original_datatype,
                                                       datatype_want = datatype_want)  
    # ZSCALE
    zscale_result <- heatmap_group_chipseq(original_name,
                                                data = data,
                                                output_folder = new_output_folder,
                                                scale_function = "zscale",
                                                cluster_number_col = cluster_number_col,
                                                cluster_number_row = cluster_number_row,
                                                custom_color = custom_color,
                                                col_annotate = col_annotate,
                                                col_color = col_color,
                                                from = s,
                                                original_datatype = original_datatype,
                                                datatype_want = datatype_want)
    
    excel_sheets <- setNames(
      list(
        maxscale_result[["col_order_table_noclustering"]],
        maxscale_result[["col_order_table"]],
        zscale_result[["col_order_table_noclustering"]],
        zscale_result[["col_order_table"]]
      ),
      paste0(
        toTitleCase(datatype_want),
        c("_Maxscale_NoCluster",
          "_Maxscale_Cluster",
          "_Zscale_NoCluster",
          "_Zscale_Cluster")
      )
    )
    
    
    write.xlsx(excel_sheets, file = file.path(new_output_folder, "Summary_SampleOrder.xlsx"), rowNames = FALSE)
    
    first_sheet <- cbind(
      data[, 1, drop = FALSE],
      setNames(
        data.frame(
          match(seq_len(nrow(data)), maxscale_result[["row_order_origin"]]),
          maxscale_result[["row_clusters"]],
          match(seq_len(nrow(data)), zscale_result[["row_order_origin"]]),
          zscale_result[["row_clusters"]]
        ),
        paste0(toTitleCase(datatype_want), c("_max_order", "_max_cluster", "_zscale_order", "_zscale_cluster"))
      ),
      data[, -1, drop = FALSE]
    )
    
    
    second_sheet <- cbind(
      data[, 1, drop = FALSE],
      setNames(
        data.frame(
          match(seq_len(nrow(data)), maxscale_result[["row_order_origin"]]),
          maxscale_result[["row_clusters"]],
          match(seq_len(nrow(data)), zscale_result[["row_order_origin"]]),
          zscale_result[["row_clusters"]]
        ),
        paste0(toTitleCase(datatype_want), c("_max_order", "_max_cluster", "_zscale_order", "_zscale_cluster"))
      ),
      maxscale_result[["data_save"]],
      zscale_result[["data_save"]]
    )
    
    wb <- createWorkbook()
    addWorksheet(wb, paste0("OriginalData_",toTitleCase(original_datatype)))
    writeData(wb, sheet = paste0("OriginalData_",toTitleCase(original_datatype)), first_sheet)
    addWorksheet(wb, "MaxScale_Zscale_Values")
    writeData(wb, sheet = "MaxScale_Zscale_Values", second_sheet)
    saveWorkbook(wb, file = file.path(new_output_folder,  "HeatmapData.xlsx"), overwrite = TRUE)
    }
  }
  return(nrow(data))
}