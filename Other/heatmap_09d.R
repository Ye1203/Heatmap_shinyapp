read_data <- function(input) {
  if (!file.exists(input)) {
    stop("Error: The file does not exist. Please check the file path.")
  }
  
  if (grepl("\\.txt$", input, ignore.case = TRUE)) {
    data <- read.table(input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl("\\.csv$", input, ignore.case = TRUE)) {
    data <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl("\\.(xls|xlsx)$", input, ignore.case = TRUE)) {
    data <- read_excel(input)
  } else {
    stop("Error: Unsupported file format. Please provide a .txt, .csv, or .xls/.xlsx file.")
  }
  
  return(data)
}

# Match column order for noncluster
match_column <- function(reference, col_need_match){
  target_order <- reference$Condition_Name
  matched_columns <- sapply(target_order, function(cond) {
    matched <- grep(cond, col_need_match, value = TRUE)
    if (length(matched) == 0) return(NA)
    matched[1]
  })
  col_need_match <- as.vector(na.omit(matched_columns))
}

# Select color
color_selected <- function(color_length) {
  color_total <- c(
    "#e6194b","#ffe119","#46f0f0","#f58231","#bcf60c", 
    "#ff00ff","#9a6324","#fffac8","#e6beff","#00bfff", 
    "#ffd8b1","#00ff7f","#f5a9bc","#1e90ff","#ffa500",
    "#98fb98","#911eb4","#afeeee","#fa8072","#9acd32",
    "#3cb44b","#000075","#808000","#cd5c5c","#dda0dd",
    "#40e0d0","#ff69b4","#8a2be2","#c71585","#5f9ea0",
    "#dc143c","#87cefa","#ff6347","#9932cc","#00ced1",
    "#ff4500","#6a5acd","#b0e0e6","#d2691e","#a9a9f5",
    "#adff2f","#8b0000","#7fffd4","#00fa9a","#ba55d3",
    "#2e8b57","#ffdab9","#b22222","#ffe4e1","#7b68ee"
  )
  
  if (color_length <= length(color_total)) {
    return(color_total[1:color_length])
  }
  
  warning("groups is larger than 50, color will randomly select")
  palette_fn <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  extra_colors <- palette_fn(color_length - length(color_total))
  colors <- c(color_total, extra_colors)
  return(colors)
}

# Get tpm in any group, return gene_to_keep
get_genes_with_high_tpm_in_any_group <- function(data, annotation_col, tpm_cutoff) {
  # Get all column names starting with "tpm_"
  tpm_cols <- grep("^tpm_", colnames(data), value = TRUE)
  
  # Extract condition names by removing "tpm_" prefix
  condition_names <- sub("^tpm_", "", tpm_cols)
  
  # Create a mapping between tpm column names and their corresponding group
  tpm_info <- data.frame(
    tpm_col = tpm_cols,
    Condition_Name = condition_names,
    stringsAsFactors = FALSE
  )
  tpm_info$Group <- annotation_col[tpm_info$Condition_Name, "Group"]
  
  # Create a TPM matrix and rename columns with group names
  tpm_matrix <- data[, tpm_info$tpm_col]
  colnames(tpm_matrix) <- tpm_info$Group
  tpm_matrix$id <- data$id
  
  # Convert wide matrix to long format for group-wise processing
  tpm_long <- tpm_matrix %>%
    pivot_longer(-id, names_to = "Group", values_to = "TPM")
  
  # For each gene and group, check if all samples in the group have TPM > cutoff
  group_check <- tpm_long %>%
    group_by(id, Group) %>%
    summarise(all_above = all(TPM >= tpm_cutoff), .groups = "drop")
  
  # Keep genes where any group satisfies the condition
  gene_to_keep <- group_check %>%
    group_by(id) %>%
    summarise(any_group_all_above = any(all_above), .groups = "drop") %>%
    filter(any_group_all_above) %>%
    pull(id)
  
  return(gene_to_keep)
}

pair_heatmap_09d <- function(input,
                             edited_gene_list,
                             condition_1,
                             condition_2,
                             output_folder,
                             save_name,
                             custom_color,
                             using_package_f,
                             scale_function,
                             log2FC,
                             FDR,
                             tpm_cutoff,
                             cluster_number_row,
                             cluster_number_col,
                             order_for_noclustering,
                             annotation_colors,
                             folder_name){

  if (is.character(input) && length(input) == 1 && input != "") {
    if (!using_package_f %in% c("EdgeR", "DESeq2")) {
      stop("Error: 'using_package_f' must be either 'EdgeR' or 'DESeq2'.")
    }
    #READ DATA
    data <- read_data(input)
    if(!is.null(edited_gene_list)){
    data <- data[tolower(data$id) %in% tolower(edited_gene_list), ]}
    
    if (nrow(data) == 0) {
      stop("No matching genes found. Please check the Custom Gene List input.")
    }

    group_data <- data %>%
      select(id,starts_with("tpm_mean"))
    data <- data %>%
      select(all_of(colnames(data)[
        which(colnames(data) == "id") : (which(grepl("^rpkm_", colnames(data)))[1] - 1)
      ]),
      starts_with("tpm"), -starts_with("tpm_mean"),
      starts_with(using_package_f),
      )
    annotation_col <- order_for_noclustering
    rownames(annotation_col) <- annotation_col$Condition_Name
    annotation_col$Condition_Name <- NULL
    if(tpm_cutoff!= 0){
      gene_to_keep <- get_genes_with_high_tpm_in_any_group(data = data, annotation_col = annotation_col, tpm_cutoff = tpm_cutoff)
      data_filtered <- data %>%
        filter(id %in% gene_to_keep) %>%
        filter(
          (abs(pull(select(., ends_with("log2FoldChange")))) > log2FC) & 
            (pull(select(., ends_with("FDR"))) < FDR)
        )}else if (tpm_cutoff == 0 & log2FC == 0 & FDR == 1){
          data_filtered <- data
        }else{
          data_filtered <- data%>%
            filter(
              (abs(pull(select(., ends_with("log2FoldChange")))) > log2FC) & 
                (pull(select(., ends_with("FDR"))) < FDR)  
            )
          }
    
    if (nrow(data_filtered) == 0) {
      msg <- paste(
        "❌ Error: No data left after filtering. Consider adjusting log2FC, FDR or tpm_cutoff thresholds:",
        file_path_sans_ext(input),
        "Please modify the number of categories in the row or adjust log2FC, FDR or tpm_cutoff thresholds.",
        sep = "<br>"
      )
      stop(msg)
          }else if(nrow(data_filtered) < cluster_number_row){
      msg <- paste(
        "❌ Error: The number of genes remaining after filtering in the following sample is less than the number of categories in the row:",
        file_path_sans_ext(input),
        "Please modify the number of categories in the row or adjust log2FC, FDR or tpm_cutoff thresholds.",
        sep = "<br>"
      )
      stop(msg)
    }
    
    
    if (scale_function == "maxscale") {
      norm_factor <- apply(select(data_filtered, starts_with("tpm")), 1, function(x) max(abs(x), na.rm = TRUE))
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("tpm"), 
                      ~ ifelse(norm_factor == 0, 0, .x / norm_factor), 
                      .names = "scale_{.col}")) %>%
        as.data.frame()  
    } 
    
    else if (scale_function == "zscale") {
      mean_x <- apply(select(data_filtered, starts_with("tpm")), 1, mean, na.rm = TRUE)
      sd_x <- apply(select(data_filtered, starts_with("tpm")), 1, sd, na.rm = TRUE)
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("tpm"), 
                      ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                      .names = "scale_{.col}")) %>%
        as.data.frame() 
      
      min_scale <- apply(select(data_filtered, starts_with("scale_")), 1, min, na.rm = TRUE)
      max_scale <- apply(select(data_filtered, starts_with("scale_")), 1, max, na.rm = TRUE)
      max_norm <- pmax(abs(min_scale), max_scale)
      data_filtered <- data_filtered %>%
        mutate(across(starts_with("scale_"), 
                      ~ case_when(
                        .x != 0  ~ (.x / max_norm),
                        .x == 0 ~ 0
                      )))
    }
    
    upregulated_count <- sum(data_filtered %>% select(starts_with(using_package_f)) %>% apply(1, function(x) any(x > 1)))
    downregulated_count <- sum(data_filtered %>% select( starts_with(using_package_f)) %>% apply(1, function(x) any(x < -1)))
    
    df <- data_filtered %>%
      select("id", starts_with("scale_"))
    rownames(df) <- df$id
    df$id <- NULL
    df <- as.matrix(df)
    df[is.na(df)] <- 0 
    colnames(df) <- gsub("^scale_tpm_", "", colnames(df))
    valid_conditions <- rownames(annotation_col) %in% colnames(df)
    annotation_col <- annotation_col[valid_conditions, , drop = FALSE]
    annotation_colors$Group <- annotation_colors$Group[unique(annotation_col$Group)]
    
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
                          annotation_col = annotation_col,
                          annotation_colors = annotation_colors,
                          treeheight_row = 200,
                          treeheight_col = 100,
                          angle_col = 90, 
                          color = custom_color,
                          silent = TRUE
    )
    
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
    
    col_order_table$Group <- annotation_col[col_order_table$COL_NAME, "Group"]
    col_order_table <- col_order_table %>%
      select(COL_NUMBER,CLUSTER_NUMBER,Group,COL_NAME)
    
    heatmap_width_ratio <- 0.75  
    heatmap_height_ratio <- 0.9 
    
    cellwidth <- (8.27 * heatmap_width_ratio * 72 - 275- 5 * cluster_number_col) / ncol(df)    
    cellheight <- if(ncol(df)>30){(11.69 * heatmap_height_ratio * 72-95- 5 * cluster_number_row) / nrow(df)}
    else{(11.69 * heatmap_height_ratio * 72-105- 5 * cluster_number_row) / nrow(df)}
    
    df <- df[, col_order_origin] 
    colnames(df) <- seq_len(ncol(df))
    annotation_col <- annotation_col[col_order_origin, , drop = FALSE]
    rownames(annotation_col) <- as.character(seq_len(ncol(df)))
    
    HM_scale <- pheatmap(df, 
                         cellwidth = cellwidth, 
                         cellheight = cellheight,
                         fontsize_col = 100/ncol(df),
                         fontsize_row = 2,
                         border_color = NA, 
                         cutree_cols = cluster_number_col, 
                         cutree_rows = cluster_number_row, 
                         annotation_col = annotation_col,
                         annotation_colors = annotation_colors,
                         cluster_rows = TRUE, 
                         cluster_cols = TRUE,
                         treeheight_row = 100,
                         treeheight_col = 100,
                         angle_col = 90, 
                         show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                         show_rownames = FALSE,
                         color = custom_color,
                         silent = TRUE)
    
    origin_pdf_width <- (ncol(df) * 20) / 72 + 6    
    origin_pdf_height <- (nrow(df) * 2) / 72 + 7 + if(!is.na(cluster_number_row)){cluster_number_row}else{0}
    
    #add row order to data
    data_filtered <- data_filtered %>%
      mutate(heatmap_order = match(1:nrow(data_filtered), row_order_origin)) 
    first_columns <- c("heatmap_order", "id")
    data_filtered$cluster_number <- row_clusters
    first_columns <- c(first_columns, "cluster_number")
    tpm_columns <- grep("^tpm_", colnames(data_filtered), value = TRUE)  
    scale_columns <- grep("^scale_", colnames(data_filtered), value = TRUE)
    
    method_columns <- grep(paste0("^", using_package_f), colnames(data_filtered), value = TRUE)
    all_other_columns <- colnames(data_filtered)
    counts_columns <- setdiff(all_other_columns, c(first_columns, scale_columns, tpm_columns, method_columns))
    data_filtered <- data_filtered %>%
      rename_with(~ paste0("counts_", .), all_of(counts_columns))
    counts_columns <- grep("^counts_", colnames(data_filtered), value = TRUE)
    
    if (scale_function == "zscale") {
      data_filtered <- data_filtered %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_filtered <- data_filtered %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    
    
    data_filtered <- data_filtered %>%
      select(any_of(first_columns), any_of(scale_columns), any_of(counts_columns), any_of(tpm_columns), any_of(method_columns))
    
    
    # Remove fc, log2fc, FDR, pvalue, ...
    data <- data %>%
      select(-starts_with(using_package_f))
    
    tmp_path <- file.path(output_folder, "tmp/")
    long_filepath <- paste0(tmp_path, file_path_sans_ext(basename(input)), "_long.pdf")
    short_filepath <- paste0(tmp_path, file_path_sans_ext(basename(input)), "_short.pdf")
    
    
    # Save a
    pdf(long_filepath, width = origin_pdf_width, height = origin_pdf_height)
    grid.newpage()
    title <- paste0(folder_name,"\n",file_path_sans_ext(basename(input)), 
                    "\nTREATMENT: ", condition_2, 
                    ",\tCONTROL: ", condition_1,
                    ",\nUPREGULATED: ", upregulated_count,
                    ",\tDOWNREGULATED: ", downregulated_count,
                    "\nSCALE FUNCTION: ", scale_function, ",\tTPM CUTOFF: ", tpm_cutoff, 
                    "\nLOG2FC: ", log2FC, ",\tFDR: ", FDR)
    grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
    pushViewport(viewport(x = 0.5, y = 0.49, width = 0.9, height = 0.7)) 
    grid.draw(HM_origin$gtable) 
    popViewport()
    
    dev.off()
    
    
    # Save b
    pdf(short_filepath, width = 8.27, height = 11.69)
    
    title <- paste0(folder_name,"\n",file_path_sans_ext(basename(input)))
    grid.newpage()
    grid.text(title, x = 0.5, y = 0.99, gp = gpar(fontsize = 10, fontface = "bold"), just = c("center","top"))
    
    pushViewport(viewport(x = 0.35, y = 0.5, width = heatmap_width_ratio, height = heatmap_height_ratio))
    print(HM_scale, newpage = FALSE)
    
    table_theme <- ttheme_default(
      core = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                  padding = unit(c(1.5, 1.5), "mm")), 
      colhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                     padding = unit(c(1.5, 1.5), "mm")), 
      rowhead = list(fg_params = list(fontsize = if(ncol(df) < 100){6}else{4}, lineheight = 0), 
                     padding = unit(c(1.5, 1.5), "mm"))
    ) 
    grid.text("Gene clusters are numbered from\ntop to bottom", x = 1, y = 0.95, gp = gpar(fontsize = 8, col = "red"), just = c("left", "top"))
    grid.text("First column: Index (->)\nSecond column: Group cluster number\nThird column: Sample Group", x = 1, y = 0.01, gp = gpar(fontsize = 8, col = "red"), just = c("left", "bottom"))
    
    info <- paste0("TREATMENT: ", condition_2, "\nCONTROL: ", condition_1,
                   "\n\nUPREGULATED: ", upregulated_count, "\nDOWNREGULATED: ",downregulated_count,
                   "\n\nSCALE FUNCTION: ", scale_function, "\nTPM CUTOFF: ", tpm_cutoff, "\nLOG2FC: ", log2FC,
                   "\nFDR: ", FDR)
    grid.text(info, x = 1, y = 0.9, gp = gpar(fontsize = 8), just = c("left","top"))
    pushViewport(viewport(x = 1, y = 0.45, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
    
    grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
    popViewport()
    
    dev.off()
    
    
    
    # Save c
    pattern <- "_(ExonCollapsed|FullGeneBody)_.*$"
    match <- str_extract(file_path_sans_ext(basename(input)), pattern)
    if (is.na(match)) {
      match <- "UnknownGeneType"
    } else {
      match <- substr(match, 2, nchar(match)) 
    }
    output_filename <- file.path(output_folder, paste0("HeatmapData_",file_path_sans_ext(basename(folder_name)), "_", match, "_log2FC_", log2FC, "_FDR_", FDR, ".xlsx"))
    write.xlsx(data_filtered, output_filename, rowNames = FALSE)
  }
  
  return(list(data = data,
              group_data = group_data,
              id_list = data_filtered$id,
              upregulated_count = upregulated_count,
              downregulated_count = downregulated_count,
              col_order_table = col_order_table,
              long_filepath = long_filepath,
              short_filepath = short_filepath))
  
}




heatmap_integrated_09d <- function(union_id_list, 
                                   all_data,
                                   save_name,
                                   log2FC,
                                   FDR,
                                   tpm_cutoff, 
                                   scale_function, 
                                   sample_type, 
                                   output_folder, 
                                   custom_color, 
                                   cluster_number_row, 
                                   cluster_number_col,
                                   order_for_noclustering,
                                   annotation_colors){
  filtered_data_list <- lapply(all_data, function(df) df[df$id %in% union_id_list, ])
  data_integrated <- filtered_data_list[[1]]
  
  for (i in 2:length(filtered_data_list)) {
    new_data <- filtered_data_list[[i]]
    new_cols <- setdiff(colnames(new_data), colnames(data_integrated))
    
    if (length(new_cols) > 0) { 
      data_integrated <- merge(data_integrated, new_data[, c("id", new_cols)], by = "id", all = TRUE)
    }
  }
  data_integrated <- data_integrated[, unique(colnames(data_integrated))]
  
  annotation_col <- order_for_noclustering
  rownames(annotation_col) <- annotation_col$Condition_Name
  annotation_col$Condition_Name <- NULL
  
  
  if (scale_function == "maxscale") {
    norm_factor <- apply(select(data_integrated, starts_with("tpm")), 1, function(x) max(abs(x), na.rm = TRUE))
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("tpm"), 
                    ~ ifelse(norm_factor == 0, 0, .x / norm_factor), 
                    .names = "scale_{.col}")) %>%
      as.data.frame()  
  }else if (scale_function == "zscale") {
    mean_x <- apply(select(data_integrated, starts_with("tpm")), 1, mean, na.rm = TRUE)
    sd_x <- apply(select(data_integrated, starts_with("tpm")), 1, sd, na.rm = TRUE)
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("tpm"), 
                    ~ ifelse(sd_x == 0, 0, (.x - mean_x) / sd_x), 
                    .names = "scale_{.col}")) %>%
      as.data.frame()  
    
    min_scale <- apply(select(data_integrated, starts_with("scale_")), 1, min, na.rm = TRUE)
    max_scale <- apply(select(data_integrated, starts_with("scale_")), 1, max, na.rm = TRUE)
    max_norm <- pmax(abs(min_scale), max_scale)
    data_integrated <- data_integrated %>%
      mutate(across(starts_with("scale_"), 
                    ~ case_when(
                      .x != 0  ~ (.x / max_norm),
                      .x == 0 ~ 0
                    )))
  }
  
  df <- data_integrated %>%
    select("id", starts_with("scale_"))
  rownames(df) <- df$id
  df$id <- NULL
  df <- as.matrix(df)
  df[is.na(df)] <- 0 
  colnames(df) <- gsub("^scale_tpm_", "", colnames(df))
  colnames(df) <- gsub("^mean_", "", colnames(df))
  #sort nonclustering col order
  if (!is.null(order_for_noclustering)) {
    target_order <- order_for_noclustering$Condition_Name
    matched_cols <- colnames(df)
    matched_index <- sapply(target_order, function(cond) {
      idx <- which(grepl(cond, matched_cols))
      if (length(idx) == 0) return(NA)
      idx[1]
    })
    matched_index <- matched_index[!is.na(matched_index)]
    df <- df[, matched_index, drop = FALSE]
  }
  
  valid_conditions <- rownames(annotation_col) %in% colnames(df)
  annotation_col <- annotation_col[valid_conditions, , drop = FALSE]
  annotation_colors$Group <- annotation_colors$Group[unique(annotation_col$Group)]
  
  
  #doing heatmap
  HM_origin <- pheatmap(df, 
                        cellwidth = 20, 
                        cellheight = 2 , 
                        fontsize_col = 10,
                        fontsize_row = 2,
                        border_color = NA, 
                        cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                        cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                        cluster_rows = TRUE, 
                        cluster_cols = TRUE,
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors,
                        treeheight_row = 200,
                        treeheight_col = 80,
                        angle_col = 90, 
                        color = custom_color,
                        silent = TRUE
  )
  
  row_order_origin <- HM_origin$tree_row$order
  col_order_origin <- HM_origin$tree_col$order

  row_clusters <- cutree(HM_origin$tree_row, k = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]})
  col_clusters <- cutree(HM_origin$tree_col, k = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]})

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
  
  col_order_table$Group <- annotation_col[col_order_table$COL_NAME, "Group"]
  col_order_table <- col_order_table %>%
    select(COL_NUMBER,CLUSTER_NUMBER,Group,COL_NAME)
  
  heatmap_width_ratio <- 0.75  
  heatmap_height_ratio <- 0.9 
  
  cellwidth <- (8.27 * heatmap_width_ratio * 72-275 - 5 *if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}) / ncol(df)    
  cellheight <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72-100 - 5*if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  else{(11.69 * heatmap_height_ratio * 72- 110 - 5 * if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  
  cellwidth_noclustering <- (8.27 * heatmap_width_ratio * 72-260) / ncol(df)    
  cellheight_noclustering <- if(ncol(df)>50){(11.69 * heatmap_height_ratio * 72 -15 - 5*if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  else{(11.69 * heatmap_height_ratio * 72- 25 - 5 * if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}) / nrow(df)}
  
  
  HM_noclustering_origin <- pheatmap(df, 
                                     cellwidth = 20, 
                                     cellheight = 2 , 
                                     fontsize_col = 10,
                                     fontsize_row = 2,
                                     border_color = NA, 
                                     cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                                     cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                                     cluster_rows = TRUE, 
                                     cluster_cols = FALSE,
                                     annotation_col = annotation_col,
                                     annotation_colors = annotation_colors,
                                     treeheight_row = 200,
                                     angle_col = 90, 
                                     color = custom_color,
                                     silent = TRUE
  )
  
  col_order_table_noclustering <- data.frame(
    COL_NUMBER = seq_along(colnames(df)),  
    COL_NAME = colnames(df),
    stringsAsFactors = FALSE
  )
  col_order_table_noclustering$Group <- annotation_col[col_order_table_noclustering$COL_NAME, "Group"]
  col_order_table_noclustering <- col_order_table_noclustering %>%
    select(COL_NUMBER,Group,COL_NAME)
  
  
  colnames(df) <- seq_len(ncol(df))
  rownames(annotation_col) <- as.character(seq_len(ncol(df)))
  
  HM_noclustering_scale <- pheatmap(df, 
                                    cellwidth = cellwidth_noclustering, 
                                    cellheight = cellheight_noclustering,
                                    fontsize_col = 100/ncol(df),
                                    fontsize_row = 2,
                                    border_color = NA, 
                                    cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                                    cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                                    cluster_rows = TRUE, 
                                    cluster_cols = FALSE,
                                    annotation_col = annotation_col,
                                    annotation_colors = annotation_colors,
                                    treeheight_row = 130,
                                    angle_col = 90, 
                                    show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                                    show_rownames = FALSE,
                                    color = custom_color,
                                    silent = TRUE)
  
  df <- df[, col_order_origin] 
  annotation_col <- annotation_col[col_order_origin, , drop = FALSE]
  #FIX HERE
  colnames(df) <- seq_len(ncol(df))
  rownames(annotation_col) <- as.character(seq_len(ncol(df)))
  
  HM_scale <- pheatmap(df, 
                       cellwidth = cellwidth, 
                       cellheight = cellheight,
                       fontsize_col = 100/ncol(df),
                       fontsize_row = 2,
                       border_color = NA, 
                       cutree_cols = if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}, 
                       cutree_rows = if(sample_type == "individuals"){cluster_number_row[length(cluster_number_row)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}, 
                       cluster_rows = TRUE, 
                       cluster_cols = TRUE,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       treeheight_row = 130,
                       treeheight_col = 80,
                       angle_col = 90, 
                       show_colnames = if(ncol(df)>50){FALSE}else{TRUE},
                       show_rownames = FALSE,
                       color = custom_color,
                       silent = TRUE)
  
  
  #PROCESSING FOR DATA_INTEGRATED
  if(sample_type == "individuals"){
    data_integrated <- data_integrated %>%
      mutate(heatmap_order = match(1:nrow(data_integrated), row_order_origin)) 
    first_columns <- c("heatmap_order", "id")  
    data_integrated$cluster_number <- row_clusters
    first_columns <- c(first_columns, "cluster_number")
    scale_columns <- grep("^scale_", colnames(data_integrated), value = TRUE)  
    tpm_columns <- grep("^tpm_", colnames(data_integrated), value = TRUE)  
    all_other_columns <- colnames(data_integrated)
    counts_columns <- setdiff(all_other_columns, c(first_columns, scale_columns, tpm_columns))
    data_integrated <- data_integrated %>%
      rename_with(~ paste0("counts_", .), all_of(counts_columns))
    counts_columns <- grep("^counts_", colnames(data_integrated), value = TRUE)
    
    if (scale_function == "zscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    if (!is.null(order_for_noclustering)) {
      scale_columns <- match_column(order_for_noclustering, scale_columns)
      tpm_columns <- match_column(order_for_noclustering, tpm_columns)
      counts_columns <- match_column(order_for_noclustering, counts_columns)}
    data_integrated <- data_integrated %>%
      select(all_of(first_columns), all_of(scale_columns), all_of(counts_columns), all_of(tpm_columns))}
  
  else if(sample_type == "groups"){
    data_integrated <- data_integrated %>%
      mutate(heatmap_order = match(1:nrow(data_integrated), row_order_origin)) 
    first_columns <- c("heatmap_order", "id") 
    data_integrated$cluster_number <- row_clusters
    first_columns <- c(first_columns, "cluster_number")
    scale_columns <- grep("^scale_", colnames(data_integrated), value = TRUE)  
    tpm_columns <- grep("^tpm_", colnames(data_integrated), value = TRUE) 
    
    if (scale_function == "zscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("z", .), all_of(scale_columns))
      scale_columns <- paste0("z", scale_columns)
    } else if (scale_function == "maxscale") {
      data_integrated <- data_integrated %>%
        rename_with(~ paste0("max", .), all_of(scale_columns))
      scale_columns <- paste0("max", scale_columns)
    }
    if (!is.null(order_for_noclustering)) {
      scale_columns <- match_column(order_for_noclustering, scale_columns)
      tpm_columns <- match_column(order_for_noclustering, tpm_columns)}
    data_integrated <- data_integrated %>%
      select(all_of(first_columns), all_of(scale_columns), all_of(tpm_columns))
    
  }
  
  
  
  origin_pdf_width <- (ncol(df) * 20) / 72 + 6 + 0.3*if(!is.na(if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]})){if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_row[length(cluster_number_row)]}}else{0}
  origin_pdf_height <- (nrow(df) * 2) / 72 + 10 + 0.3*if(!is.na(if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]})){if(sample_type == "individuals"){cluster_number_col[length(cluster_number_col)-1]}else if(sample_type == "groups"){cluster_number_col[length(cluster_number_col)]}}else{0}
  origin_pdf_width_noclustering <- (ncol(df) * 20) / 72 + 6
  origin_pdf_height_noclustering <- (nrow(df) * 2) / 72 + 8
  
  if(sample_type == "individuals"){
    cluster_long_path <- paste0("Heatmap_AllIndividuals_SamplesClustered_Long_", save_name, ".pdf")
    nocluster_long_path <- paste0("Heatmap_AllIndividuals_NoSamplesClustering_Long_", save_name, ".pdf")
    cluster_short_path <- paste0("Heatmap_AllIndividuals_SamplesClustered_Short_", save_name, ".pdf")
    nocluster_short_path <- paste0("Heatmap_AllIndividuals_NoSamplesClustering_Short_", save_name, ".pdf")
    heatmapdata_path <- paste0("HeatmapData_AllIndividuals_",save_name,"_log2FC_",log2FC,"_FDR_",FDR,".xlsx")}
  else if(sample_type == "groups"){
    cluster_long_path <- paste0("Heatmap_AllGroups_SamplesClustered_Long_", save_name, ".pdf")
    nocluster_long_path <- paste0("Heatmap_AllGroups_NoSamplesClustering_Long_", save_name, ".pdf")
    cluster_short_path <- paste0("Heatmap_AllGroups_SamplesClustered_Short_", save_name, ".pdf")
    nocluster_short_path <- paste0("Heatmap_AllGroups_NoSamplesClustering_Short_", save_name, ".pdf")
    heatmapdata_path <- paste0("HeatmapData_AllGroups_",save_name,"_log2FC_",log2FC,"_FDR_",FDR,".xlsx")}
  
  # SAVE integrated long cluster
  pdf(file.path(output_folder, cluster_long_path), width = origin_pdf_width, height = origin_pdf_height)
  grid.newpage()
  title <- paste0(save_name, "_SampleClustered_",sample_type,
                  "\nSCALE FUNCTION: ", scale_function,",\tTPM CUTOFF: ", tpm_cutoff, 
                  "\nLOG2FC: ", log2FC,",\tFDR: ", FDR,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_origin$gtable) 
  popViewport()
  dev.off()
  
  # SAVE integrated long noclustered
  pdf(file.path(output_folder, nocluster_long_path), width = origin_pdf_width_noclustering, height = origin_pdf_height_noclustering)
  grid.newpage()
  title <- paste0(save_name, "_NoSampleClustering_",sample_type,
                  "\nSCALE FUNCTION: ", scale_function,",\tTPM CUTOFF: ", tpm_cutoff, 
                  "\nLOG2FC: ", log2FC,",\tFDR: ", FDR,
                  "\nGENE NUMBER: ", nrow(df))
  grid.text(title, x = 0.5, y = 0.999, gp = gpar(fontsize = 12, fontface = "bold"), just = c("center","top"))
  pushViewport(viewport(x = 0.5, y = 0.9, width = 0.9, height = 0.8, just = c("center","top"))) 
  grid.draw(HM_noclustering_origin$gtable) 
  popViewport()
  dev.off()
  
  # SAVE integrated short noclustering
  pdf(file.path(output_folder, nocluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(save_name,"_NoSampleClustering_",sample_type)
  if (nchar(title) > 80) {
    title <- paste0(substr(title, 1, 80), "-\n", substr(title, 81, nchar(title)))
  }
  info <- paste0("\nSCALE FUNCTION: ", scale_function,",\nTPM CUTOFF: ", tpm_cutoff, 
                 "\nLOG2FC: ", log2FC,",\nFDR: ", FDR,
                 "\nGENE NUMBER: ", nrow(df))
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
  pushViewport(viewport(x = 1.025, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table_noclustering, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  
  # short clustered
  pdf(file.path(output_folder, cluster_short_path), width = 8.27, height = 11.69)
  title <- paste0(save_name,"_SampleClustered_",sample_type)
  if (nchar(title) > 80) {
    title <- paste0(substr(title, 1, 80), "-\n", substr(title, 81, nchar(title)))
  }
  info <- paste0("\nSCALE FUNCTION: ", scale_function,",\nTPM CUTOFF: ", tpm_cutoff, 
                 "\nLOG2FC: ", log2FC,",\nFDR: ", FDR,
                 "\nGENE NUMBER: ", nrow(df))
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
  pushViewport(viewport(x = 1.025, y = 0.5, width = 1-heatmap_width_ratio, height = 0.15, just = c("left","top")))
  grid.table(col_order_table, rows = NULL, cols = NULL, theme = table_theme)
  popViewport()
  dev.off()
  
  write.xlsx(data_integrated, file.path(output_folder, heatmapdata_path), rowNames = FALSE)
  #RETURN
  return(list(col_order_table_noclustering = col_order_table_noclustering, col_order_table = col_order_table))
}




heatmap_analysis_09d <- function(
    input_list,
    edited_gene_list,
    condition_1,
    condition_2,
    output_folder,
    folder_name_list,  # New parameter
    custom_color = NULL,
    using_package_f,
    scale_function,
    log2FC,
    FDR,
    tpm_cutoff,
    cluster_number_row,
    cluster_number_col,
    groups_order_for_noclustering,
    individuals_order_for_noclustering,
    annotation_colors) {
  # CHECKING INPUT
  if (missing(input_list) || is.null(input_list)) {
    stop("Error: 'input_list' is mandatory and must be a file path or a list of file paths.")
  }
  if (!is.character(input_list)) {
    stop("Error: 'input' must be a list of character vectors.")
  }
  if (any(!file.exists(input_list))) {
    stop("Error: One or more input files do not exist. Please check the file paths.")
  }
  # CHECKING OUTPUT_FOLDER
  if (missing(output_folder) || is.null(output_folder)) {
    stop("Error: 'output_folder' is mandatory and must be a valid folder path.")
  }
  if (!dir.exists(output_folder)) {
    stop("Error: The specified 'output_folder' directory does not exist.")
  }
  # CHECKING log2FC
  if (!is.numeric(log2FC) || log2FC < 0) {
    stop("Error: 'log2FC' must be a numeric value greater than or equal 0.")
  }
  # CHECKING FDR
  if (!is.numeric(FDR) || FDR <= 0 || FDR > 1) {
    stop("Error: 'FDR' must be a numeric value between 0 and 1.")
  }
  # CHECKING USING_PACKAGE_F
  if (!using_package_f %in% c("EdgeR", "DESeq2")) {
    stop("Error: 'using_package_f' must be either 'EdgeR' or 'DESeq2'.")
  }
  # CHECK cluster_number_row
  if (is.null(cluster_number_row) || length(cluster_number_row) == 0) {
    cluster_number_row <- rep(1, length(input_list) + 2)
  } else if (length(cluster_number_row) != length(input_list) + 2) {
    stop("'cluster_number_row' must be a numeric vector of length equal to the number of input_list files +2.\n e.g: cluster_number_row = c(input1_ClusterNumber, input2_ClusterNumber, ..., integrated_ClusterNumber)")
  }
  # CHECK cluster_number_col
  if (is.null(cluster_number_col) || length(cluster_number_col) == 0) {
    cluster_number_col <- rep(1, length(input_list) + 2)
  } else if (length(cluster_number_col) != length(input_list) + 2) {
    stop("'cluster_number_col' must be a numeric vector of length equal to the number of input_list files +2.\n e.g: cluster_number_col = c(input1_ClusterNumber, input2_ClusterNumber, ..., integrated_ClusterNumber)")
  }
  
  # initialize storing value
  all_col_order_tables <- list()
  summary_data <- data.frame(File = character(), Upregulated = integer(), Downregulated = integer(), stringsAsFactors = FALSE)
  long_pdfs <- c()
  short_pdfs <- c()
  all_id_lists <- list()
  all_data <- list()
  all_group_data <- list()
  # Create folder
  tmp_path <- file.path(output_folder, "tmp/")
  dir.create(tmp_path, recursive = TRUE, showWarnings = FALSE)
  for (step in seq_along(input_list)) {
    
    result <- pair_heatmap_09d(
      input = input_list[step], 
      edited_gene_list,
      condition_1 = condition_1[step], 
      condition_2 = condition_2[step], 
      output_folder = output_folder, 
      custom_color = custom_color, 
      using_package_f = using_package_f, 
      scale_function = scale_function,
      log2FC = log2FC, 
      FDR = FDR, 
      tpm_cutoff = tpm_cutoff,
      cluster_number_row = cluster_number_row[step], 
      cluster_number_col = cluster_number_col[step],
      order_for_noclustering = individuals_order_for_noclustering,
      annotation_colors = annotation_colors,
      folder_name = folder_name_list[step] 
    )
    
    all_col_order_tables[[folder_name_list[step]]] <- result[[6]]
    long_pdfs <- c(long_pdfs, result[[7]])
    short_pdfs <- c(short_pdfs, result[[8]])
    all_id_lists[[step]] <- result[[3]]
    all_data[[step]] <- result[[1]]
    all_group_data[[step]] <- result[[2]]
    
    # get upregulated and downregulated number
    summary_data <- rbind(summary_data, data.frame(
      File = folder_name_list[step],
      `Condition_1(CONTROL)` = condition_1[step],
      `Condition_2(TREATMENT)` = condition_2[step],
      Downregulated = result[[4]],
      Upregulated = result[[5]]
    ))
  }
  
  
  union_id_list <- unique(unlist(all_id_lists))
  # GET G code
  all_g_codes <- unique(unlist(lapply(all_col_order_tables, function(tbl) {
    g_codes <- gsub(".*_(G[0-9]+)_.*", "\\1", tbl$COL_NAME)  
    g_codes[grepl("^G[0-9]+$", g_codes)]  
  })))
  save_name <- paste(sort(all_g_codes), collapse = "_")
  # Combined pdf long and short
  long_pdf_output <- file.path(output_folder, paste0("Heatmap_AllPairwiseComparisons_Long_", save_name, ".pdf"))
  short_pdf_output <- file.path(output_folder, paste0("Heatmap_AllPairwiseComparisons_Short_", save_name, ".pdf"))
  system(paste("pdftk", paste(long_pdfs, collapse = " "), "cat output", long_pdf_output))
  system(paste("pdftk", paste(short_pdfs, collapse = " "), "cat output", short_pdf_output))
  unlink(tmp_path, recursive = TRUE, force = TRUE)
  
  #HERE
  individuals_result <- heatmap_integrated_09d(union_id_list = union_id_list,
                                               all_data = all_data,
                                               save_name = save_name,
                                               tpm_cutoff = tpm_cutoff,
                                               log2FC = log2FC,
                                               FDR = FDR,
                                               scale_function = scale_function,
                                               sample_type = "individuals",
                                               output_folder = output_folder,
                                               custom_color = custom_color,
                                               cluster_number_row = cluster_number_row,
                                               cluster_number_col = cluster_number_col,
                                               annotation_colors = annotation_colors,
                                               order_for_noclustering = individuals_order_for_noclustering)
  
  groups_result <- heatmap_integrated_09d(union_id_list = union_id_list, 
                                          all_data = all_group_data, 
                                          save_name = save_name,
                                          tpm_cutoff = tpm_cutoff,
                                          log2FC = log2FC,
                                          FDR = FDR,
                                          scale_function = scale_function,
                                          sample_type = "groups", 
                                          output_folder = output_folder, 
                                          custom_color = custom_color, 
                                          cluster_number_row = cluster_number_row, 
                                          cluster_number_col = cluster_number_col,
                                          annotation_colors = annotation_colors,
                                          order_for_noclustering = groups_order_for_noclustering)
  
  summary_filepath <- file.path(output_folder, paste0("Summary_DEGs_SampleOrder_", save_name, "_log2FC_", log2FC, "_FDR_", FDR, ".xlsx"))
  excel_sheets <- list(
    summary_DEG = summary_data,  
    Individuals_NoCluster_ColOrder = individuals_result[["col_order_table_noclustering"]],  
    Individuals_Cluster_ColOrder = individuals_result[["col_order_table"]],
    Groups_NoCluster_ColOrder = groups_result[["col_order_table_noclustering"]],  
    Groups_Cluster_ColOrder = groups_result[["col_order_table"]]
  )
  
  shortened_sheet_name <- function(name) {
    name <- paste0(file_path_sans_ext(basename(name)),"_Order")
    if (nchar(name) > 31) {
      return(substr(name, nchar(name) - 30, nchar(name)))  
    } else {
      return(name)
    }
  }
  shortened_names <- lapply(names(all_col_order_tables), shortened_sheet_name)
  names(all_col_order_tables) <- shortened_names
  for (sheet_name in names(all_col_order_tables)) {
    excel_sheets[[sheet_name]] <- all_col_order_tables[[sheet_name]]
  }
  write.xlsx(excel_sheets, file = summary_filepath, rowNames = FALSE)
  return(summary_data)
}

do_heatmap_09d <- function(
    input_folder,
    output_folder,
    input_folder_list,
    edited_gene_list = NULL,
    segex_transfer = FALSE,
    cluster_number_row = NULL,
    cluster_number_col = NULL,
    groups_order = NULL,
    custom_color = NULL,
    genic_region = "ExonCollapsed", # 'ExonCollapsed', 'FullGeneBody'
    using_package_f = "EdgeR",  # 'EdgeR', or 'DESeq2'
    scale_function = "both", # "maxscale", "zscale","both"
    abs_log2FC = c(1,2),
    FDR = 0.05,
    tpm_cutoff = c(0,1)) {
  # get groups_order_for_noclustering and individuals_order_for_noclustering
  label_path <- file.path(input_folder, "00_Setup_Pipeline", "Sample_Labels.txt")
  sample_label <- read.table(label_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  
  groups_order_for_noclustering <- data.frame(Condition_Name = groups_order, stringsAsFactors = FALSE)
  
  groups_order_for_noclustering <- groups_order_for_noclustering %>%
    left_join(sample_label, by = "Condition_Name")
  
  individuals_order_for_noclustering <- groups_order_for_noclustering %>%
    distinct(Group, Condition_Name, Sample_ID) %>%
    mutate(Condition_Name = paste0(Condition_Name, "_", Sample_ID)) %>%
    select(Group, Condition_Name)
  
  groups_order_for_noclustering <- groups_order_for_noclustering %>%
    distinct(Group, Condition_Name) %>%
    as.data.frame()
  
  duplicated_names <- groups_order_for_noclustering %>%
    group_by(Condition_Name) %>%
    summarise(n_groups = n_distinct(Group), .groups = "drop") %>%
    filter(n_groups > 1)
  
  if (nrow(duplicated_names) > 0) {
    full_message <- paste0(
      "❌ Error: The following Condition_Name(s) are assigned to multiple groups:\n",
      paste0(" - ", duplicated_names$Condition_Name, collapse = "\n")
    )
    stop(full_message)
  }
  
  
  input_folder <- sub("/+$", "", input_folder)
  output_folder <- sub("/+$", "", output_folder)
  input_folder_list <- file.path(input_folder, input_folder_list)
  
  # default value of custom_color
  if (is.null(custom_color) || length(custom_color) == 0 || length(custom_color) == 1) {
    custom_color <- c("red", "yellow", "#3fdc04")
  }
  custom_color <- colorRampPalette(custom_color)(200)
  
  if (!dir.exists(output_folder)) {
    stop("output_folder doesn't exist.")
  }
  
  # get color
  
  group_list <- unique(sample_label$Group)
  group_colors <- color_selected(length(group_list))
  names(group_colors) <- group_list
  annotation_colors <- list(Group = group_colors)
  
  scale_function <- if (scale_function == "both") c("maxscale", "zscale") else c(scale_function)
  
  # Get input_file_list, condition_1_list, condition_2_list (input_folders)
  input_file_list <- c()
  condition_1_list <- c()
  condition_2_list <- c()
  
  for (input in input_folder_list) {
    # Find Condition_1.txt and Condition_2.txt
    condition_1_file <- file.path(input, "Condition_1.txt")
    condition_2_file <- file.path(input, "Condition_2.txt")
    
    # Ensure both condition files exist
    if (!file.exists(condition_1_file)) {
      warning(paste("Warning: Condition_1.txt not found in", input, ". Skipping."))
      next
    }
    if (!file.exists(condition_2_file)) {
      warning(paste("Warning: Condition_2.txt not found in", input, ". Skipping."))
      next
    }
    
    # Read condition files
    condition_1_data <- read.table(condition_1_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    condition_2_data <- read.table(condition_2_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    # Ensure 'Description' column exists
    if (!"Description" %in% colnames(condition_1_data)) {
      warning(paste("Warning: 'Condition_1.txt' in", input, "is missing 'Description' column. Skipping."))
      next
    }
    if (!"Description" %in% colnames(condition_2_data)) {
      warning(paste("Warning: 'Condition_2.txt' in", input, "is missing 'Description' column. Skipping."))
      next
    }
    
    # Extract conditions
    condition_1 <- condition_1_data$Description[1]
    condition_2 <- condition_2_data$Description[1]
    
    # Locate the genic_region-specific subdirectory
    matching_dirs <- list.dirs(input, full.names = TRUE, recursive = FALSE)
    matching_dirs <- matching_dirs[grepl(paste0(genic_region, "$"), matching_dirs)]
    if (length(matching_dirs) == 0) {
      warning(paste("Warning: No subdirectories ending in", genic_region, "found in", input))
      next
    }
    
    # Collect all DiffExp files from the matching directories
    for (dir_name in matching_dirs) {
      diffexp_files <- list.files(dir_name, pattern = "^DiffExp", full.names = TRUE, recursive = FALSE)
      
      # Store input files and corresponding conditions
      for (file in diffexp_files) {
        input_file_list <- c(input_file_list, file)
        condition_1_list <- c(condition_1_list, condition_1)
        condition_2_list <- c(condition_2_list, condition_2)
      }
    }
  }
  if(!is.null(edited_gene_list)){
    edited_gene_list <- read_excel(edited_gene_list, col_names = "1")
    if(nrow(edited_gene_list) == 0 || all(is.na(edited_gene_list[[1]]))) {
      edited_gene_list <- character(0)
    } else {
      if(segex_transfer == TRUE){
        transfer_mapping <- read_excel("example_data/transfer_mapping.xlsx")
        mapping_dict <- setNames(transfer_mapping[[2]], transfer_mapping[[3]])
        edited_gene_list_converted <- mapping_dict[edited_gene_list[[1]]]
        edited_gene_list_converted <- na.omit(edited_gene_list_converted)
        edited_gene_list <- data.frame(edited_gene_list_converted)
        colnames(edited_gene_list) <- "1"
      }
      edited_gene_list <- c(na.omit(edited_gene_list[[1]]))
    }
  }
  # Record master information
  master_summary <- data.frame()
  
  for(sf in scale_function){
    for(lfc in abs_log2FC){
      for(fdr in FDR){
        for(tpm_co in tpm_cutoff){
          new_output_folder <- file.path(output_folder, paste0(genic_region, "_", sf, "_log2FC", lfc,"_FDR", fdr, "_TpmCutoff", tpm_co))
          if (!dir.exists(new_output_folder)) {
            dir.create(new_output_folder, recursive = TRUE)}
          
          summary_data <- heatmap_analysis_09d(
            input_list = input_file_list,
            edited_gene_list = edited_gene_list,
            condition_1 = condition_1_list,
            condition_2 = condition_2_list,
            output_folder = new_output_folder,
            folder_name_list = input_folder_list,  
            custom_color = custom_color,
            using_package_f = using_package_f,
            scale_function = sf,
            log2FC = lfc,
            FDR = fdr,
            tpm_cutoff = tpm_co,
            cluster_number_row = cluster_number_row,
            cluster_number_col = cluster_number_col,
            groups_order_for_noclustering = groups_order_for_noclustering,
            individuals_order_for_noclustering = individuals_order_for_noclustering,
            annotation_colors = annotation_colors
          )
            summary_data$genic_region <- genic_region
            summary_data$abs_log2FC <- lfc
            summary_data$FDR <- fdr
            summary_data$tpm_cutoff <- tpm_co
            master_summary <- rbind(master_summary, summary_data)}
      }
    }
  }
  master_summary_path <- file.path(output_folder, "Master_Summary_DEG.xlsx")
  write.xlsx(master_summary, master_summary_path, rowNames = FALSE)
}
