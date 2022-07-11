#!/usr/bin/env Rscript
# --------------
# Date:  2022-06-19 21:30:10
# Author:Dian Li
# Last update: 2022-07-06

suppressMessages(library(pheatmap))
suppressMessages(library(optparse))


# ======================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory that contains calculated expression files", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="a specific file to be plotted", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("--exclude"), type="character", default=NULL, 
              help="file that to be excluded from being combined", metavar="character"),
  make_option(c("--show_row_labels"), type="logical", default=FALSE, 
              help="logic variable indicating whether to show row labels", metavar="TRUE or FALSE")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

call_pheatmap <- function(mat, output_name, height, width, ...){
  png(paste0(output_name), height = height, width = width, units = "in", res = 300)
  pheatmap(mat, na_col = "grey", ...)
  dev.off()
}


plot_heatmap <- function(exprsn, output_name, ...){
  
  cell_line_ind = which(colnames(exprsn) == "cell_line")
  lineage_ind = which(colnames(exprsn) == "lineage")
  ind_rm = c(cell_line_ind, lineage_ind)
  ind_rm <- ind_rm[!is.na(ind_rm)]
  mat = exprsn[, -ind_rm, drop = FALSE]
  # name empty cell line names
  rownames_prep = exprsn$cell_line
  empty_index = which(rownames_prep == "")
  rownames_prep[empty_index] = paste0("filled_", 1:length(empty_index))
  
  # rename duplicated cell line names, if exists
  dup_rownames = rownames_prep[duplicated(rownames_prep)]
  if (length(dup_rownames) > 0){
    sapply(dup_rownames, FUN = function(x){
      dup_index = which(rownames_prep == x)
      rownames_prep[dup_index] <<- paste0(x, 1:length(dup_index))
    })
  }
  
  rownames(mat) = rownames_prep
  
  if (dim(mat)[1] <= 40){
    height = 4+0.1*dim(mat)[1]
    width = 5+0.2*dim(mat)[2]
  } else if (dim(mat)[1] <= 100 & dim(mat)[1] > 40){
    height = 6+0.08*dim(mat)[1]
    width = 6+0.2*dim(mat)[2]
  } else if (dim(mat)[1] <= 200 & dim(mat)[1] > 100){
    height = 6+0.001*dim(mat)[1]
    width = 4+0.2*dim(mat)[2]
  } else {
    height = 6+0.0007*dim(mat)[1]
    width = 4+0.2*dim(mat)[2]
  }
  
  # if the lineage column exists in the expression matrix
  # then generate an annotation_row dataframe
  if ("lineage" %in% colnames(exprsn)){
    annotation_row = data.frame("Lineage" = exprsn$lineage)
    rownames(annotation_row) = rownames(mat)
    call_pheatmap(mat, output_name, height, width, annotation_row = annotation_row, angle_col = 45, ...)
    
  } else{
    call_pheatmap(mat, output_name, height, width, angle_col = 45, ...)
  }
  
}


# ========================================= #
# if the input is a single file

if (!is.null(opt$f)){
  exprsn = read.csv(opt$f, check.names = F)
  output_name = file.path(opt$o, paste0('heatmap_', gsub('.csv', '', basename(opt$f)),'.png'))
  metric_cols = ifelse(dim(exprsn)[2]>4, TRUE, FALSE)
  metric_rows = ifelse(dim(exprsn)[1]>4, TRUE, FALSE)
  
  plot_heatmap(exprsn, output_name, cluster_rows = metric_rows, cluster_cols = metric_cols)
}

# ========================================= #
# if the input is a specific folder
if (!is.null(opt$d)){
  files = dir(opt$d, pattern = ".csv")
  # plot each single file
  for (eachFile in files){
    if (eachFile != basename(opt$exclude)){
      exprsn = read.csv(file.path(opt$d, eachFile), check.names = F)
      output_name = file.path(opt$o, paste0('heatmap_', gsub('.csv', '', eachFile),'.png'))
      metric_cols = ifelse(dim(exprsn)[2]>4, TRUE, FALSE)
      metric_rows = ifelse(dim(exprsn)[1]>4, TRUE, FALSE)
      
      if (dim(exprsn)[1]>100){
        show_metric = FALSE
      } else {
        show_metric = opt$show_row_labels
      }
      plot_heatmap(exprsn, output_name, show_rownames = show_metric, cluster_rows = metric_rows, cluster_cols = metric_cols)
    } else {
      exprsn = read.csv(file.path(opt$d, eachFile), check.names = F)
      output_name = file.path(opt$o, paste0('heatmap_', gsub('.csv', '', eachFile),'.png'))
      metric_cols = ifelse(dim(exprsn)[2]>4, TRUE, FALSE)
      metric_rows = ifelse(dim(exprsn)[1]>4, TRUE, FALSE)
      
      plot_heatmap(exprsn, output_name, show_rownames = T, cluster_rows = metric_rows, cluster_cols = metric_cols)
    }
  }
  # plot combined expression file, without the specified "--exclude" file
  files = files[which(files != basename(opt$exclude))]
  
  # if there is at least 1 file in the folder besides to exclude file
  if (length(files) > 0){
    exprsn = lapply(files, FUN = function(x){
      read.csv(file.path(opt$d, x), check.names = F)
    })
    exprsn = do.call(rbind, exprsn)
    output_name = file.path(opt$o, paste0('heatmap_combined_lineage.png'))
    metric_cols = ifelse(dim(exprsn)[2]>4, TRUE, FALSE)
    metric_rows = ifelse(dim(exprsn)[1]>4, TRUE, FALSE)
    
    plot_heatmap(exprsn, output_name, show_rownames = F, cluster_rows = metric_rows, cluster_cols = metric_cols)
  }
  
}