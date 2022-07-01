#!/usr/bin/env Rscript
# --------------
# Date:  2022-06-19 21:30:10
# Author:Dian Li
# Email: dianli@wustl.edu
# Last update: 2022-06-29

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
              help="file that to be excluded from being combined", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

call_pheatmap <- function(mat, output_name, height, width, annotation_row, ...){
  png(paste0(output_name), height = height, width = width, units = "in", res = 300)
  pheatmap(mat, na_col = "grey", annotation_row = annotation_row, ...)
  dev.off()
}


plot_heatmap <- function(exprsn, output_name, ...){
  mat = exprsn[, c(-1,-2)]
  
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
  
  annotation_row = data.frame("Lineage" = exprsn$lineage)
  rownames(annotation_row) = rownames(mat)
  
  if (dim(mat)[1] <= 20){
    height = 6+0.1*dim(mat)[1]
    width = 6+0.2*dim(mat)[2]
  } else if (dim(mat)[1] < 200 & dim(mat)[1] > 20){
    height = 6+0.01*dim(mat)[1]
    width = 6+0.2*dim(mat)[2]
  } else {
    height = 6+0.005*dim(mat)[1]
    width = 6+0.2*dim(mat)[2]
  }
  
  call_pheatmap(mat, output_name, height, width, annotation_row, angle_col = 45, ...)

}


# ========================================= #
# if the input is a single file

if (!is.null(opt$f)){
  exprsn = read.csv(opt$f, check.names = F)
  output_name = file.path(opt$o, paste0('heatmap_', gsub('.csv', '', basename(opt$f)),'.png'))
  plot_heatmap(exprsn, output_name)
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
      metric = ifelse(dim(exprsn)>3, TRUE, FALSE)
      plot_heatmap(exprsn, output_name, show_rownames = F, cluster_rows = metric, cluster_cols = metric)
    } else {
      exprsn = read.csv(file.path(opt$d, eachFile), check.names = F)
      output_name = file.path(opt$o, paste0('heatmap_', gsub('.csv', '', eachFile),'.png'))
      metric = ifelse(dim(exprsn)>2, TRUE, FALSE)
      plot_heatmap(exprsn, output_name, show_rownames = T, cluster_rows = metric, cluster_cols = metric)
    }
  }
  # plot combined expression file, without the specified "--exclude" file
  files = files[which(files != basename(opt$exclude))]
  exprsn = lapply(files, FUN = function(x){
    read.csv(file.path(opt$d, x), check.names = F)
  })
  exprsn = do.call(rbind, exprsn)
  output_name = file.path(opt$o, paste0('heatmap_combined_lineage.png'))
  plot_heatmap(exprsn, output_name, show_rownames = F, cluster_rows = F, cluster_cols = T)
  
}