#!/usr/bin/env Rscript
# --------------
# Date:  2022-07-08 14:30:00
# Author:Dian Li
# Last update: 2022-07-08


suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# ======================================================= #
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="a specific file to be plotted", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("--sep"), type="character", default='\t',
              help="delimiter for the input file", metavar="character"),
  make_option(c("--group"), type="character", default='group',
              help="the column used as group", metavar="character"),
  make_option(c("--value"), type="character", default='value',
              help="the column used as value", metavar="character"),
  make_option(c("--individual"), type="character", default='individual',
              help="the column used as individual", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


circular_barplot_base = function(data){
  # Get the name and the y position of each label
  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 30 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  # label_data$hjust <- ifelse( angle < -90, 1, 0)
  # label_data$angle <- ifelse(angle < -90, angle+180, angle)
  label_data$hjust <- ifelse( angle < -90 & angle >= -275, 1, 0)
  label_data$angle <- ifelse(angle < -90 & angle >= -275, angle+180, angle)
  # addition adjust for -315 to -360 region, because the origin was rotated 45
  
  
  # prepare a data frame for base lines
  base_data <- data %>% 
    group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  
  max_value = max(data$value, na.rm = T)
  min_value = min(data$value, na.rm = T)
  
  segment_y_1 = round(max_value/40) / 4 * 1 * 40
  segment_y_2 = round(max_value/40) / 4 * 2 * 40
  segment_y_3 = round(max_value/40) / 4 * 3 * 40
  segment_y_4 = round(max_value/40) / 4 * 4 * 40
  
  y_lim = max_value + 10
  
  # Make the plot
  p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = segment_y_4, xend = start, yend = segment_y_4), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = segment_y_3, xend = start, yend = segment_y_3), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = segment_y_2, xend = start, yend = segment_y_2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = segment_y_1, xend = start, yend = segment_y_1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    ggplot2::annotate("text", x = rep(max(data$id),4), y = c(segment_y_1, segment_y_2, segment_y_3, segment_y_4), 
                      label = c(as.character(segment_y_1), as.character(segment_y_2), as.character(segment_y_3), as.character(segment_y_4)) , color="grey", size=2 , angle=0, fontface="bold", hjust=1) +
    
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
    ylim(-100,y_lim) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
    coord_polar(start = 45) + 
    geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.8, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    # geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    geom_text(data=base_data, aes(x = title, y = -48, label=group), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
  
  return(p)
}

# Create dataset
# data <- data.frame(
#   individual=paste( "Mister ", seq(1,60), sep=""),
#   group=as.factor(c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6))) ,
#   value=sample( seq(10,100), 60, replace=T)
# )

data = read.csv(opt$file, sep = opt$sep, check.names = F)
# rename data column names
colnames(data)[which(colnames(data) == opt$group)] = 'group'
colnames(data)[which(colnames(data) == opt$value)] = 'value'
colnames(data)[which(colnames(data) == opt$individual)] = 'individual'

data$group = as.factor(data$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

p = circular_barplot_base(data)
ggsave(filename = file.path(opt$output), plot = p, units = "in", 
       width = 6, height = 6, dpi = 300)
