library(tidyverse)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
library(extrafont)
library(scales)
library(ggbeeswarm)

# load arial font
loadfonts(quiet = TRUE)

options(ggrepel.max.overlaps = Inf) # for crowded volcano plots

POINT_STROKE = 0.25
POINT_SIZE = 2.5
LINE_WIDTH = 0.25/2
FONT_SIZE = 6
FONT_SIZE_MM = 2 + (7/60)
FONT_FAMILY = "Arial"
BARPLOT_WIDTH = 0.75
FC_CUTOFF <- log2(1.5)
AXIS_LINE = element_line(color = "black", linetype = "solid", linewidth = LINE_WIDTH)

my_theme <- function(){
  standard_text <- element_text(
    size=FONT_SIZE, 
    color = "black", 
    family=FONT_FAMILY)
  theme_classic() + 
    theme(axis.text.x = standard_text, 
          axis.text.y = standard_text, 
          axis.title.x = standard_text, 
          axis.title.y = standard_text,
          legend.title = standard_text,
          legend.text = element_text(
            size=6, color = "black", 
            family="Arial",
            margin = margin(l = -5)),
          axis.line = element_line(colour = "black", linewidth=0),
          axis.ticks = element_line(colour="black", linewidth = LINE_WIDTH),
          panel.border = element_rect(colour = "black", fill=NA, size=LINE_WIDTH * 2),
          plot.title = element_text(
            size=6, 
            color = "black", 
            family=FONT_FAMILY,
            hjust=0.5),
    )
}

cols <- c(
  "M0" = "#AAAAAA",
  "TLR1-2"= "#EC2427",
  "TLR3" = "#F58420",
  "TLR4" = "#FFCC31",
  "LPS" = "#FFCC31",
  "TLR7" = "#4869B2",
  "TLR8" = "#5CBED7",
  "TLR9" = "#237D41",
  "STING" = "#742D16"
)

UNCHANGED_COL <- "lightgrey"

conditions <- c(
  "D4A",
  "D4C",
  "D8A",
  "D8C"
)

alphas <- c(
  "Significant Up" = 1,
  "Significant Down" = 1,
  "Not Significant" = .35,
  "Not Significant Up" = .35,
  "Not Significant Down" = .35,
  "Significant but <1.5 FC" = .35
)

alphas_binary <- c(
  "TRUE" = 1,
  "FALSE" = 0.35
)

save_plot<- function(fn, width, height){
  ggsave(paste0(fn, ".svg"), height = height, width = width)
  ggsave(paste0(fn, ".png"), height = height, width = width)
}

replot_pca <- function(pca_dir, x, y, arrow_scaling, add_grid = FALSE){
  paste0(pca_dir, "percent_explained.csv") %>%
    read_csv(show_col_types = FALSE) -> percent_explained
  pcs <- split(percent_explained$percent_explained, percent_explained$principal_component)
  
  # PCA results
  paste0(pca_dir, "pca_results.csv") %>%
    read.csv() %>%
    ggplot(aes(x=!!sym(x), y=!!sym(y),fill=condition)) + 
    geom_point(size = 2.5,
               shape=21, 
               alpha = 0.9,
               col="black", 
               stroke = POINT_STROKE,
               show.legend = FALSE) +
    my_theme() + 
    xlab(paste0(x, " (", pcs[x],"%)")) + 
    ylab(paste0(y, " (", pcs[y], "%)")) + 
    theme(
      
      aspect.ratio=1) +
    scale_fill_manual(values = cols) -> plot
  
  if(add_grid){
   plot <- plot + theme(panel.grid.major = element_line(colour = "grey90",  linewidth=LINE_WIDTH))
  }
  print(plot)
  save_plot(paste0(pca_dir, "pca_plot"), height = 4, width = 2)
  
  # PCA results with sample names
  plot + geom_text_repel(aes(label=channel_name), size=1.41111) -> plot_with_sample_names
  print(plot_with_sample_names)
  save_plot(paste0(results_dir, "pca/pca_plot_sample_names"), height = 4, width = 2)
  
  # PCA results with top loadings
  paste0(pca_dir, "loadings_results.csv") %>%
    read.csv() -> loadings_df
  num_loadings = 5
  loadings_df %>% arrange(!!sym(x)) %>% tail(num_loadings) %>% pull(variable)-> top_pc1
  loadings_df %>% arrange(!!sym(x)) %>% head(num_loadings) %>% pull(variable)-> bottom_pc1
  loadings_df %>% arrange(!!sym(y)) %>% tail(num_loadings) %>% pull(variable)-> top_pc2
  loadings_df %>% arrange(!!sym(y)) %>% head(num_loadings) %>% pull(variable)-> bottom_pc2
  
  loadings_df %>% filter(variable %in% c(top_pc1, bottom_pc1, bottom_pc2, top_pc2)) -> filtered_loadings
  plot + geom_segment(data= filtered_loadings, 
                      aes(x = 0, y = 0, 
                          xend = !!sym(x) * arrow_scaling, 
                          yend = !!sym(y) * arrow_scaling), 
                      inherit.aes = FALSE,
                      arrow = arrow(length = unit(0.1, "cm")),
                      size = LINE_WIDTH,
                      color="grey90") + 
    geom_text_repel(data= filtered_loadings,
                    mapping = aes(label=variable,x=!!sym(x)*arrow_scaling, y=!!sym(y)*arrow_scaling), 
                    inherit.aes = FALSE, 
                    segment.size = LINE_WIDTH,
                    min.segment.length = unit(0,"mm"),
                    size=1.763, 
                    color="black",
                    segment.color = "grey90") -> plot_with_loadings
  print(plot_with_loadings)
  save_plot(paste0(results_dir, "pca/pca_plot_top_loadings"), height = 4, width = 2 )
}