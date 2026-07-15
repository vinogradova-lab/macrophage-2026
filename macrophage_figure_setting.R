suppressWarnings(suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(msigdbr)
}))

library(readxl)
library(tidyverse)
library(ggrepel)
library(ggnewscale)
library(ggpubr)
library(ggh4x)
library(ggbreak)
library(circlize)
library(ComplexHeatmap)
library(extrafont)
library(lemon)
library(scales)
library(httr)
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
TEXT_ELEMENT = element_text(
  size=FONT_SIZE, 
  color = "black", 
  family=FONT_FAMILY)

my_theme <- function(){
  theme_classic() + 
    theme(axis.text.x = TEXT_ELEMENT, 
          axis.text.y = TEXT_ELEMENT, 
          axis.title.x = TEXT_ELEMENT, 
          axis.title.y = TEXT_ELEMENT,
          legend.title = TEXT_ELEMENT,
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


# Shared Fisher-exact GO enrichment via clusterProfiler::enricher() over the MSigDB C5 GO
# collection (msigdbr >= 26.1.0 = MSigDB 2026.1.Hs). Caches TERM2GENE/TERM2NAME per ontology.
.go_msigdb <- local({
  cache <- list()
  function(ont) {
    if (is.null(cache[[ont]])) {
      m <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5",
                            subcollection = paste0("GO:", ont))
      cache[[ont]] <<- list(
        t2g = dplyr::distinct(m, gs_exact_source, gene_symbol),
        t2n = dplyr::distinct(m, gs_exact_source, gs_name) |>
          dplyr::transmute(gs_exact_source,
                           name = tolower(gsub("_", " ", sub(paste0("^GO", ont, "_"), "", gs_name))))
      )
    }
    cache[[ont]]
  }
})

# Member gene SYMBOLS of one or more named MSigDB gene sets (union), from the same msigdbr
# source as go_enrich(). Lets figures pull a panel instead of hardcoding it. Caches each pull.
#   name          : one or more exact gs_names, e.g. "HALLMARK_INTERFERON_GAMMA_RESPONSE";
#                   passing several returns their union.
#   collection    : MSigDB collection, e.g. "H" (Hallmark) or "C5" (GO).
#   subcollection : sub-collection, e.g. "GO:BP"; NULL for collections without one.
msigdb_genes <- local({
  cache <- list()
  function(name, collection, subcollection = NULL) {
    key <- paste(collection, subcollection, sep = "|")
    if (is.null(cache[[key]])) {
      cache[[key]] <<- msigdbr::msigdbr(species = "Homo sapiens",
                                        collection = collection, subcollection = subcollection)
    }
    m <- cache[[key]]
    missing <- setdiff(name, m$gs_name)
    if (length(missing)) stop("No MSigDB set named '", paste(missing, collapse = "', '"),
                              "' in ", key)
    sort(unique(m$gene_symbol[m$gs_name %in% name]))
  }
})

# GOSemSim IC for simplify(). godata()'s "preparing IC data" step is slow (~30s/ontology),
# so cache it on disk (and in-memory per session). The .rds is built once and reused.
# Delete reference_dbs/go_semdata_*.rds to force a rebuild (e.g. after a GO.db update).
GO_SEMDATA_DIR <- file.path(
  "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Vinogradova Laboratory/Henry_data processing/01_data_analysis_folders/06_macrophage/macrophage",
  "reference_dbs"
)
.go_semdata <- local({
  cache <- list()
  function(ont) {
    if (is.null(cache[[ont]])) {
      rds <- file.path(GO_SEMDATA_DIR, paste0("go_semdata_", ont, ".rds"))
      if (file.exists(rds)) {
        cache[[ont]] <<- readRDS(rds)
      } else {
        cache[[ont]] <<- godata("org.Hs.eg.db", ont = ont)
        dir.create(GO_SEMDATA_DIR, showWarnings = FALSE, recursive = TRUE)
        saveRDS(cache[[ont]], rds)
      }
    }
    cache[[ont]]
  }
})

# Empty result in the standard schema (returned when enrichment finds nothing).
.go_empty <- function(ont) {
  dplyr::tibble(ontology = character(), Term = character(), GO_id = character(),
                Overlap = character(), n_proteins = integer(), p_value = double(),
                adjusted_p_value = double(), neg_log10_fdr = double(), Genes = character())
}

# Fisher-exact GO over-representation for one ontology.
#   genes       : character vector of gene SYMBOLS of interest.
#   ont         : "BP", "CC", or "MF".
#   background  : character vector of universe SYMBOLS (modality-specific), or NULL
#                 to use every gene in the MSigDB GO sets (= "all annotated genes").
#   simplify_cutoff : GOSemSim similarity cutoff for redundancy reduction; NULL disables.
go_enrich <- function(genes, ont = "BP", background = NULL,
                      simplify_cutoff = 0.5, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                      minGSSize = 10, maxGSSize = 500) {
  sets <- .go_msigdb(ont)
  res <- clusterProfiler::enricher(
    gene = genes, universe = background,
    TERM2GENE = sets$t2g, TERM2NAME = sets$t2n,
    pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize, maxGSSize = maxGSSize)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(.go_empty(ont))
  if (!is.null(simplify_cutoff)) {
    res@ontology <- ont  # enricher() leaves this "UNKNOWN"; simplify() requires BP/CC/MF
    res <- clusterProfiler::simplify(res, cutoff = simplify_cutoff, by = "p.adjust",
                                     select_fun = min, semData = .go_semdata(ont))
  }
  as.data.frame(res) |>
    dplyr::transmute(ontology = ont, Term = Description, GO_id = ID, Overlap = GeneRatio,
                     n_proteins = Count, p_value = pvalue, adjusted_p_value = p.adjust,
                     neg_log10_fdr = -log10(p.adjust), Genes = geneID)
}

# Run go_enrich() over several ontologies and bind into one tidy frame.
go_enrich_multi <- function(genes, onts = c("BP", "MF", "CC"), background = NULL, ...) {
  dplyr::bind_rows(lapply(onts, function(o) go_enrich(genes, ont = o, background = background, ...)))
}

# Shared GO-enrichment dot plot in the reactivity figure style:
# y = GO term, x = -log10(FDR), bubble size = number of proteins.
# Consumes the standard schema written by src/go_enrichment.py
# (columns: Term, adjusted_p_value, n_proteins, ontology, ...).
# `facet` (optional, a column name like "ontology" or "direction") gives one panel per
# group with per-panel y ordering (via tidytext). NULL = the single-panel attached format.
# `color_col` (optional, a column name like "omic") colors the bubbles by a discrete group,
# with `colors` a named palette for scale_fill_manual; NULL = the single fixed `fill_color`.
# `xlab` overrides the x-axis label (e.g. when plotting -log10(p) instead of FDR).
# `x_transform` maps `fdr_col` onto the x-axis (default -log10 for FDR); pass `identity`
# to plot a raw value such as an association score, and `size_name` labels the size legend.
go_dotplot <- function(df, term_col = "Term", fdr_col = "adjusted_p_value",
                       size_col = "n_proteins", facet = NULL,
                       color_col = NULL, colors = NULL,
                       fill_color = "#B3CDE3", size_range = c(1, 3), title = NULL,
                       xlab = expression("-log"[10] * "(FDR)"),
                       x_transform = function(x) -log10(x), size_name = "proteins",
                       facet_scales = "free", facet_ncol = NULL) {
  df <- df %>% mutate(.nlf = x_transform(.data[[fdr_col]]))
  if (is.null(facet)) {
    p <- ggplot(df, aes(x = .nlf, y = reorder(.data[[term_col]], .nlf), size = .data[[size_col]]))
  } else {
    p <- ggplot(df, aes(x = .nlf,
                        y = tidytext::reorder_within(.data[[term_col]], .nlf, .data[[facet]]),
                        size = .data[[size_col]]))
  }
  if (is.null(color_col)) {
    p <- p + geom_point(shape = 21, stroke = LINE_WIDTH, color = "black", fill = fill_color)
  } else {
    p <- p + geom_point(aes(fill = .data[[color_col]]), shape = 21, stroke = LINE_WIDTH,
                        color = "black")
    if (!is.null(colors)) p <- p + scale_fill_manual(values = colors, name = NULL)
  }
  p <- p +
    my_theme() +
    scale_size_continuous(range = size_range, name = size_name) +
    labs(x = xlab, y = NULL, title = title) +
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.3)))
  if (!is.null(facet)) {
    p <- p +
      facet_wrap(vars(.data[[facet]]), scales = facet_scales, ncol = facet_ncol) +
      tidytext::scale_y_reordered() +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = FONT_SIZE, family = FONT_FAMILY, color = "black"))
  }
  p
}

legend_theme <- function(plt){
  # settings for a small legend on the right side
  my_theme()  +
    theme(
      legend.title=element_blank(),
      legend.margin = margin(l = -5, b = -4.5),
      legend.box.margin = margin(l = -5, b = -5.7),
      legend.key.spacing.y = unit(0, "mm"),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(family = "Arial", color = "black", size = 6),
      legend.justification = ("center"),
      legend.position = "right",
      panel.border =element_blank(),
      axis.line.x = element_line(color = "black", linetype = "solid", linewidth = LINE_WIDTH),
      axis.line.y = element_line(color = "black", linetype = "solid", linewidth = LINE_WIDTH))
}

cols <- c(
  "M0" = "#AAAAAA",
  "TLR1-2"= "#EC2427",
  "TLR3" = "#F58420",
  "TLR4" = "#FFCC31",
  "TLR7" = "#4869B2",
  "TLR8" = "#5CBED7",
  "TLR9" = "#237D41",
  "STING" = "#742D16"
)

EXHAUSTION_COLS <- cols

SIG_UP_COL <- "indianred2"
SIG_DOWN_COL <- "steelblue2"
UNCHANGED_COL <- "lightgrey"

conditions = c("TLR1-2","TLR3","TLR4","TLR7","TLR8","TLR9","STING")
control_condition <- "M0"

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

regulation_colors <- c(
  "Higher" = SIG_UP_COL, 
  "Lower" = SIG_DOWN_COL, 
  "Unchanged" = UNCHANGED_COL, 
  "Protein expression" = "#54A868"
)

is_in_group <- function(row, fun_group){
  if (row["group"] == fun_group) {
    return(fun_group)
  } else {
    return("Other")
  }
}

save_plot<- function(fn, width, height){
  ggsave(paste0(fn, ".svg"), height = height, width = width)
  ggsave(paste0(fn, ".png"), height = height, width = width)
}

# Build a "PC1 (NN%)" axis label from the percent_explained.csv written by
# src/pca_utils.py (columns principal_component, percent_explained), so the PCA
# axis labels stay in sync with the data on every rerun.
pc_label <- function(pca_dir, pc){
  pct <- read_csv(paste0(pca_dir, "percent_explained.csv"), show_col_types = FALSE)
  paste0(pc, " (", pct$percent_explained[pct$principal_component == pc], "%)")
}

# Centralized PCA scatter used by both the whole-proteome and RNA-seq PCA panels
# (whole_proteome/wp_visualization.Rmd). Reads PCA-coordinate / loadings data frames
# produced by src/pca_utils.py::run_pca and saves a plain scatter plus a loadings overlay.
pca_plot <- function(df,loadings_df, x_column, y_column, x_lab, y_lab, fn, arrow_scaling = 250){
  df %>%
    ggplot(aes_string(
      x = x_column,
      y = y_column,
      fill = "condition"
    )) +
    geom_point(size = POINT_SIZE,
             stroke=POINT_STROKE,
             shape=21,
             alpha=0.9,
             col="black",
             show.legend = FALSE) +
    my_theme() +
    labs(x = x_lab, y = y_lab) +
    scale_fill_manual(values=cols, name="") +
    theme(
         aspect.ratio=1) -> plt
  print(plt)
  save_plot(fn, height=2, width = 2)

    # PCA results with top loadings
  num_loadings = 5
  loadings_df %>% arrange(PC1) %>% tail(num_loadings) %>% pull(variable)-> top_pc1
  loadings_df %>% arrange(PC1) %>% head(num_loadings) %>% pull(variable)-> bottom_pc1
  loadings_df %>% arrange(PC2) %>% tail(num_loadings) %>% pull(variable)-> top_pc2
  loadings_df %>% arrange(PC2) %>% head(num_loadings) %>% pull(variable)-> bottom_pc2

  loadings_df %>% filter(variable %in% c(top_pc1, bottom_pc1, bottom_pc2, top_pc2)) -> filtered_loadings
  print(filtered_loadings)
  plt + geom_segment(data= filtered_loadings,
                      aes(x = 0, y = 0,
                          xend = PC1 * arrow_scaling,
                          yend = PC2 * arrow_scaling),
                      inherit.aes = FALSE,
                    arrow = arrow(length = unit(0.1, "cm")),
                    size = LINE_WIDTH,
                    color="grey80") +
    geom_text_repel(data= filtered_loadings,
                    mapping = aes(label=variable,x=PC1*arrow_scaling, y=PC2*arrow_scaling),
                    inherit.aes = FALSE,
                    segment.size = LINE_WIDTH,
                    min.segment.length = unit(0,"mm"),
                    size=1.763,
                    color="black",
                    #bg.color = "white",
                    segment.color = "grey80") -> plot_with_loadings
  print(plot_with_loadings)
  save_plot(paste0(fn, "_with_loadings"), height=2, width = 2)

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

# Centralized grouped-volcano plot shared by the whole-proteome expression, IP-MS
# endogenous pulldown, and phosphoproteomics panels (whole_proteome/wp_visualization.Rmd
# and phosphoproteomics/visualization.Rmd). Draws one condition-vs-control volcano with
# per-group significant-count boxes and highlighted-protein labels. The differences that
# used to distinguish the three hand-copied versions are exposed as parameters:
#   color_code       - named palette mapping each `group` value to a fill colour.
#   plot_directory   - output directory; when NULL, derive a group-based subdirectory
#                      under the global `volcano_dir` (as the expression panels do).
#   id_var           - column used to flag significant points ("uniprot" or "protein").
#   fc_cutoff        - fold-change cutoff for the dashed vertical lines.
#   rasterize_points - rasterize only the point layers (keeps SVGs light for dense phospho data).
#   x_axis_word      - noun in the x-axis label ("expression" vs "pulldown").
#   y_max / x_max    - optional axis-limit overrides; computed from the data when NULL.
make_condition_plot_multi <- function(condition,
                                      volcano_data,
                                      control_condition,
                                      color_code,
                                      highlight_proteins = c(),
                                      plot_directory     = NULL,
                                      id_var             = "uniprot",
                                      fc_cutoff          = FC_CUTOFF,
                                      rasterize_points   = FALSE,
                                      x_axis_word        = "expression",
                                      y_max              = NULL,
                                      x_max              = NULL){

  if (is.null(y_max)){
    y_max <- max(abs(volcano_data$`-log10_pval`)) * 1.2
  }
  if (is.null(x_max)){
    x_max <- max(abs(volcano_data$log2_FC)) * 1.5
  }
  nudge_x <- 2

  x_var <- "log2_FC"
  y_var <- "-log10_pval"
  change_var = "Regulation"
  name_var = "protein"
  group <- "group"

  # build label dataframe
  labels <- list()
  y <- y_max * (15.5/16)
  x <- x_max * (4.7/5.3)
  for (direction in c("Significant Up", "Significant Down")){
   for (bio_group in names(color_code)){
      volcano_data[volcano_data[[change_var]] == direction &
                           volcano_data[["group"]] == bio_group,] %>%
       nrow() -> protein_count
     r <- c("direction" = direction,
            "bio_group" = bio_group,
            "count" = protein_count,
            "x" = x,
            "y" = y)
     labels[[paste0(bio_group, direction)]] <- r
     y <- y - (y_max * (1.5/16))
   }
    y <- y_max * (15.5/16)
    x <- x * -1
  }

  do.call("rbind", labels) %>% as.data.frame() -> count_label_df
  title <- paste0(control_condition," vs ", condition)

  alphas_binary <- c(0.4,0.9)

  sig <- volcano_data[ volcano_data[[change_var]] %in% c("Significant Up", "Significant Down") , ]

  volcano_data %>% mutate(significant = .data[[id_var]] %in% sig[[id_var]]) -> volcano_data
  sig$significant <- TRUE
  label_data <- sig[ sig[[name_var]] %in% highlight_proteins , ]

  set.seed(12345)
  volcano_data$group <- factor(volcano_data$group, levels = rev(names(color_code)))
  volcano_data %>%
    arrange(group) %>%
    ggplot(aes(x=log2_FC,
               y=`-log10_pval`,
                fill = group,
               alpha = significant)) +
    geom_point(data = volcano_data %>% filter(significant == FALSE),
               size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black",
               fill = 'gray88'
               ) +
    geom_point(data = volcano_data %>% filter(group == "Other", significant == TRUE),
               size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black") +
        geom_point(data = volcano_data %>% filter(group == names(color_code)[3], significant == TRUE),
              size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black") +
            geom_point(data = volcano_data %>% filter(group == names(color_code)[2], significant == TRUE),
              size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black") +

                geom_point(data = volcano_data %>% filter(group == names(color_code)[1], significant == TRUE),
              size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black") +
    geom_point(data = label_data,
              size = 1,
               shape=21,
               stroke=POINT_STROKE / 2,
               col="black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", size=LINE_WIDTH) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", size=LINE_WIDTH) +
    geom_text_repel(
      data = label_data[label_data[[change_var]] == "Significant Down",],
      aes(
        x = log2_FC,
        y = `-log10_pval`,
        color = group,
      # label with gene name instead of uniprot
      label = protein),
      direction    = "both",
      box.padding = 1,
      nudge_x = -nudge_x,
      xlim = c(NA, -1),
      size =FONT_SIZE_MM,
      #bg.color = "white",
      min.segment.length = 0,
      segment.size = LINE_WIDTH,
      inherit.aes = FALSE,
      show.legend = FALSE) +
    geom_text_repel(
      data = label_data[label_data[[change_var]] == "Significant Up",],
            aes(
        x = log2_FC,
        y = `-log10_pval`,
        color = group,
      # label with gene name instead of uniprot
      label = protein),
      direction    = "both",
      box.padding = 1,
      nudge_x = nudge_x,
      xlim = c(1,NA),
      size = FONT_SIZE_MM,
      #bg.color = "white",
      min.segment.length = 0,
      segment.size = LINE_WIDTH,
      inherit.aes = FALSE,
      show.legend = FALSE) +
  geom_label(
    data = count_label_df,
    mapping = aes(x = as.numeric(x),
                  y = as.numeric(y),
                  label = count,
                  fill = bio_group),
    color = "white",
    size=6,
    size.unit="pt",
    show.legend = FALSE,
    inherit.aes = FALSE) +
    theme_test() +
    labs(
      x = substitute(paste('protein ', w, ', log'[2], "(", n, '/', m, ")"),
                     list(w = x_axis_word, n = condition, m = control_condition)),
      y = substitute(paste('-log'[10], '(p-value)'))
    ) +
    scale_fill_manual(values=color_code,
                       breaks = rev(names(color_code))
                      ) +
        scale_color_manual(values=color_code,
                       breaks = rev(names(color_code))
                      ) +
    scale_alpha_manual(values = alphas_binary, guide="none") +
    my_theme() +
  theme(strip.text.x = element_blank(),
      strip.background = element_rect(colour="white", fill="white"),
      legend.title=element_blank(),
      legend.box.background = element_blank(),
      legend.position = "top",
      legend.justification.top = "top",
      legend.justification = "center",
      legend.box.just = "center",
      aspect.ratio=1,
      legend.key = element_blank(),
      legend.margin = margin(l = -5, b = -5.5),
      legend.box.margin = margin(l = -5, b = -0.8),
      legend.key.spacing.y = unit(-0.9, "cm"),
      legend.key.spacing.x = unit(-0.05, "cm")) +
    guides(fill = guide_legend(ncol = 2,
                               )) +
    ylim(0,y_max) +
    xlim(-x_max,x_max) -> v_plot
  if (rasterize_points){
    # Rasterize ONLY the point layers (thousands of dots) at high DPI so the SVG stays
    # light and fast; text/labels/count-boxes/axes remain crisp vector graphics.
    v_plot <- ggrastr::rasterise(v_plot, layers = "Point", dpi = 600, dev = "ragg")
  }
  print(v_plot)
  if (is.null(plot_directory)){
    plot_directory <- paste0(volcano_dir, paste(unique(volcano_data$group)[unique(volcano_data$group) != "Other"], collapse = "_"))
  }
  dir.create(plot_directory, showWarnings = FALSE)
  save_plot( paste0(plot_directory, "/", title), width=2, height=4)
}