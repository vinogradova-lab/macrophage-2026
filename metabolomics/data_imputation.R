#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(imputeLCMD)

args[1] %>% 
  read_csv(show_col_types = FALSE) %>% 
  column_to_rownames(args[3]) %>%
  as.matrix() %>% log2() -> log2_mat
log2_mat[log2_mat == -Inf]  <- NA
log2_mat %>% impute.QRILC() -> impute_results
impute_results[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(args[3]) %>%
  write_csv(args[2])