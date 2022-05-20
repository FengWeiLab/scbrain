# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May 19 18:08:35 2022
# @DESCRIPTION: integrate data and cluster

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)

# Load data ---------------------------------------------------------------

project_sct <- readr::read_rds(file = "data/rda/project_sc_sct.rds.gz")


# Function ----------------------------------------------------------------

fn_stat_cell <- function(.x, .y) {
  # .x <- project_sct$sc[[1]]
  # .y <- project_sct$sct[[1]]
  
  .xd <- .x@meta.data
  .yd <- .y@meta.data
  
  .n_x_cells <- nrow(.xd)
  .mean_reads_per_cell <- mean(.xd$nCount_RNA)
  .median_gene_per_cell <- median(.xd$nFeature_RNA)
  .n_y_cells <- nrow(.yd)
  
  tibble::tibble(
    `estimated number of cells` = .n_x_cells,
    `mean reads per cell` = .mean_reads_per_cell,
    `median genes per cell` = .median_gene_per_cell,
    `number of cells after filtering` = .n_y_cells
  )
}



# stat --------------------------------------------------------------------


project_sct %>%
  dplyr::mutate(
    stat = purrr::map2(
      .x = sc,
      .y = sct,
      .f = fn_stat_cell
    )
  ) %>% 
  dplyr::select(project, stat) %>% 
  tidyr::unnest(cols = stat) ->
  project_sct_reads

# save stat
writexl::write_xlsx(x = project_sct_reads, path = "data/result/01-basic/reads-stat.xlsx")



# Integration -------------------------------------------------------------

project_sct %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  sct_list
