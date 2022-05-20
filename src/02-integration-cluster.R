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
  sc_list

features <- Seurat::SelectIntegrationFeatures(
  object.list = sc_list,
  nfeatures = 3000
)
sc_list <- Seurat::PrepSCTIntegration(
  object.list = sc_list,
  anchor.features = features
)
anchors <- Seurat::FindIntegrationAnchors(
  object.list = sc_list,
  normalization.method = "SCT",
  anchor.features = features
)
sc_sct <- Seurat::IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

readr::write_rds(
  x = sc_sct,
  file = "data/rda/sc_sct.rds.gz"
)

# Cluster -----------------------------------------------------------------

sc_sct %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:30) %>%
  Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>%
  Seurat::FindClusters(resolution = 0.4) ->
  sc_sct_cluster

readr::write_rds(
  x = sc_sct_cluster,
  file = "data/rda/sc_sct_cluster.rds.gz"
)


# save --------------------------------------------------------------------

save.image(file = "data/rda/02-integration-cluster.rda")
