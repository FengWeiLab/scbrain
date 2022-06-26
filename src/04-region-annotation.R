# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Jun 26 16:38:21 2022
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)

# Load data ---------------------------------------------------------------

project_sct <- readr::read_rds(file = "data/rda/project_sc_sct.rds.gz")

regions <- c(
  "CB" = "Brain", 
  "CM" = "Meninge", 
  "CS" = "Skull", 
  "DB" = "Brain", 
  "DM" = "Meninge", 
  "DS" = "Skull"
)
cases <- c(
  "CB" = "Sham", 
  "CM" = "Sham", 
  "CS" = "Sham", 
  "DB" = "MCAO", 
  "DM" = "MCAO", 
  "DS" = "MCAO"
)

project_sct %>% 
  dplyr::mutate(
    region = plyr::revalue(x = project, replace = regions),
    case = plyr::revalue(x = project, replace = cases)
  ) ->
  project_sct_regions


# Split by region ---------------------------------------------------------



project_sct_regions %>% 
  dplyr::filter(region == "Brain") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  brain_sct_list


project_sct_regions %>% 
  dplyr::filter(region == "Meninge") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  meninge_sct_list


project_sct_regions %>% 
  dplyr::filter(region == "Skull") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  skull_sct_list


# Function ----------------------------------------------------------------

fn_sct_integrate <- function(.sct_list) {
  # .sct_list <- brain_sct_list
  
  .features <- Seurat::SelectIntegrationFeatures(
    object.list = .sct_list,
    nfeatures = 3000
  )
  .sct_list <- Seurat::PrepSCTIntegration(
    object.list = .sct_list,
    anchor.features = .features
  )
  .anchors <- Seurat::FindIntegrationAnchors(
    object.list = .sct_list,
    normalization.method = "SCT",
    anchor.features = .features
  )
  .sct_list <- Seurat::IntegrateData(
    anchorset = .anchors,
    normalization.method = "SCT"
  )
  
  .sct_list
}

fn_cluster <- function(.sct) {
  .sct %>%
    Seurat::RunPCA(npcs = 30) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:30) %>%
    Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>%
    Seurat::FindClusters(resolution = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) ->
    .sct_cluster
  .sct_cluster <- Seurat::PrepSCTFindMarkers(object = .sct_cluster)
  .sct_cluster
}


# Cluster -----------------------------------------------------------------


brain_meninge_skull_sct_list <- tibble::tibble(
  region = c("brain", "meninge", "skull"),
  sct_list = list(brain_sct_list, meninge_sct_list, skull_sct_list)
)


brain_meninge_skull_sct_list %>% 
  dplyr::mutate(
    sct = purrr::map2(
      .x = region,
      .y = sct_list,
      .f = function(.x, .y) {
        fn_sct_integrate(.sct_list = .y)
      }
    )
  ) %>% 
  dplyr::mutate(
    sct_cluster = purrr::map2(
      .x = region,
      .y = sct,
      .f = function(.x, .y) {
        fn_cluster(.sct = .y)
      }
    )
  ) ->
  brain_meninge_skull_sct_cluster

readr::write_rds(
  x = brain_meninge_skull_sct_cluster,
  file = "data/rda/brain_meninge_skull_sct_cluster.rds.gz"
)


# Marker genes ------------------------------------------------------------



brain_meninge_skull_sct_cluster$sct_cluster[[1]]$integrated_snn_res.0.3 %>% table()

brain_meninge_skull_sct_cluster$sct_cluster[[2]]$integrated_snn_res.0.3 %>% table() 

brain_meninge_skull_sct_cluster$sct_cluster[[3]]$integrated_snn_res.0.3 %>% table()
