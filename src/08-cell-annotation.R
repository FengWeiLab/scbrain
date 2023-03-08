# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Feb 15 20:56:40 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)
options(future.globals.maxSize = 5 * 1024^3)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

sc_sham_mcao_uv_scn_integrated <- readr::read_rds(
  file = "data/scuvrda/sc_sct.rds.gz"
)

Seurat::DefaultAssay(sc_sham_mcao_uv_scn_integrated) <- "integrated"
# body --------------------------------------------------------------------
sc_sham_mcao_uv_scn_integrated %>%
  Seurat::ScaleData(verbose = F) %>% 
  Seurat::RunPCA(npcs = 30) %>% 
  Seurat::RunUMAP(reduction = "pca", dims = 1:30) |> 
  Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>%
  Seurat::FindClusters(resolution = c(0.3)) ->
  sc_cluster

readr::write_rds(
  x = sc_cluster,
  file = "data/scuvrda/sc_cluster.rds.gz"
)



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------