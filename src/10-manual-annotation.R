# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Mar 14 16:41:18 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(HGNChelper)


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot <- readr::read_rds(
  file = "data/scuvrda/project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot.rds.gz"
) 


# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------