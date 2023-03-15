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

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------

brain_meninge_skull_sct_cluster_sctype |> 
  dplyr::select(region, sct_list, sct, sct_cluster) ->
  a

a$sct[[1]]@meta.data |> head()

a$

VariableFeaturePlot(a$sct_list[[1]]$CB)

DimPlot(a$sct_cluster[[1]])
VlnPlot(a$sct_cluster[[1]], features = c("GZMK"))
a$sct_cluster[[1]]@meta.data |> head()

PCAPlot(a$sct_cluster[[1]], split.by="case")
PCAPlot(a$sct_cluster[[1]], split.by="tissue")

DimPlot(a$sct_cluster[[3]], group.by  = "case")

VizDimLoadings(a$sct_cluster[[1]])
DimHeatmap(a$sct_cluster[[1]])
ElbowPlot(a$sct_cluster[[1]])

VlnPlot(
  a$sct_cluster[[1]],
  features = "percent.mt"
)

VlnPlot(
  a$sct_cluster[[1]],
  features = "percent.ribo"
)

FeaturePlot(
  a$sct_cluster[[1]],
  features = "G2M.Score"
)

plot_integrated_clusters(a$sct_cluster) 
# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------