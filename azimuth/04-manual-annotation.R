# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun May 21 18:33:38 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
azimuth_ref <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/project_sc_azimuth_ref_realcell_sunburst.rds"
)

# body --------------------------------------------------------------------


azimuth_ref$anno[[3]] |>
  Seurat::RunPCA() |>
  Seurat::RunTSNE() |>
  Seurat::RunUMAP(dims = 1:30) ->
  a

a$predicted.celltype.l2
Idents(a) <- a$predicted.celltype.l1


FindMarkers(
  a,
  ident.1 = "CD8 T",
  ident.2 = "B",
)


DimPlot(
  a,
  group.by = "predicted.celltype.l2",
  reduction = "ref.umap"
)

DimHeatmap(a, dims = 2)

DimPlot(
  azimuth_ref$anno_new[[3]],
  group.by = "predicted.celltype.l2"
)




# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
