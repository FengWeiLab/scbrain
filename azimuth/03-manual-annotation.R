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
library(Azimuth)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

# load data ---------------------------------------------------------------

project_sct <- readr::read_rds(
  file = "data/azimuth/project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp.rds.gz"
)


refs <- c(
  "Brain" = "mousecortexref",
  "Meninge" = "pbmcref",
  "Skull" = "bonemarrowref"
)

celllevel <- c(
  "Brain" = "predicted.cluster",
  "Meninge" = "predicted.celltype.l2",
  "Skull" = "predicted.celltype.l2"
)

# body --------------------------------------------------------------------

project_sct |>
  dplyr::select(
    project, region, case, celllevel,
    unique_celltype, unique_celltype_union,
    anno, anno_cell, cellnumber
  ) ->
  project_sct_sel

project_sct_sel |>
  dplyr::filter(case == "Sham") ->
  d

d$anno_cell[[1]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_1

umap_1 <- d$anno_cell[[1]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_1 <- dplyr::bind_cols(umap_1, .xx_1)

d$anno_cell[[2]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_2

umap_2 <- d$anno_cell[[2]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_2 <- dplyr::bind_cols(umap_2, .xx_2)

d$anno_cell[[3]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_3

umap_3 <- d$anno_cell[[3]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_3 <- dplyr::bind_cols(umap_3, .xx_3)

.xxx <- dplyr::bind_rows(.xxx_1, .xxx_2, .xxx_3)


ggplot() +
  geom_point(
    data = .xxx,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = celltype,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) ->
  p_shame



project_sct_sel |>
  dplyr::filter(case == "MCAO") ->
  d

d$anno_cell[[1]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_1

umap_1 <- d$anno_cell[[1]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_1 <- dplyr::bind_cols(umap_1, .xx_1)

d$anno_cell[[2]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_2

umap_2 <- d$anno_cell[[2]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_2 <- dplyr::bind_cols(umap_2, .xx_2)

d$anno_cell[[3]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_3

umap_3 <- d$anno_cell[[3]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_3 <- dplyr::bind_cols(umap_3, .xx_3)

.xxx <- dplyr::bind_rows(.xxx_1, .xxx_2, .xxx_3)


ggplot() +
  geom_point(
    data = .xxx,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = celltype,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) ->
  p_MCAO


project_sct_sel |>
  dplyr::filter(case == "UV") ->
  d

d$anno_cell[[1]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_1

umap_1 <- d$anno_cell[[1]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_1 <- dplyr::bind_cols(umap_1, .xx_1)

d$anno_cell[[2]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_2

umap_2 <- d$anno_cell[[2]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_2 <- dplyr::bind_cols(umap_2, .xx_2)

d$anno_cell[[3]]@meta.data |>
  dplyr::select(case, region, celltype, cluster, cluster_celltype)  ->
  .xx_3

umap_3 <- d$anno_cell[[3]]@reductions[["ref.umap"]]@cell.embeddings

.xxx_3 <- dplyr::bind_cols(umap_3, .xx_3)

.xxx <- dplyr::bind_rows(.xxx_1, .xxx_2, .xxx_3)


ggplot() +
  geom_point(
    data = .xxx,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      colour = celltype,
      shape = NULL,
      alpha = NULL
    ),
    size = 0.7
  ) ->
  p_UV

ggsave(
  filename = "SHAME_merge.pdf",
  plot = p_shame,
  device = "pdf",
  path = "/home/liuc9/github/scbrain/data/scuvresult/06-azimuth",
  width = 15,
  height = 10
)


ggsave(
  filename = "MCAO_merge.pdf",
  plot = p_MCAO,
  device = "pdf",
  path = "/home/liuc9/github/scbrain/data/scuvresult/06-azimuth",
  width = 15,
  height = 10
)


ggsave(
  filename = "UV_merge.pdf",
  plot = p_UV,
  device = "pdf",
  path = "/home/liuc9/github/scbrain/data/scuvresult/06-azimuth",
  width = 15,
  height = 10
)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
