# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun May  7 21:39:54 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
project_sc_azimuth_ref_realcell_sunburst <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/project_sc_azimuth_ref_realcell_sunburst.rds"
)

recell_color_final <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/recell_color.rds"
)
azimuth_ref_sunburst_cell <-  readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
)


# body --------------------------------------------------------------------
synap_SIG <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/synap_SIG.rds"
) |>
  dplyr::mutate(geneID = stringr::str_to_upper(geneID))


a3 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[2]]
a6 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[5]]
a9 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[8]]

a_gene <- union(union(rownames(a3), rownames(a6)),rownames(a9))

gg <- a_gene[grepl("GAB", a_gene)]

synap_SIG |>
  dplyr::filter(geneID %in% a_gene) ->
  a_diff_gene

a_diff_gene |>
  dplyr::filter(grepl("GAB", geneID ))

FeaturePlot(
  object = a3,
  features = "GLS",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)

FeaturePlot(
  object = a6,
  features = "GLS",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)

FeaturePlot(
  object = a6,
  features = "AHSP",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)


Idents(a6) <- "predicted.celltype.l2"
#
VlnPlot(a6, features = "HBM")


a3_diff_gene_m <- AverageExpression(
  a3,
  genes  = a_diff_gene$geneID
)

a3_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a33

a6_diff_gene_m <- AverageExpression(
  a6,
  genes  = a_diff_gene$geneID
)

a6_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a66

a9_diff_gene_m <- AverageExpression(
  a9,
  genes  = a_diff_gene$geneID
)

a9_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a99

head(a33, 20)
head(a66, 20)
head(a99, 20)

FeaturePlot(
  object = a3,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
FeaturePlot(
  object = a6,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
FeaturePlot(
  object = a9,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
