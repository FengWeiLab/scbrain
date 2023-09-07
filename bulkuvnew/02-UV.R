#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Sep  6 21:43:55 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------

source(file = "utils/utils.R")


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------



UV_data <- readr::read_rds(
  file = file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew",
    "all_expr_count_matrix_uv.rds.gz"
  )
)

tibble::tibble(
  barcode = colnames(UV_data)[-1]
) |>
  dplyr::mutate(
    group = gsub(pattern = "_\\d", replacement = "", x = barcode)
  ) |>
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "B", x = group), "Brain", "Skull")
  ) |>
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "M", x = group), "Meninge", seq)
  ) |>
  dplyr::mutate(
    case = ifelse(grepl(pattern = "N", x = group), "Control", "UVB 0.5h")
  ) |>
  dplyr::mutate(
    case = ifelse(grepl(pattern = "1", x = group), "UVB 24h", case)
  ) ->
  barcode_group

barcode_group |> print(n = Inf)

readr::write_tsv(
  barcode_group,
  file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew",
    "barcode_group.tsv"
  )
)
# body --------------------------------------------------------------------

UV_data |>
  tibble::column_to_rownames("GeneName") |>
  as.matrix() ->
  UV_data.matrix

barcode_group %>%
  dplyr::mutate(new_barcode = barcode) %>%
  tibble::column_to_rownames("new_barcode") ->
  UV_data.matrix.meta

brain_skull_pca <- PCAtools::pca(
  mat = UV_data.matrix,
  metadata = UV_data.matrix.meta,
  removeVar = 0.2
)

PCAtools::biplot(
  pcaobj = brain_skull_pca,
  colby = "seq",
  colLegendTitle = "Seq",
  shape = "case",
  shapeLegendTitle = "Case",
  legendPosition = "right",
  # encircle = TRUE,
  # encircleFill = TRUE,
  labSize = 5,
  # showLoadings = T
  pointSize = 5,
  sizeLoadingsNames = 5,
  ellipse = TRUE,
  ellipseType = 't',
  ellipseLevel = 0.95,
  ellipseFill = TRUE,
  ellipseAlpha = 1/4,
  ellipseLineSize = 1.0
) ->
  brain_skull_pca_plot;brain_skull_pca_plot

PCAtools::biplot(
  pcaobj = brain_skull_pca,
  colby = "seq",
  colLegendTitle = "Seq",
  shape = "case",
  shapeLegendTitle = "Case",
  legendPosition = "right",
  encircle = TRUE,
  encircleFill = TRUE,
  pointSize = 5,
  sizeLoadingsNames = 5,
) ->
  brain_skull_pca_plot

ggsave(
  filename = "pca-brain-skull.pdf",
  plot = brain_skull_pca_plot,
  path = "data/uvresultnew/",
  device = "pdf",
  width = 7, height = 5
)


# make se -----------------------------------------------------------------


se <- SummarizedExperiment::SummarizedExperiment(
  assays = UV_data.matrix,
  colData = UV_data.matrix.meta
)


readr::write_rds(
  x = se,
  file = file.path(
    "data/uvrdanew",
    "UV_count_matrix.se.rds.gz"
  )
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
# save.image(file = "data/uvrdanew/02-UV.rda")
load(file = "data/uvrdanew/02-UV.rda")
