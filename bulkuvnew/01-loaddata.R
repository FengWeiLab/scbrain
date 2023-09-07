#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Sep  6 20:54:26 2023
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
dedir <- "/home/liuc9/data/projnet/2022-02-08-single-cell/bulknew/HGC20230803001-0001_report/Result/03.DiffExpresion"

depaths <- tibble::tibble(
  vsfilepath = list.files(
    path = dedir,
    pattern = "*.All.xls",
    recursive = T,
    full.names = T
  )
)

# body --------------------------------------------------------------------

depaths |>
  dplyr::mutate(
    loaddata = purrr::map(
      .x = vsfilepath,
      .f = \(.vsfilepath) {
        data.table::fread(input = .vsfilepath) |>
          dplyr::select(1, 5:10) |>
          dplyr::rename(GeneName = AccID)
      }
    )
  ) |>
  dplyr::select(loaddata) ->
  depaths_loadeddata

depaths_loadeddata |>
  dplyr::pull(loaddata) |>
  purrr::map(.f = \(.x) {.x$GeneName}) |>
  purrr::reduce(.f = intersect) ->
  allgenenames

depaths_loadeddata |>
  dplyr::pull(loaddata) |>
  purrr::map(.f = \(.x) {
    .x |> dplyr::slice(match(allgenenames, GeneName))
  }) |>
  purrr::reduce(.f = dplyr::inner_join, by = "GeneName") ->
  depaths_loadeddata_join



cols <- tibble::tibble(
  rawname = colnames(depaths_loadeddata_join)
) |>
  tibble::rowid_to_column(var = "i") |>
  dplyr::mutate(
    sucname = gsub(
      pattern = "\\..*",
      replacement = "",
      rawname
    )
  ) |>
  dplyr::distinct(sucname, .keep_all = T)



depaths_loadeddata_join |>
  dplyr::select(cols$i) ->
  depaths_loadeddata_join_sel

# relabel UVB1_1 -> UVM0_3
newname <- cols$sucname
stringr::str_match(newname, pattern = "UVB1_1")
stringr::str_match(newname, pattern = "UVM0_3")

newname[14] <- "UVM0_3"
newname[19] <- "UVB1_1"


colnames(depaths_loadeddata_join_sel) <- newname

depaths_loadeddata_join_sel |>
  readr::write_rds(
    file.path(
      "/home/liuc9/github/scbrain/data/uvrdanew",
      "all_expr_count_matrix_uv.rds.gz"
    )
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(
  file = file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew",
    "01-loaddata.rda"
  )
)
