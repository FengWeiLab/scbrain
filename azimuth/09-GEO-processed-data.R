#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Dec 27 16:41:55 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
library(Seurat)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
project_sc <- readr::read_rds(
  file = "data/azimuth/project_sc.rds.gz"
)

outdir <- "data/azimuth/GEO-scRNAseq"

dir.create(outdir)
# body --------------------------------------------------------------------
project_sc$sc[[1]]
project_sc$dir_path[[1]]

project_sc |>
  dplyr::select(-sc, -ref) |>
  dplyr::mutate(targetdir = glue::glue("{region}_{case}_{project}")) |>
  dplyr::select(dir_path, project, targetdir) ->
  project_path


project_path |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .x = dir_path,
        .y = targetdir,
        .z = project
      ),
      .f = \(.x, .y, .z) {
        # .x <- project_path$dir_path[[1]]
        # .y <- project_path$targetdir[[1]]

        .td <- file.path(outdir, .y)
        # dir.create(.td, showWarnings = F, recursive = T)

        file.copy(
          from = .x,
          to = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/scbrian-geo",
          recursive = T,
        )

        file.rename(
          from = file.path("/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/scbrian-geo", .z),
          to = .td
        )
      }
    )
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
