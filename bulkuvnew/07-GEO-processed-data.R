#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Dec 27 16:32:54 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


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

readr::write_tsv(
  x = UV_data,
  file = file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew/GEO-Bulk-RNAseq",
    "Bulk-RNA-seq-read-count.tsv"
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

readr::write_tsv(
  x = barcode_group,
  file = file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew/GEO-Bulk-RNAseq",
    "Bulk-RNA-seq-group.tsv"
  )
)


# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
