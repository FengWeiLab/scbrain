#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Sep  6 22:33:19 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)
library(DESeq2)
# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
UV.se <- readr::read_rds(
  file = file.path(
    "data/uvrdanew",
    "UV_count_matrix.se.rds.gz"
  )
)

colData(UV.se) %>% as.data.frame() -> se_coldata




# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
