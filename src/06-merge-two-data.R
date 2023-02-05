# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Feb  5 17:36:13 2023
# @DESCRIPTION: 06-merge-two-data.R

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

sc_sham_mcao <- readr::read_rds(file = "data/rda/project_sc.rds.gz")

sc_uv <- readr::read_rds(file = "scuvdata/rda/project_uvsc.rds.gz")

# body --------------------------------------------------------------------

project_sc <- dplyr::bind_rows(
  sc_sham_mcao,
  sc_uv
)





# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------