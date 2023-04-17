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


# body --------------------------------------------------------------------

project_sct |> 
  dplyr::select(
    project, region, case, 
    unique_celltype, unique_celltype_union,
    anno_cell, cellnumber
  ) ->
  project_sct_sel


project_sct_sel

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------