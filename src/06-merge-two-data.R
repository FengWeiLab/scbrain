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


sc_sham_mcao_uv <- readr::read_rds(
  file = "scuvdata/rda/project_sc_raw_all.rds.gz"
)
# body --------------------------------------------------------------------



# Normalization -----------------------------------------------------------


sc_sham_mcao_uv %>% 
  dplyr::mutate(
    scn = purrr::map(
      .x = sc,
      .f = function(.sc) {
        Seurat::FindVariableFeatures(
          object = Seurat::NormalizeData(
            object = .sc
          ),
          selection.method = "vst", 
          nfeatures = 2000
        )
      }
    )
  ) ->
  sc_sham_mcao_uv_scn

# select features taht repeatedly variable across datasets from integration

features <- Seurat::SelectIntegrationFeatures(
  object.list = sc_sham_mcao_uv_scn$scn
)

# Integration -------------------------------------------------------------

anchors <- Seurat::FindIntegrationAnchors(
  object.list = sc_sham_mcao_uv_scn$scn,
  anchor.features = features
)

sc_sham_mcao_uv_scn_integrated <- Seurat::IntegrateData(
  anchorset = anchors
)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(
  file = "scuvdata/rda/project_uvsc.rda"
)
