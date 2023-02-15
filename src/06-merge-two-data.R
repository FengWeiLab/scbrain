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
  file = "data/scuvrda/project_sc_raw_all.rds.gz"
)
# body --------------------------------------------------------------------



# Normalization -----------------------------------------------------------


sc_sham_mcao_uv %>% 
  dplyr::mutate(
    scn = purrr::map(
      .x = sc,
      .f = function(.sc) {
        
        .sc_sub <- subset(
          x = .sc, 
          subset = nFeature_RNA > 500 & 
            nFeature_RNA < 6000 &
            percent.mt < 25 &
            percent.ribo < 50 &
            Percent.Largest.Gene < 40
        )
        
        Seurat::FindVariableFeatures(
          object = Seurat::NormalizeData(
            object = .sc_sub
          ),
          selection.method = "vst",
          nfeatures = 3000
        )
        #
      }
    )
  ) ->
  sc_sham_mcao_uv_scn


readr::write_rds(
  x = sc_sham_mcao_uv_scn,
  file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
)


# save stat ---------------------------------------------------------------


sc_sham_mcao_uv_scn %>% 
  dplyr::mutate(
    stat = purrr::map2(
      .x = sc,
      .y = scn,
      .f = function(.x, .y) {
        .xd <- .x@meta.data
        .yd <- .y@meta.data
        
        .n_x_cells <- nrow(.xd)
        .mean_reads_per_cell <- mean(.xd$nCount_RNA)
        .median_gene_per_cell <- median(.xd$nFeature_RNA)
        .n_y_cells <- nrow(.yd)
        
        tibble::tibble(
          `estimated number of cells` = .n_x_cells,
          `mean reads per cell` = .mean_reads_per_cell,
          `median genes per cell` = .median_gene_per_cell,
          `number of cells after filtering` = .n_y_cells
        )
      }
    )
  ) %>% 
  dplyr::select(project, stat) %>% 
  tidyr::unnest(cols = stat) ->
  project_stat

writexl::write_xlsx(
  x = project_stat, 
  path = "data/scuvresult/01-basic/reads-stat.xlsx"
  )

# to lsit -----------------------------------------------------------------



sc_sham_mcao_uv_scn %>% 
  dplyr::select(project, scn) %>% 
  tibble::deframe() ->
  sc_sham_mcao_uv_scn_list

# select features that repeatedly variable across data set from integration

features <- Seurat::SelectIntegrationFeatures(
  object.list = sc_sham_mcao_uv_scn_list,
  nfeatures = 3000
)

# Integration -------------------------------------------------------------

anchors <- Seurat::FindIntegrationAnchors(
  object.list = sc_sham_mcao_uv_scn_list,
  anchor.features = features,
  k.filter = NA
)


sc_sham_mcao_uv_scn_integrated <- Seurat::IntegrateData(
  anchorset = anchors,
  features = features,
  dims = 1:20
)

readr::write_rds(
  sc_sham_mcao_uv_scn_integrated,
  file = "data/scuvrda/sc_sham_mcao_uv_scn_integrated.rds.gz"
)

Seurat::DefaultAssay(sc_sham_mcao_uv_scn_integrated) <- "integrated"

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(
  file = "data/scuvrda/06-merge-two-data.rda"
)

load(file = "data/scuvrda/06-merge-two-data.rda")
