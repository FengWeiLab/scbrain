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
library(Seurat)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


fn_filter_sct <- function(.sc) {
  .sc_sub <- subset(
    x = .sc, 
    subset = nFeature_RNA > 500 & 
      nFeature_RNA < 6000 &
      percent.mt < 25 &
      percent.ribo < 50 &
      Percent.Largest.Gene < 40
  )
  
  .sc_sub_sct <- Seurat::SCTransform(
    object = .sc_sub,
    do.scale = FALSE,
    do.center = FALSE
  )
  
  .sc_sub_sct <- Seurat::CellCycleScoring(
    object = .sc_sub_sct,
    g2m.features = Seurat::cc.genes$g2m.genes,
    s.features = Seurat::cc.genes$s.genes
  )
  
  .sc_sub_sct$CC.Difference <- .sc_sub_sct$S.Score - .sc_sub_sct$G2M.Score
  
  Seurat::DefaultAssay(.sc_sub_sct) <- "RNA"
  
  .sc_sub_sct_sct <- Seurat::SCTransform(
    object = .sc_sub_sct,
    method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"),
    do.scale = TRUE,
    do.center = TRUE
  )
  .sc_sub_sct_sct
  
}

fn_filter_sc <- function(.sc) {
  
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
    nfeatures = 2000
  )
  #
}

fn_filter_sct_v2 <- function(.sc) {
  .sc_sub <- subset(
    x = .sc, 
    subset = nFeature_RNA > 500 & 
      nFeature_RNA < 6000 &
      percent.mt < 25 &
      percent.ribo < 50 &
      Percent.Largest.Gene < 40
  )
  
  
  Seurat::SCTransform(
    .sc,
    vst.flavor = "v2", 
    verbose = FALSE
  ) %>% 
    Seurat::RunPCA(
      npcs = 30, 
      verbose = FALSE
    )
}

# load data ---------------------------------------------------------------


sc_sham_mcao_uv <- readr::read_rds(
  file = "data/scuvrda/project_sc_raw_all.rds.gz"
)

# body --------------------------------------------------------------------



# Normalization -----------------------------------------------------------
sc_sham_mcao_uv$sc 
  purrr::map(nrow)

sc_sham_mcao_uv %>% 
  dplyr::mutate(
    scn = purrr::map(
      .x = sc,
      .f = fn_filter_sct_v2
    )
  ) ->
  sc_sham_mcao_uv_scn

sc_sham_mcao_uv_scn$scn %>% 
  purrr::map(nrow)

# readr::write_rds(
#   x = sc_sham_mcao_uv_scn,
#   file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
# )
# 
# 
# sc_sham_mcao_uv_scn <- readr::read_rds(
#   file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
# )

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
  tibble::deframe() %>% 
  purrr::map(
    .f = function(.x) {
      .t <- .x$tissue %>% unique()
      
      Seurat::RenameCells(object = .x, add.cell.id = .t)
    }
  ) ->
  sc_sham_mcao_uv_scn_list

sc_sham_mcao_uv_scn_list <- sc_sham_mcao_uv_scn_list[c(1, 4)]


features <- Seurat::SelectIntegrationFeatures(
  object.list = sc_sham_mcao_uv_scn_list,
  nfeatures = 3000
)

sc_list <- Seurat::PrepSCTIntegration(
  object.list = sc_sham_mcao_uv_scn_list,
  anchor.features = features
)

anchors <- Seurat::FindIntegrationAnchors(
  object.list = sc_list,
  normalization.method = "SCT",
  anchor.features = features
)

sc_sct <- Seurat::IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

readr::write_rds(
  x = sc_sct,
  file = "data/scuvrda/sc_sct.rds.gz"
)


# Integration -------------------------------------------------------------


a <- sc_sham_mcao_uv_scn_list[1:6] %>% 
  purrr::map(
    .f = function(.x) {
      .t <- .x$tissue %>% unique()
      
      Seurat::RenameCells(object = .x, add.cell.id = .t)
    }
  )

# a_genename <- a %>% 
#   purrr::map(rownames)

a_feature <- Seurat::SelectIntegrationFeatures(
  object.list = a,
  # nfeatures = 3000
 
)


# ggvenn::ggvenn(
#   data = c(
#     list(feature = a_feature),
#     a_genename[4:6]
#   )
# )


a_anchors <- Seurat::FindIntegrationAnchors(
  object.list = a,
  anchor.features = features,
  # k.filter = NA
)
  

anchors <- Seurat::FindIntegrationAnchors(
  object.list = a,
  anchor.features = features,
  # k.filter = NA
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
