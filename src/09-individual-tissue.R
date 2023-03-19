# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Mar 19 14:39:34 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)
library(HGNChelper)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_cluster <- function(.sct) {
  .sct %>%
    Seurat::RunPCA(npcs = 30) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:30) %>%
    Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>%
    Seurat::FindClusters(resolution = c(0.3)) ->
    .sct_cluster
  .sct_cluster <- Seurat::PrepSCTFindMarkers(object = .sct_cluster)
  .sct_cluster
}

fn_find_all_markers <- function(.sct_cluster) {
  Seurat::FindAllMarkers(
    object = .sct_cluster,
    assay = "SCT",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

fn_heatmap10 <- function(.sct_cluster,.region,.case,.allmarkers, .outdir){
  # .sct_cluster <- project_sct_region_cluster_allmarkers$sct_cluster[[1]]
  # .region <- project_sct_region_cluster_allmarkers$region[[1]]
  # .case <- project_sct_region_cluster_allmarkers$case[[1]]
  # .allmarkers <- project_sct_region_cluster_allmarkers$allmarkers[[1]]
  # .outdir <-  "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
  # 
  .filename <- glue::glue("{.region}_{.case}_heatmap_top10.pdf")
  
  .allmarkers |> 
    dplyr::group_by(cluster) |> 
    dplyr::top_n(10, wt = avg_log2FC) ->
    .top10
  
  p <- DoHeatmap(.sct_cluster, features = .top10$gene) + NoLegend()
  
  ggsave(
    filename = .filename,
    plot = p,
    device = "pdf",
    path = .outdir,
    width = 12,
    height = 15
  )

}

# load data ---------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


regions <- c(
  "UVB" = "Brain",
  "UVM" = "Meninge",
  "UVS" = "Skull",
  "CB" = "Brain", 
  "CM" = "Meninge", 
  "CS" = "Skull", 
  "DB" = "Brain", 
  "DM" = "Meninge", 
  "DS" = "Skull"
)
cases <- c(
  "UVB" = "UV",
  "UVM" = "UV",
  "UVS" = "UV",
  "CB" = "Sham", 
  "CM" = "Sham", 
  "CS" = "Sham", 
  "DB" = "MCAO", 
  "DM" = "MCAO", 
  "DS" = "MCAO"
)

# individuals <- c(
#   "UVB" = "Brain UV",
#   "UVM" = "Meninge UV",
#   "UVS" = "Skull UV",
#   "CB" = "Brain Sham ", 
#   "CM" = "Meninge Sham", 
#   "CS" = "Skull Sham", 
#   "DB" = "Brain MCAO", 
#   "DM" = "Meninge MCAO", 
#   "DS" = "Skull MCAO"
# )

sc_sct_raw<- readr::read_rds(
  file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
)

sc_sct_raw |> 
  dplyr::mutate(
    region = plyr::revalue(x = project, replace = regions),
    case = plyr::revalue(x = project, replace = cases),
    # individual = plyr::revalue(x = project, replace = individuals)
  ) ->
  project_sct_region


# body --------------------------------------------------------------------


# cluster -----------------------------------------------------------------


project_sct_region |> 
  dplyr::mutate(
    sct_cluster = purrr::map(
      .x = sct,
      .f = fn_cluster
    )
  ) ->
  project_sct_region_cluster


project_sct_region_cluster |> 
  dplyr::mutate(
    allmarkers = purrr::map(
      .x = sct_cluster,
      .f = fn_find_all_markers
    )
  ) ->
  project_sct_region_cluster_allmarkers

readr::write_rds(
  x = project_sct_region_cluster_allmarkers,
  file = "data/scuvrda/project_sct_region_cluster_allmarkers.rds.gz"
)



# Heatmap top10 -----------------------------------------------------------

project_sct_region_cluster_allmarkers |> 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .sct_cluster = sct_cluster,
        .region = region,
        .case = case,
        .allmarkers = allmarkers
      ),
      .f = fn_heatmap10,
      .outdir = "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
    )
  )



# annotation --------------------------------------------------------------

a <- project_sct_region_cluster_allmarkers$sct_cluster[[1]]


a


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/scuvrda/09-individual-tissue.rda.gz")
