# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar 29 14:45:00 2023
# @DESCRIPTION: filename

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

fn_load_sc_10x <- function(.x) {
  # .x <- "scuvdata/UVB"
  .regions <- c(
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
  .cases <- c(
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
  
  .project <- basename(.x)
  .case <- plyr::revalue(x = .project, replace = .cases)
  .region <- plyr::revalue(x = .project, replace = .regions)
  
  .dir <- file.path(.x, "filtered_feature_bc_matrix")
  .counts <- Seurat::Read10X(data.dir = .dir)
  
  .sc <- Seurat::CreateSeuratObject(
    counts = .counts,
    project = .project,
    min.cells = 3, 
    min.features = 200
  )
  
  .sc$tissue <- .project
  .sc$case <- .case
  .sc$region <- .region
  
  .sc <- Seurat::PercentageFeatureSet(
    object = .sc, 
    pattern = "^mt-", 
    col.name = "percent.mt"
  )
  .sc <- Seurat::PercentageFeatureSet(
    object = .sc, 
    # pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", 
    pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",
    col.name = "percent.ribo"
  )
  apply(
    .sc@assays$RNA@counts,
    2,
    function(x) (100 * max(x)) / sum(x)
  ) ->
    .sc$Percent.Largest.Gene
  
  .sc
}


# load data ---------------------------------------------------------------


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

refs <- c(
  "Brain" = "mousecortexref",
  "Meninge" = "pbmcref",
  "Skull" = "bonemarrowref"
)

dplyr::bind_rows(
  tibble::tibble(
    dir_path = list.dirs(
      file.path(
        "data/raw/singlecell/report", 
        "1_Cellranger"
      ),
      recursive = FALSE
    )
  ),
  tibble::tibble(
    dir_path = list.dirs(
      file.path("scuvdata"),
      recursive = FALSE
    )
  ) 
) %>% 
  dplyr::filter(!grepl(pattern = "rda", x = dir_path)) %>% 
  dplyr::mutate(project = basename(dir_path)) |> 
  dplyr::mutate(
    region = plyr::revalue(
      x = project,
      replace = regions
    )
  ) |> 
  dplyr::mutate(
    case = plyr::revalue(
      x = project,
      replace = cases
    )
  ) |> 
  dplyr::mutate(ref = plyr::revalue(
    x = region,
    replace = refs
  )) ->
  project_path


project_sc <- project_path %>% 
  dplyr::mutate(
    sc = furrr::future_map(
      .x = dir_path,
      .f = fn_load_sc_10x
    )
  ) 




# body --------------------------------------------------------------------



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------