# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Oct 16 16:20:07 2022
# @DESCRIPTION: UV

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# srouce ----------------------------------------------------------------
source(file = "utils/utils.R")


# load data ---------------------------------------------------------------
UV_data <- readr::read_rds(
  file = file.path(
    "data/uvrda",
    "all_expr_count_matrix_uv.rds.gz"
  )
) %>% 
  dplyr::select(
    1, 
    dplyr::contains(match = "BC"),
    dplyr::contains(match = "UVB0"),
    dplyr::contains(match = "UVB1"),
    dplyr::contains(match = "SC"),
    dplyr::contains(match = "UVS0"),
    dplyr::contains(match = "UVS1")
  )


tibble::tibble(
  barcode = colnames(UV_data)[-1]
) %>% 
  dplyr::mutate(
    group = gsub(pattern = "_\\d", replacement = "", x = barcode)
  ) %>% 
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "B", x = group), "Brain", "Skull")
  ) %>% 
  dplyr::mutate(
    case = ifelse(group == "BC", "Brain Control", group)
  ) %>% 
  dplyr::mutate(case = ifelse(case == "SC", "Skull Control", case)) ->
  barcode_group


UV_data %>% 
  tibble::column_to_rownames("GeneName") %>% 
  as.matrix() ->
  UV_data.matrix

barcode_group %>% 
  dplyr::mutate(new_barcode = barcode) %>% 
  tibble::column_to_rownames("new_barcode") ->
  UV_data.matrix.meta


brain_skull_pca <- PCAtools::pca(
  mat = UV_data.matrix,
  metadata = UV_data.matrix.meta,
  removeVar = 0.2
)

PCAtools::biplot(
  pcaobj = brain_skull_pca,
  colby = "seq",
  colLegendTitle = "Seq",
  shape = "case",
  shapeLegendTitle = "Case",
  legendPosition = "right",
  encircle = TRUE, 
  encircleFill = TRUE,
  lab = NULL
) ->
  brain_skull_pca_plot

ggsave(
  filename = "pca-brain-skull.pdf",
  plot = brain_skull_pca_plot,
  path = "data/uvresult/",
  device = "pdf",
  width = 7, height = 5
)


se <- SummarizedExperiment::SummarizedExperiment(
  assays = UV_data.matrix,
  colData = UV_data.matrix.meta
)

readr::write_rds(
  x = se,
  file = file.path(
    "data/uvrda",
    "UV_count_matrix.se.rds.gz"
  )
)

# seq brain
bb <- c("BC", "UVB0", "UVB1")
ss <- c("SC", "UVS0", "UVS1")


# save image --------------------------------------------------------------

save.image(file = "data/uvrda/13-uv-count-se.rda")
