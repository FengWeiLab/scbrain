# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun May 21 18:33:38 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
azimuth_ref <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/project_sc_azimuth_ref_realcell_sunburst.rds"
)


azimuth_ref$anno_new |>
  purrr::map(
    .f = function(.x) {
      .x@meta.data |>
        dplyr::select(celltype.l1, celltype.l2) |>
        dplyr::distinct() |>
        tibble::as_tibble()
    }
  ) |>
  dplyr::bind_rows() |>
  dplyr::distinct() ->
  celltypes

recell3 <- c(
  # Myeloid

  "Neutrophils" = "Neutrophils",
  "Proliferating neutrophils" = "Neutrophils",

  "Basophils" = "Neutrophils",
  "Mast cells" = "Neutrophils",

  "cl. Monocytes" = "CD14 Monocytes",
  "Non-cl. monocytes" = "CD16 Monocytes",
  "MC" = "CD14 Monocytes",
  "MdC" = "CD14 Monocytes",
  "CD14 Mono" = "CD14 Monocytes",
  "CD16 Mono" = "CD16 Monocytes",

  "Monocytes" = "CD14 Monocytes",

  "D-BAM" = "Macrophage",
  "Macrophage" = "Macrophage",

  "Macrophages" = "Macrophage",

  "pDC" = "pDC",
  "cDC2" = "mDC",
  "cDC1" = "mDC",
  "migDC" = "mDC",
  "pre-mDC" = "mDC",
  "pre-pDC" = "pDC",
  "ASDC" = "mDC",

  "DCs cells" = "mDC",

  # Lymphoctyes
  "B cells" = "Mature B cells",
  "Memory B" = "Memory B cells",
  "transitional B" = "Naive B cells",
  "pre B" = "Pre-B cells",
  "Naive B" = "Naive B cells",
  "Plasma" = "Mature B cells",
  "pro B" = "Pro-B cells",

  "Mature B cells" = "Mature B cells",
  "Immature B cells" = "Immature B cells",
  "Pro-B cells" = "Pro-B cells",
  "Pre-B cells" = "Pre-B cells",
  "Plasma cells" = "Mature B cells",


  "CD4 Naive" = "CD4 naive T cells",
  "CD4 Memory" = "CD4 memory T cells",
  "CD8 Memory" = "CD8 memory T cells",
  "CD8 Effector_1" = "CD8 effector T cells",
  "CD8 Effector_3" = "CD8 effector T cells",
  "CD8 Effector_2" = "CD8 effector T cells",
  "T Proliferating" = "CD8 T cells",

  "T cells" = "T cells",


  "NK" = "NK cells",
  "NK cells" =  "NK cells",
  "T/NKT cells" = "T/NKT cells",

  "Stromal" = "other",
  "LMPP" = "HSPC",
  "GMP" = "HSPC",
  "BaEoMa" = "HSPC",
  "Late Eryth" = "Erythroid",
  "Early Eryth" = "Erythroid",

  "Erythroblasts" = "Erythroid",

  "Prog Mk" = "HSPC",
  "CLP" = "HSPC",
  "Platelet" = "Platelet",
  "Endo" = "Endothelial",

  "HSCs" = "HSPC",

  "Micro" = "Microglia",
  "PVM_1" = "Microglia",

  "Astro Aqp4_Slc7a10" = "Astrocyte Aqp4_Slc7a10",
  "Astro Aqp4_Gfap" = "Astrocyte Aqp4_Gfap",


  "VLMC_6" = "VLMC",
  "VLMC_2" = "VLMC",
  "VLMC_3" = "VLMC",

  "Oligo Opalin_1" = "OLG",
  "Oligo Opalin_2" = "OLG",
  "Oligo Enpp6_2" = "OLG",
  "Oligo Opalin_4" = "OLG",

  "OPC Pdgfra" = "OPC",

  "Sensory neurons" = "Sensory neurons",
  "Peri" = "Pericytes",
  "Meis2" = "Neuron",
  "L6 IT" = "Neuron",
  "Meis2_Top2a" = "Neuron",
  "L6 IT_1" = "Neuron"
  )

recell2 <- c(
  # Myeloid

  "Neutrophils" = "Neutrophils",
  "Proliferating neutrophils" = "Neutrophils",

  "Basophils" = "Neutrophils",
  "Mast cells" = "Neutrophils",

  "cl. Monocytes" = "Monocytes",
  "Non-cl. monocytes" = "Monocytes",
  "MC" = "Monocytes",
  "MdC" = "Monocytes",
  "CD14 Mono" = "Monocytes",
  "CD16 Mono" = "Monocytes",

  "Monocytes" = "Monocytes",

  "D-BAM" = "Macrophage",
  "Macrophage" = "Macrophage",

  "Macrophages" = "Macrophage",

  "pDC" = "DC",
  "cDC2" = "DC",
  "cDC1" = "DC",
  "migDC" = "DC",
  "pre-mDC" = "DC",
  "pre-pDC" = "DC",
  "ASDC" = "DC",

  "DCs cells" = "DC",

  # Lymphoctyes
  "B cells" = "B cells",
  "Memory B" = "B cells",
  "transitional B" = "B cells",
  "pre B" = "B cells",
  "Naive B" = "B cells",
  "Plasma" = "B cells",
  "pro B" = "B cells",

  "Mature B cells" = "B cells",
  "Immature B cells" = "B cells",
  "Pro-B cells" = "B cells",
  "Pre-B cells" = "B cells",
  "Plasma cells" = "B cells",

  "CD4 Naive" = "T cells",
  "CD4 Memory" = "T cells",
  "CD8 Memory" = "T cells",
  "CD8 Effector_1" = "T cells",
  "CD8 Effector_3" = "T cells",
  "CD8 Effector_2" = "T cells",
  "T Proliferating" = "T cells",

  "T cells" = "T cells",

  "NK" = "NK cells",
  "NK cells" =  "NK cells",
  "T/NKT cells" = "T/NKT cells",

  "Stromal" = "other",
  "LMPP" = "HSPC",
  "GMP" = "HSPC",
  "BaEoMa" = "HSPC",
  "Late Eryth" = "Erythroid",
  "Early Eryth" = "Erythroid",

  "Erythroblasts" = "Erythroid",

  "Prog Mk" = "HSPC",
  "CLP" = "HSPC",
  "Platelet" = "Platelet",
  "Endo" = "Endothelial",

  "HSCs" = "HSPC",

  "Micro" = "Microglia",
  "PVM_1" = "Microglia",

  "Astro Aqp4_Slc7a10" = "Astrocyte",
  "Astro Aqp4_Gfap" = "Astrocyte",


  "VLMC_6" = "VLMC",
  "VLMC_2" = "VLMC",
  "VLMC_3" = "VLMC",

  "Oligo Opalin_1" = "OLG",
  "Oligo Opalin_2" = "OLG",
  "Oligo Enpp6_2" = "OLG",
  "Oligo Opalin_4" = "OLG",

  "OPC Pdgfra" = "OPC",


  "Sensory neurons" = "Sensory neurons",
  "Peri" = "Pericytes",
  "Meis2" = "Neuron",
  "L6 IT" = "Neuron",
  "Meis2_Top2a" = "Neuron",
  "L6 IT_1" = "Neuron"
)


recell1 <- c(
  "Astrocyte" = "Astrocyte",
  "B cells" = "Lymphoctyes",
  "DC" = "Myeloid cells",
  "Endothelial" = "Endothelial",
  "Erythroid" = "Erythroid",
  "HSPC" = "HSPC",
  "Macrophage" = "Myeloid cells",
  "Microglia" = "Microglia",
  "Monocytes" = "Myeloid cells",
  "Neuron" = "Neuron",
  "Neutrophils" = "Myeloid cells",
  "NK cells" = "Lymphoctyes",
  "OLG" = "OLG",
  "OPC" = "OPC",
  "other" = "other",
  "Pericytes" = "Pericytes",
  "Platelet" = "Erythroid",
  "Sensory neurons" = "Neuron",
  "T cells" = "Lymphoctyes",
  "T/NKT cells" = "Lymphoctyes",
  "VLMC" = "VLMC"
)


celltypes |>
  dplyr::mutate(
    cell3 = plyr::revalue(
      x = celltype.l2,
      replace = recell3
    )
  ) |>
  dplyr::mutate(
    cell2 = plyr::revalue(
      x = celltype.l2,
      replace = recell2
    )
  ) |>
  dplyr::mutate(
    cell1 = plyr::revalue(
      x = cell2,
      replace = recell1
    )
  ) |>
  dplyr::select(-celltype.l1) |>
  dplyr::distinct() ->
  celltypes_recell

celltypes_recell |>
  dplyr::select(cell1, cell2, cell3) |>
  dplyr::arrange(cell1, cell2, cell3) |>
  dplyr::distinct() ->
  recells



paletteer::palettes_d_names |>
  dplyr::filter(grepl("material", palette)) |>
  dplyr::mutate(
    color_name = glue::glue("{package}::{palette}")
  ) |>
  dplyr::arrange(length) |>
  dplyr::select(color_name) |>
  dplyr::mutate(
    color = purrr::map(
      .x = color_name,
      .f = paletteer::paletteer_d
    )
  ) ->
  d

d |>
  tibble::deframe()


recells |>
  dplyr::count(cell1, cell2) |>
  dplyr::arrange(cell1, -n) |>
  dplyr::mutate(
    color_name = c(
      "ggsci::purple_material",
      "ggsci::brown_material",
      "ggsci::grey_material",
      "ggsci::blue_grey_material",
      "ggsci::yellow_material",
      "ggsci::pink_material",
      "ggsci::red_material", # T cell
      "ggsci::red_material", # T cell,
      # "ggsci::purple_material",
      "ggsci::red_material",
      "ggsci::light_blue_material",
      "ggsci::blue_material",
      "ggsci::cyan_material",
      "ggsci::teal_material",
      "ggsci::green_material", # "Neuron"
      "ggsci::green_material", # "Sensory neurons"
      "ggsci::deep_purple_material",
      "ggsci::deep_purple_material",
      "ggsci::indigo_material",
      "ggsci::blue_grey_material"
    )
  ) |>
  dplyr::mutate(
    color = purrr::map(
      .x = color_name,
      .f = paletteer::paletteer_d
    )
  ) |>
  # dplyr::mutate(
  #   m = glue::glue("{cell2}, {n}")
  # ) |>
  # dplyr::select(m, color) |>
  # tibble::deframe()
  dplyr::mutate(
    ns = c(
      list(c(10, 9, 5)),
      list(c(3)),
      list(c(3)),
      list(c(3)),
      # list(c(3)),
      list(c(7, 9, 7, 5, 3)),
      list(c(10)),
      list(c(10)), #T,
      list(c(10)), #T
      list(c(7)),
      # list(c(10)),
      list(c(10, 7, 5)),
      list(c(10, 9, 5)),
      list(c(10)),
      list(c(10)), # Neutrophil
      # list(c(6)),
      list(c(3)), # Neurons
      list(c(4)), # Sensory neurons
      list(c(4)),
      list(c(5)),
      list(c(2)),
      list(c(1))
    )
  ) |>
  dplyr::mutate(
    color_select = purrr::map2(
      .x = color,
      .y = ns,
      .f = function(.x, .y) {
        .x[.y]
      }
    )
  ) |>
  # dplyr::mutate(
  #   g = glue::glue("{cell2} {n}")
  # ) |>
  # dplyr::select(g, color_select)
  dplyr::select(cell2, color_select) |>
  dplyr::left_join(
    recells,
    by = "cell2"
  ) |>
  dplyr::mutate(
    cell2_color = purrr::map(
      .x = color_select,
      .f = function(.x) {
        .x[1]
      }
    )
  ) |>
  tidyr::unnest(cell2_color) ->
  recell_color

recell_color |>
  dplyr::select(cell2, cell3, color_select) |>
  dplyr::group_by(cell2) |>
  dplyr::filter(dplyr::n() > 1) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    cell3_color = purrr::map(
      .x = data,
      .f = function(.x) {
        .color <- .x$color_select[[1]][-1]
        .x |>
          dplyr::select(-color_select) |>
          dplyr::mutate(cell3_color = .color)
      }
    )
  ) |>
  dplyr::select(cell3_color) |>
  tidyr::unnest(cols = cell3_color) ->
  recell_color_m


recell_color |>
  dplyr::left_join(
    recell_color_m,
    by = "cell3"
  ) |>
  dplyr::mutate(
    cell3_color = ifelse(
      cell3_color == "#FFFFFF00",
      cell2_color,
      cell3_color
    )
  ) |>
  dplyr::mutate(
    cell1_color = cell2_color
  ) |>
  dplyr::mutate(
    cell1_color = ifelse(
      cell1 == "Erythroid",
      "#EEEEEEFF",
      cell1_color
    )
  ) |>
  dplyr::mutate(
    cell1_color = ifelse(
      cell1 == "Lymphoctyes",
      "#B71C1CFF",
      cell1_color
    )
  ) |>
  dplyr::mutate(
    cell1_color = ifelse(
      cell1 == "Myeloid cells",
      "#006064FF",
      cell1_color
    )
  ) |>
  dplyr::mutate(
    cell1_color = prismatic::clr_lighten(
      prismatic::clr_saturate(
        cell1_color,
        shift = 0.7
      ),
      shift = 0.4
    ),
    cell2_color = prismatic::clr_lighten(
      prismatic::clr_saturate(
        cell2_color,
        shift = 0.8
      ),
      shift = 0.3
    ),
    cell3_color = prismatic::clr_lighten(
      prismatic::clr_saturate(
        cell3_color,
        shift = 0.8
      ),
      shift = 0.3
    )
  ) |>
  dplyr::select(
    cell1, cell2, cell3, cell1_color, cell2_color, cell3_color
  ) |>
  dplyr::mutate(
    cell1 = ifelse(
      cell2 == "Platelet",
      "Platelet",
      cell1
    )
  ) |>
  dplyr::mutate(
    cell1_color = ifelse(
      cell2 == "Platelet",
      "#A0DDFFFF",
      cell1_color
    )
  ) ->
  recell_color_final

recell_color_final |>
  ggplot(aes(
    x = cell1,
    y = 1,
    fill = as.character(cell1_color)
  )) +
  geom_tile() +
  scale_fill_identity() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) ->
  pc1;pc1

recell_color_final |>
  ggplot(aes(
    x = cell2,
    y = 1,
    fill = as.character(cell2_color)
  )) +
  geom_tile() +
  scale_fill_identity() +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) ->
  pc2

recell_color_final |>
  ggplot(aes(
    x = cell3,
    y = 1,
    fill = as.character(cell3_color)
  )) +
  geom_tile() +
  scale_fill_identity() +
  # theme_bw() +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) ->
  pc3

pc1 / pc2 / pc3 ->
  pc123;pc123

ggsave(
  filename = "colors.pdf",
  plot = pc123,
  device = "pdf",
  width = 10,
  height = 8,
  path = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13"
)

readr::write_rds(
  x = recell_color_final,
  file = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13/recell_color.rds"
)

# body --------------------------------------------------------------------


# Sunburst ----------------------------------------------------------------


azimuth_ref |>
  dplyr::mutate(
    sunburst = purrr::pmap(
      .l = list(
        .region = project,
        .anno = anno_new
      ),
      .f = function(.region, .anno) {
        # .anno <- azimuth_ref$anno_new[[1]]
        # .region <- azimuth_ref$region[[1]]
        # print(.region)

        .anno@meta.data |>
          dplyr::left_join(
            celltypes_recell |>
              dplyr::mutate(
                cell1 = ifelse(
                  cell2 == "Platelet",
                  "Platelet",
                  cell1
                )
              ) |>
              dplyr::left_join(
                recell_color_final,
                by = c("cell3", "cell2", "cell1")
              ),
            by = "celltype.l2"
          ) ->
          .d

        .d |>
          dplyr::count(cell1, cell2, cell3) |>
          dplyr::mutate(cell3_r = n / sum(n)) |>
          dplyr::group_by(cell2) |>
          dplyr::mutate(cell2_r = sum(cell3_r)) |>
          dplyr::ungroup() |>
          dplyr::group_by(cell1) |>
          dplyr::mutate(cell1_r = sum(cell3_r)) |>
          dplyr::ungroup() |>
          dplyr::mutate(
            cell1 =  glue::glue("{cell1} {round(cell1_r * 100, 2)}%"),
            cell2 = glue::glue("{cell2} {round(cell2_r * 100, 2)}%"),
            cell3 = glue::glue("{cell3} {round(cell3_r * 100, 2)}%"),
          ) ->
          .dd

        .dd |>
          dplyr::select(1, 2,3,4) |>
          plotme::count_to_sunburst() ->
          .p

        .new_anno <- .anno
        .new_anno@meta.data$cell1 <- .d$cell1
        .new_anno@meta.data$cell2 <- .d$cell2
        .new_anno@meta.data$cell3 <- .d$cell3
        .new_anno@meta.data$cell1_color <- .d$cell1_color
        .new_anno@meta.data$cell2_color <- .d$cell2_color
        .new_anno@meta.data$cell3_color <- .d$cell3_color

        tibble::tibble(
          sunburst = list(.p),
          anno_new_new = list(.new_anno),
          cellratio = list(.dd)
        )

      }


    )
  ) |>
  tidyr::unnest(cols = sunburst) ->
  azimuth_ref_sunburst



azimuth_ref_sunburst |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .p = sunburst
      ),
      .f = function(.region, .case, .p, .outdir) {
        dir.create(
          path = .outdir,
          showWarnings = F,
          recursive = T
        )
        reticulate::py_run_string("import sys")

        .filename <- glue::glue("Sunburst_{.region}_{.case}")
        plotly::save_image(
          p = .p,
          file = file.path(
            .outdir,
            glue::glue("{.filename}.pdf")
          ),
          width = 800,
          height = 800,
          device = "pdf"
        )

        htmlwidgets::saveWidget(
          .p,
          file = file.path(
            .outdir,
            glue::glue("{.filename}.html")
          )
        )

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13"
    )
  )

# Cell level --------------------------------------------------------------


azimuth_ref_sunburst |>
  dplyr::select(region, case, cellratio) |>
  dplyr::group_by(region) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    cell_factor = purrr::map2(
      .x = data,
      .y = region,
      .f = function(.x, .y) {
        # .x <- ddd$data[[3]]
        # .y <- ddd$region[[3]]

        .x$cellratio |>
          dplyr::bind_rows() |>
          dplyr::select(cell1, cell2, cell3) |>
          dplyr::mutate_all(
            .funs = function(.s) {
              gsub(
                pattern = " [[0-9]]*.*%$",
                replacement = "",
                x = .s
              )
            }
          ) |>
          dplyr::distinct() ->
          .xx

        .xx$cell1 |> unique()
        .xx$cell2 |> unique()
        .xx$cell3 |> unique()
        # .y
        if(.y == "Meninge") {
          list(
            cell1 = c("Lymphoctyes", "Myeloid cells"),
            cell2 = c("T/NKT cells", "B cells", "NK cells", "Monocytes", "Macrophage", "DC", "Neutrophils"),
            cell3 = c("T/NKT cells", "Mature B cells", "NK cells", "CD14 Monocytes", "CD16 Monocytes", "Macrophage", "mDC", "pDC", "Neutrophils")
          )
        } else if(.y == "Skull") {
          list(
            cell1 = c("Lymphoctyes", "Myeloid cells", "HSPC", "Erythroid", "Neuron"),
            cell2 = c("T cells", "B cells", "NK cells", "Monocytes", "Macrophage", "DC", "Neutrophils",  "HSPC", "Erythroid", "Sensory neurons"),
            cell3 = c("T cells", "Mature B cells", "Immature B cells", "Pre-B cells", "Pro-B cells", "NK cells", "CD14 Monocytes", "Macrophage", "mDC", "Neutrophils", "HSPC", "Erythroid", "Sensory neurons")
          )
        } else if(.y == "Brain") {
          list(
            cell1 = c("Lymphoctyes", "Myeloid cells", "Astrocyte", "Microglia", "OLG", "OPC", "VLMC", "Endothelial", "Pericytes", "Neuron"),
            cell2 = c("T/NKT cells", "B cells", "NK cells", "Monocytes", "Macrophage", "DC", "Neutrophils", "Astrocyte", "Microglia", "OLG", "OPC", "VLMC", "Endothelial", "Pericytes", "Neuron"),
            cell3 = c("T/NKT cells", "Mature B cells", "NK cells", "CD14 Monocytes", "CD16 Monocytes", "Macrophage", "mDC", "Neutrophils", "Astrocyte Aqp4_Gfap", "Astrocyte Aqp4_Slc7a10", "Microglia", "OLG", "OPC", "VLMC", "Endothelial", "Pericytes", "Neuron" )
          )
        }

      }
    )
  ) ->
  azimuth_ref_sunburst_sel

azimuth_ref_sunburst |>
  dplyr::left_join(
    azimuth_ref_sunburst_sel |>
      tidyr::unnest(cols = data) |>
      dplyr::select(-cellratio),
    by = c("region", "case")
  ) ->
  azimuth_ref_sunburst_cell

# readr::write_rds(
#   x = azimuth_ref_sunburst_cell,
#   file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
# )

azimuth_ref_sunburst_cell <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
)

# Col plot ----------------------------------------------------------------



azimuth_ref_sunburst_sel |>
  dplyr::mutate(
    ratiop = purrr::pmap(
      .l = list(
        .x = region,
        .y = data,
        .cell_factor = cell_factor
      ),
      .f = function(.x, .y, .cell_factor) {
        # .x <- azimuth_ref_sunburst_sel$region[[1]]
        # .y <- azimuth_ref_sunburst_sel$data[[1]]
        # .cell_factor <- azimuth_ref_sunburst_sel$cell_factor[[1]]


        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell1, ratio = cell1_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) |>
          dplyr::mutate(
            cluster = factor(x = cluster, levels = .cell_factor$cell1)
          ) ->
          .d

        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell2, ratio = cell2_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) |>
          dplyr::mutate(
            cluster = factor(x = cluster, levels = .cell_factor$cell2)
          ) ->
          .dd

        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell3, ratio = cell3_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) |>
          dplyr::mutate(
            cluster = factor(x = cluster, levels = .cell_factor$cell3)
          ) ->
          .ddd

        # .scale_fill <- if (.x == "Brain") {
        #   ggthemes::scale_fill_tableau(
        #     palette = "Tableau 20",
        #     name = "Cell types",
        #     direction = 1
        #   )
        # } else if(.x == "Meninge") {
        #   ggsci::scale_fill_npg(
        #     name = "Cell types"
        #   )
        # } else {
        #   ggsci::scale_fill_npg(
        #     name = "Cell types"
        #   )
        # }
        recell_color_final |>
          dplyr::select(cell1, cell1_color) |>
          dplyr::distinct() |>
          dplyr::filter(cell1 %in% .d$cluster) |>
          dplyr::mutate(
            cell1 = factor(
              x = cell1,
              levels = .cell_factor$cell1
            )
          ) |>
          dplyr::arrange(cell1) ->
          .cell1_color

        .d |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 0.9,
            # color = 1,
            linewidth = 0,
            # alpha = 0.8
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            labels = c("Sham", "tMCAO", "tMCAO+UVB"),
            expand = expansion(mult= c(0.3, 0.3), add = 0)
          ) +
          scale_y_continuous(
            labels = scales::percent_format(),
            expand = c(0, 0.01)
          ) +
          scale_fill_manual(
            name = "Cell types",
            limits = .cell1_color$cell1,
            values =alpha(
              .cell1_color$cell1_color,
              # 0.7
            )
          ) +
          # .scale_fill +
          # ggthemes::scale_fill_tableau(
          #   palette = "Tableau 20",
          #   name = "Cell types",
          #   direction = 1
          # ) +
          theme(
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_line(),
            axis.text = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.title = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.text = element_text(
              size = 14,
              color = "black",
              face = "bold"
            ),
            axis.text.x = element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            )
          ) ->
          .p;.p

        recell_color_final |>
          dplyr::select(cell2, cell2_color) |>
          dplyr::distinct() |>
          dplyr::filter(cell2 %in% .dd$cluster) |>
          dplyr::mutate(
            cell2 = factor(
              x = cell2,
              levels = .cell_factor$cell2
            )
          ) |>
          dplyr::arrange(cell2) ->
          .cell2_color

        .dd |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 0.9,
            # color = 1,
            linewidth = 0,
            # alpha = 0.8
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            labels = c("Sham", "tMCAO", "tMCAO+UVB"),
            expand = expansion(mult= c(0.3, 0.3), add = 0)
          ) +
          scale_y_continuous(
            labels = scales::percent_format(),
            expand = c(0, 0.01)
          ) +
          scale_fill_manual(
            name = "Cell types",
            limits = .cell2_color$cell2,
            values =.cell2_color$cell2_color
          ) +
          theme(
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_line(),
            axis.text = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.title = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.text = element_text(
              size = 14,
              color = "black",
              face = "bold"
            ),
            axis.text.x = element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            )
          ) ->
          .pp;.pp

        recell_color_final |>
          dplyr::select(cell3, cell3_color) |>
          dplyr::distinct() |>
          dplyr::filter(cell3 %in% .ddd$cluster) |>
          dplyr::mutate(
            cell3 = factor(
              x = cell3,
              levels = .cell_factor$cell3
            )
          ) |>
          dplyr::arrange(cell3) ->
          .cell3_color

        .ddd |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 0.9,
            # color = 1,
            linewidth = 0,
            # alpha = 0.8
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            labels = c("Sham", "tMCAO", "tMCAO+UVB"),
            expand = expansion(mult= c(0.3, 0.3), add = 0)
          ) +
          scale_y_continuous(
            labels = scales::percent_format(),
            expand = c(0, 0.01)
          ) +
          scale_fill_manual(
            name = "Cell types",
            limits = .cell3_color$cell3,
            values =.cell3_color$cell3_color
          ) +
          # ggthemes::scale_fill_tableau(
          #   palette = "Tableau 20",
          #   name = "Cell types",
          #   direction = 1
          # ) +
          theme(
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.line = element_line(),
            axis.text = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.title = element_text(
              size = 16,
              color = "black",
              face = "bold"
            ),
            legend.text = element_text(
              size = 14,
              color = "black",
              face = "bold"
            ),
            axis.text.x = element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            )
          ) ->
          .ppp;.ppp


        list(
          p_cell1 = .p,
          p_cell2 = .pp,
          p_cell3 = .ppp,
          d_cell1 = .d,
          d_cell2 = .dd,
          d_cell3 = .ddd
        )
      }
    )
  ) ->
  azimuth_ref_sunburst_sel_ratiop

azimuth_ref_sunburst_sel_ratiop |>
  dplyr::mutate(
    a = purrr::map2(
      .x = region,
      .y = ratiop,
      .f = function(.region, .ratiop, .outdir) {
        ggsave(
          filename = glue::glue("Propertion_{.region}_cell1.2.pdf"),
          plot = .ratiop$p_cell1,
          device = "pdf",
          path = .outdir,
          width = 4,
          height = 7
        )

        ggsave(
          filename = glue::glue("Propertion_{.region}_cell2.2.pdf"),
          plot = .ratiop$p_cell2,
          device = "pdf",
          path = .outdir,
          width = 4,
          height = 7
        )

        wid <- ifelse(
          .region == "Brain",
          5,
          4
        )
        ggsave(
          filename = glue::glue("Propertion_{.region}_cell3.2.pdf"),
          plot = .ratiop$p_cell3,
          device = "pdf",
          path = .outdir,
          width = wid,
          height = 7
        )

        tox <- list(
          Cell1 = .ratiop$d_cell1,
          Cell3 = .ratiop$d_cell2,
          Cell3 = .ratiop$d_cell3
        )

        writexl::write_xlsx(
          x = tox,
          path = file.path(
            .outdir,
            glue::glue("Propertion_{.region}_cell_ratio.xlsx")
          )
        )


      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13"
    )
  )


# Merge -------------------------------------------------------------------


azimuth_ref_sunburst_sel_ratiop |>
  dplyr::mutate(
    a = purrr::map(
      .x = ratiop,
      .f = function(.x) {
        .x$d_cell1
      }
    )
  ) |>
  dplyr::select(region, a) |>
  tidyr::unnest(cols = a) |>
  dplyr::mutate(
    case = factor(x = case, levels = c("Sham", "MCAO", "UV"))
  ) ->
  azimuth_ref_sunburst_sel_ratiop_cell1_forplot

recell_color_final |>
  dplyr::select(cell1, cell1_color) |>
  dplyr::distinct() |>
  dplyr::mutate(
    cell1 = factor(
      x = cell1,
      levels = levels(azimuth_ref_sunburst_sel_ratiop_cell1_forplot$cluster)
    )
  ) |>
  dplyr::arrange(cell1) ->
  cell1_color

azimuth_ref_sunburst_sel_ratiop_cell1_forplot |>
  dplyr::group_by(region, case) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |>
          dplyr::ungroup() %>%
          dplyr::arrange(cluster) |>
          dplyr::mutate(csum = rev(cumsum(rev(ratio)))) %>%
          dplyr::mutate(pos = ratio/2 + dplyr::lead(csum, 1)) %>%
          dplyr::mutate(pos = dplyr::if_else(is.na(pos), ratio/2, pos)) %>%
          dplyr::mutate(percentage = ratio/sum(ratio))
      }
    )
  ) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = 1,
    y = percentage
  )) +
  geom_col(
    aes(
      fill = cluster
    ),
    width = 1.5,
    color = "white"
  ) +
  geom_col(aes(x = 0, y = 0)) +
  coord_polar(
    theta = "y",
    start = 0
  ) +
  ggrepel::geom_text_repel(
    aes(
      y = pos,
      label = glue::glue("{scales::percent(percentage, accuracy = 0.01)}"),
    ),
    size = 8,
    nudge_x = 0.2,
    show.legend = FALSE,
  ) +
  # ggsci::scale_fill_simpsons() +
  # paletteer::scale_fill_paletteer_d(
  #   palette = "ggsci::category20c_d3",
  #   direction = 1,
  #   name = "Cell types"
  # ) +
  scale_fill_manual(
    name = "Cell types",
    limits = cell1_color$cell1,
    values = cell1_color$cell1_color
  ) +
  facet_grid(
    rows = dplyr::vars(region),
    cols = dplyr::vars(case),
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      color = "black",
      size = 20,
      face = "bold"
    ),
    panel.margin = unit(0, "lines"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
  ) ->
  p1;p1


ggsave(
  filename = glue::glue("Pie_plot_cell1.pdf"),
  plot = p1,
  device = "pdf",
  path =  "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13",
  width = 15,
  height = 15
)

azimuth_ref_sunburst_sel_ratiop |>
  dplyr::mutate(
    a = purrr::map(
      .x = ratiop,
      .f = function(.x) {
        .x$d_cell2
      }
    )
  ) |>
  dplyr::select(region, a) |>
  tidyr::unnest(cols = a) |>
  dplyr::mutate(
    case = factor(x = case, levels = c("Sham", "MCAO", "UV"))
  ) ->
  azimuth_ref_sunburst_sel_ratiop_cell2_forplot

recell_color_final |>
  dplyr::select(cell2, cell2_color) |>
  dplyr::distinct() |>
  dplyr::mutate(
    cell2 = factor(
      x = cell2,
      levels = levels(azimuth_ref_sunburst_sel_ratiop_cell2_forplot$cluster)
    )
  ) |>
  dplyr::arrange(cell2) ->
  cell2_color


azimuth_ref_sunburst_sel_ratiop_cell2_forplot |>
  dplyr::group_by(region, case) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |>
          dplyr::ungroup() %>%
          dplyr::arrange(cluster) |>
          dplyr::mutate(csum = rev(cumsum(rev(ratio)))) %>%
          dplyr::mutate(pos = ratio/2 + dplyr::lead(csum, 1)) %>%
          dplyr::mutate(pos = dplyr::if_else(is.na(pos), ratio/2, pos)) %>%
          dplyr::mutate(percentage = ratio/sum(ratio))
      }
    )
  ) |>
  tidyr::unnest(cols = data) |>
  ggplot(aes(
    x = 1,
    y = percentage
  )) +
  geom_col(
    aes(
      fill = cluster
    ),
    width = 1.5,
    color = "white"
  ) +
  geom_col(aes(x = 0, y = 0)) +
  coord_polar(
    theta = "y",
    start = 0
  ) +
  ggrepel::geom_text_repel(
    aes(
      y = pos,
      label = glue::glue("{scales::percent(percentage, accuracy = 0.01)}"),
    ),
    size = 6,
    nudge_x = 0.2,
    show.legend = FALSE,
  ) +
  # ggsci::scale_fill_simpsons() +
  # paletteer::scale_fill_paletteer_d(
  #   palette = "ggsci::category20c_d3",
  #   direction = 1,
  #   name = "Cell types"
  # ) +
  scale_fill_manual(
    name = "Cell types",
    limits = cell2_color$cell2,
    values = cell2_color$cell2_color
  ) +
  facet_grid(
    rows = dplyr::vars(region),
    cols = dplyr::vars(case),
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      color = "black",
      size = 20,
      face = "bold"
    ),
    panel.margin = unit(0, "lines"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
  ) ->
  p2;p2


ggsave(
  filename = glue::glue("Pie_plot_cell2.pdf"),
  plot = p2,
  device = "pdf",
  path =  "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel13",
  width = 20,
  height = 20
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
# save.image(file = "data/azimuth/04-manual-annotation.rda")
load(file = "data/azimuth/04-manual-annotation.rda")
#
