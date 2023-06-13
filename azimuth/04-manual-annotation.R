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

  "cl. Monocytes" = "CD14 Monocytes",
  "Non-cl. monocytes" = "CD16 Monocytes",
  "MC" = "CD14 Monocytes",
  "MdC" = "CD14 Monocytes",
  "CD14 Mono" = "CD14 Monocytes",
  "CD16 Mono" = "CD16 Monocytes",

  "D-BAM" = "Macrophage",
  "Macrophage" = "Macrophage",

  "pDC" = "pDC",
  "cDC2" = "mDC",
  "cDC1" = "mDC",
  "migDC" = "mDC",
  "pre-mDC" = "mDC",
  "pre-pDC" = "pDC",
  "ASDC" = "mDC",

  # Lymphoctyes
  "B cells" = "Mature B cells",
  "Memory B" = "Memory B cells",
  "transitional B" = "Naive B cells",
  "pre B" = "Pre-B cells",
  "Naive B" = "Naive B cells",
  "Plasma" = "Mature B cells",
  "pro B" = "Pro-B cells",

  "CD4 Naive" = "CD4 naive T cells",
  "CD4 Memory" = "CD4 memory T cells",
  "CD8 Memory" = "CD8 memory T cells",
  "CD8 Effector_1" = "CD8 effector T cells",
  "CD8 Effector_3" = "CD8 effector T cells",
  "CD8 Effector_2" = "CD8 effector T cells",
  "T Proliferating" = "CD8 T cells",

  "NK" = "NK cells",
  "NK cells" =  "NK cells",
  "T/NKT cells" = "CD8 T cells",

  "Stromal" = "other",
  "LMPP" = "HSPC",
  "GMP" = "HSPC",
  "BaEoMa" = "HSPC",
  "Late Eryth" = "Erythroid",
  "Early Eryth" = "Erythroid",
  "Prog Mk" = "HSPC",
  "CLP" = "HSPC",
  "Platelet" = "Platelet",
  "Endo" = "Endothelial",

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

  "OPC Pdgfra" = "OPC"
  )

recell2 <- c(
  # Myeloid

  "Neutrophils" = "Neutrophils",

  "cl. Monocytes" = "Monocytes",
  "Non-cl. monocytes" = "Monocytes",
  "MC" = "Monocytes",
  "MdC" = "Monocytes",
  "CD14 Mono" = "Monocytes",
  "CD16 Mono" = "Monocytes",

  "D-BAM" = "Macrophage",
  "Macrophage" = "Macrophage",

  "pDC" = "DC",
  "cDC2" = "DC",
  "cDC1" = "DC",
  "migDC" = "DC",
  "pre-mDC" = "DC",
  "pre-pDC" = "DC",
  "ASDC" = "DC",

  # Lymphoctyes
  "B cells" = "B cells",
  "Memory B" = "B cells",
  "transitional B" = "B cells",
  "pre B" = "B cells",
  "Naive B" = "B cells",
  "Plasma" = "B cells",
  "pro B" = "B cells",

  "CD4 Naive" = "CD4 T cells",
  "CD4 Memory" = "CD4 T cells",
  "CD8 Memory" = "CD8 T cells",
  "CD8 Effector_1" = "CD8 T cells",
  "CD8 Effector_3" = "CD8 T cells",
  "CD8 Effector_2" = "CD8 T cells",
  "T Proliferating" = "CD8 T cells",

  "NK" = "NK cells",
  "NK cells" =  "NK cells",
  "T/NKT cells" = "CD8 T cells",

  "Stromal" = "other",
  "LMPP" = "HSPC",
  "GMP" = "HSPC",
  "BaEoMa" = "HSPC",
  "Late Eryth" = "Erythroid",
  "Early Eryth" = "Erythroid",
  "Prog Mk" = "HSPC",
  "CLP" = "HSPC",
  "Platelet" = "Platelet",
  "Endo" = "Endothelial",

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

  "OPC Pdgfra" = "OPC"
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
  "Neutrophils" = "Myeloid cells",
  "OLG" = "OLG",
  "OPC" = "OPC",
  "Platelet" = "Erythroid",
  "CD8 T cells" = "Lymphoctyes",
  "CD4 T cells" = "Lymphoctyes",
  "NK cells" = "Lymphoctyes",
  "VLMC" = "VLMC",
  "other" = "other"
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
      "ggsci::deep_orange_material",
      "ggsci::yellow_material",
      "ggsci::red_material",
      "ggsci::pink_material",
      "ggsci::purple_material",
      "ggsci::light_blue_material",
      "ggsci::blue_material",
      "ggsci::cyan_material",
      "ggsci::teal_material",
      "ggsci::green_material",
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
  dplyr::mutate(
    ns = c(
      list(c(10, 9, 5)),
      list(c(3)),
      list(c(3)),
      list(c(3)),
      list(c(3)),
      list(c(10, 9, 7, 5, 4,3)),
      list(c(10, 7, 4, 2)),
      list(c(10, 7, 4)),
      list(c(10)),
      list(c(10)),
      list(c(10, 7, 5)),
      list(c(10, 9, 5)),
      list(c(10)),
      list(c(10)),
      list(c(6)),
      list(c(5)),
      list(c(5)),
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
  # dplyr::mutate(
  #   cell2_color = as.color
  # )
  dplyr::select(
    cell1, cell2, cell3, cell1_color, cell2_color, cell3_color
  ) ->
  recell_color_final

readr::write_rds(
  x = recell_color_final,
  file = "/home/liuc9/github/scbrain/data/azimuth/recell_color.rds"
)

recell_color_final |>
  dplyr::select(cell3, cell3_color) |>
  dplyr::mutate(cell3_color = prismatic::color(cell3_color)) |>
  tibble::deframe() |>
  prismatic::clr_saturate(shift = 0.8) |>
  # prismatic::clr_desaturate(shift = 0.4) |> |>
  prismatic::clr_lighten(shift = 0.3) |>
  # prismatic::clr_grayscale(method = "blue_channel") |>
  # prismatic::clr_negate() |>
  # prismatic::clr_tritan() |>
  plot()
  prismatic::color()
  alpha(0.8) |>
  scales::show_col()

# body --------------------------------------------------------------------


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
            celltypes_recell,
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

readr::write_rds(
  x = azimuth_ref_sunburst,
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst.rds"
)

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
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel12"
    )
  )


azimuth_ref_sunburst |>
  dplyr::select(region, case, cellratio) |>
  dplyr::group_by(region) |>
  tidyr::nest() |>
  dplyr::ungroup() ->
  azimuth_ref_sunburst_sel


azimuth_ref_sunburst_sel |>
  dplyr::mutate(
    ratiop = purrr::map2(
      .x = region,
      .y = data,
      .f = function(.x, .y) {
        # .x <- azimuth_ref_sunburst_sel$region[[3]]
        # .y <- azimuth_ref_sunburst_sel$data[[3]]


        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell1, ratio = cell1_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) ->
          .d

        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell2, ratio = cell2_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) ->
          .dd

        .y |>
          tidyr::unnest(cols = cellratio) |>
          dplyr::select(case, cluster = cell3, ratio = cell3_r) |>
          dplyr::distinct() |>
          dplyr::mutate(cluster = gsub(
            pattern = " [[0-9]]*.*%$",
            replacement = "",
            x = cluster
          )) ->
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
          dplyr::arrange(cell1) ->
          .cell1_color

        .d |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 1,
            color = 1,
            size = 0.05,
            alpha = 0.8
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            expand = c(0, 0)
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
            )
          ) ->
          .p;.p

        recell_color_final |>
          dplyr::select(cell2, cell2_color) |>
          dplyr::distinct() |>
          dplyr::filter(cell2 %in% .dd$cluster) |>
          dplyr::arrange(cell2) ->
          .cell2_color

        .dd |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 1,
            color = 1,
            size = 0.05
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            expand = c(0, 0)
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
            )
          ) ->
          .pp

        recell_color_final |>
          dplyr::select(cell3, cell3_color) |>
          dplyr::distinct() |>
          dplyr::arrange(cell3) |>
          dplyr::filter(cell3 %in% .ddd$cluster) ->
          .cell3_color

        .ddd |>
          ggplot(aes(
            x = case,
            y = ratio,
            fill = cluster
          )) +
          geom_col(
            width = 1,
            color = 1,
            size = 0.05
          ) +
          scale_x_discrete(
            limits = c("Sham", "MCAO", "UV"),
            expand = c(0, 0)
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
            )
          ) ->
          .ppp


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
          filename = glue::glue("Propertion_{.region}_cell1.pdf"),
          plot = .ratiop$p_cell1,
          device = "pdf",
          path = .outdir,
          width = 10,
          height = 8
        )

        ggsave(
          filename = glue::glue("Propertion_{.region}_cell2.pdf"),
          plot = .ratiop$p_cell2,
          device = "pdf",
          path = .outdir,
          width = 10,
          height = 8
        )
        ggsave(
          filename = glue::glue("Propertion_{.region}_cell3.pdf"),
          plot = .ratiop$p_cell3,
          device = "pdf",
          path = .outdir,
          width = 10,
          height = 8
        )

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel12"
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
  azimuth_ref_sunburst_sel_ratiop_forplot

recell_color_final |>
  dplyr::select(cell1, cell1_color) |>
  dplyr::distinct() |>
  dplyr::arrange(cell1) ->
  cell1_color

azimuth_ref_sunburst_sel_ratiop_forplot |>
  dplyr::group_by(region, case) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |>
          dplyr::ungroup() %>%
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
    width = 1,
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
      fill = NA
    ),
  ) ->
  p;p


ggsave(
  filename = glue::glue("Pie_plot_cell1.pdf"),
  plot = p,
  device = "pdf",
  path =  "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel12",
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
  azimuth_ref_sunburst_sel_ratiop_forplot

recell_color_final |>
  dplyr::select(cell2, cell2_color) |>
  dplyr::distinct() |>
  dplyr::arrange(cell2) ->
  cell2_color


azimuth_ref_sunburst_sel_ratiop_forplot |>
  dplyr::group_by(region, case) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |>
          dplyr::ungroup() %>%
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
    width = 1,
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
  p;p


ggsave(
  filename = glue::glue("Pie_plot_cell2.pdf"),
  plot = p,
  device = "pdf",
  path =  "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel12",
  width = 15,
  height = 15
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/04-manual-annotation.rda")
