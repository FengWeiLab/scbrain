#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Nov 21 16:29:34 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)
library(Seurat)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich.rds.gz"
  )

angiogenesis <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231121_163027_angiogenesis.txt"
)

angiogenesis |>
  dplyr::select(Symbol, term = `Annotated Term`) |>
  dplyr::group_by(term) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::filter(
    purrr::map_lgl(
      .x = data,
      .f = \(.x) {
        nrow(.x) > 3
      }
    )
  ) |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = \(.x) {
        .x$Symbol
      }
    )
  ) |>
  tibble::deframe() ->
  angiogenesis_nest

inflammatory <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231121_163530_inflammatory.txt"
)

inflammatory |>
  dplyr::select(Symbol, term = `Annotated Term`) |>
  dplyr::group_by(term) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::filter(
    purrr::map_lgl(
      .x = data,
      .f = \(.x) {
        nrow(.x) > 3
      }
    )
  ) |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = \(.x) {
        .x$Symbol
      }
    )
  ) |>
  tibble::deframe() ->
  inflammatory_nest

negative_immune <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231121_163818_negative_immune.txt"
)

negative_immune |>
  dplyr::select(Symbol, term = `Annotated Term`) |>
  dplyr::group_by(term) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::filter(
    purrr::map_lgl(
      .x = data,
      .f = \(.x) {
        nrow(.x) > 3
      }
    )
  ) |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = \(.x) {
        .x$Symbol
      }
    )
  ) |>
  tibble::deframe() ->
  negative_immune_nest

positive_immune <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231121_163914_positive_immune.txt"
)

positive_immune |>
  dplyr::select(Symbol, term = `Annotated Term`) |>
  dplyr::group_by(term) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::filter(
    purrr::map_lgl(
      .x = data,
      .f = \(.x) {
        nrow(.x) > 3
      }
    )
  ) |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = \(.x) {
        .x$Symbol
      }
    )
  ) |>
  tibble::deframe() ->
  positive_immune_nest

neuron_regeneration <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231121_164313_neuron_regeneration.txt"
)

neuron_regeneration |>
  dplyr::select(Symbol, term = `Annotated Term`) |>
  dplyr::group_by(term) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::filter(
    purrr::map_lgl(
      .x = data,
      .f = \(.x) {
        nrow(.x) > 3
      }
    )
  ) |>
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = \(.x) {
        .x$Symbol
      }
    )
  ) |>
  tibble::deframe() ->
  neuron_regeneration_nest

# body --------------------------------------------------------------------
the_case_color <- tibble::tibble(
  case = c("Sham", "MCAO", "UV"),
  label = c("Sham", "tMCAO", "tMCAO+UVB"),
  color = c("#1B1919FF", "#0099B4FF", "#FDAF91FF") |> prismatic::color()
)

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  dplyr::select(region, norm) |>
  dplyr::mutate(
    ucell = purrr::map(
      .x = norm,
      .f = \(.norm) {

        a1 <- UCell::AddModuleScore_UCell(
          .norm,
          features = neuron_regeneration_nest
        )

        a1@meta.data |>
          tibble::rownames_to_column(
            var = "barcode"
          ) |>
          dplyr::select(
            case, barcode,
            cell1, cell2, cell3,
            dplyr::contains("_color"),
            dplyr::contains("_UCell")
          ) |>
          tibble::as_tibble() |>
          tibble::add_column(
            geneset = "neuron_regeneration",
            .before = 1
          ) ->
          a1_sel

        a2 <- UCell::AddModuleScore_UCell(
          .norm,
          features = positive_immune_nest
        )

        a2@meta.data |>
          tibble::rownames_to_column(
            var = "barcode"
          ) |>
          dplyr::select(
            case, barcode,
            cell1, cell2, cell3,
            dplyr::contains("_color"),
            dplyr::contains("_UCell")
          ) |>
          tibble::as_tibble() |>
          tibble::add_column(
            geneset = "positive_immune",
            .before = 1
          ) ->
          a2_sel

        a3 <- UCell::AddModuleScore_UCell(
          .norm,
          features = negative_immune_nest
        )

        a3@meta.data |>
          tibble::rownames_to_column(
            var = "barcode"
          ) |>
          dplyr::select(
            case, barcode,
            cell1, cell2, cell3,
            dplyr::contains("_color"),
            dplyr::contains("_UCell")
          ) |>
          tibble::as_tibble() |>
          tibble::add_column(
            geneset = "negative_immune",
            .before = 1
          ) ->
          a3_sel

        a4 <- UCell::AddModuleScore_UCell(
          .norm,
          features = inflammatory_nest
        )

        a4@meta.data |>
          tibble::rownames_to_column(
            var = "barcode"
          ) |>
          dplyr::select(
            case, barcode,
            cell1, cell2, cell3,
            dplyr::contains("_color"),
            dplyr::contains("_UCell")
          ) |>
          tibble::as_tibble() |>
          tibble::add_column(
            geneset = "inflammatory",
            .before = 1
          ) ->
          a4_sel

        a5 <- UCell::AddModuleScore_UCell(
          .norm,
          features = angiogenesis_nest
        )

        a5@meta.data |>
          tibble::rownames_to_column(
            var = "barcode"
          ) |>
          dplyr::select(
            case, barcode,
            cell1, cell2, cell3,
            dplyr::contains("_color"),
            dplyr::contains("_UCell")
          ) |>
          tibble::as_tibble() |>
          tibble::add_column(
            geneset = "angiogenesis",
            .before = 1
          ) ->
          a5_sel

        dplyr::bind_rows(
          a1_sel, a2_sel, a3_sel,
          a4_sel, a5_sel
        )
      }
    )
  ) ->
  azimuth_ref_ucell

outdir <- "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/11-geneset"

azimuth_ref_ucell |>
  dplyr::mutate(
    p = purrr::map2(
      .x = region,
      .y = ucell,
      .f = \(.x, .y) {
        # .x <- azimuth_ref_ucell$region[[1]]
        # .y <- azimuth_ref_ucell$ucell[[1]]

        # .x

        outdir1 <- file.path(
          outdir, .x
        )
        dir.create(outdir1, showWarnings = F, recursive = F)

        the_case_color |>
          dplyr::left_join(
            .y,
            by = "case"
          ) |>
          dplyr::group_by(
            geneset
          ) |>
          tidyr::nest() |>
          dplyr::ungroup() ->
          .d

        .d |>
          dplyr::mutate(
            a = purrr::map2(
              .x = data,
              .y = geneset,
              .f = \(.m, .gs) {
                # .d$data[[1]] -> .m
                # .d$geneset[[1]] -> .gs

                .regnames <- colnames(
                  .m |>
                    dplyr::select(-c(
                      case, label, color, barcode, cell1, cell2, cell3,
                      cell1_color, cell2_color, cell3_color
                    ))
                )

                .m$cell1 |> as.character() |>  unique() -> .m_cell1
                .m$cell2 |> as.character() |>  unique() -> .m_cell2
                .m$cell3 |> as.character() |>  unique() -> .m_cell3


                .m_cell1 |>
                  purrr::map(
                    .f = \(.mc) {
                      # .mc <- .m_cell1[[1]]
                      outdir2 <- file.path(
                        outdir1, .mc
                      )
                      dir.create(outdir2, showWarnings = F, recursive = F)

                      outdir3 <- file.path(
                        outdir2, .gs
                      )

                      dir.create(outdir3, showWarnings = F, recursive = F)


                      .m |>
                        dplyr::filter(cell1 == .mc) ->
                        .mcm

                      .regnames |>
                        purrr::map(
                          .f = \(.regname) {

                            # .regname <- .regnames[[1]]

                            .ylabel <- gsub(
                              pattern = "\\.|_UCell",
                              replacement = " ",
                              x = .regname
                            ) |>
                              stringr::str_to_title()

                            ggpubr::ggviolin(
                              .mcm,
                              x = "label", y = .regname,
                              fill = "label",
                              palette = the_case_color$color,
                              add = "boxplot",
                              add.params = list(fill = "white"),
                              order = the_case_color$label,
                            ) +
                              ggpubr::stat_compare_means(
                                comparisons = list(
                                  c("tMCAO", "Sham"),
                                  c("tMCAO+UVB", "tMCAO")
                                ),
                              ) +
                              ggpubr::stat_compare_means(label.y = 0.4) +
                              theme(
                                legend.position = "none"
                              ) +
                              labs(
                                x = "",
                                y = .ylabel
                              ) ->
                              p

                            .filename <- "{.ylabel}.pdf" |> glue::glue()

                            ggsave(
                              filename = .filename,
                              plot = p,
                              path = outdir3,
                              width = 4,
                              height = 4.5
                            )
                          }
                        )
                    }
                  )
              }
            )
          )
      }
    )
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/08-geneset.rda")
load(file = "data/azimuth/08-geneset.rda")
