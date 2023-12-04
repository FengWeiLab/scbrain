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

establishment_BBB <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231126_144257_establishement_of_blood_brain_barrier.txt"
)

establishment_BBB |>
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
  establishment_BBB_nest

maintenance_BBB <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231126_144602_maintenance_of_blood_brain_barrier.txt"
)


maintenance_BBB |>
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
  maintenance_BBB_nest




reg_nervous_system <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231126_144957_regulation_of_nervous_system_development.txt"
)


reg_nervous_system |>
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
  reg_nervous_system_nest



tissue_regeneration <- readr::read_tsv(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/GO_term_summary_20231126_145119_tissue_regeneration.txt"
)


tissue_regeneration |>
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
  tissue_regeneration_nest

combined_gene_set_list <- c(
  angiogenesis_nest,
  inflammatory_nest,
  negative_immune_nest,
  positive_immune_nest,
  neuron_regeneration_nest,
  establishment_BBB_nest,
  maintenance_BBB_nest,
  reg_nervous_system_nest,
  tissue_regeneration_nest
)

the_case_color <- tibble::tibble(
  case = c("Sham", "MCAO", "UV"),
  label = c("Sham", "tMCAO", "tMCAO+UVB"),
  color = c("#1B1919FF", "#0099B4FF", "#FDAF91FF") |> prismatic::color()
)
# body --------------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  dplyr::select(region, norm) |>
  dplyr::mutate(
    ucell = purrr::map(
      .x = norm,
      .f = \(.norm) {
        # .norm <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$norm[[1]]

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

        list(
          neuron_regeneration = a1_sel,
          positive_immune = a2_sel,
          negative_immune = a3_sel,
          inflammatory = a4_sel,
          angiogenesis = a5_sel
        ) |>
          tibble::enframe(
            name = "gs",
            value = "data"
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
        # .x <- azimuth_ref_ucell$region[[3]]
        # .y <- azimuth_ref_ucell$ucell[[3]]

        # .x

        outdir1 <- file.path(
          outdir, .x
        )
        dir.create(outdir1, showWarnings = F, recursive = F)

        # the_case_color |>
        #   dplyr::left_join(
        #     .y,
        #     by = "case"
        #   ) |>
        #   dplyr::group_by(
        #     geneset
        #   ) |>
        #   tidyr::nest() |>
        #   dplyr::ungroup() ->
        #   .d

        .y |>
          dplyr::mutate(
            data = purrr::map(
              .x = data,
              .f = \(.xx) {
                the_case_color |>
                  dplyr::left_join(
                    .xx,
                    by = "case"
                  )
              }
            )
          ) ->
          .d

        .d |>
          dplyr::mutate(
            a = purrr::map2(
              .x = data,
              .y = gs,
              .f = \(.m, .gs) {
                # .d$data[[1]] -> .m
                # .d$gs[[1]] -> .gs

                .regnames <- colnames(
                  .m |>
                    dplyr::select(-c(
                      case, label, color, geneset, barcode, cell1, cell2, cell3,
                      cell1_color, cell2_color, cell3_color
                    ))
                )

                .m$cell1 |> as.character() |>  unique() -> .m_cell1
                .m$cell2 |> as.character() |>  unique() |> gsub("/", "_", x = _) -> .m_cell2
                .m$cell3 |> as.character() |>  unique() |> gsub("/", "_", x = _) -> .m_cell3


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

                .m_cell2 |>
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
                        dplyr::filter(cell2 == .mc) ->
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

                .m_cell3 |>
                  purrr::map(
                    .f = \(.mc) {
                      # .mc <- .m_cell3[[1]]
                      outdir2 <- file.path(
                        outdir1, .mc
                      )
                      dir.create(outdir2, showWarnings = F, recursive = F)

                      outdir3 <- file.path(
                        outdir2, .gs
                      )

                      dir.create(outdir3, showWarnings = F, recursive = F)


                      .m |>
                        dplyr::filter(cell3 == .mc) ->
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

# Update ------------------------------------------------------------------

future::plan(future::multisession, workers = 3)
azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  dplyr::select(region, norm) |>
  dplyr::mutate(
    ucell = furrr::future_map(
      .x = norm,
      .f = \(.norm) {
        # .norm <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$norm[[1]]

        UCell::AddModuleScore_UCell(
          .norm,
          features = combined_gene_set_list
        )

      }
    )
  ) ->
  azimuth_ref_ucell
future::plan(future::sequential)

geneset_name <- names(combined_gene_set_list)

azimuth_ref_ucell |>
  dplyr::mutate(
    selected_metadata = purrr::map(
      .x = ucell,
      .f = \(.x) {
        # .x <- azimuth_ref_ucell$ucell[[3]]
        .x@meta.data |>
          dplyr::select(
            case, cell3,
            dplyr::contains("_UCell")
          ) |>
          dplyr::select(-dplyr::contains("signature_")) |>
          dplyr::mutate(case = factor(case, levels =  the_case_color$case)) ->
          .xx

        .xx |>
          tidyr::gather(
            -case, -cell3,
            key = geneset,
            value = score
          ) |>
          dplyr::group_by(
            cell3, geneset
          ) |>
          tidyr::nest() |>
          dplyr::ungroup() ->
          .d

        future::plan(future::multisession, workers = 10)
        .d |>
          dplyr::mutate(
            t = purrr::map(
              .x = data,
              .f = \(.m) {
                vs1 <- c("Sham", "MCAO")
                vs2 <- c("MCAO", "UV")

                .m |>
                  dplyr::filter(case %in% vs1) ->
                  .m_vs1

                .m |>
                  dplyr::filter(case %in% vs2) ->
                  .m_vs2

                .mm <- tryCatch(
                  expr = {
                    t.test(score ~ case, .m_vs1) |>
                      broom::tidy() |>
                      dplyr::select(Sham = estimate1, MCAO = estimate2, MCAO_vs_Sham_pval = p.value) ->
                      t_vs1

                    t.test(score ~ case, .m_vs2) |>
                      broom::tidy() |>
                      dplyr::select(MCAO = estimate1, UV = estimate2, UV_vs_MCAO_pval = p.value) ->
                      t_vs2

                    t_vs1 |>
                      dplyr::bind_cols(t_vs2 |> dplyr::select(-MCAO)) |>
                      dplyr::relocate(UV, .before = 3) |>
                      dplyr::mutate(MCAO_Sham = MCAO - Sham, UV_MCAO = UV - MCAO, UV_Sham = UV - Sham)
                  },
                  error = \(e) {
                    tibble::tibble(
                      Sham = NULL,
                      MCAO = NULL,
                      UV = NULL,
                      MOCA_vs_Sham_pval = NULL,
                      UV_vs_MCAO_pval = NULL,
                      MCAO_Sham = NULL,
                      UV_MCAO = NULL,
                      UV_Sham = NULL
                    )
                  }
                )
                .mm
              }
            )
          ) |>
          tidyr::unnest(cols = t) ->
          .dd
        future::plan(future::sequential)

        .dd
      }
    )
  ) ->
  azimuth_ref_ucell_ttest


azimuth_ref_ucell_ttest |>
  dplyr::mutate(
    sigs = purrr::map(
      .x = selected_metadata,
      .f = \(.x) {
        .x |>
          dplyr::select(-data) |>
          dplyr::filter(MCAO_vs_Sham_pval < 0.05) |>
          dplyr::filter(UV_vs_MCAO_pval < 0.05)
      }
    )
  ) |>
  dplyr::select(region, sigs) ->
  azimuth_ref_ucell_ttest_sigs


azimuth_ref_ucell_ttest_sigs |>
  tibble::deframe() |>
  writexl::write_xlsx(
    "/home/liuc9/github/scbrain/scuvresult/11-geneset/geneset_ttest.xlsx"
  )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/08-geneset.rda")
load(file = "data/azimuth/08-geneset.rda")
