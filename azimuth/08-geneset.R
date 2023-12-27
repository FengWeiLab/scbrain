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

combined_gene_set_list |>
  tibble::enframe() |>
  dplyr::filter()


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


# Inde --------------------------------------------------------------------

fn_tile <- function(a) {
  a |>
    dplyr::mutate(
      score = scale(score)
    ) |>
    dplyr::mutate(
      type = factor(type, c("Sham", "MCAO", "UV"))
    ) |>
    ggplot(aes(
      x = type,
      y = cell3,
      fill = score
    )) +
    geom_tile(
      linewidth = 0.2,
      colour = "grey",
    ) +
    geom_tile() +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    # geom_text(aes(label = n), size = 5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradient2(
      # breaks = round(seq(-0.4, 0.4,length.out = 5), digits = 2),
      # labels = format(seq(-0.4, 0.4,length.out = 5), digits = 2),
      low = "#33cbfe",
      mid = "#000000",
      high = "#fdfe00"
    ) +
    theme(
      text = element_text(family = 'Times'),
      title = element_text(family = 'Times'),

      # axis text
      axis.title = element_blank(),
      axis.text = element_text(size = 16, colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),

      # ticks
      axis.ticks = element_blank(),

      # legend
      legend.background = element_rect(colour = NA, fill = NA),
      legend.key = element_rect(colour = NA),
      legend.position = "right",

      #panel
      panel.background = element_blank(),
      # panel.spacing = unit(0, units = 'pt'),
      panel.grid = element_blank(),

      # plot
      plot.background = element_rect(fill = NA, colour = NA),
      # plot.margin = unit(c(0,0,0,0), units = 'pt'),
      plot.title = element_text(hjust = 0.5),

      validate = TRUE
    )
}



# Brain
azimuth_ref_ucell_ttest_sigs$sigs[[3]] |>
  dplyr::filter(cell3 %in% c("Astrocyte Aqp4_Gfap", "Astrocyte Aqp4_Slc7a10", "Microglia", "OLG", "Endothelial", "Pericytes" )) |>
  dplyr::mutate(
    cell3 = factor(cell3, levels = c("Astrocyte Aqp4_Gfap", "Astrocyte Aqp4_Slc7a10", "Microglia", "OLG", "Endothelial", "Pericytes" ) |> rev())
  ) ->
  brain_geneset


brain_geneset |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, MCAO_vs_Sham_pval, UV_vs_MCAO_pval
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    geneset
  ) |>
  tidyr::nest() |>
  dplyr::ungroup()->
  brain_geneset_nest

brain_geneset_nest |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_tile
    )
  ) ->
  brain_geneset_nest_p

brain_geneset_nest_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::map2(
      .x = p,
      .y = geneset,
      .f = \(.p, .geneset) {
        .filename <- glue::glue("{.geneset}.pdf")
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 5,
          height = 6,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/brain"
        )
      }
    )
  )


azimuth_ref_ucell_ttest_sigs$sigs[[1]] |>
  dplyr::filter(cell3 %in% c("Mature B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils")) |>
  dplyr::mutate(cell3 = factor(cell3, c("Mature B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils") |> rev())) ->
  meninge_geneset

meninge_geneset |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, MCAO_vs_Sham_pval, UV_vs_MCAO_pval
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    geneset
  ) |>
  tidyr::nest() |>
  dplyr::ungroup()->
  meninge_geneset_nest

meninge_geneset_nest |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_tile
    )
  ) ->
  meninge_geneset_nest_p


meninge_geneset_nest_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::map2(
      .x = p,
      .y = geneset,
      .f = \(.p, .geneset) {
        .filename <- glue::glue("{.geneset}.pdf")
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 5,
          height = 6,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/meninge"
        )
      }
    )
  )


azimuth_ref_ucell_ttest_sigs$sigs[[2]] |>
  dplyr::filter(cell3 %in% c("Mature B cells",  "Pre-B cells", "Pro-B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils")) |>
  dplyr::mutate(cell3 = factor(cell3, c("Mature B cells",  "Pre-B cells", "Pro-B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils") |> rev())) ->
  skull_geneset



skull_geneset |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, MCAO_vs_Sham_pval, UV_vs_MCAO_pval
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    geneset
  ) |>
  tidyr::nest() |>
  dplyr::ungroup()->
  skull_geneset_nest

skull_geneset_nest |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_tile
    )
  ) ->
  skull_geneset_p


skull_geneset_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::map2(
      .x = p,
      .y = geneset,
      .f = \(.p, .geneset) {
        .filename <- glue::glue("{.geneset}.pdf")
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 5,
          height = 6,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/skull"
        )
      }
    )
  )






# Load gene set -----------------------------------------------------------

fn_new_tile <- function(.a) {

  .a |>
    dplyr::mutate(
      geneset = gsub("_UCell", "", geneset)
    ) |>
    dplyr::mutate(
      geneset = gsub("\\.", " ", geneset)
    ) |>
    dplyr::mutate(
      geneset = stringr::str_to_sentence(
        string = geneset
      )
    ) |>
    dplyr::group_by(geneset) |>
    dplyr::mutate(score = scale(score)) |>
    dplyr::ungroup() ->
    for_plot


  for_plot |>
    dplyr::select(Terms) |>
    dplyr::distinct() |>
    dplyr::arrange(Terms) |>
    dplyr::mutate(
      color = ggsci::pal_aaas()(length(unique(for_plot$Terms)))
    ) ->
    parent_term_color

  for_plot |>
    dplyr::group_by(Terms, geneset) |>
    dplyr::summarise(s = sum(score)) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(Terms), s
    ) |>
    dplyr::left_join(
      parent_term_color,
      by = "Terms"
    ) ->
    for_plot_e_rank


  for_plot_e_rank |>
    dplyr::mutate(
      Terms = factor(
        Terms,
        levels = parent_term_color$Terms
      )
    ) |>
    dplyr::mutate(
      geneset = factor(
        geneset,
        levels = for_plot_e_rank$geneset
      )
    ) |>
    dplyr::arrange(Terms) |>
    tibble::rowid_to_column() |>
    dplyr::group_by(Terms) |>
    dplyr::mutate(
      n = mean(rowid)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      n = as.integer(floor(n))
    ) |>
    dplyr::mutate(
      label = ifelse(
        rowid == n,
        stringr::str_wrap(
          string = Terms,
          width = 10
        ),
        ""
      )
    ) |>
    ggplot(aes(
      x = 1,
      y = geneset,
      fill = Terms,
      color = Terms,
      label = label
    )) +
    geom_tile(
      width = 0.3
    ) +
    geom_text(
      color = "white",
      angle = 90,
      size = 6
    ) +
    scale_fill_manual(
      values = parent_term_color$color
    ) +
    scale_color_manual(
      values = parent_term_color$color
    ) +
    scale_y_discrete(
      expand = c(0,0),
      position = "left"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        color = for_plot_e_rank$color,
        size = 8
      ),
      axis.ticks.y = element_line(
        color = for_plot_e_rank$color,
      ),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) ->
    p_bar

  CPCOLS <- c("#191970", "#F8F8FF", "#FF4040")

  for_plot |>
    dplyr::mutate(
      geneset = factor(geneset, for_plot_e_rank$geneset),
      type = factor(type, c("Sham", "MCAO", "UV"))
    ) |>
    ggplot(aes(
      x = type,
      y = geneset,
    )) +
    geom_tile(aes(fill = score)) +
    scale_x_discrete(position = "top") +
    coord_fixed(ratio = 1) +
    # scale_fill_gradient2(
    #   high = "#B71C1CFF",
    #   mid = "#EEEEEEFF",
    #   low = "#0D47A1FF"
    # ) +
    scale_fill_gradient2(
      name = "Normalized score",
      low = CPCOLS[1],
      mid = CPCOLS[2],
      high = CPCOLS[3],
    ) +
    theme(
      # panel.grid = element_blank(),
      # panel.background = element_blank(),
      panel.background = element_rect(colour = NA, fill = NA),
      # panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid = element_blank(),
      # panel.grid.major = element_line(
      #   colour = "grey",
      #   linetype = "dashed",
      #   size = 0.2
      # ),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 0,
        vjust = 0.5,
        face = "bold",
        size = 12
      ),
      # axis.text.y = element_text(
      #   face = "bold",
      #   color = for_plot_e_rank$color,
      #   size = 12
      # ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      # complete = 1,
      legend.position = "right",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) ->
    p_tile

  cowplot::plot_grid(
    plotlist = list(p_bar, p_tile),
    align = 'h',
    rel_widths = c(.7, 1)
  ) ->
    p
  p

}

fn_linear <- function(.a) {

  terms_color <- c("Angiogenesis&BBB" = "#E41A1CFF", "Immune response" = "#FFFF33FF", "Inflammatory response" = "#984EA3FF", "Neurogenesis" = "#66C2A5FF")

  .a$Terms |> sort() |> unique() -> .terms

  tibble::tibble(
    Terms = .terms
  ) |>
    dplyr::mutate(
      termcolor = terms_color[Terms]
    ) |>
    dplyr::mutate(
      Terms = factor(Terms, levels = .terms)
    ) |>
    dplyr::arrange(Terms) ->
    .terms_color

  .a |>
    dplyr::mutate(
      geneset = gsub("_UCell", "", geneset)
    ) |>
    dplyr::mutate(
      geneset = gsub("\\.", " ", geneset)
    ) |>
    dplyr::mutate(
      geneset = stringr::str_to_sentence(
        string = geneset
      )
    ) |>
    dplyr::group_by(geneset) |>
    dplyr::mutate(score = scale(score)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Terms = factor(Terms, levels = .terms)
    ) ->
    for_plot

  for_plot |>
    dplyr::mutate(
      type = factor(type, c("Sham", "MCAO", "UV"))
    ) |>
    dplyr::mutate(group = as.numeric(type)) |>
    ggplot(aes(x = type, y = score, group = geneset, color = Terms)) +
    geom_point() +
    geom_line() +
    scale_color_manual(
      values = .terms_color$termcolor,
      guide = guide_legend(
        nrow = 2
      )
    ) +
    scale_x_discrete(
      labels = c("Sham", "tMCAO", "tMCAO+UVB"),
      expand = expansion(mult = 0, add = 0.3),
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title = element_text(
        color = "black",
        size = 16,
        face = "bold"
      ),
      axis.text = element_text(
        color = "black",
        size = 14,
        face = "bold"
      )
    ) +
    labs(
      y = "Gene signature score"
    ) ->
    p;p
}

geneset_filename <- "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/scuvresult/geneset_ttest-updata.xlsx"

brain_geneset <- readxl::read_xlsx(
  path = geneset_filename,
  sheet = "Brain"
)

brain_geneset |>
  dplyr::filter(cell3 %in% c("Astrocyte Aqp4_Gfap", "Astrocyte Aqp4_Slc7a10", "Microglia", "OLG", "Endothelial", "Pericytes", "VLMC" )) |>
  dplyr::mutate(
    cell3 = factor(cell3, levels = c("Astrocyte Aqp4_Gfap", "Astrocyte Aqp4_Slc7a10", "Microglia", "OLG", "Endothelial", "Pericytes", "VLMC" ) |> rev())
  ) |>
  dplyr::filter(!is.na(Terms)) ->
  brain_geneset_terms

brain_geneset_terms |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, Terms
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    cell3
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_new_tile
    )
  ) |>
  dplyr::mutate(
    plinear = purrr::map(
      .x = data,
      .f = fn_linear
    )
  ) ->
  brain_geneset_terms_nest_p

brain_geneset_terms_nest_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .p = p,
        .cell3  = cell3,
        .plinear = plinear
      ),
      .f = \(.p, .cell3, .plinear) {
        .filename <- glue::glue("{.cell3}.pdf")
        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/brain",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 14,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/brain"
        )

        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/brain",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .plinear,
          device = "pdf",
          width = 7,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/brain"
        )
      }
    )
  )


meninge_geneset <- readxl::read_xlsx(
  path = geneset_filename,
  sheet = "Meninge"
)

meninge_geneset |>
  dplyr::filter(cell3 %in% c("Mature B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils")) |>
  dplyr::mutate(cell3 = factor(cell3, c("Mature B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils") |> rev())) |>
  dplyr::filter(!is.na(Terms)) ->
  meninge_geneset_terms



meninge_geneset_terms |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, Terms
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    cell3
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_new_tile
    )
  ) |>
  dplyr::mutate(
    plinear = purrr::map(
      .x = data,
      .f = fn_linear
    )
  ) ->
  meninge_geneset_terms_nest_p


meninge_geneset_terms_nest_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .p = p,
        .cell3  = cell3,
        .plinear = plinear
      ),
      .f = \(.p, .cell3, .plinear)  {
        .filename <- glue::glue("{.cell3}.pdf")
        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/meninge",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 14,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/meninge"
        )


        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/meninge",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .plinear,
          device = "pdf",
          width = 7,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/meninge"
        )
      }
    )
  )


skull_geneset <- readxl::read_xlsx(
  path = geneset_filename,
  sheet = "Skull"
)

skull_geneset |>
  dplyr::filter(cell3 %in% c("Immature B cells", "Mature B cells",  "Pre-B cells", "Pro-B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils")) |>
  dplyr::mutate(cell3 = factor(cell3, c("Immature B cells", "Mature B cells",  "Pre-B cells", "Pro-B cells", "CD14 Monocytes", "Macrophage",  "mDC", "Neutrophils") |> rev()))|>
  dplyr::filter(!is.na(Terms)) ->
  skull_geneset_terms


skull_geneset_terms |>
  dplyr::select(
    cell3, geneset, Sham, MCAO, UV, Terms
  ) |>
  tidyr::pivot_longer(
    cols = c(Sham, MCAO, UV),
    names_to = "type",
    values_to = "score"
  ) |>
  dplyr::group_by(
    cell3
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p = purrr::map(
      .x = data,
      .f = fn_new_tile
    )
  ) |>
  dplyr::mutate(
    plinear = purrr::map(
      .x = data,
      .f = fn_linear
    )
  ) ->
  skull_geneset_terms_nest_p

skull_geneset_terms_nest_p |>
  # head(1)  |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .p = p,
        .cell3  = cell3,
        .plinear = plinear
      ),
      .f = \(.p, .cell3, .plinear)  {
        .filename <- glue::glue("{.cell3}.pdf")
        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/skull",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .p,
          device = "pdf",
          width = 14,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/12-geneset-tileplot/skull"
        )


        dir.create(
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/skull",
          recursive = T
        )
        ggsave(
          filename = .filename,
          plot = .plinear,
          device = "pdf",
          width = 7,
          height = 7,
          path = "/home/liuc9/github/scbrain/scuvresult/13-geneset-linearplot/skull"
        )
      }
    )
  )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/08-geneset.rda")
load(file = "data/azimuth/08-geneset.rda")
