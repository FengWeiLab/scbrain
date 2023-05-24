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

# body --------------------------------------------------------------------


azimuth_ref |>
  dplyr::mutate(
    sunburst = purrr::pmap(
      .l = list(
        .region = region,
        .anno = anno_new
      ),
      .f = function(.region, .anno) {
        # .anno <- azimuth_ref$anno_new[[2]]
        # .region <- azimuth_ref$region[[2]]

        .d <- if(.region == "Brain"){
          .anno@meta.data
        } else {
          .anno@meta.data |>
            dplyr::mutate(
              celltype = dplyr::case_when(
                celltype.l1 %in% c("DC", "Mono", "cDC1", "cDC2", "cl. Monocytes", "MC", "MdC", "migDC", "Neutrophils", "Non-cl. monocytes", "pDC") ~ "Myeloid",
                celltype.l1 %in% c("B", "CD4 T", "CD 8T", "NK", "other T", "B cells", "NK cells", "T/NKT cells") ~ "Lymphoid",
                celltype.l1 %in% c("HSPC") ~ "HSPC",
                TRUE ~ "other"
              )
            )
        }


        .dd <- if(.region == "Brain") {
          .d |>
            dplyr::count(class, subclass, cluster) |>
            dplyr::mutate(cluster_r = n / sum(n)) |>
            dplyr::group_by(subclass) |>
            dplyr::mutate(subclass_r = sum(cluster_r)) |>
            dplyr::ungroup() |>
            dplyr::group_by(class) |>
            dplyr::mutate(class_r = sum(cluster_r)) |>
            dplyr::ungroup() |>
            dplyr::mutate(
              class = glue::glue("{class} {round(class_r * 100, 2)}%"),
              subclass = glue::glue("{subclass} {round(subclass_r * 100, 2)}%"),
              cluster = glue::glue("{cluster} {round(cluster_r * 100, 2)}%")
            )
        } else {
          .d |>
            dplyr::count(celltype, celltype.l1, celltype.l2) |>
            dplyr::mutate(celltype.l2_r = n / sum(n)) |>
            dplyr::group_by(celltype.l1) |>
            dplyr::mutate(celltype.l1_r = sum(celltype.l2_r)) |>
            dplyr::ungroup() |>
            dplyr::group_by(celltype) |>
            dplyr::mutate(celltype_r = sum(celltype.l2_r)) |>
            dplyr::ungroup() |>
            dplyr::mutate(
              celltype =  glue::glue("{celltype} {round(celltype_r * 100, 2)}%"),
              celltype.l1 = glue::glue("{celltype.l1} {round(celltype.l1_r * 100, 2)}%"),
              celltype.l2 = glue::glue("{celltype.l2} {round(celltype.l2_r * 100, 2)}%"),
            )
        }
        .dd |>
          dplyr::select(1, 2,3,4) |>
          plotme::count_to_sunburst() ->
          .p

        .new_anno <- .anno
        .new_anno@meta.data <- .d

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
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel6"
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
        # .x <- azimuth_ref_sunburst_sel$region[[1]]
        # .y <- azimuth_ref_sunburst_sel$data[[1]]

        .d <- if (.x == "Brain") {
          .y |>
            tidyr::unnest(cols = cellratio) |>
            dplyr::select(case, cluster = subclass, ratio = subclass_r ) |>
            dplyr::distinct() |>
            dplyr::mutate(cluster = gsub(
              pattern = " .*$",
              replacement = "",
              x = cluster
            ))
        } else {
          .y |>
            tidyr::unnest(cols = cellratio) |>
            dplyr::select(case, cluster = celltype, ratio = celltype_r) |>
            dplyr::distinct() |>
            dplyr::mutate(cluster = gsub(
              pattern = " .*$",
              replacement = "",
              x = cluster
            ))
        }


        .d |>
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
          # scale_fill_manual(
          #   name = "Cell type",
          #   limits = p_d_fill$cluster,
          #   labels = p_d_fill$celltype,
          #   values = pcc$color
          # )
          ggsci::scale_fill_npg(
            name = "Cell type"
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
          .p
        .p
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
          filename = glue::glue("Propertion_{.region}.pdf"),
          plot = .ratiop,
          device = "pdf",
          path = .outdir,
          width = 10,
          height = 8
        )

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel6"
    )
  )

azimuth_ref_sunburst_sel_ratiop$ratiop |>
  wrap_plots()



azimuth_ref_sunburst_sel |>
  dplyr::mutate(
    ratiop = purrr::map2(
      .x = region,
      .y = data,
      .f = function(.x, .y) {
        # .x <- azimuth_ref_sunburst_sel$region[[2]]
        # .y <- azimuth_ref_sunburst_sel$data[[2]]

        .d <- if (.x == "Brain") {
          .y |>
            tidyr::unnest(cols = cellratio) |>
            dplyr::select(case, cluster = subclass, ratio = subclass_r ) |>
            dplyr::distinct() |>
            dplyr::mutate(cluster = gsub(
              pattern = " .*$",
              replacement = "",
              x = cluster
            ))
        } else {
          .y |>
            tidyr::unnest(cols = cellratio) |>
            dplyr::select(case, cluster = celltype.l1, ratio = celltype.l1_r) |>
            dplyr::distinct() |>
            dplyr::mutate(cluster = gsub(
              pattern = " [[0-9]]⁠.*$",
              replacement = "",
              x = cluster
            ))
        }


        .d |>
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
          # scale_fill_manual(
          #   name = "Cell type",
          #   limits = p_d_fill$cluster,
          #   labels = p_d_fill$celltype,
          #   values = pcc$color
          # )
          ggsci::scale_fill_npg(
            name = "Cell type"
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
          .p
        .p
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
          filename = glue::glue("Propertion_{.region}.pdf"),
          plot = .ratiop,
          device = "pdf",
          path = .outdir,
          width = 10,
          height = 8
        )

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel6"
    )
  )
# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
