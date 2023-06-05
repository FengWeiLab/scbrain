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

recell1 <- c("NK cells" =  "NK cells", "D-BAM" = "Macrophage", "B cells" = "B cells", "T/NKT cells" = "T/NKT cells", "Neutrophils" = "Neutrophils", "cl. Monocytes" = "Monocytes", "MC" = "Monocytes", "Non-cl. monocytes" = "Monocytes", "pDC" = "DC", "cDC2" = "DC", "cDC1" = "DC", "MdC" = "Monocytes", "migDC" = "DC", "Memory B" = "B cells", "transitional B" = "B cells", "CD4 Naive" = "T cells", "CD14 Mono" = "Monocytes", "pre B" = "B cells", "CD16 Mono" = "Monocytes", "pre-mDC" = "DC", "CD8 Memory" = "T cells", "Stromal" = "other", "Late Eryth" = "Erythroid", "Naive B" = "B cells", "CD8 Effector_1" = "T cells", "Plasma" = "B cells", "CD4 Memory" = "T cells", "GMP" = "HSPC", "LMPP" = "HSPC", "BaEoMa" = "HSPC", "NK" = "NK cells", "Early Eryth" = "Erythroid", "pro B" = "B cells", "Prog Mk" = "HSPC", "Macrophage" = "Macrophage", "pre-pDC" = "DC", "ASDC" = "DC", "CLP" = "HSPC", "T Proliferating" = "T cells", "CD8 Effector_3" = "T cells", "Platelet" = "Platelet", "CD8 Effector_2" = "T cells", "Endo" = "Endothelial", "Micro" = "Microglia", "Astro Aqp4_Slc7a10" = "Astrocyte", "VLMC_6" = "VLMC", "Oligo Opalin_2" = "OLG", "Astro Aqp4_Gfap" = "Astrocyte", "PVM_1" = "Microglia", "VLMC_2" = "VLMC", "VLMC_3" = "VLMC", "Oligo Opalin_1" = "OLG", "Oligo Enpp6_2" = "OLG", "Oligo Opalin_4" = "OLG", "OPC Pdgfra" = "OPC")

recell2 <- c(
  "Astrocyte" = "Astrocyte",
  "B cells" = "Lymphoctyes",
  "DC" = "Myeloid cells",
  "Endothelial" = "Endothelial",
  "Erythroid" = "Erythroid",
  "HSPC" = "HSPC",
  "Macrophage" = "Myeloid cells",
  "Microglia" = "Microglia",
  "Monocytes" = "Myeloid cells",
  "NK cells" = "Lymphoctyes",
  "Neutrophils" = "Myeloid cells",
  "OLG" = "OLG",
  "OPC" = "OPC",
  "Platelet" = "Erythroid",
  "T cells" = "Lymphoctyes",
  "T/NKT cells" = "Lymphoctyes",
  "VLMC" = "VLMC",
  "other" = "other"
)


celltypes |>
  dplyr::mutate(
    cell1 = plyr::revalue(
      x = celltype.l2,
      replace = recell1
    )
  ) |>
  dplyr::mutate(
    cell2 = plyr::revalue(
      x = cell1,
      replace = recell2
    )
  ) |>
  dplyr::select(-celltype.l1) |>
  dplyr::distinct() ->
  celltypes_recell



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
          dplyr::count(cell2, cell1, celltype.l2) |>
          dplyr::mutate(celltype.l2_r = n / sum(n)) |>
          dplyr::group_by(cell1) |>
          dplyr::mutate(cell1_r = sum(celltype.l2_r)) |>
          dplyr::ungroup() |>
          dplyr::group_by(cell2) |>
          dplyr::mutate(cell2_r = sum(celltype.l2_r)) |>
          dplyr::ungroup() |>
          dplyr::mutate(
            cell2 =  glue::glue("{cell2} {round(cell2_r * 100, 2)}%"),
            cell1 = glue::glue("{cell1} {round(cell1_r * 100, 2)}%"),
            celltype.l2 = glue::glue("{celltype.l2} {round(celltype.l2_r * 100, 2)}%"),
          ) ->
          .dd

        .dd |>
          dplyr::select(1, 2,3,4) |>
          plotme::count_to_sunburst() ->
          .p

        .new_anno <- .anno
        .new_anno@meta.data$cell1 <- .d$cell1
        .new_anno@meta.data$cell2 <- .d$cell2

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
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel11"
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

        .scale_fill <- if (.x == "Brain") {
          ggthemes::scale_fill_tableau(
            palette = "Tableau 20",
            name = "Cell types",
            direction = 1
          )
        } else if(.x == "Meninge") {
          ggsci::scale_fill_npg(
            name = "Cell types"
          )
        } else {
          ggsci::scale_fill_npg(
            name = "Cell types"
          )
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
          .scale_fill +
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
          .p
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
          ggsci::scale_fill_npg(
            name = "Cell types"
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

        list(
          p_cell1 = .p,
          p_cell2 = .pp,
          d_cell1 = .d,
          d_cell2 = .dd
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

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel11"
    )
  )


# Merge -------------------------------------------------------------------


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
  ) +
  coord_polar(
    theta = "y",
    start = 0
  ) +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      label = glue::glue("{scales::percent(percentage)}"),
    ),
    size = 6,
    nudge_x = 0.2,
    show.legend = FALSE,
  ) +
  ggsci::scale_fill_simpsons() +
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
  p


ggsave(
  filename = glue::glue("Pie_plot_cell2.pdf"),
  plot = p,
  device = "pdf",
  path =  "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel11",
  width = 15,
  height = 15
)



# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
