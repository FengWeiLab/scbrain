#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Sep  7 16:21:24 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)
library(DESeq2)
library(ComplexHeatmap)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------

colors <- tibble::tibble(
  group = unique(se$group) |> sort(),
  color = c(
    rep("grey",3),
    rep(
      c("#3B4992FF", "#EE0000FF"),
      3
    )
  )
)
# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_heatmap_se <- function(.se, .g = NULL) {
  # .se <- se_brain_deg


  .matrix <- assay(.se)
  .meta <- colData(.se)

  .matrix_scale <- .matrix %>%
    apply(1, scale) %>%
    t()

  colnames(.matrix_scale) <- .meta$barcode

  colors %>%
    dplyr::filter(group %in% .meta$group) %>%
    tibble::deframe() ->
    .cluster_col

  hma_top = ComplexHeatmap::HeatmapAnnotation(
    df = as.data.frame(.meta) %>%
      dplyr::select(Group = group),
    gap = unit(c(2,2), "mm"),
    col = list(Group = .cluster_col),
    which = "column"
  )

  if(!is.null(.g)){
    hma_index_right <- match(.g, rownames(.matrix_scale))
    hma_right <- ComplexHeatmap::rowAnnotation(
      link = anno_mark(
        at = hma_index_right,
        labels = rownames(.matrix_scale)[hma_index_right],
        which = "row",
        side = "left",
        lines_gp = gpar(
          lwd = 0.5
          # col = .g$color
        ),
        labels_gp = gpar(
          fontsize = 7
          # col = .g$color
        ),
        padding = unit(0.5, "mm"),
        link_width = unit(5, "mm"),
      )
    )
  } else {
    hma_right <- NULL
  }



  ComplexHeatmap::Heatmap(
    # data and color
    matrix = .matrix_scale,
    col =  circlize::colorRamp2(
      breaks = c(-1.1, 0, 1.1),
      colors = c("blue", "white", "red"),
      space = "RGB"
    ),
    name = "Normalized counts",
    na_col = 'grey',
    color_space = 'LAB',
    rect_gp = gpar(col = NA),
    border = NA,
    cell_fun = NULL,
    layer_fun = NULL,
    jitter = FALSE,

    # title
    # row_title = 'Selected genes', # OC44
    row_title = 'Gene name',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 10),
    column_title_rot = 0,

    # clustering of row
    cluster_rows = T,
    cluster_row_slices = T,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    row_dend_side = 'left',
    row_dend_width = unit(10, 'mm'),
    show_row_dend = F,
    row_dend_reorder = T,
    row_dend_gp = gpar(),

    # clustering of column
    cluster_columns = F,
    # cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    # column_dend_side = 'top',
    # column_dend_height = unit(10, 'mm'),
    # show_column_dend = F, # OC521
    # column_dend_gp = gpar(),
    # column_dend_reorder = F,

    row_order = NULL,
    # column_order = c(7,8,9,4,5,6,1,2,3),

    # row labels
    row_labels = rownames(.matrix_scale),
    row_names_side = 'right',
    show_row_names = F,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 3),
    row_names_rot = 0,
    row_names_centered = FALSE,

    # column labels
    # column_labels = colnames(.matrix_scale),
    column_names_side = 'bottom',
    show_column_names = T,
    column_names_max_height = unit(6, 'cm'),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_centered = FALSE,

    # annotation
    top_annotation = hma_top,
    bottom_annotation = NULL,
    left_annotation = hma_right,
    # right_annotation = hma_left,
    # left_annotation = NULL,
    right_annotation = NULL,


    # kmeans cluster number
    # row cluster is 1
    # column cluster is 2 with 10 repeats
    # km = 1,
    # split = NULL,
    row_km = 3,
    # row_km_repeats = 3,
    row_split = NULL,

    # column_km = 3
    # column_km_repeats = 5,
    column_split = rep(c("A", "B", "C"), each = 3),

    row_gap = unit(1, 'mm'),
    column_gap = unit(1, 'mm'),

    show_heatmap_legend = T,
    heatmap_legend_param = list(title = 'Z-score'),

    raster_quality = 2,
    raster_device_param = list(),

    raster_resize = F,

    post_fun = NULL
  ) ->
    .heatmap #; .heatmap

  .row_order <- ComplexHeatmap::row_order(.heatmap)

  .row_order |>
    purrr::map(
      .f = \(.rr) {
        rownames(.matrix_scale)[.rr]
      }
    ) ->
    row_names_cluster

  list(
    row_names_cluster = row_names_cluster,
    heatmap = .heatmap
  )

}

# load data ---------------------------------------------------------------

reactome_immune <- readr::read_lines(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/uvrdanew/REACTOME_IMMUNE_SYSTEM.v2023.1.Mm.gmt"
) |>
  strsplit(split = "\t")
reactome_immune <- reactome_immune[[1]][-c(1:3)]

reactome_inflam <- readr::read_lines(
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/uvrdanew/REACTOME_INFLAMMASOMES.v2023.1.Mm.gmt"
) |>
  strsplit(split = "\t")
reactome_inflam <- reactome_inflam[[1]][-c(1:3)]

gs_immune_inflam <- c(reactome_immune, reactome_inflam)


se_group_de <- readr::read_rds(
  file = "data/uvrdanew/se_group_de_volcano.rds.gz"
) %>%
  dplyr::select(bs, group, seq, se, vs, des_color)

se <- readr::read_rds(
  file = "data/uvrdanew/UV_count_matrix.se.rds.gz"
)

barcode_group <- readr::read_tsv(
  file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew",
    "barcode_group.tsv"
  )
)


# body --------------------------------------------------------------------


se_group_de |>
  dplyr::mutate(
    deg = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>%
          dplyr::filter(color != "grey") %>%
          dplyr::pull(GeneName)
      }
    ),
    deg_up = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>%
          dplyr::filter(color == "red") %>%
          dplyr::pull(GeneName)
      }
    ),
    deg_down = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>%
          dplyr::filter(color == "green") %>%
          dplyr::pull(GeneName)
      }
    )
  ) ->
  se_group_de_deg

se_group_de_deg |>
  dplyr::select(seq, deg) |>
  dplyr::group_by(seq) |>
  tidyr::nest() |>
  dplyr::ungroup() ->
  se_group_de_deg_nest

se_group_de_deg_nest |>
  dplyr::mutate(
    heatmap = purrr::map2(
      .x = data,
      .y = seq,
      .f = \(.x, .y, se = se, gs_immune_inflam = gs_immune_inflam) {
        # .x <- se_group_de_deg_nest$data[[1]]
        # .y <-se_group_de_deg_nest$seq[[1]]
        .x |>
          dplyr::pull(deg) |>
          purrr::reduce(.f = union) ->
          .x_deg_union

        barcode_group |>
          dplyr::filter(seq == .y) |>
          dplyr::pull(barcode) |>
          sort() ->
          .y_barcode

        .g <- intersect(gs_immune_inflam, .x_deg_union)

        .se <- se[.x_deg_union, .y_barcode]

        # chm <- fn_heatmap_se(.se = .se, .g = .g)
        chm <- fn_heatmap_se(.se = .se, .g = NULL)


      },
      se = se,
      gs_immune_inflam = gs_immune_inflam
    )
  ) |>
  dplyr::mutate(
    heatmap2 = purrr::map2(
      .x = data,
      .y = seq,
      .f = \(.x, .y, se = se, gs_immune_inflam = gs_immune_inflam) {
        # .x <- se_group_de_deg_nest$data[[1]]
        # .y <-se_group_de_deg_nest$seq[[1]]
        .x |>
          dplyr::pull(deg) |>
          purrr::reduce(.f = union) ->
          .x_deg_union

        barcode_group |>
          dplyr::filter(seq == .y) |>
          dplyr::pull(barcode) |>
          sort() ->
          .y_barcode

        .g <- intersect(gs_immune_inflam, .x_deg_union)

        .se <- se[.x_deg_union, .y_barcode]

        chm <- fn_heatmap_se(.se = .se, .g = .g)
        # chm <- fn_heatmap_se(.se = .se, .g = NULL)


      },
      se = se,
      gs_immune_inflam = gs_immune_inflam
    )
  ) ->
  se_group_de_deg_nest_heatmap


se_group_de_deg_nest_heatmap |>
  dplyr::mutate(
    a = purrr::map2(
      .x = seq,
      .y = heatmap,
      .f = \(.x, .y) {
        .filename <- glue::glue("{.x}-heatmap.pdf")
          pdf(
            file = file.path(
              "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
              .filename
            ),
            width = 5,
            height = 6
          )
          ComplexHeatmap::draw(
            object = .y$heatmap
          )
          dev.off()
      }
    )
  )

se_group_de_deg_nest_heatmap |>
  dplyr::mutate(
    a = purrr::map2(
      .x = seq,
      .y = heatmap2,
      .f = \(.x, .y) {
        .filename <- glue::glue("{.x}-heatmap-mark.pdf")
        pdf(
          file = file.path(
            "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
            .filename
          ),
          width = 5,
          height = 10
        )
        ComplexHeatmap::draw(
          object = .y$heatmap
        )
        dev.off()
      }
    )
  )


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/uvrdanew/06-heatmap.rda")
