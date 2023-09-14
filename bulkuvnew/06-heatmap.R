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
library(simplifyEnrichment)

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


fn_gobp <- function(.genename, .color = "red") {
  .cc <- c("red" = "#AE1700", "green" = "#112a13")

  .go_bp <- clusterProfiler::enrichGO(
    gene = .genename,
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
  )

  if(is.null(.go_bp)) {
    return(
      tibble::tibble(
        gobp = list(NULL),
        goplot = list(NULL)
      )
    )
  }

  .go_bp %>%
    tibble::as_tibble()  %>%
    dplyr::mutate(
      Description = stringr::str_wrap(
        stringr::str_to_sentence(string = Description),
        width = 60
      )
    ) %>%
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count) %>%
    head(20) %>%
    dplyr::arrange(adjp, Count) %>%
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    )  ->
    .go_bp_for_plot

  .go_bp_for_plot %>%
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .cc[.color], color = NA, width = 0.7) +
    geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
    labs(y = "-log10(Adj. P value)") +
    scale_y_continuous(expand = c(0, 0.02)) +
    coord_flip() +
    theme(
      panel.background = element_rect(fill = NA),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 13, hjust = 1),
      axis.ticks.length.y = unit(3, units = "mm"),
      axis.text.x = element_text(color = "black", size = 12)
    ) ->
    .go_bp_plot

  tibble::tibble(
    gobp = list(.go_bp),
    goplot = list(.go_bp_plot)
  )
}

fn_heatmap_se <- function(.se) {

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
    row_title = NULL,
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    column_title = NULL,
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
    left_annotation = NULL,
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

  .heatmap


}

fn_heatmap_go <- function(.heatmap) {
  .row_order <- ComplexHeatmap::row_order(.heatmap)

  .row_order |>
    tibble::enframe(name = "cluster", value = "index") |>
    tidyr::unnest(index) |>
    dplyr::arrange(index) |>
    dplyr::mutate(genename = .heatmap@row_names_param$labels[index]) ->
    .row_order_index


  .row_order_index |>
    dplyr::select(cluster, genename) |>
    dplyr::group_by(cluster) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    tibble::deframe() |>
    purrr::map(.f = \(.x) {.x$genename}) |>
    purrr::map(.f = fn_gobp) |>
    tibble::enframe(name = "cluster", value = "gobp") |>
    tidyr::unnest(cols = gobp) ->
    .row_order_gobp

  .row_order_gobp |>
    dplyr::select(cluster, gobp) |>
    dplyr::mutate(goid = purrr::map(
      .x = gobp,
      .f = \(.x) {.x$ID}
    )) |>
    dplyr::select(cluster, goid) |>
    dplyr::arrange(cluster) |>
    dplyr::filter(purrr::map_lgl(.x = goid, .f = \(.x) {length(.x) != 0})) |>
    tibble::deframe() ->
    .goidlist


  .row_order_index |>
    dplyr::pull(cluster) ->
    .align_to_id_cluster

  .anno_go <-  anno_word_cloud_from_GO(
    align_to = .align_to_id_cluster,
    go_id = .goidlist,
    max_words = 30,
    bg_gp = gpar(fill = "#F0F0F0", col = "#AAAAAAFF")
  )

  .heatmap +  rowAnnotation(go = .anno_go) ->
    .heatmap_go

  tibble::tibble(
    row_order_index = list(.row_order_index),
    row_order_gobp = list(.row_order_gobp),
    heatmap_go = list(.heatmap_go)
  )


}


fn_heatmap_genelabel <- function(.se, .g = NULL) {
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
      .f = \(.x, .y, se = se) {
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

        # .g <- intersect(gs_immune_inflam, .x_deg_union)

        .se <- se[.x_deg_union, .y_barcode]

        # chm <- fn_heatmap_se(.se = .se, .g = .g)
        chm <- fn_heatmap_se(.se = .se)

      },
      se = se
    )
  ) ->
  se_group_de_deg_nest_heatmap



future::plan(future::multisession, workers = 4)
se_group_de_deg_nest_heatmap |>
  dplyr::mutate(
    heatmap_gobp = furrr::future_map(
      .x = heatmap,
      .f = \(.hm) {
        fn_heatmap_go(.hm)
      }
    )
  ) ->
  se_group_de_deg_nest_heatmap_heatmap_gobp
future::plan(future::sequential)

se_group_de_deg_nest_heatmap_heatmap_gobp |>
  tidyr::unnest(cols = heatmap_gobp) ->
  se_group_de_deg_nest_heatmap_heatmap_gobp_unnest




# save gopb plot ----------------------------------------------------------

se_group_de_deg_nest_heatmap_heatmap_gobp_unnest |>
  dplyr::select(seq, row_order_gobp) |>
  tidyr::unnest(cols = row_order_gobp) |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .x = seq,
        .y = cluster,
        .z = goplot,
        .zz = gobp
      ),
      .f = \(.x, .y, .z, .zz) {
        .xlsx_filename <- glue::glue("{.x}-{.y}-gobp.xlsx")
        .pdf_filename <- glue::glue("{.x}-{.y}-gobp.pdf")

        ggsave(
          plot = .z,
          filename = .pdf_filename,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
          width = 10,
          height = 6.5
        )

        writexl::write_xlsx(
          x = as.data.frame(.zz),
          path = file.path(
            "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
            .xlsx_filename
          )
        )
      }
    )
  )


# save heatmap ------------------------------------------------------------




se_group_de_deg_nest_heatmap_heatmap_gobp_unnest |>
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
            width = 8,
            height = 8
          )
          ComplexHeatmap::draw(object = .y)
          dev.off()
      }
    )
  )

se_group_de_deg_nest_heatmap_heatmap_gobp_unnest |>
  dplyr::mutate(
    a = purrr::map2(
      .x = seq,
      .y = heatmap_go,
      .f = \(.x, .y) {
        .filename <- glue::glue("{.x}-heatmap-gobp.pdf")
        pdf(
          file = file.path(
            "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
            .filename
          ),
          width = 8,
          height = 8
        )
        ComplexHeatmap::draw(object = .y)
        dev.off()
      }
    )
  )



# Heatmap gene label ------------------------------------------------------

# Brain -------------------------------------------------------------------


se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$row_order_gobp[[1]] |>
  dplyr::arrange(cluster) ->
  brain_cluster


brain_cluster$gobp[[1]] |> as.data.frame() |>
  dplyr::filter(grepl(
    pattern = "immune",
    x = Description,
    ignore.case = T
  )) |>
  dplyr::pull(geneID) |>
  strsplit("/") |>
  unlist() |>
  unique() ->
  brain_gene_1

brain_cluster$gobp[[2]] |> as.data.frame() |>
  dplyr::filter(grepl(
    pattern = "behavior|momory|cognition",
    x = Description,
    ignore.case = T
  )) |>
  dplyr::pull(geneID) |>
  strsplit("/") |>
  unlist() |>
  unique() ->
  brain_gene_2

se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$row_order_index[[1]] |>
  dplyr::slice(match(c(brain_gene_1, brain_gene_2) |> unique(), genename)) ->
  brain_genes

ComplexHeatmap::rowAnnotation(
  link = anno_mark(
    at = brain_genes$index,
    labels = brain_genes$genename,
    which = "row",
    side = "left",
    lines_gp = gpar(lwd = 0.5),
    labels_gp = gpar(fontsize = 10),
    padding = unit(0.5, "mm"),
    link_width = unit(5, "mm"),
  )
) +
  se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$heatmap_go[[1]] ->
  brain_heatmap_mark;brain_heatmap_mark
{
  pdf(
    file = file.path(
      "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
      "Brain-heatmap-mark.pdf"
    ),
    width = 10,
    height = 8
  )
  ComplexHeatmap::draw(object = brain_heatmap_mark)
  dev.off()
}




# Skull -----------------------------------------------------------------

se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$row_order_gobp[[2]] |>
  dplyr::arrange(cluster) ->
  skull_cluster

skull_cluster$gobp[[1]] |>
  as.data.frame() |>
  dplyr::filter(grepl(
    pattern = "Leukocyte",
    x = Description,
    ignore.case = T
  )) |>
  dplyr::pull(geneID) |>
  strsplit("/") |>
  unlist() |>
  unique() |>
  sort() ->
  skull_gene_1

skull_gene_1_1 <- skull_gene_1[stringr::str_detect(skull_gene_1, "Cxcl|Il|B2m|Gzm|Cd|Gata|Jun|Lck|Tnf")]

skull_cluster$gobp[[2]] |>
  as.data.frame() |>
  # dplyr::select(2) |> head(10)
  dplyr::filter(grepl(
    pattern = "immune",
    x = Description,
    ignore.case = T
  )) |>
  dplyr::pull(geneID) |>
  strsplit("/") |>
  unlist() |>
  unique() |>
  sort() ->
  skull_gene_2

skull_gene_2_2 <- skull_gene_2[stringr::str_detect(skull_gene_2, "Cxcl|Il|B2m|Gzm|Cd")]

skull_cluster$gobp[[3]] |>
  as.data.frame() |>
  # dplyr::select(2) |> head(10)
  dplyr::filter(grepl(
    pattern = "chemotaxis|would healing",
    x = Description,
    ignore.case = T
  )) |>
  dplyr::pull(geneID) |>
  strsplit("/") |>
  unlist() |>
  unique() ->
  skull_gene_3

se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$row_order_index[[2]] |>
  dplyr::slice(match(c(skull_gene_1_1, skull_gene_2_2, skull_gene_3) |> unique(), genename)) ->
  skull_genes

ComplexHeatmap::rowAnnotation(
  link = anno_mark(
    at = skull_genes$index,
    labels = skull_genes$genename,
    which = "row",
    side = "left",
    lines_gp = gpar(lwd = 0.5),
    labels_gp = gpar(fontsize = 10),
    padding = unit(0.5, "mm"),
    link_width = unit(5, "mm"),
  )
) +
  se_group_de_deg_nest_heatmap_heatmap_gobp_unnest$heatmap_go[[2]] ->
  skull_heatmap_mark;skull_heatmap_mark

{
  pdf(
    file = file.path(
      "/home/liuc9/github/scbrain/data/uvresultnew/05-heatmap",
      "Skull-heatmap-mark.pdf"
    ),
    width = 10,
    height = 8
  )
  ComplexHeatmap::draw(object = skull_heatmap_mark)
  dev.off()
}

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/uvrdanew/06-heatmap.rda")

load(file = "data/uvrdanew/06-heatmap.rda")
