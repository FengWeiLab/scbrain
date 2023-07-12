# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun May  7 21:39:54 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

recell_color_final <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/recell_color.rds"
)
azimuth_ref_sunburst_cell |>
  dplyr::select(project, region, case, cell_factor) ->
  azimuth_ref_sunburst_cell_cell_factor

azimuth_ref_sunburst_cell_cell_factor |>
  readr::write_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/recell_cell_factor.rds"
  )

azimuth_ref_sunburst_cell <-  readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
)


azimuth_ref_sunburst_cell |>
  dplyr::select(region, case, anno_new_new) |>
  dplyr::group_by(region) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    merge_object = purrr::map2(
      .x = region,
      .y = data,
      .f = function(.x, .y) {
        # .x <- .d$region[[1]]
        # .y <- .d$data[[1]]

        purrr::map(
          .y$anno_new_new,
          VariableFeatures
        ) |>
          purrr::reduce(.f = intersect) ->
          variable_features

        glue::glue("{.x}_{.y$case}")

        .m <- merge(
          .y$anno_new_new[[1]],
          y = c(
            .y$anno_new_new[[2]],
            .y$anno_new_new[[3]]
          ),
          add.cell.id = glue::glue("{.x}_{.y$case}"),
          project = .x
        )
        tibble::tibble(
          merge_object = list(.m),
          variable_features = list(variable_features)
        )
      }
    )
  ) |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = merge_object) ->
  azimuth_ref_sunburst_cell_merge




# body --------------------------------------------------------------------

azimuth_ref_sunburst_cell_merge |>
  dplyr::mutate(
    norm = purrr::pmap(
      .l = list(
        .r = region,
        .m = merge_object,
        .v = variable_features
      ),
      .f = function(.r, .m, .v) {
        # .r <- azimuth_ref_sunburst_cell_merge$region[[1]]
        # .m <- azimuth_ref_sunburst_cell_merge$merge_object[[1]]
        # .v <- azimuth_ref_sunburst_cell_merge$variable_features[[1]]

        Seurat::DefaultAssay(.m) <- "RNA"
        .m <- Seurat::NormalizeData(.m)
        .all.genes <- rownames(.m)
        .m <- Seurat::ScaleData(.m, features = .all.genes)
        .m <- Seurat::RunPCA(.m, features = .v)
        .m <- FindNeighbors(.m, dims = 1:10)
        .m <- RunUMAP(.m, dims = 1:10)
        .m <- RunTSNE(.m, dims = 1:10)
        .m
      }
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm


fn_plot_umap_tsne_new <- function(.x, .y) {

  azimuth_ref_sunburst_cell_cell_factor |>
    dplyr::filter(region == .x) ->
    .d
  .d$cell_factor[[1]] -> .cell_factor

  .cellemb <- .y@reductions[['tsne']]@cell.embeddings
  dplyr::bind_cols(
    .y@meta.data,
    .cellemb
  ) |>
    dplyr::mutate(celltype = cell3) |>
    dplyr::mutate(celltype_f = factor(celltype, levels = .cell_factor$cell3)) |>
    dplyr::mutate(cluster = as.factor(as.numeric(celltype_f))) |>
    dplyr::mutate(UMAP_1 = tSNE_1, UMAP_2 = tSNE_2) ->
    .xxx

  .xxx %>%
    dplyr::select(cluster, celltype) %>%
    dplyr::group_by(cluster, celltype) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}")) ->
    .xxx_celltype

  .xxx %>%
    dplyr::group_by(cluster) %>%
    tidyr::nest() %>%
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d %>%
      #   dplyr::filter(cluster == 18) %>%
      #   dplyr::pull(data) %>%
      #   .[[1]] ->
      #   .m

      .m %>%
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m %>%
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd

      .mmd %>%
        dplyr::group_by(u1, u2) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-n) ->
        .mmm

      .fc <- tryCatch(
        expr = {
          .mmm$n[[1]] / .mmm$n[[2]]
        }, # 1.1
        error = function(e) {
          1
        }
      )
      .mmm
      .fc

      if(.fc > 1.1) {
        .mmd %>%
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd %>%
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }

    })) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(cols = u) %>%
    dplyr::select(-data) %>%
    dplyr::left_join(.xxx_celltype, by = "cluster") %>%
    dplyr::arrange(cluster) ->
    .xxx_label

  recell_color_final |>
    dplyr::select(cell3, cell3_color) |>
    dplyr::distinct() |>
    dplyr::filter(cell3 %in% .xxx$cell3) |>
    dplyr::mutate(
      cell3 = factor(
        x = cell3,
        levels = .cell_factor$cell3
      )
    ) |>
    dplyr::arrange(cell3) ->
    .cell3_color

  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = cluster,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    )+
    # geom_text(
    #   data = .xxx_label,
    #   aes(
    #     label = cluster,
    #     x = UMAP_1,
    #     y = UMAP_2,
    #   ),
    #   size = 6
    # ) +
  scale_color_manual(
    name = "",
    values = .cell3_color$cell3_color,
    labels = .xxx_label$celltype,
    guide = guide_legend(
      ncol = 1,
      override.aes = list(size=2)
    )
  ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        size = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) ->
    p1

  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = case,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    scale_color_manual(
      name = "",
      limits = c("Sham", "MCAO", "UV"),
      labels = c("Sham", "tMCAO", "tMCAO+UVB"),
      guide = guide_legend(
        ncol = 1,
        override.aes = list(size=2)
      ),
      values = c("#1B1919FF", "#0099B4FF", "#FDAF91FF")
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        size = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) ->
    p2

  tibble::tibble(
    p_celltype = list(p1),
    p_case = list(p2)
  )
}

azimuth_ref_sunburst_cell_merge_norm |>
  dplyr::mutate(
    p = purrr::map2(
      .x = region,
      .y = norm,
      .f = fn_plot_umap_tsne_new
    )
  ) |>
  tidyr::unnest(cols = p) ->
  azimuth_ref_sunburst_cell_merge_norm_p

azimuth_ref_sunburst_cell_merge_norm_p |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .r = region,
        .p1 = p_celltype,
        .p2 = p_case
      ),
      .f = function(.r, .p1, .p2) {
        ggsave(
          filename = "{.r}-case-cluster-tsne.pdf" |> glue::glue(),
          plot = .p2,
          width = 10,
          height = 8,
          path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
        )
        ggsave(
          filename = "{.r}-celltype-cluster-tsne.pdf" |> glue::glue(),
          plot = .p1,
          width = 10,
          height = 8,
          path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
        )
      }
    )
  )

# tmp ---------------------------------------------------------------------






fn_plot_umap_tsne <- function(.x) {


  .x |>
    purrr::map(
      .f = function(.s) {

        dplyr::bind_cols(
          .s@meta.data,

          .s@reductions$ref.umap@cell.embeddings
        )
      }
    ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(celltype = cell3) |>
    dplyr::mutate(celltype_f = factor(celltype)) |>
    dplyr::mutate(cluster = as.factor(as.numeric(celltype_f))) ->
    .xxx

  .xxx %>%
    dplyr::select(cluster, celltype) %>%
    dplyr::group_by(cluster, celltype) %>%
    dplyr::count() %>%
    dplyr::arrange(-n) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}")) ->
    .xxx_celltype

  .xxx %>%
    dplyr::group_by(cluster) %>%
    tidyr::nest() %>%
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d %>%
      #   dplyr::filter(cluster == 18) %>%
      #   dplyr::pull(data) %>%
      #   .[[1]] ->
      #   .m

      .m %>%
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m %>%
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd

      .mmd %>%
        dplyr::group_by(u1, u2) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-n) ->
        .mmm

      .fc <- tryCatch(
        expr = {
          .mmm$n[[1]] / .mmm$n[[2]]
        }, # 1.1
        error = function(e) {
          1
        }
      )
      .mmm
      .fc

      if(.fc > 1.1) {
        .mmd %>%
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd %>%
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }

    })) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(cols = u) %>%
    dplyr::select(-data) %>%
    dplyr::left_join(.xxx_celltype, by = "cluster") %>%
    dplyr::arrange(cluster) ->
    .xxx_label


  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = cluster,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    geom_text(
      data = .xxx_label,
      aes(
        label = cluster,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 6
    ) +
    scale_colour_manual(
      name = NULL,
      values = pcc$color,
      labels = .xxx_label$celltype,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(size=4)
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        size = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) ->
    p1

  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = case,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    ggsci::scale_color_aaas(
      name = ""
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        size = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) ->
    p2

  p1 | p2
}
ggsave(
  filename = "Brain-cluster.pdf",
  plot = fn_plot_umap_tsne(
    azimuth_ref_sunburst_cell |>
      dplyr::filter(region == "Brain") |>
      dplyr::pull(anno_new_new)
  ),
  width = 16,
  height = 8,
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
)

ggsave(
  filename = "Skull-cluster.pdf",
  plot = fn_plot_umap_tsne(
      azimuth_ref_sunburst_cell |>
        dplyr::filter(region == "Skull") |>
        dplyr::pull(anno_new_new)
    )
  ,
  width = 16,
  height = 8,
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
)


ggsave(
  filename = "Meninge-cluster.pdf",
  plot = fn_plot_umap_tsne(
    azimuth_ref_sunburst_cell |>
      dplyr::filter(region == "Meninge") |>
      dplyr::pull(anno_new_new)
  ),
  width = 16,
  height = 8,
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
)









# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
