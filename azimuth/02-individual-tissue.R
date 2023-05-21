# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar 29 15:02:36 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)

# src ---------------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(dplyr::desc(cancer_types))

# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


fn_azimuth <- function(.sc, .ref) {
  # .sc <- d$sc[[5]]
  # .ref <- d$refs[[5]]

  .sct <- RunAzimuth(
    query = .sc,
    reference = .ref
  )

  # p1 <- DimPlot(.sct, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3)
  # p1

  .sct
}

fn_plot_umap_tsne <- function(.x, .celltype = "", .reduction = "ref.umap", .facet = FALSE) {
  # .x = project_sc_azimuth$anno[[3]]
  # .celltype = "predicted.celltype.l2"
  # .reduction = "ref.umap"


  .x@meta.data[[.celltype]] |>
    unique() |>
    as.factor() ->
    .celltype_u
  .replace <- as.numeric(.celltype_u) - 1
  names(.replace) <- .celltype_u


  .umap <- as.data.frame(.x@reductions[[.reduction]]@cell.embeddings)
  colnames(.umap) <- c("UMAP_1", "UMAP_2")

  .xx <- .x@meta.data[, c("case", "region", .celltype)] |>
    dplyr::rename(celltype = .celltype) |>
    dplyr::mutate(cluster = plyr::revalue(
      x = celltype,
      replace = .replace
    )) |>
    dplyr::mutate(cluster = as.numeric(cluster)) |>
    dplyr::mutate(cluster = factor(cluster))

  .xxx <- dplyr::bind_cols(.umap, .xx)

  .xxx |>
    dplyr::select(cluster, celltype) |>
    dplyr::group_by(cluster, celltype) |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup() |>
    dplyr::group_by(cluster) |>
    dplyr::top_n(1) |>
    dplyr::ungroup() |>
    dplyr::select(-n) |>
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}")) ->
    .xxx_celltype

  .xxx |>
    dplyr::group_by(cluster) |>
    tidyr::nest() |>
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d |>
      #   dplyr::filter(cluster == 14) |>
      #   dplyr::pull(data) |>
      #   .[[1]] ->
      #   .m

      .m |>
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m |>
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd

      .mmd |>
        dplyr::group_by(u1, u2) |>
        dplyr::count() |>
        dplyr::ungroup() |>
        dplyr::arrange(-n) ->
        .mmm

      if (nrow(.mmm) == 1) {
        return(
          .mmd |>
            # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
            dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        )
      }

      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1

      if (.fc > 1.1) {
        .mmd |>
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd |>
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }
    })) |>
    dplyr::ungroup() |>
    tidyr::unnest(cols = u) |>
    dplyr::select(-data) |>
    dplyr::left_join(.xxx_celltype, by = "cluster") |>
    dplyr::arrange(cluster) ->
    .xxx_label

  .labs <- switch(
    EXPR = .reduction,
    "umap" = {
      labs(
        x = "UMAP1",
        y = "UMAP2"
      )
    },
    "ref.umap" = {
      labs(
        x = "UMAP1",
        y = "UMAP2"
      )
    },
    "tsne" = {
      labs(
        x = "tSNE1",
        y = "tSNE2"
      )
    }
  )

  .split <- if (.facet) {
    facet_wrap(~case, nrow = 1)
  } else {
    theme()
  }


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
        override.aes = list(size = 4)
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5,
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
    ) +
    .split +
    .labs
}

fn_plot_refumap <- function(.x, .celltype = "celltype", .reduction = "ref.umap", .facet = FALSE) {
  # .x = project_sc_azimuth_refumap_unique_celltype_union_anno_cell$anno_cell[[1]]
  # .celltype = "celltype"
  # .reduction = "ref.umap"
  # .facet = FALSE


  # .x@

  # .x@meta.data[[.celltype]] |>
  #   unique() |>
  #   as.factor() ->
  #   .celltype_u
  # .replace <- as.numeric(.celltype_u) - 1
  # names(.replace) <- .celltype_u


  .umap <- as.data.frame(.x@reductions[[.reduction]]@cell.embeddings)
  colnames(.umap) <- c("UMAP_1", "UMAP_2")

  .x@meta.data |>
    dplyr::select(case, region, tidyselect::all_of(.celltype), cluster, cluster_celltype) |>
    dplyr::rename(celltype = .celltype) ->
    .xx

  .xxx <- dplyr::bind_cols(.umap, .xx)

  .xxx |>
    dplyr::select(cluster, celltype = cluster_celltype) |>
    dplyr::group_by(cluster, celltype) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::mutate(ratio = n / sum(n)) |>
    dplyr::mutate(celltype_ratio = glue::glue("{celltype}, {round(ratio * 100, 2)}%")) ->
    .xxx_celltype

  .xxx |>
    dplyr::group_by(cluster) |>
    tidyr::nest() |>
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d |>
      #   dplyr::filter(cluster == 14) |>
      #   dplyr::pull(data) |>
      #   .[[1]] ->
      #   .m

      .m |>
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m |>
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd

      .mmd |>
        dplyr::group_by(u1, u2) |>
        dplyr::count() |>
        dplyr::ungroup() |>
        dplyr::arrange(-n) ->
        .mmm

      if (nrow(.mmm) == 1) {
        return(
          .mmd |>
            # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
            dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        )
      }

      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1

      if (.fc > 1.1) {
        .mmd |>
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd |>
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) |>
          dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }
    })) |>
    dplyr::ungroup() |>
    tidyr::unnest(cols = u) |>
    dplyr::select(-data) |>
    dplyr::left_join(.xxx_celltype, by = "cluster") |>
    dplyr::arrange(cluster) ->
    .xxx_label

  .labs <- switch(
    EXPR = .reduction,
    "umap" = {
      labs(
        x = "UMAP1",
        y = "UMAP2"
      )
    },
    "ref.umap" = {
      labs(
        x = "UMAP1",
        y = "UMAP2"
      )
    },
    "tsne" = {
      labs(
        x = "tSNE1",
        y = "tSNE2"
      )
    }
  )

  .split <- if (.facet) {
    facet_wrap(~case, nrow = 1)
  } else {
    theme()
  }


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
      labels = .xxx_label$celltype_ratio,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(size = 4)
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5,
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
    ) +
    .split +
    .labs ->
    .p

  tibble::tibble(
    new_p = list(.p),
    cellnumber = list(.xxx_label)
  )
}


# load data ---------------------------------------------------------------

project_sc <- readr::read_rds(file = "data/azimuth/project_sc.rds.gz")


refs <- c(
  "Brain" = "mousecortexref",
  # "Meninge" = "pbmcref",
  "Meninge" = "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_dura",
  "Skull" = "bonemarrowref"
)

celllevel <- c(
  "Brain" = "predicted.cluster",
  "Meninge" = "predicted.annotation.l1",
  "Skull" = "predicted.celltype.l2"
)

refcelllevel <- tibble::tibble(
  region = c("Brain", "Meninge", "Skull"),
  refs = c("mousecortexref", "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_dura", "bonemarrowref"),
  celllevel = c("predicted.cluster", "predicted.annotation.l1", "predicted.celltype.l2"),
  supercelllevel = c("predicted.subclass", "predicted.annotation.l1", "predicted.celltype.l1")
)

# body --------------------------------------------------------------------


project_sc |>
  dplyr::left_join(
    refcelllevel,
    by = "region"
  ) |>
  dplyr::mutate(
    anno = purrr::map2(
      .x = sc,
      .y = refs,
      .f = fn_azimuth
    )
  ) ->
  project_sc_azimuth

readr::write_rds(
  x = project_sc_azimuth,
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/project_sc_azimuth_newref.rds"
)

# project_sc_azimuth |>
#   dplyr::mutate(
#     celllevel = plyr::revalue(
#       x = region,
#       replace = celllevel
#     )
#   ) ->
#   project_sc_azimuth_refumap
#
#
#
# readr::write_rds(
#   x = project_sc_azimuth_refumap,
#   file = "data/azimuth/project_sc_azimuth_refumap.rds.gz"
# )



# merge annotation --------------------------------------------------------

project_sc_azimuth_refumap |>
  dplyr::mutate(
    unique_celltype = purrr::map2(
      .x = anno,
      .y = celllevel,
      .f = function(.anno, .celllevel) {
        .anno@meta.data[[.celllevel]] |>
          tibble::enframe() |>
          dplyr::select(celltype = value) |>
          dplyr::group_by(celltype) |>
          dplyr::count() |>
          dplyr::ungroup() |>
          dplyr::arrange(celltype) |>
          dplyr::mutate(ratio = n / sum(n))
      }
    )
  ) ->
  project_sc_azimuth_refumap_unique_celltype


project_sc_azimuth_refumap_unique_celltype |>
  dplyr::filter(region == "Brain") |>
  dplyr::pull(unique_celltype) |>
  dplyr::bind_rows() |>
  dplyr::pull(celltype) |>
  unique() |>
  sort() |>
  tibble::enframe(name = "cluster", value = "celltype") |>
  dplyr::mutate(cluster = factor(cluster)) |>
  dplyr::mutate(cluster_celltype = glue::glue("{cluster}, {celltype}")) ->
  brain_celltype

project_sc_azimuth_refumap_unique_celltype |>
  dplyr::filter(region == "Meninge") |>
  dplyr::pull(unique_celltype) |>
  dplyr::bind_rows() |>
  dplyr::pull(celltype) |>
  unique() |>
  sort() |>
  tibble::enframe(name = "cluster", value = "celltype") |>
  dplyr::mutate(cluster = factor(cluster)) |>
  dplyr::mutate(cluster_celltype = glue::glue("{cluster}, {celltype}")) ->
  meninge_celltype

project_sc_azimuth_refumap_unique_celltype |>
  dplyr::filter(region == "Skull") |>
  dplyr::pull(unique_celltype) |>
  dplyr::bind_rows() |>
  dplyr::pull(celltype) |>
  unique() |>
  sort() |>
  tibble::enframe(name = "cluster", value = "celltype") |>
  dplyr::mutate(cluster = factor(cluster)) |>
  dplyr::mutate(cluster_celltype = glue::glue("{cluster}, {celltype}")) ->
  skull_celltype


project_sc_azimuth_refumap_unique_celltype |>
  dplyr::mutate(
    unique_celltype_union = purrr::map(
      .x = region,
      .f = function(.region) {
        switch(
          EXPR = .region,
          "Brain" = brain_celltype,
          "Meninge" = meninge_celltype,
          "Skull" = skull_celltype
        )
      }
    )
  ) ->
  project_sc_azimuth_refumap_unique_celltype_union

project_sc_azimuth_refumap_unique_celltype_union |>
  dplyr::mutate(
    anno_cell = purrr::pmap(
      .l = list(
        .anno = anno,
        .celllevel = celllevel,
        .unique_celltype_union = unique_celltype_union
      ),
      .f = function(.anno, .celllevel, .unique_celltype_union) {
        .anno_cell <- .anno


        .anno_cell@meta.data |>
          dplyr::mutate(
            celltype = .anno@meta.data[, .celllevel]
          ) |>
          dplyr::left_join(
            .unique_celltype_union,
            by = "celltype"
          ) ->
          .metadata

        .anno_cell@meta.data <- .metadata
        .anno_cell
      }
    )
  ) ->
  project_sc_azimuth_refumap_unique_celltype_union_anno_cell


# plot --------------------------------------------------------------------


project_sc_azimuth_refumap_unique_celltype_union_anno_cell |>
  dplyr::mutate(
    a = purrr::map(
      .x = anno_cell,
      .f = fn_plot_refumap
    )
  ) |>
  tidyr::unnest(cols = a) ->
  project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp


readr::write_rds(
  x = project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp,
  file = "data/azimuth/project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp.rds.gz"
)


# save plot ---------------------------------------------------------------


project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .p = new_p
      ),
      .f = function(.region, .case, .p, .outdir) {
        .filename <- glue::glue("{.region}_{.case}")


        ggsave(
          filename = glue::glue("{.filename}_umap_new.pdf"),
          plot = .p,
          device = "pdf",
          path = .outdir,
          width = 14,
          height = 8
        )
      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth"
    )
  )



# col plot ----------------------------------------------------------------

celllevel |>
  names() |>
  purrr::map(
    .f = function(.region) {
      project_sc_azimuth_refumap_unique_celltype_union_anno_cell_newp |>
        dplyr::select(region, case, cellnumber) |>
        dplyr::filter(region == .region) |>
        tidyr::unnest(cols = cellnumber) ->
        p_d

      p_d |>
        dplyr::select(cluster, celltype) |>
        dplyr::distinct() |>
        dplyr::arrange(cluster) ->
        p_d_fill

      p_d |>
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
          name = "Cell type",
          limits = p_d_fill$cluster,
          labels = p_d_fill$celltype,
          values = pcc$color
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

      ggsave(
        filename = glue::glue("{.region}_celltype_proportion.pdf"),
        plot = .p,
        device = "pdf",
        path = "/home/liuc9/github/scbrain/scuvresult/06-azimuth",
        width = 10,
        height = 8
      )

      p_d |>
        dplyr::select(
          -c(UMAP_1, UMAP_2, celltype_ratio)
        ) |>
        dplyr::mutate(ratio = round(ratio * 100, 2))
    }
  ) ->
  celltype_ratio

names(celltype_ratio) <- c("Brain", "Meninge", "Skull")
writexl::write_xlsx(
  celltype_ratio,
  path = "/home/liuc9/github/scbrain/scuvresult/06-azimuth/celltype_ratio.xlsx"
)



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/02-individual-tissue.rda")


# load(file = "data/azimuth/02-individual-tissue.rda")
