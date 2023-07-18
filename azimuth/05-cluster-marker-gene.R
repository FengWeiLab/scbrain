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

azimuth_ref_sunburst_cell_cell_factor <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/recell_cell_factor.rds"
  )

# azimuth_ref_sunburst_cell <-  readr::read_rds(
#   file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
# )


# merge object ------------------------------------------------------------


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




# cluster plot --------------------------------------------------------------------


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
        .m <- Seurat::FindNeighbors(.m, dims = 1:10)
        .m <- Seurat::RunUMAP(.m, dims = 1:10)
        .m <- Seurat::RunTSNE(.m, dims = 1:10)

        azimuth_ref_sunburst_cell_cell_factor |>
          dplyr::filter(region == .r) ->
          .d

        .d$cell_factor[[1]] -> .cell_factor



        .m$cell1 <- factor(x = .m$cell1, levels = .cell_factor$cell1)
        .m$cell2 <- factor(x = .m$cell2, levels = .cell_factor$cell2)
        .m$cell3 <- factor(x = .m$cell3, levels = .cell_factor$cell3)

        .m$cell1_cluster <- as.factor(as.numeric(.m$cell1))
        .m$cell2_cluster <- as.factor(as.numeric(.m$cell2))
        .m$cell3_cluster <- as.factor(as.numeric(.m$cell3))

        .m
      }
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm


fn_plot_umap_tsne_new <- function(.x, .y) {
  # .x <- azimuth_ref_sunburst_cell_merge_norm$region[[1]]
  # .y <- azimuth_ref_sunburst_cell_merge_norm$norm[[1]]

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
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}"))  ->
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
    geom_text(
      data = .xxx_label,
      aes(
        label = cluster,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 6
    ) +
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
    ) +
    labs(
      x = "tSNE_1",
      y = "tSNE_2"
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
    ) +
    labs(
      x = "tSNE_1",
      y = "tSNE_2"
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


# Marker genes ------------------------------------------------------------


fn_find_all_markers <- function(.norm) {

  Idents(.norm) <- "cell3_cluster"
  future::plan(future::multisession, workers = 10)
  .allmarkers <- Seurat::FindAllMarkers(
    object = .norm,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  future::plan(future::sequential)
  .allmarkers

}

azimuth_ref_sunburst_cell_merge_norm |>
  dplyr::mutate(
    allmarkers = purrr::map(
      .x = norm,
      .f = fn_find_all_markers
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_allmarkers

# azimuth_ref_sunburst_cell_merge_norm_allmarkers |>
#   readr::write_rds(
#     file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_allmarkers.rds.gz"
#   )

azimuth_ref_sunburst_cell_merge_norm_allmarkers <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_allmarkers.rds.gz"
)

# Heatmap -----------------------------------------------------------------


fn_heatmap <- function(.norm,.region, .allmarkers, .topn=5){
  # .norm <- azimuth_ref_sunburst_cell_merge_norm_allmarkers$norm[[1]]
  # .region <- azimuth_ref_sunburst_cell_merge_norm_allmarkers$region[[1]]
  # .allmarkers <- azimuth_ref_sunburst_cell_merge_norm_allmarkers$allmarkers[[1]]
  # .topn <- 10

  Idents(.norm) <- "cell3_cluster"

  .allmarkers |>
    dplyr::group_by(cluster) |>
    dplyr::top_n(.topn, wt = avg_log2FC) ->
    .top

  p <- Seurat::DoHeatmap(.norm, features = .top$gene) + NoLegend()
  p
}

azimuth_ref_sunburst_cell_merge_norm_allmarkers |>
  dplyr::mutate(
    heatmap10 = purrr::pmap(
      .l = list(
        .norm = norm,
        .region = region,
        .allmarkers = allmarkers
      ),
      .f = fn_heatmap
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap

azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .r = region,
        .p = heatmap10
      ),
      .f = function(.r, .p) {
        ggsave(
          filename = "{.r}-markergenes-heatmap.pdf" |> glue::glue(),
          plot = .p,
          width = 20,
          height = 15,
          path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot"
        )
      }
    )
  )


fn_marker_gene_dotplot <- function(object, assay = NULL, features, cols = c(
  "lightgrey",
  "blue"
), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 8,
idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) {

  # object <- .norm
  # assay = NULL
  # features = unique(.marker_head$gene)
  # cols = c("blue", "red")
  # col.min = -2.5
  # col.max = 2.5
  # dot.min = 0
  # dot.scale = 8
  # idents = NULL
  # group.by = NULL
  # split.by = NULL
  # cluster.idents = FALSE
  # scale = TRUE
  # scale.by = "radius"
  # scale.min = NA
  # scale.max = NA

  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in%
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by,
                       size = scale_size,
                       radius = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.",
              call. = FALSE,
              immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(
    object = object, vars = features,
    cells = cells
  )
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(
      rep(x = id.levels, each = length(x = unique.splits)),
      "_", rep(x = unique(x = splits), times = length(x = id.levels))
    )
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1),
                              drop = FALSE
    ]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(
      X = data.use, MARGIN = 2, FUN = PercentAbove,
      threshold = 0
    )
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(
      X = data.plot,
      FUN = unlist
    ))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results",
            call. = FALSE, immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot ==
                              x, "avg.exp"]
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(
          data = data.use, min = col.min,
          max = col.max
        )
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(
      x = avg.exp.scaled,
      breaks = 20
    ))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0(
        "^((",
        paste(sort(x = levels(x = object), decreasing = TRUE),
              collapse = "|"
        ), ")_)"
      ), replacement = "",
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  # object@meta.data |> dplyr::glimpse()
  object@meta.data |>
    dplyr::select(cell3_cluster, cell3_color) |>
    dplyr::distinct() |>
    dplyr::arrange(cell3_cluster) |>
    ggplot(aes(
      x = 1,
      y = cell3_cluster,
      label = cell3_cluster,
      fill = cell3_color
    )) +
    geom_tile(
      aes(
        width = 0.8,
        height = 0.8
      )
    ) +
    scale_fill_identity() +
    geom_text(
      color = "black",
      size = 4
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    coord_fixed() ->
    p1

  data.plot |>
    dplyr::mutate(
      avg.exp.scaled = ifelse(
        avg.exp.scaled > 2,
        2,
        avg.exp.scaled
      )
    ) |>
    dplyr::mutate(
      pct.exp = ifelse(
        pct.exp > 75,
        75,
        pct.exp
      )
    ) |>
    ggplot(mapping = aes_string(
      x = "features.plot",
      y = "id"
    )) +
    geom_point(mapping = aes_string(
      size = "pct.exp",
      color = color.by
    )) +
    scale.func(
      range = c(0, dot.scale),
      limits = c(0, 75),
      breaks = c(0, 25, 50, 75),
      labels = c(0, 25, 50, 75)
    ) +
    scale_color_gradient(
      low = cols[1],
      high = cols[2],
      limits = c(-1, 2),
      name = "Average Expression"
    ) +
    scale_x_discrete(
      position = "top"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.top = element_text(
        angle = 90,
        hjust = 0,
        vjust = 0.5,
        size = 10,
        color = "black",
        face = "italic"
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    guides(size = guide_legend(title = "Percent Expressed (%)")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                    yes = "Identity", no = "Split Identity"
    ))  ->
    p2

  cowplot::plot_grid(
    plotlist = list( p1, p2),
    align = 'h',
    rel_widths = c(0.2, 2)
  ) ->
    plot

  return(plot)
}


fn_gene_dotplot <- function(.norm, .allmarkers, .n = 2) {
  # .norm <- azimuth_ref_sunburst_cell_merge_norm_allmarkers$norm[[1]]
  # .allmarkers <- azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap$allmarkers[[1]]

  .allmarkers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = .n, order_by = avg_log2FC)  ->
    .marker_head

  DefaultAssay(.norm) <- "RNA"
  Idents(.norm) <- "cell3_cluster"


  fn_marker_gene_dotplot(
    .norm,
    features = unique(.marker_head$gene),
    cols = c("blue", "red"),
    dot.scale = 6
  )
}


azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap |>
  dplyr::mutate(
    markerdot = purrr::map2(
      .x = norm,
      .y = allmarkers,
      .f = fn_gene_dotplot
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot

ggsave(
  filename = glue::glue("{azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$region[[1]]}-markergenes-dot.pdf"),
  plot = azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$markerdot[[1]],
  device = "pdf",
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot",
  width = 7,
  height = 3
)


ggsave(
  filename = glue::glue("{azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$region[[2]]}-markergenes-dot.pdf"),
  plot = azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$markerdot[[2]],
  device = "pdf",
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot",
  width = 10,
  height = 5
)

ggsave(
  filename = glue::glue("{azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$region[[3]]}-markergenes-dot.pdf"),
  plot = azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot$markerdot[[3]],
  device = "pdf",
  path = "/home/liuc9/github/scbrain/scuvresult/07-cluster-dot",
  width = 9,
  height = 4
)


# marker gene -------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/05-cluster-marker-gene.rda")
