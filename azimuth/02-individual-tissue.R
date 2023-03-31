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
pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_azimuth <- function(.sc, .ref) {
  # .sc <- project_sc$sc[[1]]
  # .ref <- project_sc$ref[[1]]
  
  .sct <- RunAzimuth(
    query = .sc,
    reference = .ref
    # reference = "mousecortexref"
  )
  
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
    dplyr::select(cluster, celltype)  |> 
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
  
  .xxx %>%
    dplyr::group_by(cluster) %>%
    tidyr::nest() %>%
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d %>%
      #   dplyr::filter(cluster == 14) %>%
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
      
      if(nrow(.mmm) == 1) {
        return(
          .mmd %>%
            # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
            dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        )
        
      }
      
      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1
      
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
  
  .split <- if(.facet) {
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
        override.aes = list(size=4)
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


# load data ---------------------------------------------------------------

project_sc <- readr::read_rds(file = "data/azimuth/project_sc.rds.gz")


# body --------------------------------------------------------------------


refs <- c(
  "Brain" = "mousecortexref",
  "Meninge" = "pbmcref",
  "Skull" = "bonemarrowref"
)


project_sc |> 
  dplyr::mutate(
    anno = purrr::map2(
      .x = sc,
      .y = ref,
      .f = fn_azimuth
    )
  ) ->
  project_sc_azimuth


celllevel <- c(
  "Brain" = "predicted.cluster",
  "Meninge" = "predicted.celltype.l2",
  "Skull" = "predicted.celltype.l2"
)

project_sc_azimuth |> 
  dplyr::mutate(
    celllevel = plyr::revalue(
      x = region,
      replace = celllevel
    )
  ) |> 
  dplyr::mutate(
    p = purrr::map2(
      .x = anno,
      .y = celllevel,
      .f = fn_plot_umap_tsne,
      .reduction = "ref.umap"
    )
  ) ->
  project_sc_azimuth_refumap


readr::write_rds(
  x = project_sc_azimuth_refumap,
  file = "data/azimuth/project_sc_azimuth_refumap.rds.gz"
)


project_sc_azimuth_refumap |> 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .p = p
      ),
      .f = function(.region, .case, .p, .outdir) {
        .filename <- glue::glue("{.region}_{.case}")
        
        
        ggsave(
          filename = glue::glue("{.filename}_umap.pdf"),
          plot = .p,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
        
      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth"
    )
  )

# fn_plot_umap_tsne(
#   .x = project_sc_azimuth$anno[[1]],
#   .celltype = "predicted.cluster",
#   .reduction = "ref.umap"
# )
# 
# fn_plot_umap_tsne(
#   .x = project_sc_azimuth$anno[[2]],
#   .celltype = "predicted.celltype.l2",
#   .reduction = "ref.umap"
# )
# 
# fn_plot_umap_tsne(
#   .x = project_sc_azimuth$anno[[3]],
#   .celltype = "predicted.celltype.l2",
#   .reduction = "ref.umap"
# )





# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/azimuth/02-individual-tissue.rda")
