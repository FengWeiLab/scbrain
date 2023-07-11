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
azimuth_ref_sunburst_cell <-  readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell.rds"
)



# body --------------------------------------------------------------------

azimuth_ref_sunburst_cell |>
  dplyr::filter(region == "Brain") ->
  testa


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










testa$anno_new_new

testam <- merge(
  testa$anno_new_new[[1]],
  y = c(testa$anno_new_new[[2]], testa$anno_new_new[[3]])
)


DimPlot(testam, reduction = "ref.umap")
DefaultAssay(testam) <- "RNA"
testam <- Seurat::NormalizeData(testam)



synap_SIG <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/synap_SIG.rds"
) |>
  dplyr::mutate(geneID = stringr::str_to_upper(geneID))


a3 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[2]]
a6 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[5]]
a9 <- project_sc_azimuth_ref_realcell_sunburst$anno_new[[8]]

a_gene <- union(union(rownames(a3), rownames(a6)),rownames(a9))

gg <- a_gene[grepl("GAB", a_gene)]

synap_SIG |>
  dplyr::filter(geneID %in% a_gene) ->
  a_diff_gene

a_diff_gene |>
  dplyr::filter(grepl("GAB", geneID ))

FeaturePlot(
  object = a3,
  features = "GLS",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)

FeaturePlot(
  object = a6,
  features = "GLS",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)

FeaturePlot(
  object = a6,
  features = "AHSP",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)


Idents(a6) <- "predicted.celltype.l2"
#
VlnPlot(a6, features = "HBM")


a3_diff_gene_m <- AverageExpression(
  a3,
  genes  = a_diff_gene$geneID
)

a3_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a33

a6_diff_gene_m <- AverageExpression(
  a6,
  genes  = a_diff_gene$geneID
)

a6_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a66

a9_diff_gene_m <- AverageExpression(
  a9,
  genes  = a_diff_gene$geneID
)

a9_diff_gene_m$RNA |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "geneID") |>
  dplyr::arrange(-all) |>
  dplyr::filter(geneID %in% a_diff_gene$geneID) ->
  a99

head(a33, 20)
head(a66, 20)
head(a99, 20)

FeaturePlot(
  object = a3,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
FeaturePlot(
  object = a6,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
FeaturePlot(
  object = a9,
  features = "SNAP25",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "ref.umap",
  # max.cutoff = 3,
)
# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
