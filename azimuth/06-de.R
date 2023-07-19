# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Jul 18 17:27:43 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_allmarkers_heatmap_markerdot_feature_gene_integration.rds.gz"
  ) |>
  dplyr::select(
    region, norm, allmarkers
  )


# body --------------------------------------------------------------------

fn_case_de <- function(.x) {
  DefaultAssay(.x) <- "RNA"
  Idents(.x) <- "case"

  .cells <- .x$cell3 |> levels()

  tibble::tibble(
    cell = .cells
  ) |>
    dplyr::mutate(
      a = purrr::map(
        .x = cell,
        .f = function(.c) {
          .cc <- subset(
            .x,
            cell3 == .c
          )
          .cc_avg <- AverageExpression(.cc, verbose = FALSE)$RNA |> as.data.frame()
          .mcao_vs_sham <- FindMarkers(.cc, ident.1 = "MCAO", ident.2 = "Sham")
          .uv_vs_sham <-  FindMarkers(.cc, ident.1 = "UV", ident.2 = "Sham")
          .uv_vs_mcao <- FindMarkers(.cc, ident.1 = "UV", ident.2 = "MCAO")
          tibble::tibble(
            mcao_vs_sham = list(.mcao_vs_sham),
            uv_vs_sham = list(.uv_vs_sham),
            uv_vs_mcao = list(.uv_vs_mcao),
            case_avg = list(.cc_avg)
          )
        }
      )
    ) |>
    tidyr::unnest(cols = a) ->
    .cell_diff

  .mcao_vs_sham <- FindMarkers(.x, ident.1 = "MCAO", ident.2 = "Sham")
  .uv_vs_sham <-  FindMarkers(.x, ident.1 = "UV", ident.2 = "Sham")
  .uv_vs_mcao <- FindMarkers(.x, ident.1 = "UV", ident.2 = "MCAO")
  .cc_avg <- AverageExpression(.x, verbose = FALSE)$RNA |> as.data.frame()

  .bulk_diff <- tibble::tibble(
    cell = "Pseudo bulk",
    mcao_vs_sham = list(.mcao_vs_sham),
    uv_vs_sham = list(.uv_vs_sham),
    uv_vs_mcao = list(.uv_vs_mcao),
    case_avg = list(.cc_avg)
  )

  dplyr::bind_rows(
  .cell_diff,
  .bulk_diff
  )
}

azimuth_ref_sunburst_cell_merge_norm |>
  dplyr::mutate(
    de = purrr::map(
      .x = norm,
      .f = fn_case_de
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de

azimuth_ref_sunburst_cell_merge_norm_de |>
  readr::write_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de.rds.gz"
  )

# tmp ---------------------------------------------------------------------

DefaultAssay(a)

a <- azimuth_ref_sunburst_cell_merge_norm$norm[[3]]

a@meta.data |> dplyr::glimpse()

Idents(a) <- "case"
Idents(a)
future::plan(future::multisession, workers = 10)
m <- FindMarkers(a, ident.1 = "MCAO", ident.2 = "Sham")
future::plan(future::sequential)

m |>
  dplyr::arrange(
    p_val_adj,
    -avg_log2FC
  ) |> head()

FeaturePlot(
  a,
  features = c("Lyz"),
  split.by = "case"
)




azimuth_ref_sunburst_cell_merge_norm

a <- azimuth_ref_sunburst_cell_merge_norm$norm[[3]]

a@meta.data |> dplyr::glimpse()
b <- subset(a, subset = cell2 == "Astrocyte")
Idents(b) <- "case"
avg.b <- log1p(AverageExpression(b, verbose = FALSE)$RNA) |> as.data.frame()
avg.b$gene <- rownames(avg.b)

ggplot(avg.b, aes(MCAO, UV)) + geom_point()


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/06-de.rda")
