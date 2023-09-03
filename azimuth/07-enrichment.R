# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Jul 24 15:51:46 2023
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


azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano.rds.gz"
  )

msig_df <- msigdbr::msigdbr(species = "Mus musculus")


msig_df %>%
  dplyr::mutate(gs_name = glue::glue("{gs_cat}#{gs_name}")) %>%
  dplyr::select(gs_name, gene_symbol) ->
  msig_df_s

# body --------------------------------------------------------------------

fn_gobp <- function(.dde, .change="up") {
  .d <- .dde %>% dplyr::filter(change == .change)
  .color <- c("up" = "#AE1700", "down" = "#112a13")

  .go_bp <- clusterProfiler::enrichGO(
    gene = .d$gene,
    # universe = .dde$gene,
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
  )
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
    ) ->
    .go_bp_for_plot

  .go_bp_for_plot |>
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .color[.change], color = NA, width = 0.7) +
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

fn_go_enrichment <- function(.de, .nn) {
  # .de <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano$de_change[[1]]
  # .nn <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano$nn[[1]]

  .nnd <- .nn$d |>
    dplyr::filter(change != "nosig") |>
    tidyr::spread(key = change, value = n) |>
    tidyr::replace_na(
      replace = list(
        up = 0,
        down = 0
      )
    )

  .de |>
    dplyr::select(-case_avg) |>
    tidyr::pivot_longer(
      cols = -cell,
      names_to = "type",
      values_to = "de"
    ) |>
    dplyr::mutate(
      type = dplyr::case_when(
        type == "mcao_vs_sham" ~ "ms",
        type == "uv_vs_mcao" ~ "um"
      )
    ) ->
    .xx

  .nnd |>
    dplyr::left_join(
      .xx, by = c("cell", "type")
    ) ->
    .nndxx

  .nndxx |>
    dplyr::mutate(
      down_en = furrr::future_map2(
        .x = down,
        .y = de,
        .f = function(.down, .dde) {
          if(.down == 0){return(NULL)}
          tryCatch(
            expr = {
              fn_gobp(.dde, .change = "down")
            },
            error = function(e) {return(NULL)}
          )
        }
      )
    ) |>
    dplyr::mutate(
      up_en = furrr::future_map2(
        .x = up,
        .y = de,
        .f = function(.up, .dde) {
          if(.up == 0) {return(NULL)}
          tryCatch(
            expr = {
              fn_gobp(.dde, .change = "up")
            },
            error = function(e) {return(NULL)}
          )
        }
      )
    )
}

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano |>
  dplyr::mutate(
    enrichment = purrr::map2(
      .x = de_change,
      .y = nn,
      .f = fn_go_enrichment
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich


azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  readr::write_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich.rds.gz"
  )


azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  dplyr::mutate(
    a = purrr::map2(
      .x = region,
      .y = enrichment,
      .f = function(.r, .e) {
        # .r <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$region[[1]]
        # .e <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$enrichment[[1]]

        dir.create(
          path = file.path(
            "/home/liuc9/github/scbrain/scuvresult/10-enrichment-2",
            .r
          )
        )

        typeconv <- c("ms" = "tMCAO vs. Sham", "um" = "UVB+tMCAO vs. tMCAO")

        .e |>
          dplyr::mutate(
            a = purrr::pmap(
              .l = list(
                .cell = cell,
                .type = type,
                .down = down,
                .up = up,
                .down_en = down_en,
                .up_en = up_en
              ),
              .f = function(.cell, .type, .down, .up, .down_en, .up_en) {
                .filename <- "{.r}-{.cell}-{typeconv[.type]}-down{.down}-up{.up}" |> glue::glue()

                tryCatch(
                  expr = {
                    ggsave(
                      plot = .down_en$goplot[[1]],
                      filename = glue::glue("{.filename}-down-enrichment.pdf"),
                      device = "pdf",
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-2",
                        .r
                      ),
                      width = 10,
                      height = 6.5
                    )

                    writexl::write_xlsx(
                      x = as.data.frame(.down_en$gobp[[1]]),
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-2",
                        .r,
                        glue::glue("{.filename}-down-enrichment.xlsx")
                      )
                    )

                  },
                  error = function(e) {NULL}
                )

                tryCatch(
                  expr = {
                    ggsave(
                      plot = .up_en$goplot[[1]],
                      filename = glue::glue("{.filename}-up-enrichment.pdf"),
                      device = "pdf",
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-2",
                        .r
                      ),
                      width = 10,
                      height = 6.5
                    )

                    writexl::write_xlsx(
                      x = as.data.frame(.up_en$gobp[[1]]),
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-2",
                        .r,
                        glue::glue("{.filename}-up-enrichment.xlsx")
                      )
                    )

                  },
                  error = function(e) {NULL}
                )

              }
            )
          )

      }
    )
  )

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/07-enrichment.rda")
