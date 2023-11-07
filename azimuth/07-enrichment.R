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

#
# azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
#   readr::write_rds(
#     file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich.rds.gz"
#   )


azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich.rds.gz"
  )


fn_plot_enrich <- function(.m, .change="up") {
  # .d <- .dde %>% dplyr::filter(change == .change)
  .go_bp <- .m$gobp[[1]]
  .color <- c("up" = "#EE0000FF", "down" = "#008280FF")

  if(is.null(.m)) {
    return(NULL)
  }

  .go_bp %>%
    tibble::as_tibble()  %>%
    # dplyr::mutate(
    #   Description = stringr::str_wrap(
    #     stringr::str_to_sentence(string = Description),
    #     width = 60
    #   )
    # ) %>%
    dplyr::mutate(Description = stringr::str_to_sentence(string = Description)) |>
    dplyr::mutate(adjp = -log10(p.adjust)) %>%
    dplyr::select(ID, Description, adjp, Count) %>%
    head(20) %>%
    dplyr::arrange(adjp, Count) %>%
    dplyr::mutate(
      Description = factor(Description, levels = Description)
    ) |>
    dplyr::mutate(lb = glue::glue("{Description} (count = {Count})")) ->
    .go_bp_for_plot

  if (.change == "down") {
    .go_bp_for_plot |>
      ggplot(aes(x = Description, y = adjp)) +
      geom_col(fill = .color[.change], color = NA, width = 0.9, alpha = 0.7) +
      scale_y_reverse(expand = expansion(mult=0, add = c(0.5, 0.02))) +
      scale_x_discrete(name = "", position = "top") +
      coord_flip() +
      geom_text(aes(label = lb, y = 0.05), hjust = 1, color = "black") +
      labs(y = "-log10(Adj. P value)") +
      theme(
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_blank(),
        # axis.text.y = element_text(color = "black", size = 13, hjust = 1),
        axis.text.y = element_blank(),
        # axis.ticks.length.y = unit(3, units = "mm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.line.y.right = element_blank()

      ) ->
      .go_bp_plot;.go_bp_plot
  } else {
    .go_bp_for_plot |>
      ggplot(aes(x = adjp, y = Description)) +
      geom_col(fill = .color[.change], color = NA, width = 0.9, alpha = 0.7) +
      geom_text(aes(label = lb, x = 0.05), hjust = 0, color = "black") +
      labs(x = "-log10(Adj. P value)") +
      scale_x_continuous(expand = expansion(mult=0, add = c(0.02, 0.5))) +
      theme(
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_blank(),
        # axis.text.y = element_text(color = "black", size = 13, hjust = 1),
        axis.text.y = element_blank(),
        # axis.ticks.length.y = unit(3, units = "mm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 12, face = "bold"),
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.line.y.left = element_blank()

      ) ->
      .go_bp_plot;.go_bp_plot
  }

  tibble::tibble(
    gobp = list(.go_bp),
    goplot = list(.go_bp_plot)
  )
}

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
  dplyr::mutate(
    b = purrr::map(
      .x = enrichment,
      .f = \(.x) {
        # .x <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$enrichment[[2]]

        .x |>
          dplyr::mutate(
            up_en = purrr::map(
              .x = up_en,
              .f = fn_plot_enrich,
              .change = "up"
            ),
            down_en = purrr::map(
              .x = down_en,
              .f = fn_plot_enrich,
              .change = "down"
            )
          )
      }
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich_newenrichplot


# azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich |>
azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich_newenrichplot |>
  dplyr::mutate(
    a = purrr::map(
      .x = region,
      .y = b,
      .f = function(.r, .e) {
        # .r <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich_newenrichplot$region[[2]]
        # .e <- azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich_newenrichplot$b[[2]]

        dir.create(
          path = file.path(
            "/home/liuc9/github/scbrain/scuvresult/10-enrichment-5",
            .r
          ),
          recursive = T
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
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-5",
                        .r
                      ),
                      width = 3,
                      height = 5
                    )

                    writexl::write_xlsx(
                      x = as.data.frame(.down_en$gobp[[1]]),
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-5",
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
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-5",
                        .r
                      ),
                      width = 3,
                      height = 5
                    )

                    writexl::write_xlsx(
                      x = as.data.frame(.up_en$gobp[[1]]),
                      path = file.path(
                        "/home/liuc9/github/scbrain/scuvresult/10-enrichment-5",
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
