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
# psudo bulk
# https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

# https://academic.oup.com/bib/article/23/5/bbac286/6649780

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
          .cc_avg <-
            tryCatch(
              expr = {
                  AverageExpression(.cc, verbose = FALSE)$RNA |> as.data.frame()
                },
              error = function(e) {
                NULL
              }
            )
          .mcao_vs_sham <- tryCatch(
            expr = {
              FindMarkers(.cc, ident.1 = "MCAO", ident.2 = "Sham")
            },
            error = function(e) {
              NULL
            }
          )
          .uv_vs_sham <-  tryCatch(
            expr = {
              FindMarkers(.cc, ident.1 = "UV", ident.2 = "Sham")
            },
            error = function(e) {
              NULL
            }
          )
          .uv_vs_mcao <- tryCatch(
            expr = {
              FindMarkers(.cc, ident.1 = "UV", ident.2 = "MCAO")
            },
            error = function(e) {
              NULL
            }
          )
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

# azimuth_ref_sunburst_cell_merge_norm_de |>
#   readr::write_rds(
#     file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de.rds.gz"
#   )

# New load ----------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm_de <- readr::read_rds(
  file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de.rds.gz"
)


# Prop changes ------------------------------------------------------------

fn_plot_prop_bulb <- function(.norm) {
  .norm@meta.data |>
    dplyr::select(case, cell3, cell3_color) |>
    dplyr::filter(cell3 != "Platelet") |>
    dplyr::count(case, cell3, cell3_color) |>
    dplyr::group_by(case) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup() ->
    .d_n_prop

  .d_n_prop |>
    dplyr::select(-prop) |>
    tidyr::spread(key = case, value = n) |>
    dplyr::select(-cell3_color) |>
    tidyr::replace_na(
      replace = list(
        MCAO = 0,
        Sham = 0,
        UV = 0
      )
    ) ->
    .d_n_prop_n

  .d_n_prop |>
    dplyr::select(-n) |>
    tidyr::spread(key = case, value = prop) |>
    tidyr::replace_na(
      replace = list(
        MCAO = 0,
        Sham = 0,
        UV = 0
      )
    ) |>
    dplyr::mutate(mcao_vs_sham = MCAO / Sham) |>
    dplyr::mutate(
      mcao_vs_sham = dplyr::case_when(
        mcao_vs_sham >= 16 ~ 16,
        mcao_vs_sham <= 1/16 ~ 1/16,
        is.nan(mcao_vs_sham) ~ 1,
        TRUE ~ mcao_vs_sham
      )
    ) |>
    dplyr::mutate(mcao_vs_sham_log2 = log2(mcao_vs_sham)) |>
    dplyr::mutate(mcao_vs_sham_m = (MCAO + Sham) / 2) |>
    dplyr::mutate(uv_vs_mcao = UV / MCAO) |>
    dplyr::mutate(
      uv_vs_mcao = dplyr::case_when(
        uv_vs_mcao >= 16 ~ 16,
        uv_vs_mcao <= 1/16 ~ 1/16,
        is.nan(uv_vs_mcao) ~ 1,
        TRUE ~ uv_vs_mcao
      )
    ) |>
    dplyr::mutate(uv_vs_mcao_log2 = log2(uv_vs_mcao)) |>
    dplyr::mutate(uv_vs_mcao_m = (UV + MCAO) / 2) ->
    .prop


  .prop |>
    ggplot(
      aes(
        x = mcao_vs_sham_log2,
        y = forcats::fct_reorder(cell3, mcao_vs_sham_log2) ,
        size = mcao_vs_sham_m,
        color = cell3_color
      )
    ) +
    geom_point() +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_identity() +
    scale_size(name = "Mean Prop of tMCAO and Sham") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_line(
        colour = "grey",
        linetype = "dashed"
      ),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        linewidth = 0.2
      ),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(
        color = "black",
        size = 12,
      ),
      axis.title = element_text(
        color = "black",
        size = 14,
      ),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.key = element_blank()
    ) +
    labs(
      x = "log2 fold change of tMCAO vs. Sham",
    ) ->
    p_tmca_vs_sham


  .prop |>
    ggplot(
      aes(
        x = uv_vs_mcao_log2,
        y = forcats::fct_reorder(cell3, uv_vs_mcao_log2) ,
        size = uv_vs_mcao_m,
        color = cell3_color
      )
    ) +
    geom_point() +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_identity() +
    scale_size(name = "Mean Prop of UVB+tMCAO and tMCAO") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_line(
        colour = "grey",
        linetype = "dashed"
      ),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        linewidth = 0.2
      ),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(
        color = "black",
        size = 12,
      ),
      axis.title = element_text(
        color = "black",
        size = 14,
      ),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.key = element_blank()
    ) +
    labs(
      x = "log2 fold change of UVB+tMCAO and tMCAO",
    ) ->
    p_uv_vs_mcao

  .prop |>
    dplyr::select(cell3, cell3_color, mcao_vs_sham_log2, uv_vs_mcao_log2, mcao_vs_sham_m, uv_vs_mcao_m) |>
    dplyr::mutate(cell3 = forcats::fct_reorder(cell3, uv_vs_mcao_log2)) ->
    .prop_forp

  .prop_forp |>
    dplyr::select(-c(mcao_vs_sham_m, uv_vs_mcao_m)) |>
    tidyr::pivot_longer(
      cols = -c(cell3, cell3_color),
      names_to = "vs",
      values_to = "v"
    ) |>
    dplyr::mutate(
      vs_l = gsub(
        pattern = "_log2$",
        replacement = "",
        x = vs
      )
    ) |>
    dplyr::left_join(
      .prop_forp |>
        dplyr::select(1, 5, 6) |>
        tidyr::pivot_longer(
          cols = -cell3,
          names_to = "vs_l",
          values_to = "m"
        ) |>
        dplyr::mutate(
          vs_l = gsub(
            pattern = "_m$",
            replacement = "",
            x = vs_l
          )
        ),
      by = c("cell3", "vs_l")
    ) |>
    ggplot(aes(
      x = v,
      y = cell3,
    )) +
    geom_segment(
      data = .prop_forp,
      aes(
        x = mcao_vs_sham_log2,
        y = cell3,
        xend = uv_vs_mcao_log2,
        yend = cell3,
      ),
      arrow = arrow(length = unit(10, "points")),
      # linetype = "dotted"
      linewidth = 0.3
    ) +
    geom_point(
      aes(
        color = vs,
        size = m
      )
    ) +
    scale_color_manual(
      name = "VS",
      limits = c("mcao_vs_sham_log2", "uv_vs_mcao_log2"),
      label = c("tMCAO vs Sham", "UVB+tMCAO vs tMCAO"),
      values = c("blue", "red")
    ) +
    geom_vline(xintercept = 0, color = "red", linetype = "dotted") +
    scale_size(name = "Mean prop") +
    scale_x_continuous(
      limits = c(-4, 4),
      breaks = c(-4, -3,  -2, -1, 0, 1, 2, 3, 4),
      labels = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_line(
        colour = "grey",
        linetype = "dashed"
      ),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        linewidth = 0.2
      ),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(
        color = "black",
        size = 12,
      ),
      axis.title = element_text(
        color = "black",
        size = 14,
      ),
      axis.title.y = element_blank(),
      legend.position = "top",
      legend.key = element_blank()
    ) +
    labs(
      x = "log2 Fold Change",
    ) ->
    p_merge


  list(
    p_tmca_vs_sham = p_tmca_vs_sham,
    p_uv_vs_mcao = p_uv_vs_mcao,
    p_merge = p_merge,
    prop = .prop,
    prop_n = .d_n_prop_n
  )

}

azimuth_ref_sunburst_cell_merge_norm_de |>
  dplyr::mutate(
    prop_p = purrr::map(
      .x = norm,
      .f = fn_plot_prop_bulb
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de_prop

azimuth_ref_sunburst_cell_merge_norm_de_prop |>
  dplyr::mutate(
    a = purrr::map2(
      .x = region,
      .y = prop_p,
      .f = function(.x, .y) {

        ggsave(
          filename = glue::glue("{.x}-prop-ratio-tmcao-vs-sham.pdf"),
          plot = .y$p_tmca_vs_sham,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/scuvresult/08-prop-ratio",
          width = 8,
          height = 5
        )

        ggsave(
          filename = glue::glue("{.x}-prop-ratio-uvb-vs-tmcao.pdf"),
          plot = .y$p_uv_vs_mcao,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/scuvresult/08-prop-ratio",
          width = 8,
          height = 5
        )

        ggsave(
          filename = glue::glue("{.x}-prop-ratio-uvb-merge.pdf"),
          plot = .y$p_merge,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/scuvresult/08-prop-ratio",
          width = 8,
          height = 5
        )

        .y$prop |>
          dplyr::mutate(
            `Sham (%)` = Sham * 100,
            `tMCAO (%)` = MCAO * 100,
            `tMCAO+UVB (%)` = UV * 100,
          ) |>
          dplyr::select(
            cell3,
            `Sham (%)`,
            `tMCAO (%)`,
            `tMCAO+UVB (%)`,
            `log2(tMCAO vs. Sham)` = mcao_vs_sham_log2,
            `log2(UVB+tMCAO vs. tMCAO)` = uv_vs_mcao_log2
          ) ->
          .yy1
        .y$prop_n |>
          dplyr::mutate(
            `Sham (n)` = Sham,
            `tMCAO (n)` = MCAO,
            `tMCAO+UVB (n)` = UV
          ) |>
          dplyr::select(-c(2,3,4)) ->
          .yy2

        .yy1 |>
          dplyr::left_join(.yy2, by = "cell3") |>
          dplyr::rename(
            `Cell type` = cell3
          ) ->
          .yyd

        writexl::write_xlsx(
          x = .yyd |>
            dplyr::select(1, 7, 8, 9, 2, 3, 4, 5, 6),
          "/home/liuc9/github/scbrain/scuvresult/08-prop-ratio/{.x}-celltype-number-prop.xlsx" |> glue::glue(),
        )

      }
    )
  )


# n DEG -------------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm_de_prop$de[[1]] ->.de

.de$mcao_vs_sham[[1]] |>
  dplyr::filter()

azimuth_ref_sunburst_cell_merge_norm_de_prop |>
  dplyr::select(-prop_p) |>
  dplyr::mutate(
    de_change = purrr::map(
      .x = de,
      .f = function(.de) {
        .de |>
          dplyr::select(-uv_vs_sham) |>
          dplyr::mutate(
            mcao_vs_sham = purrr::map(
              .x = mcao_vs_sham,
              .f = function(.m) {
                if(is.null(.m)) {return(NULL)}
                .m |>
                  tibble::rownames_to_column(
                    var = "gene"
                  ) |>
                  dplyr::mutate(
                    change = dplyr::case_when(
                      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "up",
                      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "down",
                      TRUE ~ "nosig"
                    )
                  )
              }
            )
          ) |>
          dplyr::mutate(
            uv_vs_mcao = purrr::map(
              .x = uv_vs_mcao,
              .f = function(.m) {
                if(is.null(.m)) {return(NULL)}
                .m |>
                  tibble::rownames_to_column(
                    var = "gene"
                  ) |>
                  dplyr::mutate(
                    change = dplyr::case_when(
                      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "up",
                      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "down",
                      TRUE ~ "nosig"
                    )
                  )
              }
            )
          )
      }
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de_change



azimuth_ref_sunburst_cell_merge_norm_de_change |>
  dplyr::mutate(
    nn = purrr::map(
      .x = de_change,
      .f = function(.x) {
        .x |>
          dplyr::select(-case_avg) |>
          dplyr::mutate(
            mcao_vs_sham = purrr::map(
              .x = mcao_vs_sham,
              .f = function(.m) {
                if(is.null(.m)) {return(NULL)}
                .m |>
                  dplyr::count(change) |>
                  dplyr::mutate(type = "ms")
              }
            )
          ) |>
          dplyr::mutate(
            uv_vs_mcao = purrr::map(
              .x = uv_vs_mcao,
              .f = function(.m) {
                if(is.null(.m)) {return(NULL)}
                .m |>
                  dplyr::count(change) |>
                  dplyr::mutate(type = "um")
              }
            )
          ) |>
          dplyr::mutate(
            n = purrr::map2(
              .x = mcao_vs_sham,
              .y = uv_vs_mcao,
              .f = function(.x, .y) {
                dplyr::bind_rows(.x, .y)
              }
            )
          ) |>
          dplyr::select(cell, n) |>
          tidyr::unnest(cols = n) |>
          dplyr::filter(cell != "Platelet") ->
          .d

        .d |>
          dplyr::filter(type == "ms") ->
          .d_ms
        .d |>
          dplyr::filter(type == "um") ->
          .d_um

        .d_ms |>
          dplyr::filter(change != "nosig") |>
          dplyr::filter(cell != "Pseudo bulk") |>
          dplyr::mutate(cell = forcats::fct_reorder(cell, n)) |>
          dplyr::mutate(n = ifelse(change == "down", -n, n)) |>
          dplyr::mutate(hjust = ifelse(change == "down", 1, 0)) |>
          ggplot(aes(
            x = n,
            y = cell,
            fill = change
          )) +
          geom_col() +
          geom_text(
            aes(
              label = abs(n),
              hjust = hjust
            ),
            size = 5,
            fontface = 1,
            show.legend = F
          ) +
          scale_x_continuous(
            expand = expansion(
              mult = c(0.1, 0.1)
            )
          ) +
          geom_vline(xintercept = 0, size = 1) +
          scale_fill_brewer(
            palette = "Set1",
            direction = -1,
            name = "Regulations",
            label = c("Down", "Up")
          ) +
          theme(
            plot.background = element_blank(),
            panel.background = element_blank(),
            # panel.grid =  element_blank(),
            panel.grid.major.y = element_line(
              colour = "black",
              linetype = "dotted",
              linewidth = 0.2,
            ),
            # axis.line = element_line(size = 1),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            # axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1, face = "bold"),
            axis.text.x = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            # axis.title = element_text(size =14, color = "black", face = "bold"),
            legend.position = "top",
            plot.title = element_text(size =14, color = "black", face = "bold", hjust = 0.5)
          ) +
          labs(
            title = "tMCAO vs. Sham"
          ) ->
          p_ms


        .d_um |>
          dplyr::filter(change != "nosig") |>
          dplyr::filter(cell != "Pseudo bulk") |>
          dplyr::mutate(cell = forcats::fct_reorder(cell, n)) |>
          dplyr::mutate(n = ifelse(change == "down", -n, n)) |>
          dplyr::mutate(hjust = ifelse(change == "down", 1, 0)) |>
          ggplot(aes(
            x = n,
            y = cell,
            fill = change
          )) +
          geom_col() +
          geom_text(
            aes(
              label = abs(n),
              hjust = hjust
            ),
            size = 5,
            fontface = 1,
            show.legend = F
          ) +
          scale_x_continuous(
            expand = expansion(
              mult = c(0.1, 0.1)
            )
          ) +
          geom_vline(xintercept = 0, size = 1) +
          scale_fill_brewer(
            palette = "Set1",
            direction = -1,
            name = "Regulations",
            label = c("Down", "Up")
          ) +
          theme(
            plot.background = element_blank(),
            panel.background = element_blank(),
            # panel.grid =  element_blank(),
            panel.grid.major.y = element_line(
              colour = "black",
              linetype = "dotted",
              linewidth = 0.2,
            ),
            # axis.line = element_line(size = 1),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            # axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1, face = "bold"),
            axis.text.x = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            # axis.title = element_text(size =14, color = "black", face = "bold"),
            legend.position = "top",
            plot.title = element_text(size =14, color = "black", face = "bold", hjust = 0.5)
          ) +
          labs(
            title = "UVB+tMCAO vs. tMCAO"
          ) ->
          p_um

        list(
          p_ms = p_ms,
          p_um = p_um,
          d = .d
        )
      }
    )
  ) ->
  azimuth_ref_sunburst_cell_merge_norm_de_change_nn

azimuth_ref_sunburst_cell_merge_norm_de_change_nn$nn[[2]] -> .nn
azimuth_ref_sunburst_cell_merge_norm_de_change_nn$region[[2]] -> .r

azimuth_ref_sunburst_cell_merge_norm_de_change_nn |>
  dplyr::mutate(
    a = purrr::map2(
      .x = region,
      .y = nn,
      .f = function(.r, .nn) {
        ggsave(
          filename = glue::glue("{.r}-DEG-n-tMCAO_vs_Sham.pdf"),
          plot = .nn$p_ms,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/scuvresult/09-de",
          width = 8,
          height = 5
        )

        ggsave(
          filename = glue::glue("{.r}-DEG-n-UVB_vs_tMCAO.pdf"),
          plot = .nn$p_um,
          device = "pdf",
          path = "/home/liuc9/github/scbrain/scuvresult/09-de",
          width = 8,
          height = 5
        )

        .nn$d |>
          dplyr::filter(change != "nosig") ->
          .dd

        .dd |>
          dplyr::filter(type == "ms") |>
          dplyr::select(-type) |>
          tidyr::spread(key = change, value = n) |>
          tidyr::replace_na(
            replace = list(
              up = 0,
              down = 0
            )
          ) |>
          dplyr::select(`Cell type` = cell, Down = down, Up = up) |>
          dplyr::mutate(DEG = Down + Up) |>
          dplyr::arrange(-DEG) ->
          .dd_ms

        .dd |>
          dplyr::filter(type == "um") |>
          dplyr::select(-type) |>
          tidyr::replace_na(
            replace = list(
              up = 0,
              down = 0
            )
          ) |>
          tidyr::spread(key = change, value = n) |>
          dplyr::select(`Cell type` = cell, Down = down, Up = up) |>
          dplyr::mutate(DEG = Down + Up) |>
          dplyr::arrange(-DEG) ->
          .dd_um

        list(
          "tMCAO vs. Sham" = .dd_ms,
          "UVB+tMCAO vs. tMCAO" = .dd_um
        ) |>
          writexl::write_xlsx(
            "/home/liuc9/github/scbrain/scuvresult/09-de/{.r}-n-deg.xlsx" |> glue::glue(),
          )

      }
    )
  )


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/06-de.rda")
load(file = "data/azimuth/06-de.rda")
