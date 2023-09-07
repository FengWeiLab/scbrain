#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Sep  7 14:49:19 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_boxplot3 <- function(.uv, .levels = c("BC", "UVB0", "UVB1"), .celltype = NULL) {

  .uv %>%
    dplyr::group_by(cell_type) %>%
    tidyr::nest() %>%
    # dplyr::mutate(
    #   f = purrr::map_lgl(
    #     .x = data,
    #     .f = function(.x) {
    #       all(.x$score < 0.01)
    #     }
    #   )
    # ) %>%
    # dplyr::filter(!f) %>%
    # dplyr::select(-f) %>%
    dplyr::mutate(
      test = purrr::map(
        .x = data,
        .f = function(.x) {

          .xx <- kruskal.test(
            score ~ group,
            data = .x
          )

          .p <- broom::tidy(.xx)$p.value
          .m <- mean(.x$score)

          tibble::tibble(
            pval = .p,
            m = .m
          )

        }
      )
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(col = test) %>%
    dplyr::arrange(-m) %>%
    dplyr::mutate(
      cell_type = factor(
        x = cell_type,
        levels = cell_type
      )
    ) %>%
    dplyr::mutate(
      sig = ifelse(
        test = pval < 0.05,
        yes = "sig",
        no = "non-sig"
      )
    ) %>%
    tidyr::unnest(cols = data) %>%
    dplyr::mutate(
      group = factor(
        x = group,
        levels = .levels
      )
    ) ->
    .uv_forplot

  if(!is.null(.celltype)) {
    .uv_forplot %>%
      dplyr::mutate(
        cell_type = factor(
          x = cell_type,
          levels = .celltype
        )
      ) ->
      .uv_forplot
  }

  .uv_forplot %>%
    dplyr::select(cell_type, pval, m, sig) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      mark = dplyr::case_when(
        pval < 0.05 ~ "**",
        pval < 0.1 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::mutate(score = m + 0.2) ->
    .uv_forplot_mark



  .uv_forplot %>%
    ggplot(
      aes(
        x = cell_type,
        y = score,
        color = group
      )
    ) +
    geom_boxplot() +
    scale_color_manual(
      values = c("#00F5FF", "#CD0000",  "#191970"),
      name = "Group"
    ) +
    labs(
      x = "Immune cell types",
      y = "Immune infiltration score"
    ) +
    annotate(
      geom = "text",
      x = .uv_forplot_mark$cell_type,
      y = .uv_forplot_mark$score,
      label = .uv_forplot_mark$mark
    ) +
    theme(
      panel.background = element_rect(fill = NA),

      panel.grid = element_blank(),

      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 16),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),

      legend.position = c(0.9, 0.7),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) ->
    .uv_boxplot

  list(
    UV_forplot = .uv_forplot,
    mark = .uv_forplot_mark,
    UV_boxplot = .uv_boxplot
  )
}

fn_boxplot2 <- function(.uv, .levels = c("SC", "UVS0", "UVS1")) {
  .uv %>%
    dplyr::group_by(cell_type) %>%
    tidyr::nest() %>%
    dplyr::ungroup() ->
    .uvd

  .uvd %>%
    dplyr::mutate(
      a = purrr::map(
        .x = data,
        .f = function(.x) {
          .x %>%
            dplyr::mutate(
              group = factor(
                x = group,
                levels = c(.levels[1], .levels[2])
              )
            ) %>%
            dplyr::filter(!is.na(group)) ->
            .x0

          .x0p <- t.test(score ~ group, data = .x0)

          .p0 <- broom::tidy(.x0p)$p.value
          .m0 <- mean(.x0$score)

          .x %>%
            dplyr::mutate(
              group = factor(
                x = group,
                levels = c(.levels[1], .levels[3])
              )
            ) %>%
            dplyr::filter(!is.na(group)) ->
            .x1

          .x1p <- t.test(score ~ group, data = .x1)

          .p1 <- broom::tidy(.x1p)$p.value
          .m1 <- mean(.x1$score)

          tibble::tibble(
            pval0 = .p0,
            m0 = .m0,
            d0 = list(.x0),
            d0p = list(
              .x0p |>
                broom::tidy() |>
                dplyr::select(
                  df = estimate,
                  case = estimate1,
                  control = estimate2,
                  p_val = p.value
                )
            ),
            pval1 = .p1,
            m1 = .m1,
            d1 = list(.x1),
            d1p = list(
              .x1p |>
                broom::tidy() |>
                dplyr::select(
                  df = estimate,
                  case = estimate1,
                  control = estimate2,
                  p_val = p.value
                )
            )
          )
        }
      )
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = a) ->
    .uvt

  .uvt %>%
    dplyr::mutate(
      sig0 = ifelse(
        test = pval0 < 0.05,
        yes = "sig",
        no = "non-sig"
      )
    ) %>%
    dplyr::mutate(
      sig1 = ifelse(
        test = pval1 < 0.05,
        yes = "sig",
        no = "non-sig"
      )
    ) ->
    .uvt_sig

  .uvt_sig %>%
    dplyr::arrange(-m0) %>%
    dplyr::mutate(
      cell_type = factor(
        x = cell_type,
        levels = cell_type
      )
    ) %>%
    tidyr::unnest(cols = c(d0)) ->
    .uvt_sig0_forplot

  .uvt_sig %>%
    dplyr::select(cell_type, pval0, m0, sig0) %>%
    dplyr::mutate(
      mark = dplyr::case_when(
        pval0 < 0.05 ~ "**",
        pval0 < 0.1 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::mutate(score = m0 + 0.2) ->
    .uvt_sig0_forplot_mark

  .uvt_sig0_forplot %>%
    ggplot(
      aes(
        x = cell_type,
        y = score,
        color = group
      )
    ) +
    geom_boxplot() +
    scale_color_manual(
      values = c("#00F5FF", "#CD0000",  "#191970"),
      name = "Group"
    ) +
    labs(
      x = "Immune cell types",
      y = "Immune infiltration score"
    ) +
    annotate(
      geom = "text",
      x = .uvt_sig0_forplot_mark$cell_type,
      y = .uvt_sig0_forplot_mark$score,
      label = .uvt_sig0_forplot_mark$mark
    ) +
    theme(
      panel.background = element_rect(fill = NA),

      panel.grid = element_blank(),

      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 16),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),

      legend.position = c(0.9, 0.7),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) ->
    .uv0_boxplot

  .uvt_sig %>%
    dplyr::arrange(-m1) %>%
    dplyr::mutate(
      cell_type = factor(
        x = cell_type,
        levels = cell_type
      )
    ) %>%
    tidyr::unnest(cols = c(d1)) ->
    .uvt_sig1_forplot

  .uvt_sig %>%
    dplyr::select(cell_type, pval1, m1, sig1) %>%
    dplyr::mutate(
      mark = dplyr::case_when(
        pval1 < 0.05 ~ "**",
        pval1 < 0.1 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::mutate(score = m1 + 0.2) ->
    .uvt_sig1_forplot_mark

  .uvt_sig1_forplot %>%
    ggplot(
      aes(
        x = cell_type,
        y = score,
        color = group
      )
    ) +
    geom_boxplot() +
    scale_color_manual(
      values = c("#00F5FF", "#CD0000",  "#191970"),
      name = "Group"
    ) +
    labs(
      x = "Immune cell types",
      y = "Immune infiltration score"
    ) +
    annotate(
      geom = "text",
      x = .uvt_sig1_forplot_mark$cell_type,
      y = .uvt_sig1_forplot_mark$score,
      label = .uvt_sig1_forplot_mark$mark
    ) +
    theme(
      panel.background = element_rect(fill = NA),

      panel.grid = element_blank(),

      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 16),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),

      legend.position = c(0.9, 0.7),
      legend.background = element_blank(),
      legend.key = element_blank(),
    ) ->
    .uv1_boxplot

  list(
    p0 = .uv0_boxplot,
    p0_forplot = .uvt_sig0_forplot,
    p0_mark = .uvt_sig0_forplot_mark,
    p1 = .uv1_boxplot,
    p1_forplot = .uvt_sig1_forplot,
    p1_mark = .uvt_sig1_forplot_mark
  )

}


# load data ---------------------------------------------------------------

barcode_group <- readr::read_tsv(
  file.path(
    "/home/liuc9/github/scbrain/data/uvrdanew",
    "barcode_group.tsv"
  )
)

list(
  "Brain", "Meninge", "Skull"
) |>
  purrr::map(
    .f = \(.x) {
      readr::read_tsv(
        file = file.path(
          "/home/liuc9/github/scbrain/data/uvrdanew",
          "ImmuCellAI_mouse_abundance_result_{.x}.txt" |> glue::glue()
        )
      ) |>
        dplyr::rename(barcode = `...1`) %>%
        tidyr::gather(key = "cell_type", value = "score", - barcode) |>
        dplyr::filter(cell_type != "Infiltration_score")
    }
  ) |>
  dplyr::bind_rows() |>
  dplyr::left_join(barcode_group, by = "barcode") ->
  immunes


# body --------------------------------------------------------------------

immunes |>
  dplyr::group_by(seq) |>
  tidyr::nest() |>
  dplyr::ungroup() ->
  immunes_nest

immunes_nest |>
  dplyr::mutate(
    box3 = purrr::map2(
      .x = seq,
      .y = data,
      .f = \(.seq, .data) {
        fn_boxplot3(
          .uv = .data ,
          .levels = .data |>
            dplyr::pull(group) |>
            unique()
        )
      }
    )
  ) |>
  dplyr::mutate(
    box2 = purrr::map2(
      .x = seq,
      .y = data,
      .f = \(.seq, .data) {
        fn_boxplot2(
          .uv = .data ,
          .levels = .data |>
            dplyr::pull(group) |>
            unique()
        )
      }
    )
  ) ->
  immunes_boxplot

immunes_boxplot |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        seq = seq,
        box3 = box3,
        box2 = box2
      ),
      .f = \(seq, box3, box2, outdir) {
        ggsave(
          plot = box3$UV_boxplot,
          path = outdir,
          width = 9, height = 4,
          filename = "{seq}-immune-cell-distribution.pdf" |> glue::glue()
        )
        ggsave(
          plot = box2$p0,
          path = outdir,
          width = 9, height = 4,
          filename = "{seq}-immune-cell-distribution-t0.pdf" |> glue::glue()
        )
        ggsave(
          plot = box2$p1,
          path = outdir,
          width = 9, height = 4,
          filename = "{seq}-immune-cell-distribution-t1.pdf" |> glue::glue()
        )
      },
      outdir = "/home/liuc9/github/scbrain/data/uvresultnew/04-immune"
    )
  )

# Bubble ------------------------------------------------------------------


immunes_boxplot |>
  dplyr::mutate(
    a = purrr::map2(
      .x = seq,
      .y = box2,
      .f = \(.x, .y) {
        .y$p0_forplot |>
          dplyr::select(cell_type, d0p) |>
          dplyr::distinct() |>
          dplyr::mutate(type = glue::glue("{.x} t0")) |>
          tidyr::unnest(cols = d0p)  ->
          t0

        .y$p0_forplot |>
          dplyr::select(cell_type, d1p) |>
          dplyr::distinct() |>
          dplyr::mutate(type = glue::glue("{.x} t1")) |>
          tidyr::unnest(cols = d1p) ->
          t1

        dplyr::bind_rows(
          t0, t1
        ) |>
          dplyr::filter(
            !(case < 0.0001 | control < 0.0001)
          ) |>
          dplyr::arrange(cell_type) |>
          dplyr::mutate(
            cell_type = gsub(
              pattern = "_",
              replacement = " ",
              x = cell_type
            )
          ) |>
          dplyr::mutate(p_val_log = -log10(p_val)) ->
          forplot
      }
    )
  ) |>
  dplyr::select(seq, a) |>
  tidyr::unnest(cols = a) |>
  dplyr::filter(
    !(case < 0.0001 | control < 0.0001)
  ) |>
  dplyr::arrange(cell_type) |>
  dplyr::mutate(
    cell_type = gsub(
      pattern = "_",
      replacement = " ",
      x = cell_type
    )
  ) |>
  dplyr::mutate(p_val_log = -log10(p_val)) ->
  forplot
#
# immunes_boxplot$box2[[1]]$p0_forplot |>
#   dplyr::select(cell_type, d0p) |>
#   dplyr::distinct() |>
#   dplyr::mutate(type = "t0") |>
#   tidyr::unnest(cols = d0p) ->
#   a0
#
# immunes_boxplot$box2[[1]]$p0_forplot |>
#   dplyr::select(cell_type, d1p) |>
#   dplyr::distinct() |>
#   dplyr::mutate(type = "t1") |>
#   tidyr::unnest(cols = d1p) ->
#   a1
# dplyr::bind_rows(
#   a0, a1
# ) |>
#   dplyr::filter(
#     !(case < 0.0001 | control < 0.0001)
#   ) |>
#   dplyr::arrange(cell_type) |>
#   dplyr::mutate(
#     cell_type = gsub(
#       pattern = "_",
#       replacement = " ",
#       x = cell_type
#     )
#   ) |>
#   dplyr::mutate(p_val_log = -log10(p_val)) ->
#   forplot

forplot |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(r = sum(df)) |>
  dplyr::arrange(-r) ->
  rank_celltype

forplot |>
  ggplot(aes(
    x = cell_type,
    y = type,
    fill = df,
  )) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = c(-0.08, -0.05, -0.03, 0, 0.03, 0.05, 0.08),
    labels = c(-0.08, -0.05, -0.03, 0, 0.03, 0.05, 0.08) |> as.character(),
    name = "Diff"
  ) +
  scale_x_discrete(
    limits = rank_celltype$cell_type,
    position = "top",
    expand = expansion(add = c(0, 2))
  ) +
  scale_y_discrete(
    limits = c("Skull t1", "Skull t0", "", "Meninge t1", "Meninge t0", "", "Brain t1", "Brain t0"),
    labels = c("Skull 24h", "Skull 0.5h", "", "Meninge 24h", "Meninge 0.5h", "", "Brain 24h", "Brain 0.5h")
  ) +
  geom_point(
    data = forplot |> dplyr::filter(p_val < 0.1),
    shape = "asterisk"
  ) +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_text(vjust = 0.8),
    legend.key.width = unit(30, units = "points"),
    legend.key.height = unit(12, units = "points"),
    axis.title = element_blank(),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.ticks = element_blank(),
  ) ->
  p_tile;p_tile

ggsave(
  filename = "immune-tile.pdf",
  plot = p_tile,
  path = "/home/liuc9/github/scbrain/data/uvresultnew/04-immune",
  width = 12, height = 8
)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/uvrdanew/05-uv-immune.rda")
