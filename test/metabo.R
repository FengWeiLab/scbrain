filename <- "/mnt/isilon/xing_lab/liuc9/projnet/154.08596-Dopamine.xlsx"

d514 <- readxl::read_xlsx(
  path = filename
)

d514 |>
  dplyr::group_by(
    area
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    test = purrr::map(
      .x = data,
      .f = \(.d) {

        .levels <- c("C", "I", "A")

        .d |>
          dplyr::mutate(
            group = factor(group, .levels)
          ) ->
          .dd
        combn(.levels, 2) |>
          as.data.frame() |>
          as.list() |>
          purrr::map(
            .f = \(.xx) {
              t.test(
                value ~ group,
                data = .dd |> dplyr::filter(group %in% .xx)
              ) |>
                broom::tidy() |>
                dplyr::select(1, 2, p.value) |>
                dplyr::mutate(vs = glue::glue("{.xx[[1]]} vs. {.xx[[2]]}"))
            }
          ) |>
          dplyr::bind_rows()
      }
    )
  ) |>
  dplyr::select(area, test) |>
  tidyr::unnest(cols = test) ->
  d514_test

writexl::write_xlsx(
  x = d514_test,
  path = glue::glue("{filename}.compare.xlsx"),
)

d514 |>
  dplyr::group_by(area) |>
  dplyr::summarise(m = mean(value)) |>
  dplyr::arrange(-m) ->
  d514_rank_area

d514 |>
  dplyr::mutate(
    group = factor(group, c("C", "I", "A")),
    area = factor(area, d514_rank_area$area)
  ) |>
  # dplyr::filter(mode == "neg") |>
  ggplot(aes(
    x = area,
    y = value,
    fill = group,
  )) +
  stat_boxplot(
    position = position_dodge(width = 0.75),
    geom = "errorbar",
    width = 0.5,
  ) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.0, 0.1)),
    # limits = c(0, 50)
  ) +
  scale_fill_manual(
    name = "Group",
    values = c("#AEB5B5", "#CD2722", "#2E4D87")
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.text = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(
      colour = "black",
      size = 14,
      face = "bold"
    ),
    axis.title.x = element_blank(),
    legend.position = c(0.9, 0.75),
    legend.key = element_blank()
  ) +
  labs(
    y = "Expression"
  ) ->p1;p1

ggsave(
  filename = glue::glue("{filename}.pdf"),
  plot = p1,
  device = "pdf",
  width = 12,
  height = 7
)



filename <- "/mnt/isilon/xing_lab/liuc9/projnet/215.03273-D-Glucose.xlsx"

d514 <- readxl::read_xlsx(
  path = filename
)


d514 |>
  dplyr::filter(mode == "neg") |>
  dplyr::group_by(
    area
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    test = purrr::map(
      .x = data,
      .f = \(.d) {

        .levels <- c("C", "I", "A")

        .d |>
          dplyr::mutate(
            group = factor(group, .levels)
          ) ->
          .dd
        combn(.levels, 2) |>
          as.data.frame() |>
          as.list() |>
          purrr::map(
            .f = \(.xx) {
              t.test(
                value ~ group,
                data = .dd |> dplyr::filter(group %in% .xx)
              ) |>
                broom::tidy() |>
                dplyr::select(1, 2, p.value) |>
                dplyr::mutate(vs = glue::glue("{.xx[[1]]} vs. {.xx[[2]]}"))
            }
          ) |>
          dplyr::bind_rows()
      }
    )
  ) |>
  dplyr::select(area, test) |>
  tidyr::unnest(cols = test) ->
  d514_test

writexl::write_xlsx(
  x = d514_test,
  path = glue::glue("{filename}.compare.xlsx"),
)

d514 |>
  dplyr::group_by(area) |>
  dplyr::summarise(m = mean(value)) |>
  dplyr::arrange(-m) ->
  d514_rank_area

d514 |>
  dplyr::mutate(
    group = factor(group, c("C", "I", "A")),
    area = factor(area, d514_rank_area$area)
  ) |>
  dplyr::filter(mode == "neg") |>
  ggplot(aes(
    x = area,
    y = value,
    fill = group,
  )) +
  stat_boxplot(
    position = position_dodge(width = 0.75),
    geom = "errorbar",
    width = 0.5,
  ) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.0, 0.1)),
    limits = c(0, 50)
  ) +
  scale_fill_manual(
    name = "Group",
    values = c("#AEB5B5", "#CD2722", "#2E4D87")
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.text = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(
      colour = "black",
      size = 14,
      face = "bold"
    ),
    axis.title.x = element_blank(),
    legend.position = c(0.9, 0.75),
    legend.key = element_blank()
  ) +
  labs(
    y = "Expression"
  ) ->p2;p2


ggsave(
  filename = glue::glue("{filename}.pdf"),
  plot = p2,
  device = "pdf",
  width = 12,
  height = 7
)

