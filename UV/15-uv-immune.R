# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Oct 16 19:20:00 2022
# @DESCRIPTION: immune

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# Load data ---------------------------------------------------------------


UVB <- readr::read_tsv(file = "data/uvrda/ImmuCellAI_mouse_abundance_result_UVB.txt") %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::gather(key = "cell_type", value = "score", - barcode) %>% 
  dplyr::mutate(
    group = gsub(pattern = "_\\d", replacement = "", x = barcode)
  ) %>% 
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "B", x = group), "Brain", "Skull")
  ) %>% 
  dplyr::mutate(
    case = ifelse(group == "BC", "Brain Control", group)
  ) %>% 
  dplyr::mutate(case = ifelse(case == "SC", "Skull Control", case)) %>% 
  dplyr::filter(cell_type != "Infiltration_score")

UVS <- readr::read_tsv(file = "data/uvrda/ImmuCellAI_mouse_abundance_result_UVS.txt") %>% 
  dplyr::rename(barcode = `...1`) %>% 
  tidyr::gather(key = "cell_type", value = "score", - barcode) %>% 
  dplyr::mutate(
    group = gsub(pattern = "_\\d", replacement = "", x = barcode)
  ) %>% 
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "B", x = group), "Brain", "Skull")
  ) %>% 
  dplyr::mutate(
    case = ifelse(group == "BC", "Brain Control", group)
  ) %>% 
  dplyr::mutate(case = ifelse(case == "SC", "Skull Control", case)) %>% 
  dplyr::filter(cell_type != "Infiltration_score")


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
            pval1 = .p1,
            m1 = .m1,
            d1 = list(.x1)
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
    p1 = .uv1_boxplot
  )
  
}

# Plot --------------------------------------------------------------------

UVS_boxplot <- fn_boxplot3(
  .uv = UVS,
  .levels = c("SC", "UVS0", "UVS1")
)

ggsave(
  filename = "UVS-immune-cell-distribution.pdf",
  plot = UVS_boxplot$UV_boxplot,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)


UVB_boxplot <- fn_boxplot3(
  .uv = UVB,
  .levels = c("BC", "UVB0", "UVB1"),
  .celltype = levels(UVS_boxplot$mark$cell_type)
)
ggsave(
  filename = "UVB-immune-cell-distribution.pdf",
  plot = UVB_boxplot$UV_boxplot,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)



# two two -----------------------------------------------------------------

UVS_boxplot2 <- fn_boxplot2(
  .uv = UVS,
  .levels = c("SC", "UVS0", "UVS1")
)

ggsave(
  filename = "UVS-immune-cell-distribution-t0.pdf",
  plot = UVS_boxplot2$p0,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)

ggsave(
  filename = "UVS-immune-cell-distribution-t1.pdf",
  plot = UVS_boxplot2$p1,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)

UVB_boxplot2 <- fn_boxplot2(
  .uv = UVB,
  .levels = c("BC", "UVB0", "UVB1")
)

ggsave(
  filename = "UVB-immune-cell-distribution-t0.pdf",
  plot = UVB_boxplot2$p0,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)

ggsave(
  filename = "UVB-immune-cell-distribution-t1.pdf",
  plot = UVB_boxplot2$p1,
  path = "data/uvresult/02-immune",
  width = 9, height = 4
)

# ggvenn::ggvenn(
#   data = list(
#     UVB_celltype = UVB_boxplot$mark$cell_type,
#     UVS_celltype = UVS_boxplot$mark$cell_type
#   )
# )

# save image --------------------------------------------------------------


save.image(file = "data/uvrda/15-uv-immune.rda")
