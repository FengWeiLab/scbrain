# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Oct 16 16:41:10 2022
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

library(DESeq2)

# Load data ---------------------------------------------------------------

UV.se <- readr::read_rds(
  file = file.path(
    "data/uvrda",
    "UV_count_matrix.se.rds.gz"
  )
)

bc_uvb0 <- c("BC", "UVB0")
bc_uvb1 <- c("BC", "UVB1")
uvb0_uvb1 <- c("UVB0", "UVB1")

sc_uvs0 <- c("SC", "UVS0")
sc_uvs1 <- c("SC", "UVS1")
uvs0_uvs1 <- c("UVS0", "UVS1")

colData(UV.se) %>% as.data.frame() -> se_coldata

list(
  bc_uvb0 = bc_uvb0,
  bc_uvb1 =  bc_uvb1,
  uvb0_uvb1 = uvb0_uvb1,
  
  sc_uvs0 = sc_uvs0,
  sc_uvs1 = sc_uvs1,
  uvs0_uvs1 = uvs0_uvs1
) %>% 
  tibble::enframe(name = "bs", value = "group") %>% 
  dplyr::mutate(
    seq = ifelse(grepl(pattern = "b", x = bs), "Brain", "Skull")
  ) %>% 
  dplyr::mutate(
    se = purrr::map(
      .x = group,
      .f = function(.x) {
        se_coldata %>% 
          dplyr::filter(group %in% .x) ->
          .xx
        
        .xxx <- assay(UV.se)[, .xx$barcode]
        .se <- SummarizedExperiment::SummarizedExperiment(
          assays = .xxx,
          colData = .xx
        )
      }
    )
  ) ->
  se_group

readr::write_rds(
  x = se_group,
  file = file.path(
    "data/uvrda",
    "count_matrix_UV.se.group.rds.gz"
  )
)


# Function ----------------------------------------------------------------

fn_se_de <- function(.x, .y) {
  # .x <- se_group$group[[1]]
  # .y <- se_group$se[[1]]
  
  .x1 <- .x[c(1, 2)]
  # .x2 <- .x[c(1, 3)]
  # .x3 <- .x[c(2, 3)]
  
  .y1 <- .y[,.y@colData[.y$group %in% .x1, ]$barcode]
  .y1$group <- factor(.y1$group, .x1)
  .y1_dds <- DESeqDataSet(.y1, design = ~ group)
  .y1_res <- results(DESeq(.y1_dds))
  
  # .y2 <- .y[,.y@colData[.y$group %in% .x2, ]$barcode]
  # .y2$group <- factor(.y2$group, .x2)
  # .y2_dds <- DESeqDataSet(.y2, design = ~ group)
  # .y2_res <- results(DESeq(.y2_dds))
  # 
  # .y3 <- .y[,.y@colData[.y$group %in% .x3, ]$barcode]
  # .y3$group <- factor(.y3$group, .x3)
  # .y3_dds <- DESeqDataSet(.y3, design = ~ group)
  # .y3_res <- results(DESeq(.y3_dds))
  
  .n1 <- paste0(rev(.x1), collapse = "_vs_")
  # .n2 <- paste0(rev(.x2), collapse = "_vs_")
  # .n3 <- paste0(rev(.x3), collapse = "_vs_")
  
  tibble::tibble(
    n1 = list(.y1_res),
    # n2 = list(.y2_res),
    # n3 = list(.y3_res)
  ) ->
    .d
  # colnames(.d) <- c(.n1, .n2, .n3)
  colnames(.d) <- c(.n1)
  
  .d
}

fn_volcano <- function(.x, .y) {
  # .x <- se_group_de$seq[[2]]
  # .y <- se_group_de$de[[2]]
  .prefix <- .x
  
  .y %>% 
    tidyr::gather(key = "vs", value = "des") ->
    .yy
  
  .yy %>% 
    dplyr::mutate(
      des_color = purrr::map(
        .x = des,
        .f = function(.des) {
          .des %>% 
            as.data.frame() %>%
            tibble::rownames_to_column(var = "GeneName") %>% 
            dplyr::filter(!is.na(!padj)) %>% 
            dplyr::filter(!is.na(log2FoldChange)) %>% 
            # dplyr::mutate(log2FC = ifelse(log2FoldChange > 4, 4, log2FoldChange)) %>% 
            # dplyr::mutate(log2FC = ifelse(log2FC < -4, -4, log2FC)) %>% 
            dplyr::mutate(log2FC = log2FoldChange) %>% 
            dplyr::mutate(FDR = -log10(padj)) %>% 
            dplyr::mutate(
              color = dplyr::case_when(
                log2FC > log2(2) & padj < 0.05 ~ "red",
                log2FC < log2(1/2) & padj < 0.05 ~ "green",
                TRUE ~ "grey"
              )
            )
        }
      )
    ) ->
    .yyy
  
  .yyy %>% 
    dplyr::mutate(
      vp = purrr::map2(
        .x = vs,
        .y = des_color,
        .f = function(.x, .y) {
          # .x <- .yyy$vs[[1]]
          # .y <- .yyy$des_color[[1]]
          .xd <- .y
          
          .xd %>% 
            dplyr::filter(color != "grey") %>% 
            dplyr::group_by(color) %>% 
            dplyr::count() %>% 
            dplyr::ungroup() %>% 
            tibble::deframe() ->
            .xxx
          
          .xd %>% 
            ggplot(aes(x = log2FC, y = FDR, color = color)) +
            geom_point(alpha = 0.8) +
            scale_color_manual(values =  c("#6BAB62", "grey", "#9A84B2")) +
            geom_segment(
              aes(x = x1, y = y1, xend = x2, yend = y2),
              color = "black",
              linetype = 88,
              data = tibble::tibble(
                x1 = c(-Inf, 1, -1, 1),
                y1 = -log10(0.05),
                x2 = c(-1, Inf, -1, 1 ),
                y2 = c(-log10(0.05), -log10(0.05), Inf, Inf))
            ) +
            ggrepel::geom_text_repel(
              aes(label = GeneName),
              data = subset(.xd, color == "red") %>% 
                dplyr::arrange(-FDR, -abs(log2FC)) %>% 
                dplyr::slice(1:10),
              box.padding = 0.5,
              max.overlaps = Inf,
              # size = 6
            ) +
            ggrepel::geom_text_repel(
              aes(label = GeneName),
              data = subset(.xd, color == "green") %>% 
                dplyr::arrange(-FDR, -abs(log2FC)) %>% 
                dplyr::slice(1:10),
              box.padding = 0.5,
              max.overlaps = Inf,
              # size = 6
            ) +
            scale_x_continuous(
              # limits = c(-4, 4),
              expand = c(0.02, 0)
            ) +
            scale_y_continuous(
              expand = c(0.01, 0),
              limits = c(
                0,
                ceiling(
                  max(
                    .xd %>% 
                      dplyr::filter(!is.infinite(FDR)) %>% 
                      dplyr::pull(FDR)
                  ) / 10
                ) * 10
              )
            ) +
            theme(
              panel.background = element_rect(fill = NA, color = NA),
              axis.line.x.bottom = element_line(color = "black"),
              axis.line.y.left = element_line(color = "black"),
              axis.text = element_text(color = "black", size = 16),
              axis.title = element_text(color = "black", size = 18),
              
              legend.position = "none",
              
            ) +
            labs(
              x = "log2FC",
              y = "-log10(FDR)",
              title = glue::glue("{.x}; Up (n={.xxx[2]}), Down (n={.xxx[1]})")
            ) ->
            .p
          .plotfilename <- glue::glue("volcano_plot_{.prefix}_{.x}.pdf")
          ggsave(
            filename = .plotfilename,
            plot = .p,
            device = "pdf",
            path = "data/uvresult/01-de",
            width = 7,
            height = 5
          )
          .p
        }
      )
    )
}

fn_enrichment <- function(.gs, .color = "red") {
  .cc <- c("red" = "#AE1700", "green" = "#112a13")
  .go_bp <- clusterProfiler::enrichGO(
    gene = .gs,
    # universe = .de$GeneName,
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
  
  .go_bp_for_plot %>%
    ggplot(aes(x = Description, y = adjp)) +
    geom_col(fill = .cc[.color], color = NA, width = 0.7) +
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


fn_save_enrichment_xlsx_pdf <- function(.d, .p, .filename, .path = "data/uvresult/01-de") {
  .xlsx_filename <- glue::glue("{.filename}.xlsx")
  .pdf_filename <- glue::glue("{.filename}.pdf")
  
  writexl::write_xlsx(
    x = as.data.frame(.d),
    path = file.path(
      .path,
      .xlsx_filename
    )
  )
  
  ggsave(
    plot = .p,
    filename = .pdf_filename,
    device = "pdf",
    path = .path,
    width = 10,
    height = 6.5
  )
}

# DE ----------------------------------------------------------------------

se_group %>% 
  dplyr::mutate(
    de = purrr::map2(
      .x = group,
      .y = se,
      .f = fn_se_de
    )
  ) ->
  se_group_de


readr::write_rds(
  x = se_group_de,
  file = file.path(
    "data/uvrda",
    "UV_se_group_de.rds.gz"
  )
)


unlist(se_group_de$de) %>% 
  purrr::map(.f = function(.x) {
    .x %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "GeneName") %>% 
      dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1) %>% 
      dplyr::arrange(padj, -abs(log2FoldChange))
  }) %>% 
  writexl::write_xlsx(
    path = file.path(
      "data/uvresult/01-de",
      "01-basic_de.xlsx"
    )
  )

# volcano plot ------------------------------------------------------------


se_group_de %>% 
  dplyr::mutate(
    volcano_plot = purrr::map2(
      .x = seq,
      .y = de,
      .f = fn_volcano
    )
  ) %>% 
  tidyr::unnest(cols = volcano_plot) ->
  se_group_de_volcano


readr::write_rds(
  x = se_group_de_volcano,
  file = "data/uvrda/se_group_de_volcano.rds.gz"
)

(se_group_de_volcano$vp[[1]] | se_group_de_volcano$vp[[2]] | se_group_de_volcano$vp[[3]]) / (se_group_de_volcano$vp[[4]] | se_group_de_volcano$vp[[5]] | se_group_de_volcano$vp[[6]])  +
  plot_annotation(
    title = glue::glue("Differential expression"),  
    tag_levels = "A"
  ) ->
  .p

ggsave(
  filename = "Brain-skull-differential.pdf",
  plot = .p,
  device = "pdf",
  path = "data/uvresult/01-de",
  width = 17,
  height = 12
)


# Barplot -----------------------------------------------------------------
se_group_de_volcano %>% 
  dplyr::select(vs, des_color) %>% 
  dplyr::mutate(n = purrr::map(
    .x = des_color,
    .f = function(.x) {
      .x %>% 
        dplyr::group_by(color) %>% 
        dplyr::count() %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(color != "grey") %>% 
        tidyr::spread(key = color, value = n)
    }
  )) %>% 
  dplyr::select(-des_color) %>% 
  tidyr::unnest(cols = n) %>% 
  tidyr::gather(
    key = color,
    value = n,
    -vs
  ) ->
  se_group_de_volcano_n

xlimits <- c(se_group_de_volcano$vs[c(1:3)], "", se_group_de_volcano$vs[c(4:6)])
xlabels <- gsub(pattern = "_vs_", replacement = " vs. ", x = xlimits)


se_group_de_volcano_n %>% 
  dplyr::mutate(n = ifelse(color == "green", -n, n)) %>% 
  dplyr::mutate(vjust = ifelse(color == "green", 1.5, -0.5)) %>% 
  ggplot(aes(
    x = vs,
    y = n,
    fill = color
  )) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  geom_text(
    aes(
      label = abs(n), 
      vjust = vjust
    ), 
    size = 5, 
    fontface = 2, 
    show.legend = F
  ) +
  scale_x_discrete(
    limits = xlimits,
    labels = xlabels
  ) +
  geom_hline(yintercept = 0, size = 1) +
  scale_fill_brewer(
    palette = "Set1", 
    direction = -1,
    name = "Regulations",
    label = c("Down", "Up")
  ) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid =  element_blank(),
    
    axis.line = element_line(size = 1),
    # axis.line.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_text(size =14, color = "black", face = "bold"),
    
    legend.position = "top"
  ) +
  labs(
    y = "Number of DEGs"
  ) ->
  p;p
  
ggsave(
  filename = "Number-of-DEG-barplot.pdf",
  plot = p,
  device = "pdf",
  path = "data/uvresult/01-de",
  width = 7,
  height = 7
)




# Intersection ------------------------------------------------------------


se_group_de_volcano %>% 
  dplyr::select(vs, des_color) %>% 
  dplyr::mutate(n = purrr::map(
    .x = des_color,
    .f = function(.x) {
      .x %>% 
        dplyr::select(GeneName, color) %>% 
        dplyr::filter(color != "grey") %>% 
        dplyr::group_by(color) %>% 
        tidyr::nest() %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(data = purrr::map(.x = data, .f = function(.xx) {
          .xx$GeneName
        })) %>% 
        tidyr::spread(key = color, value = data)
    }
  )) %>% 
  dplyr::select(-des_color) %>% 
  tidyr::unnest(cols = n) %>% 
  dplyr::slice(1,2,4,5) ->
  se_group_de_volcano_red_green


up_inter = se_group_de_volcano_red_green$red %>% 
  purrr::reduce(.f = intersect) 

down_inter <- se_group_de_volcano_red_green$green %>% 
  purrr::reduce(.f = intersect) 


# ggplot() + 
#   annotate(
#     "text",
#     x = 1,
#     y = 1,
#     size = 6,
#     label = stringr::str_wrap(glue::glue("Up: {paste0(up_inter, collapse = ',')}\n\nDown:{paste0(down_inter, collapse = ',')}"), width = 10)
#   ) + 
#   theme_void()

ggplot() + 
  annotate(
  "text",
  x = 1,
  y = 1,
  size = 6,
  label = glue::glue("Up:{paste0(up_inter[c(1:8)], collapse = ',')}\n{paste0(up_inter[c(9:15)], collapse = ',')}\n{paste0(up_inter[c(16:20)], collapse = ',')}\n\nDown:{paste0(down_inter[c(1:8)], collapse = ',')}\n{paste0(down_inter[c(9:15)], collapse = ',')}\n{paste0(down_inter[c(16:21)], collapse = ',')}")
  ) + 
  theme_void() ->
  p3


ggvenn::ggvenn(
  data = list(
    "UVB0 vs. BC" = se_group_de_volcano_red_green$red[[1]],
    "UVB1 vs. BC" = se_group_de_volcano_red_green$red[[2]],
    "UVS0 vs. SC" = se_group_de_volcano_red_green$red[[3]],
    "UVS1 vs. SC" = se_group_de_volcano_red_green$red[[4]]
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_jco()(5),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 6
) ->
  p1



ggvenn::ggvenn(
  data = list(
    "UVB0 vs. BC" = se_group_de_volcano_red_green$green[[1]],
    "UVB1 vs. BC" = se_group_de_volcano_red_green$green[[2]],
    "UVS0 vs. SC" = se_group_de_volcano_red_green$green[[3]],
    "UVS1 vs. SC" = se_group_de_volcano_red_green$green[[4]]
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_lancet()(5),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 6
) ->
  p2

(p1 | p2) / p3 +
  plot_annotation(
    title = glue::glue("Up and down DEGs intersection"),  
    tag_levels = "A"
  ) ->
  venn_plot

  

ggsave(
  filename = "Number-of-DEG-venn.pdf",
  plot = venn_plot,
  device = "pdf",
  path = "data/uvresult/01-de",
  width = 12,
  height = 12
)


ggvenn::ggvenn(
  data = list(
    "UVB0/1 vs. BC" = union(
      se_group_de_volcano_red_green$red[[1]],
      se_group_de_volcano_red_green$red[[2]]
    ),
    "UVS0/1 vs. SC" = union(
      se_group_de_volcano_red_green$red[[3]],
      se_group_de_volcano_red_green$red[[4]]
    )
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_jco()(5),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 6
) ->
  p4

ggvenn::ggvenn(
  data = list(
    "UVB0/1 vs. BC" = union(
      se_group_de_volcano_red_green$green[[1]],
      se_group_de_volcano_red_green$green[[2]]
    ),
    "UVS0/1 vs. SC" = union(
      se_group_de_volcano_red_green$green[[3]],
      se_group_de_volcano_red_green$green[[4]]
    )
  ),
  show_percentage = FALSE,
  fill_color = ggsci::pal_lancet()(5),
  fill_alpha = 0.9,
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 6
) ->
  p5

(p4 | p5) +
  plot_annotation(
    title = glue::glue("Up and down DEGs intersection"),  
    tag_levels = "A"
  ) ->
  venn_plot2

ggsave(
  filename = "Number-of-DEG-venn2.pdf",
  plot = venn_plot2,
  device = "pdf",
  path = "data/uvresult/01-de",
  width = 8,
  height = 6
)


# Intersected genes enrichment --------------------------------------------

up_list <- list(
  "UVB0/1 vs. BC" = union(
    se_group_de_volcano_red_green$red[[1]],
    se_group_de_volcano_red_green$red[[2]]
  ),
  "UVS0/1 vs. SC" = union(
    se_group_de_volcano_red_green$red[[3]],
    se_group_de_volcano_red_green$red[[4]]
  )
)

up_inter <- intersect(
  up_list$`UVB0/1 vs. BC`,
  up_list$`UVS0/1 vs. SC`
)

up_inter_enrichment <- fn_enrichment(
  .gs = up_inter,
  .color = "red"
)

down_list <- list(
  "UVB0/1 vs. BC" = union(
    se_group_de_volcano_red_green$green[[1]],
    se_group_de_volcano_red_green$green[[2]]
  ),
  "UVS0/1 vs. SC" = union(
    se_group_de_volcano_red_green$green[[3]],
    se_group_de_volcano_red_green$green[[4]]
  )
)

down_inter <- intersect(
  down_list$`UVB0/1 vs. BC`,
  down_list$`UVS0/1 vs. SC`
)

down_inter_enrichment <- fn_enrichment(
  .gs = down_inter,
  .color = "green"
)

dplyr::bind_rows(
  up_inter_enrichment,
  down_inter_enrichment
) %>% 
  dplyr::mutate(reg = c("go_up", "go_down")) ->
  up_down_gobp


up_down_gobp %>% 
  purrr::pmap(
    .f = function(gobp, goplot, reg) {
      .filename <- glue::glue("Inter-GOBP-{reg}") %>%
        toupper()
      
      fn_save_enrichment_xlsx_pdf(
        .d = gobp,
        .p = goplot,
        .filename = .filename
      )
    }
  )


# save image --------------------------------------------------------------


save.image(file = "data/uvrda/14-uv-de.rda")
load(file = "data/uvrda/14-uv-de.rda")
