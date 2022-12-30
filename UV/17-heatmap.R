# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Dec 27 23:26:46 2022
# @DESCRIPTION: heatmap

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(DESeq2)
library(ComplexHeatmap)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_heatmap_se <- function(.se) {
  # .se <- se_brain_deg
  
  
  .matrix <- assay(.se)
  .meta <- colData(.se)
  
  .matrix_scale <- .matrix %>% 
    apply(1, scale) %>%
    t()
  
  colnames(.matrix_scale) <- .meta$barcode
  
  colors %>% 
    dplyr::filter(group %in% .meta$group) %>% 
    tibble::deframe() ->
    .cluster_col
  
  hma_top = ComplexHeatmap::HeatmapAnnotation(
    df = as.data.frame(.meta) %>% 
      dplyr::select(Group = group),
    gap = unit(c(2,2), "mm"),
    col = list(Group = .cluster_col),
    which = "column"
  )
  
  ComplexHeatmap::Heatmap(
    # data and color
    matrix = .matrix_scale,
    col =  circlize::colorRamp2(
      breaks = c(-1.1, 0, 1.1), 
      colors = c("blue", "white", "red"), 
      space = "RGB"
    ),
    name = "Normalized counts",
    na_col = 'grey', 
    color_space = 'LAB', 
    rect_gp = gpar(col = NA),
    border = NA, 
    cell_fun = NULL, 
    layer_fun = NULL,
    jitter = FALSE,
    
    # title
    # row_title = 'Selected genes', # OC44
    row_title = 'Gene name',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    column_title = 'Samples',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 10),
    column_title_rot = 0,
    
    # clustering of row
    cluster_rows = T,
    cluster_row_slices = T,
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    row_dend_side = 'left',
    row_dend_width = unit(10, 'mm'),
    show_row_dend = F,
    row_dend_reorder = T,
    row_dend_gp = gpar(),
    
    # clustering of column
    cluster_columns = F,
    # cluster_column_slices = T,
    # clustering_distance_columns = "pearson",
    # clustering_method_columns = "ward.D",
    # column_dend_side = 'top',
    # column_dend_height = unit(10, 'mm'),
    # show_column_dend = F, # OC521
    # column_dend_gp = gpar(),
    # column_dend_reorder = F,
    
    row_order = NULL,
    # column_order = c(7,8,9,4,5,6,1,2,3),
    
    # row labels
    row_labels = rownames(.matrix_scale),
    row_names_side = 'right',
    show_row_names = F,
    row_names_max_width = unit(6, 'cm'),
    row_names_gp = gpar(fontsize = 3),
    row_names_rot = 0,
    row_names_centered = FALSE,
    
    # column labels
    # column_labels = colnames(.matrix_scale),
    column_names_side = 'bottom',
    show_column_names = T,
    column_names_max_height = unit(6, 'cm'),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_centered = FALSE,
    
    # annotation
    top_annotation = hma_top,
    bottom_annotation = NULL,
    # left_annotation = hma_right,
    # right_annotation = hma_left,
    left_annotation = NULL,
    right_annotation = NULL,
    
    
    # kmeans cluster number
    # row cluster is 1
    # column cluster is 2 with 10 repeats
    # km = 1,
    split = NULL,
    row_km = 2,
    row_km_repeats = 5,
    row_split = NULL,
    column_km = 2,
    column_km_repeats = 2,
    # column_split = rep(c("A", "B", "C", "D", "E"), each = 3), ###AAAAA
    # gap = unit(1, 'mm'),
    row_gap = unit(1, 'mm'),
    column_gap = unit(1, 'mm'),
    
    show_heatmap_legend = T,
    heatmap_legend_param = list(title = 'Z-score'),
    
    # raster_device = 'tiff',
    raster_quality = 2,
    raster_device_param = list(),
    
    raster_resize = F,
    
    post_fun = NULL
  ) -> 
    .heatmap
  
  .row_order <- ComplexHeatmap::row_order(.heatmap)
  
  .row_order %>% 
    tibble::enframe() %>% 
    dplyr::mutate(
      p = purrr::map2(
        .x = name,
        .y = value,
        .f = function(.rn, .rr) {
          # .rn <- .u$name[[1]]
          # .rr <- .u$value[[1]]
          
          .rr_genename <- rownames(.matrix_scale)[.rr]
          
          .go_bp <- clusterProfiler::enrichGO(
            gene = .rr_genename,
            universe = rownames(se),
            keyType = "SYMBOL",
            OrgDb = org.Mm.eg.db::org.Mm.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
          )
          
          .go_bp %>%
            tibble::as_tibble()  %>% 
            dplyr::mutate(
              Description = stringr::str_wrap(
                stringr::str_to_sentence(string = Description), width = 60
              )
            ) %>%
            dplyr::mutate(adjp = -log10(p.adjust)) %>%
            dplyr::select(ID, Description, adjp, Count) %>%
            head(20) %>%
            dplyr::arrange(adjp, Count) %>%
            dplyr::mutate(Description = factor(Description, levels = Description)) ->
            .go_bp_for_plot
          
          .go_bp_for_plot %>%
            ggplot(aes(x = Description, y = adjp)) +
            geom_col(fill = "#AE1700", color = NA, width = 0.7) +
            geom_text(aes(label = Count), hjust = 4, color = "white", size = 5) +
            labs(
              title = glue::glue("{.type} {.rn}; n=({length(.rr_genename)})"),
              y = "-log10(Adj. P value)"
            ) +
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
          
          .filename_prefix <- glue::glue("heatmap_{.rn}_n=({length(.rr_genename)})")
          
          writexl::write_xlsx(
            x = as.data.frame(.go_bp), 
            path = file.path(
              "data/uvresult/04-heatmap",
              glue::glue("{.filename_prefix}.xlsx")
            )
          )
          ggsave(
            plot = .go_bp_plot,
            filename = glue::glue("{.filename_prefix}.pdf"),
            device = "pdf",
            path = "data/uvresult/04-heatmap",
            width = 10, 
            height = 6.5
          )
          
        }
      )
    )
  
  
  .heatmap
}

# load data ---------------------------------------------------------------
se_group_de_enrichment <- readr::read_rds(
  file = "data/uvrda/se_group_de_enrichment.rds.gz"
)

se <- readr::read_rds(
  file = file.path(
    "data/uvrda",
    "UV_count_matrix.se.rds.gz"
  )
)

colors <- tibble::tibble(
  group = unique(se$group),
  color = c("black", RColorBrewer::brewer.pal(4, name = "Paired")[c(1,2)], "grey", RColorBrewer::brewer.pal(4, name = "Paired")[c(3,4)] )
)

# body --------------------------------------------------------------------

se_group_de_enrichment %>% 
  dplyr::filter(seq == "Skull") %>% 
  dplyr::mutate(
    deg = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color != "grey") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_up = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "red") %>% 
          dplyr::pull(GeneName)
      }
    ),
    deg_down = purrr::map(
      .x = des_color,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(color == "green") %>% 
          dplyr::pull(GeneName)
      }
    )
  ) ->
  se_group_de_enrichment_skull



{
  pdf(
    file = "data/uvresult/04-heatmap/Heatmap.pdf", 
    width = 8, 
    height = 10
  )
  ComplexHeatmap::draw(
    object = fn_heatmap_se(
      .se = se[se_group_de_enrichment_skull$deg %>% purrr::reduce(.f = union),]
    )
  )
  dev.off()
}

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/uvrda/14-heatmap.rda")
