library(ComplexHeatmap)

d <- readxl::read_excel(
  "/scr1/users/liuc9/tmp/CD47 Heatmap.xlsx",
  sheet = 1
)

d |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "Gene") |>
  as.matrix() ->
  d.mat

d.mat |>
  apply(1, scale) |>
  t() ->
  d.mat.scale

colnames(d.mat.scale) <- colnames(d.mat)

cluster_col <- ggsci::pal_aaas()(3)[c(1, 3, 2)]
names(cluster_col) <- c("SC/LNP", "Unmodified/LNP", "CD47/LNP")

hma_top <- ComplexHeatmap::HeatmapAnnotation(
  df = data.frame(
    Group = c("SC/LNP", "Unmodified/LNP", "CD47/LNP")
  ),
  gap = unit(c(2, 2), "mm"),
  col = list(
    Group = cluster_col
  ),
  which = "column"
)

ComplexHeatmap::Heatmap(
  matrix = d.mat.scale,
  col =  circlize::colorRamp2(
    breaks = c(-1.1, 0, 1.1),
    colors = c("blue", "white", "red"),
    space = "RGB"
  ),
  name = "Normalized",
  na_col = "grey",
  color_space = "LAB",
  rect_gp = gpar(col = NA),
  border = NA,
  cell_fun = NULL,
  layer_fun = NULL,
  jitter = FALSE,
  top_annotation = hma_top,

  # column
  cluster_columns = F,

  # row
  cluster_rows = T,
  cluster_row_slices = T,
  row_dend_side = "left",
  show_row_dend = F,
  row_dend_reorder = T
) ->
  p;p

{
  pdf(
    file = file.path(
      "/home/liuc9/tmp/CD47_Heatmap.pdf"
    ),
    width = 5,
    height = 7
  )
  ComplexHeatmap::draw(object = p)
  dev.off()
}

d.mat.scale |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "Gene") |>
  writexl::write_xlsx(
    "/scr1/users/liuc9/tmp/CD47_Heatmap_scale.xlsx"
  )
