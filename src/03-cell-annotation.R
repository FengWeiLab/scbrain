# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri May 20 23:08:28 2022
# @DESCRIPTION: cell annotation

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(HGNChelper)
library(Seurat)
library(patchwork)

# Load data ---------------------------------------------------------------

sc_sct_cluster <- readr::read_rds(
  file = "data/rda/sc_sct_cluster.rds.gz"
)


# sc-Type -----------------------------------------------------------------

# load gene set preparation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list_immune_system = gene_sets_prepare(db_, "Immune system")
gs_list_brain = gene_sets_prepare(db_, "Brain")

gs_list_immune_system %>%
  purrr::map(.f = function(.x) {
    names(.x) <- paste(names(.x), "(Immune cell)")
    .x
  }) ->
  gs_list_immune_system
gs_list_brain %>%
  purrr::map(.f = function(.x) {
    names(.x) <- paste(names(.x), "(Brain)")
    .x
  }) ->
  gs_list_brain

gs_list <- list(
  gs_positive = c(
    gs_list_immune_system$gs_positive,
    gs_list_brain$gs_positive
  ),
  gs_negative = c(
    gs_list_immune_system$gs_negative,
    gs_list_brain$gs_negative
  )
)

es.max <- sctype_score(
  scRNAseqData = sc_sct_cluster[["SCT"]]@scale.data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)

cL_results <- do.call(
  "rbind",
  lapply(
    unique(sc_sct_cluster@meta.data$seurat_clusters),
    function(cl) {
      es.max.cl <- sort(
        rowSums(es.max[, rownames(sc_sct_cluster@meta.data[sc_sct_cluster@meta.data$seurat_clusters == cl, ])]),
        decreasing = TRUE
      )
      head(
        data.frame(
          cluster = cl,
          type = names(es.max.cl),
          scores = es.max.cl,
          ncells = sum(sc_sct_cluster@meta.data$seurat_clusters == cl)
        ),
        10
      )
    }
  )
)

cL_results %>%
  dplyr::group_by(cluster) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mm = purrr::map(
    .x = data,
    .f = function(.x) {
      .x %>%
        dplyr::mutate(
          cell = gsub(
            pattern = " \\(.*\\)",
            replacement = "",
            x = type
          )
        ) %>%
        dplyr::mutate(
          tissue = gsub(
            pattern = ".*\\(|\\)",
            replacement = "",
            x = type
          )
        ) %>%
        dplyr::arrange(tissue, -scores)
    }
  ))
  

sctype_scores <-  cL_results %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 1, wt = scores)  

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sc_sct_cluster@meta.data$customclassif <- ""
for(j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  sc_sct_cluster@meta.data$customclassif[sc_sct_cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

sc_sct_cluster$seurat_clusters <- sc_sct_cluster$integrated_snn_res.0.25

pcc <- readr::read_tsv(file = "https://chunjie-sam-liu.life/data/pcc.tsv")
sc_sct_cluster <- SetIdent(sc_sct_cluster, value = "seurat_clusters")

theme_umap <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(face = "bold"),
  # legend.position = "none"
)

fn_plot_umap <- function(.x) {
  # .x <- sc_sct_cluster
  .umap <- as.data.frame(.x@reductions$umap@cell.embeddings)
  .xx <- .x@meta.data[, c("seurat_clusters", "customclassif")] %>% 
    dplyr::rename(cluster = seurat_clusters, celltype =  customclassif)
  .xxx <- dplyr::bind_cols(.umap, .xx)
  .xxx %>% 
    dplyr::select(cluster, celltype) %>% 
    dplyr::group_by(cluster, celltype) %>%
    dplyr::count() %>% 
    dplyr::arrange(-n) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-n) %>% 
    dplyr::mutate(celltype = glue::glue("{cluster} {celltype}")) ->
    .xxx_celltype
  .xxx %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>% 
    dplyr::left_join(.xxx_celltype, by = "cluster") ->
    .xxx_label
  
  ggplot(data = .xxx) +
    geom_point(
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = cluster,
        shape = NULL,
        alpha = NULL
      )
    ) +
    geom_text(
      data = .xxx_label,
      aes(
        label = cluster,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 6
    ) +
    scale_colour_manual(
      name = NULL,
      values = pcc$color,
      labels = .xxx_label$celltype,
      guide = guide_legend(
        ncol = 1, 
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 14)
    ) +
    labs(
      x = "UMAP1",
      y = "UMAP2"
    )
}

p <- fn_plot_umap(.x = sc_sct_cluster)


ggsave(
  filename = "cluster-plot-umap.pdf",
  plot = p,
  device = "pdf",
  path = "data/result/",
  width = 20,
  height = 10
)

Seurat::DimPlot(
  sc_sct_cluster,
  reduction = "umap",
  label = TRUE,
  label.size = 6,
  cols = pcc$color,
  split.by = "region"
) +
  labs(x = "UAMP1", y = "UMAP2") +
  theme_umap ->
  p1;p1

ggsave(
  filename = "cluster-plot-umap-region.pdf",
  plot = p1,
  device = "pdf",
  path = "data/result/",
  width = 20,
  height = 8
)

Seurat::DimPlot(
  sc_sct_cluster,
  reduction = "umap",
  label = TRUE,
  label.size = 6,
  cols = pcc$color,
  split.by = "tissue",
  ncol = 3
) +
  labs(x = "UAMP1", y = "UMAP2") +
  theme_umap ->
  p2;p2

ggsave(
  filename = "cluster-plot-umap-tisue.pdf",
  plot = p2,
  device = "pdf",
  path = "data/result/",
  width = 20,
  height = 15
)


# Seurat::DimPlot(
#   sc_sct_cluster,
#   reduction = "umap",
#   label = TRUE,
#   label.size = 4,
#   cols = pcc$color,
#   group.by = 'customclassif',
# ) +
#   theme_umap +
#   labs(x = "UAMP1", y = "UMAP2")  ->
#   p2;p2

# Seurat::DimPlot(
#   sc_sct_cluster,
#   reduction = "umap",
#   label = TRUE,
#   # repel = TRUE,
#   group.by = 'customclassif',
#   split.by = "tissue",
#   ncol = 3
# ) +
#   theme(
#     legend.position = "bottom"
#   ) ->
#   p

ggsave(
  filename = "test-annotation-plot-map.pdf",
  plot = p,
  device = "pdf",
  path = "data/result/",
  width = 25,
  height = 20
)

# scMCA -------------------------------------------------------------------

mca <- scMCA::scMCA(
  scdata = exp(sc_sct_cluster[["integrated"]]@scale.data),
  numbers_plot = 3
)

mca$cors_matrix[1:4,1:3]
mca$scMCA %>%
  tibble::enframe() ->
  dd
dd %>%
  dplyr::group_by(value) %>%
  dplyr::count()


ref.expr %>%
  colnames() %>%
  tibble::enframe() %>%
  dplyr::mutate(cell = gsub(pattern = "\\(.*\\)", replacement = "", x = value)) %>%
  dplyr::mutate(tissue = gsub(pattern = ".*\\(|\\)", replacement = "", x = value)) ->
  cells

cells %>%
  dplyr::filter(grepl(pattern = "Brain", x = tissue)) %>%
  print(n = Inf)



# save image --------------------------------------------------------------

save.image(file = "data/rda/03-cell-annotation.rda")
