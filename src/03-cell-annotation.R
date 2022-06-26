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
library(SingleR)

# Load data ---------------------------------------------------------------

sc_sct_cluster <- readr::read_rds(
  file = "data/rda/sc_sct_cluster.rds.gz"
)


sc_sct_cluster$seurat_clusters <- sc_sct_cluster$integrated_snn_res.0.2
sc_sct_cluster <- SetIdent(sc_sct_cluster, value = "seurat_clusters")



# sc-Type -----------------------------------------------------------------

# load gene set preparation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/sctype_score_.R")
# DB file
db_full = "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_full.xlsx";
db_short <- "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_short.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
gs_list_immune_system = gene_sets_prepare(db_full, "Immune system")
gs_list_brain = gene_sets_prepare(db_full, "Brain")

gs_list_immune_system %>%
  purrr::map(.f = function(.x) {
    names(.x) <- paste(names(.x), "(Immune cell)")
    .x$`Cancer cells (Immune cell)` <- NULL
    .x
  }) ->
  gs_list_immune_system_mod

gs_list_brain %>%
  purrr::map(.f = function(.x) {
    names(.x) <- paste(names(.x), "(Brain)")
    .x$`Cancer cells (Brain)` <- NULL
    .x
  }) ->
  gs_list_brain_mod

gs_list <- list(
  gs_positive = c(
    gs_list_immune_system_mod$gs_positive,
    gs_list_brain_mod$gs_positive
  ),
  gs_negative = c(
    gs_list_immune_system_mod$gs_negative,
    gs_list_brain_mod$gs_negative
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

cL_results %>% 
  dplyr::filter(cluster == 3)

sctype_scores <-  cL_results %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 1, wt = scores)  

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sc_sct_cluster@meta.data$sctype <- ""
for(j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  sc_sct_cluster@meta.data$sctype[sc_sct_cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}



# SingleR -----------------------------------------------------------------

# sc_sct_cluster.sce <- as.SingleCellExperiment(sc_sct_cluster)
# 
# mouse.se <- celldex::MouseRNAseqData()
# 
# sc_sct_cluster.sce.pred <- SingleR::SingleR(
#   test = sc_sct_cluster.sce,
#   ref = mouse.se,
#   labels = mouse.se$label.main
# )
# 
# 
# sc_sct_cluster$singler <- sc_sct_cluster.sce.pred$labels


# scMCA -------------------------------------------------------------------

# mca <- scMCA::scMCA(
#   scdata = exp(sc_sct_cluster[["integrated"]]@scale.data),
#   numbers_plot = 3
# )
# 
# ref.expr %>%
#   colnames() %>%
#   tibble::enframe() %>%
#   dplyr::mutate(cell = gsub(pattern = "\\(.*\\)", replacement = "", x = value)) %>%
#   dplyr::mutate(tissue = gsub(pattern = ".*\\(|\\)", replacement = "", x = value)) ->
#   cells
# 
# cells %>%
#   dplyr::filter(grepl(pattern = "Brain|Bone-Marrow$|Peripheral_Blood", x = tissue)) ->
#   cells_cand
# 
# cells_cand %>% 
#   dplyr::filter(grepl(pattern = "Endo", x = cell))
# 
# sc_sct_cluster@meta.data %>% head()
# mca$cors_matrix[, 1] %>% sort(decreasing = T) %>% head()
# mca$scMCA[1:4]
# 
# sc_sct_cluster$scmca <- mca$scMCA

# Umap plot ---------------------------------------------------------------


pcc <- readr::read_tsv(file = "https://chunjie-sam-liu.life/data/pcc.tsv")

theme_umap <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(face = "bold"),
  # legend.position = "none"
)

fn_plot_umap <- function(.x, .celltype="sctype", .reduction="umap") {
  # .x <- sc_sct_cluster
  # .celltype="sctype"
  # .reduction="tsne"
  
  .umap <- as.data.frame(.x@reductions[[.reduction]]@cell.embeddings)
  colnames(.umap) <- c("UMAP_1", "UMAP_2")
  .xx <- .x@meta.data[, c("seurat_clusters", .celltype)] %>% 
    dplyr::rename(cluster = seurat_clusters, celltype =  .celltype)
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
    tidyr::nest() %>%
    dplyr::mutate(u = purrr::map(.x = data, .f = function(.m) {
      # d %>% 
      #   dplyr::filter(cluster == 14) %>% 
      #   dplyr::pull(data) %>% 
      #   .[[1]] ->
      #   .m

      .m %>%
        dplyr::summarise(u1 = mean(UMAP_1), u2 = mean(UMAP_2)) ->
        .mm

      .m %>%
        dplyr::mutate(u1 = UMAP_1 > .mm$u1, u2 = UMAP_2 > .mm$u2) ->
        .mmd
      
      .mmd %>%
        dplyr::group_by(u1, u2) %>%
        dplyr::count() %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-n) ->
        .mmm
      
      .fc <- .mmm$n[[1]] / .mmm$n[[2]] # 1.1
      .mmm
      .fc
      
      if(.fc > 1.1) {
        .mmd %>%
          dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      } else {
        .mmd %>%
          # dplyr::filter(u1 == .mmm$u1[[1]], u2 == .mmm$u2[[1]]) %>%
          dplyr::summarise(UMAP_1  = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
      }

    })) %>% 
    dplyr::ungroup() %>% 
    tidyr::unnest(cols = u) %>% 
    dplyr::select(-data) %>% 
    dplyr::left_join(.xxx_celltype, by = "cluster") %>% 
    dplyr::arrange(cluster) ->
    .xxx_label
  
  # .xxx %>% 
  #   dplyr::group_by(cluster) %>% 
  #   dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>% 
  #   dplyr::left_join(.xxx_celltype, by = "cluster")
  #   
    
  
  .labs <- if (.reduction == "umap") {
    labs(
      x = "UMAP1",
      y = "UMAP2"
    )
  } else {
    labs(
      x = "tSNE1",
      y = "tSNE2"
    )
  }
  
  ggplot() +
    geom_point(
      data = .xxx,
      aes(
        x = UMAP_1,
        y = UMAP_2,
        colour = cluster,
        shape = NULL,
        alpha = NULL
      ),
      size = 0.7
    ) +
    geom_text(
      data = .xxx_label,
      aes(
        label = cluster,
        x = UMAP_1,
        y = UMAP_2,
      ),
      size = 7
    ) +
    scale_colour_manual(
      name = NULL,
      values = pcc$color,
      labels = .xxx_label$celltype,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(size=4)
      )
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(
        colour = "black",
        size = 0.5,
        arrow = grid::arrow(
          angle = 5,
          length = unit(5, "npc"),
          type = "closed"
        )
      ),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(
        size = 12, 
        face = "bold", 
        hjust = 0.05
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(
        face = "bold",
        color = "black",
        size = 12
      )
    ) +
    coord_fixed(
      ratio = 1,
    ) +
    .labs
}

p <- fn_plot_umap(.x = sc_sct_cluster, .celltype = "sctype", .reduction = "tsne")
p


Seurat::DimPlot(
  sc_sct_cluster,
  reduction = "tsne",
  label = TRUE,
  label.size = 6,
  cols = pcc$color,
)
# 
# ggsave(
#   filename = "cluster-plot-umap.pdf",
#   plot = p,
#   device = "pdf",
#   path = "data/result/",
#   width = 20,
#   height = 10
# )




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

# ggsave(
#   filename = "cluster-plot-umap-region.pdf",
#   plot = p1,
#   device = "pdf",
#   path = "data/result/",
#   width = 20,
#   height = 8
# )

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

# ggsave(
#   filename = "cluster-plot-umap-tisue.pdf",
#   plot = p2,
#   device = "pdf",
#   path = "data/result/",
#   width = 20,
#   height = 15
# )


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
# 
# ggsave(
#   filename = "test-annotation-plot-map.pdf",
#   plot = p,
#   device = "pdf",
#   path = "data/result/",
#   width = 25,
#   height = 20
# )



# Marker genes ------------------------------------------------------------

DefaultAssay(sc_sct_cluster) <- "integrated"
sc_sct_cluster <- PrepSCTFindMarkers(object = sc_sct_cluster)


all.markers <- FindAllMarkers(
  object = sc_sct_cluster,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


readr::write_rds(
  x = all.markers,
  file = "data/rda/sc_sct_cluster_marker_genes.rds.gz"
)


# Gene plot ---------------------------------------------------------------
all.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_max(n = 5, order_by = avg_log2FC) %>% 
  print(n = Inf)

DefaultAssay(sc_sct_cluster) <- "RNA"

VlnPlot(
  object = sc_sct_cluster,
  features = c("Vpreb2"),
)

FeaturePlot(
  object = sc_sct_cluster,
  features = "Vpreb2",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "tsne"
) 

# Features ----------------------------------------------------------------

all_marker_genes <- c("Cd3e", "Cd8a", "Nkg7", "Hexb", "Cx3cr1", "Ccl3", "Ccl4", "Il1b", "Lyz2", "Clec10a")

# Lymphoid lineage clusters
Bcells <- c("Cd79a")
Mature_Bcells <- c(Bcells, "Ms4a1")
CD4T <- c("Cd3e")
CD8T <- c("Cd3e", "Cd8a")
NKcells <- c("Nkg7", "Gzmk")

# myeloid lineage clusters
myeloid <- c("Lyz2", "Cd68")
microglia <- c("Hexb", "Cx3cr1")
cns_macrophage <- c("C1qc", "Fcna", "Pf4", "Cd14", "Cd163")
macrophages <- c("Cd68", "RT1-Db1")
monocytes <- c("Fcnb", "Fcgr3a")
granulocytes <- c("S100a9", "Hp")
myeloid_dendritc_cells <- c("Clec9a", "Xcr1")
plasmacytoid_dendritic_cells <- c("Siglech", "Runx2")

SMC_vascular_smooth_muscle_cells <- c("Acta2")
FB_perivascular_fibroblast_like <- c("Dcn")
CAM_cns_associated_macrophages <- ("Pf4")
MdC_monocyte_derived <- c("Ccr2")
vEC_venous_endothelial <- c("Itm2a")
capEC_capillary_endothelial <- c("Itm2a")
aEC_endothelial_cells <- c("Itm2a")
PC_pericytes <- c("Kcnj8")
CPC_choroid_plexus_capillary_endothelia <- c("Plvap")
EPC_ependymocytes <- c("Ttr", "Tmem212")
MG_microglia <- c("Hexb")
NEUT_neutrophils <- c("S100a8")
ASC_astrocytes <- c("Aldoc")
OLG_oligodendrocytes <- c("Plp1")
NPC_neural_progenitor <- c()
LYM_lymphocytes <- c()
RBC_red_blood_cell <- c()

cL_results %>% 
  dplyr::filter(cluster == 18)

DefaultAssay(sc_sct_cluster) <- "SCT"

FeaturePlot(
  object = sc_sct_cluster,
  features = "Gzmk",
  cols = c("lightgrey", "#CD0000"),
  order = TRUE,
  reduction = "tsne"
) 

# PAX6,HES5,GFAP,SLC1A3

# GFAP,SLC1A3,SLC1A2,GLUL,S100B,ALDH1L1,AQP4,IGFBP3,ATP13A4,CBS,SOX9,CD40,CD80,CD86,C5AR1



# save image --------------------------------------------------------------

save.image(file = "data/rda/03-cell-annotation.rda")
load(file = "data/rda/03-cell-annotation.rda")
