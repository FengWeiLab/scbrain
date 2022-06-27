# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Jun 26 16:38:21 2022
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)
library(HGNChelper)


# Load data ---------------------------------------------------------------
pcc <- readr::read_tsv(file = "https://chunjie-sam-liu.life/data/pcc.tsv")


project_sct <- readr::read_rds(file = "data/rda/project_sc_sct.rds.gz")

regions <- c(
  "CB" = "Brain", 
  "CM" = "Meninge", 
  "CS" = "Skull", 
  "DB" = "Brain", 
  "DM" = "Meninge", 
  "DS" = "Skull"
)
cases <- c(
  "CB" = "Sham", 
  "CM" = "Sham", 
  "CS" = "Sham", 
  "DB" = "MCAO", 
  "DM" = "MCAO", 
  "DS" = "MCAO"
)

project_sct %>% 
  dplyr::mutate(
    region = plyr::revalue(x = project, replace = regions),
    case = plyr::revalue(x = project, replace = cases)
  ) ->
  project_sct_regions


# Split by region ---------------------------------------------------------



project_sct_regions %>% 
  dplyr::filter(region == "Brain") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  brain_sct_list


project_sct_regions %>% 
  dplyr::filter(region == "Meninge") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  meninge_sct_list


project_sct_regions %>% 
  dplyr::filter(region == "Skull") %>% 
  dplyr::select(project, sct) %>% 
  tibble::deframe() ->
  skull_sct_list


# Function ----------------------------------------------------------------

fn_sct_integrate <- function(.sct_list) {
  # .sct_list <- brain_sct_list
  
  .features <- Seurat::SelectIntegrationFeatures(
    object.list = .sct_list,
    nfeatures = 3000
  )
  .sct_list <- Seurat::PrepSCTIntegration(
    object.list = .sct_list,
    anchor.features = .features
  )
  .anchors <- Seurat::FindIntegrationAnchors(
    object.list = .sct_list,
    normalization.method = "SCT",
    anchor.features = .features
  )
  .sct_list <- Seurat::IntegrateData(
    anchorset = .anchors,
    normalization.method = "SCT"
  )
  
  .sct_list
}

fn_cluster <- function(.sct) {
  .sct %>%
    Seurat::RunPCA(npcs = 30) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:30) %>%
    Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>%
    Seurat::FindClusters(resolution = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) ->
    .sct_cluster
  .sct_cluster <- Seurat::PrepSCTFindMarkers(object = .sct_cluster)
  .sct_cluster
}

fn_gs_list <- function() {
  source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/sctype_score_.R")
  
  db_full = "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_full.xlsx";
  db_short <- "https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/ScTypeDB_short.xlsx"
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
  
  gs_list
}

fn_sctype <- function(.sct_cluster, .res = 0.3, gs_list) {
  .snn_res <- glue::glue("integrated_snn_res.{.res}")
  
  .sct_cluster$seurat_clusters <- .sct_cluster[[.snn_res]]
  .sct_cluster <- SetIdent(.sct_cluster, value = "seurat_clusters")
  
  es.max <- sctype_score(
    scRNAseqData = .sct_cluster[["SCT"]]@scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  
  cL_results <- do.call(
    "rbind",
    lapply(
      unique(.sct_cluster@meta.data$seurat_clusters),
      function(cl) {
        es.max.cl <- sort(
          rowSums(es.max[, rownames(.sct_cluster@meta.data[.sct_cluster@meta.data$seurat_clusters == cl, ])]),
          decreasing = TRUE
        )
        head(
          data.frame(
            cluster = cl,
            type = names(es.max.cl),
            scores = es.max.cl,
            ncells = sum(.sct_cluster@meta.data$seurat_clusters == cl)
          ),
          10
        )
      }
    )
  )
  
  sctype_scores <-  cL_results %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n = 1, wt = scores)  
  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  .sct_cluster@meta.data$sctype <- ""
  for(j in unique(sctype_scores$cluster)) {
    cl_type = sctype_scores[sctype_scores$cluster == j, ]
    .sct_cluster@meta.data$sctype[.sct_cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  .sct_cluster
}

fn_plot_umap_tsne <- function(.x, .celltype="sctype", .reduction="umap") {
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
      size = 6
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

fn_find_all_markers <- function(.sct_cluster, .res = 0.3) {
  .snn_res <- glue::glue("integrated_snn_res.{.res}")
  .sct_cluster$seurat_clusters <- .sct_cluster[[.snn_res]]
  .sct_cluster <- SetIdent(.sct_cluster, value = "seurat_clusters")
  
  Seurat::FindAllMarkers(
    object = .sct_cluster,
    assay = "SCT",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  

}

fn_gene_dotplot <- function(.sct_cluster, .marker, .n = 3) {
  
  .marker %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::slice_max(n = .n, order_by = avg_log2FC) %>% 
    print(n = Inf) ->
    .marker_head
  
  DefaultAssay(.sct_cluster) <- "SCT"
  
  DotPlot(
    .sct_cluster, 
    features = unique(.marker_head$gene), 
    cols = c("blue", "red"), 
    dot.scale = 8
  ) +
    RotatedAxis()
}


fn_qc <- function(.region, .sct) {
  # .sct
  Seurat::VlnPlot(
    object = .sct, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
    cols = viridis::viridis_pal()(2),
    ncol = 2,
    group.by = "case",
    pt.size = 0
  )  ->
    .vlnplot #;.vlnplot
  
  .plot1 <- Seurat::FeatureScatter(
        object = .sct,
        feature1 = "nCount_RNA",
        feature2 = "percent.mt",
        group.by = "case",
        cols = viridis::viridis_pal()(2),
      )
  .plot2 <- Seurat::FeatureScatter(
        object = .sct,
        feature1 = "nCount_RNA",
        feature2 = "percent.ribo",
        group.by = "case",
        cols = viridis::viridis_pal()(2),
      )
  .plot3 <- Seurat::FeatureScatter(
        object = .sct,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
        group.by = "case",
        cols = viridis::viridis_pal()(2),
      )
  # .plot1 + .plot2 + .plot3
  # # .plot1 + .plot2 - .plot3
  # .plot1 | .plot2 | .plot3
  # .plot1 | .plot2 / .plot3
  # (.plot1 | .plot2) / .plot3
  # (.plot1 + (.plot2 + .plot3) + .plot3  +
  #   plot_layout(ncol = 1)) *
  #   theme_bw()
  # (.plot1 + (.plot2 + .plot3) + .plot3  +
  #     plot_layout(ncol = 1)) &
  #   theme_bw()

  
  .vlnplot / .plot3 / (.plot1 | .plot2)  +
    plot_annotation(
      title = glue::glue("Quality control {.region}"),  
      tag_levels = "A"
    ) +
    plot_layout(guides='collect') &
    theme(legend.position='bottom') -> .qc; .qc
  
  ggsave(
    filename = glue::glue("quality_control_{.region}.pdf"),
    plot = .qc,
    device = "pdf",
    path = "data/result/01-qc",
    width = 9,
    height = 15
    
  )

}


#
# Cluster -----------------------------------------------------------------


brain_meninge_skull_sct_list <- tibble::tibble(
  region = c("brain", "meninge", "skull"),
  sct_list = list(brain_sct_list, meninge_sct_list, skull_sct_list)
)


brain_meninge_skull_sct_list %>% 
  dplyr::mutate(
    sct = purrr::map2(
      .x = region,
      .y = sct_list,
      .f = function(.x, .y) {
        fn_sct_integrate(.sct_list = .y)
      }
    )
  ) %>% 
  dplyr::mutate(
    sct_cluster = purrr::map2(
      .x = region,
      .y = sct,
      .f = function(.x, .y) {
        fn_cluster(.sct = .y)
      }
    )
  ) ->
  brain_meninge_skull_sct_cluster

readr::write_rds(
  x = brain_meninge_skull_sct_cluster,
  file = "data/rda/brain_meninge_skull_sct_cluster.rds.gz"
)






# Metrics -----------------------------------------------------------------


brain_meninge_skull_sct_cluster %>% 
  dplyr::mutate(
    qc = purrr::map2(
      .x = region,
      .y = sct,
      .f = fn_qc
    )
  )
  

#
# Annotation --------------------------------------------------------------


gs_list <- fn_gs_list()

brain_meninge_skull_sct_cluster %>% 
  dplyr::mutate(
    sct_cluster_sctype = purrr::map(
      .x = sct_cluster,
      .f = fn_sctype,
      .res = 0.1,
      gs_list = gs_list
    )
  ) %>% 
  dplyr::mutate(
    sct_cluster_sctype_umap = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "umap"
    )
  ) %>% 
  dplyr::mutate(
    sct_cluster_sctype_tsne = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "tsne"
    )
  ) ->
  brain_meninge_skull_sct_cluster_sctype


# Marker genes ------------------------------------------------------------

brain_meninge_skull_sct_cluster_sctype %>% 
  dplyr::mutate(
    marker_genes = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_find_all_markers,
      .res = 0.1
    )
  ) %>% 
  dplyr::mutate(
    gene_dotplot = purrr::map2(
      .x = sct_cluster_sctype,
      .y = marker_genes,
      .f = fn_gene_dotplot,
      .n = 2
    )
  ) ->
  brain_meninge_skull_sct_cluster_sctype_marker

readr::write_rds(
  x = brain_meninge_skull_sct_cluster_sctype_marker,
  file = "data/rda/brain_meninge_skull_sct_cluster_sctype_marker.rds.gz"
)


# Manual marker gene ------------------------------------------------------
list(
  `0` = list(
    markers = c("Cldn5", "Itm2a"),
    fullname = "endothelial cells",
    shortname = "EC"
  ),
  `1` = list(
    markers = c("Hexb", "Ctss"),
    fullname = "microglia",
    shortname = "Microglia"
  ),
  
  `2` = list(
    markers = c("Acat2", "Myl9"),
    fullname = "vascular smooth muscle cells",
    shortname = "SMC"
  ),
  `3` = list(
    markers = c("Aldoc", "Clu"),
    fullname = "astrocytes",
    shortname = "ASC"
  ),
  `4` = list(
    markers = c("Cd74", "Lyz2"),
    fullname = "Monocyte",
    shortname = "Lyz2+ Mono"
  ),
  `5` = list(
    markers = c("Ttr", "Enpp2"),
    fullname = "Ependymal cell",
    shortname = "Ttr+ Ependymal"
  ),
  `6` = list(
    markers = c("Apoe", "Pf4"),
    fullname = "central nervous system (CNS)-associated macrophages",
    shortname = "Macrophages"
  ),
  `7` = list(
    markers = c("Plp1", "Ptgds"),
    fullname = "Oligodendrocyte",
    shortname = "Oligodendrocyte"
  ),
  `8` = list(
    markers = c("mt-Atp8", "mt-Co3"),
    fullname = "remove",
    shortname = "remove"
  ),
  `9` = list(
    markers = c("S100a8", "S100a9"),
    fullname = "Neutrophil",
    shortname = "Neutrophil"
  ),
  
  `10` = list(
    markers = c("Dcn", "Apod"),
    fullname = "perivascular fibroblast-like cells",
    shortname = "FB"
  ),
  `11` = list(
    markers = c("Plvap", "Plpp1"),
    fullname = "choroid plexus capillary endothelial cells",
    shortname = "CPC"
  ),
  `12` = list(
    markers = c("Rarres2", "Tmem212"),
    fullname = "Ciliated ependymal cell",
    shortname = "CEC"
  )
) %>% 
  tibble::enframe() %>% 
  dplyr::mutate(
    a = purrr::map(
      .x = value,
      .f = function(.x) {
        .x$markers <- paste0(.x$markers, collapse = ",")
        .x %>% 
          tibble::enframe() %>%  
          tidyr::unnest(cols = value) %>% 
          tidyr::spread(key = name, value = value)
      }
    )
  ) %>% 
  dplyr::select(-value) %>% 
  tidyr::unnest(cols = a) %>% 
  dplyr::rename(cluster = name) ->
  brain_marker_celltype 

list(
  `0` = list(
    markers = c("Apoe", "C1qb", "C1qa", "Ctss"),
    fullname = "Macrophage",
    shortname = "Macrophage"
  ),
  `1` = list(
    markers = c("Retnlg", "Il1b", "Ifitm1"),
    fullname = "Neutrophil",
    shortname = "Neutrophil"
  ),
  `2` = list(
    markers = c("Igkc", "Cd79a", "Ly6d", "Ighm", "Ms4a1"),
    fullname = "mature B Cell",
    shortname = "Mature B Cell"
  ),
  `3` = list(
    markers = c("Ccl5", "Ms4a4b", "Cd3d", "Ighm", "Ms4a1"),
    fullname = "CD8+ T cells",
    shortname = "CD8+ T cells"
  )
)


{
  b <- brain_meninge_skull_sct_cluster_sctype_marker$sct_cluster_sctype[[2]]
  d <- brain_meninge_skull_sct_cluster_sctype_marker$marker_genes[[2]] 
  
  DefaultAssay(b) <- "SCT"
  
  d %>% 
    dplyr::filter(!grepl("^mt-", gene)) %>% 
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 6, order_by = avg_log2FC) %>%
    dplyr::filter(cluster == 3)
  
  FeaturePlot(
    object = b,
    features = "Cd8a",
    cols = c("lightgrey", "#CD0000"),
    order = TRUE,
    reduction = "tsne",
    max.cutoff = 3,
  )
}

  




# Save image --------------------------------------------------------------

save.image(file = "data/rda/04-region-annotation.rda")
