# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed Mar  8 19:06:14 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)
library(HGNChelper)

# src ---------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

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
    Seurat::FindClusters(resolution = c(0.3)) ->
    .sct_cluster
  .sct_cluster <- Seurat::PrepSCTFindMarkers(object = .sct_cluster)
  .sct_cluster
}

fn_qc <- function(.region, .sct) {
  # .region <- brain_meninge_skull_sct_cluster$region[[1]]
  # .sct <- brain_meninge_skull_sct_cluster $sct[[1]]
  
  ncolors <- length(unique(.sct$case))
  vcolors <- viridis::viridis_pal()(ncolors)
  # # .sct
  Seurat::VlnPlot(
    object = .sct,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    cols = vcolors,
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
    cols = vcolors
  )
  .plot2 <- Seurat::FeatureScatter(
    object = .sct,
    feature1 = "nCount_RNA",
    feature2 = "percent.ribo",
    group.by = "case",
    cols = vcolors,
  )
  .plot3 <- Seurat::FeatureScatter(
    object = .sct,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "case",
    cols = vcolors
  )

  
  
  .vlnplot / .plot3 / (.plot1 | .plot2)  +
    plot_annotation(
      title = glue::glue("Quality control {.region}"),
      tag_levels = "A"
    ) +
    plot_layout(guides='collect') &
    theme(legend.position='bottom') -> 
    .qc#; .qc
  
  ggsave(
    filename = glue::glue("quality_control_{.region}.pdf"),
    plot = .qc,
    device = "pdf",
    path = "data/scuvresult/03-qc",
    width = 9,
    height = 15
    
  )
  
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
  # .sct_cluster <- brain_meninge_skull_sct_cluster$sct_cluster[[1]]
  .snn_res <- glue::glue("integrated_snn_res.{.res}")
  source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/chunjie-sam-liu/sc-type/master/R/sctype_score_.R")
  
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

fn_plot_umap_tsne <- function(.x, .celltype="sctype", .reduction="umap", .facet = FALSE) {
  # .x <- brain_meninge_skull_sct_cluster_sctype$sct_cluster_sctype[[1]]
  # .celltype="sctype"
  # .reduction="umap"
  
  # DimPlot(.x, reduction = "umap", split.by = "case")
  
  .umap <- as.data.frame(.x@reductions[[.reduction]]@cell.embeddings)
  colnames(.umap) <- c("UMAP_1", "UMAP_2")
  .xx <- .x@meta.data[, c("seurat_clusters", "case", .celltype)] %>%
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
  
  .split <- if(.facet) {
    facet_wrap(~case, nrow = 1)
  } else {
    theme()
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
    .split +
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

fn_gene_dotplot_new <- function(.sct_cluster, .marker) {
  .marker %>%
    dplyr::pull(markers) %>%
    stringr::str_split(pattern = ",") %>%
    unlist() ->
    .m
  .mm <- intersect(.m, rownames(.sct_cluster))
  DefaultAssay(.sct_cluster) <- "SCT"
  
  Idents(.sct_cluster) <-  "celltype"
  
  DotPlot(
    .sct_cluster,
    features = unique(.mm),
    cols = c("blue", "red"),
    dot.scale = 8
  ) +
    RotatedAxis() +
    labs(
      x = "Genes",
      y = "Cell types"
    )
}


# load data ---------------------------------------------------------------
regions <- c(
  "UVB" = "Brain",
  "UVM" = "Meninge",
  "UVS" = "Skull",
  "CB" = "Brain", 
  "CM" = "Meninge", 
  "CS" = "Skull", 
  "DB" = "Brain", 
  "DM" = "Meninge", 
  "DS" = "Skull"
)
cases <- c(
  "UVB" = "UV",
  "UVM" = "UV",
  "UVS" = "UV",
  "CB" = "Sham", 
  "CM" = "Sham", 
  "CS" = "Sham", 
  "DB" = "MCAO", 
  "DM" = "MCAO", 
  "DS" = "MCAO"
)

sc_sct_raw<- readr::read_rds(
  file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
)

sc_sct_raw |> 
  dplyr::mutate(
    region = plyr::revalue(x = project, replace = regions),
    case = plyr::revalue(x = project, replace = cases)
  ) ->
  project_sct_region


# split by region ---------------------------------------------------------

project_sct_region |> 
  dplyr::filter(region == "Brain") |> 
  dplyr::select(project, sct) |> 
  tibble::deframe() ->
  brain_sct_list

project_sct_region %>%
  dplyr::filter(region == "Meninge") %>%
  dplyr::select(project, sct) %>%
  tibble::deframe() ->
  meninge_sct_list


project_sct_region %>%
  dplyr::filter(region == "Skull") %>%
  dplyr::select(project, sct) %>%
  tibble::deframe() ->
  skull_sct_list


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
  ) |> 
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
  file = "data/scuvrda/brain_meninge_skull_sct_cluster.rds.gz"
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

# Annotation --------------------------------------------------------------


gs_list <- fn_gs_list()

brain_meninge_skull_sct_cluster %>%
  dplyr::mutate(
    sct_cluster_sctype = purrr::map(
      .x = sct_cluster,
      .f = fn_sctype,
      .res = 0.3,
      gs_list = gs_list
    )
  ) ->
  brain_meninge_skull_sct_cluster_annotation

readr::write_rds(
  x = brain_meninge_skull_sct_cluster_annotation,
  file = "data/scuvrda/brain_meninge_skull_sct_cluster_annotation.rds.gz"
)


# tsne and umap plot ------------------------------------------------------



brain_meninge_skull_sct_cluster_annotation |> 
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

brain_meninge_skull_sct_cluster_sctype |> 
  dplyr::mutate(
    sct_cluster_sctype_umap_split = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "umap",
      .facet = TRUE
    )
  ) |> 
  dplyr::mutate(
    sct_cluster_sctype_tsne_split = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "tsne",
      .facet = TRUE
    )
  ) ->
  brain_meninge_skull_sct_cluster_sctype_split

# save image

brain_meninge_skull_sct_cluster_sctype_split |> 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .umap = sct_cluster_sctype_umap,
        .tsne = sct_cluster_sctype_umap,
        .umap_split = sct_cluster_sctype_umap_split,
        .tsne_split = sct_cluster_sctype_tsne_split
      ),
      .f = function(.region, .umap, .tsne, .umap_split, .tsne_split) {
        # .filename <- glue::glue("{.region}")
        
        # .region <- brain_meninge_skull_sct_cluster_sctype_split$region[[1]]
        # .umap <- brain_meninge_skull_sct_cluster_sctype_split$sct_cluster_sctype_umap[[1]]
        # .tsne <- brain_meninge_skull_sct_cluster_sctype_split$sct_cluster_sctype_tsne[[1]]
        # .umap_split <- brain_meninge_skull_sct_cluster_sctype_split$sct_cluster_sctype_umap_split[[1]]
        # .tsne_split <- brain_meninge_skull_sct_cluster_sctype_split$sct_cluster_sctype_tsne_split[[1]]
        ggsave(
          filename = glue::glue("{.region}_umap.pdf"),
          plot = .umap,
          device = "pdf",
          path = "data/scuvresult/04-cluster",
          width = 12,
          height = 6
        )
        ggsave(
          filename = glue::glue("{.region}_tsne.pdf"),
          plot = .tsne,
          device = "pdf",
          path = "data/scuvresult/04-cluster",
          width = 12,
          height = 6
        )
        ggsave(
          filename = glue::glue("{.region}_umap_split.pdf"),
          plot = .umap_split,
          device = "pdf",
          path = "data/scuvresult/04-cluster",
          width = 18,
          height = 8
        )
        ggsave(
          filename = glue::glue("{.region}_tsne_split.pdf"),
          plot = .tsne_split,
          device = "pdf",
          path = "data/scuvresult/04-cluster",
          width = 18,
          height = 8
        )
      }
    )
  )
  

# Marker genes ------------------------------------------------------------

brain_meninge_skull_sct_cluster_sctype %>%
  dplyr::mutate(
    marker_genes = purrr::map(
      .x = sct_cluster_sctype,
      .f = fn_find_all_markers,
      .res = 0.3
    )
  ) ->
  brain_meninge_skull_sct_cluster_sctype_allmarkers

brain_meninge_skull_sct_cluster_sctype_allmarkers |> 
  dplyr::mutate(
    gene_dotplot = purrr::map2(
      .x = sct_cluster_sctype,
      .y = marker_genes,
      .f = fn_gene_dotplot,
      .n = 2
    )
  ) ->
  brain_meninge_skull_sct_cluster_sctype_allmarkers_dotplot

readr::write_rds(
  x = brain_meninge_skull_sct_cluster_sctype_allmarkers,
  file = "data/scuvrda/brain_meninge_skull_sct_cluster_sctype_allmarkers.rds.gz"
)

# body --------------------------------------------------------------------





# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/scuvrda/09-region-annotation.rda.gz")
