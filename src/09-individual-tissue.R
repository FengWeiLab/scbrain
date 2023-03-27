# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Mar 19 14:39:34 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)
library(Seurat)
library(HGNChelper)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

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

fn_find_all_markers <- function(.sct_cluster) {
  Seurat::FindAllMarkers(
    object = .sct_cluster,
    assay = "SCT",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
}

fn_heatmap10 <- function(.sct_cluster,.region,.case,.allmarkers, .outdir){
  # .sct_cluster <- project_sct_region_cluster_allmarkers$sct_cluster[[1]]
  # .region <- project_sct_region_cluster_allmarkers$region[[1]]
  # .case <- project_sct_region_cluster_allmarkers$case[[1]]
  # .allmarkers <- project_sct_region_cluster_allmarkers$allmarkers[[1]]
  # .outdir <-  "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
  # 
  .filename <- glue::glue("{.region}_{.case}_heatmap_top10.pdf")
  
  .allmarkers |> 
    dplyr::group_by(cluster) |> 
    dplyr::top_n(10, wt = avg_log2FC) ->
    .top10
  
  p <- DoHeatmap(.sct_cluster, features = .top10$gene) + NoLegend()
  p
  # ggsave(
  #   filename = .filename,
  #   plot = p,
  #   device = "pdf",
  #   path = .outdir,
  #   width = 12,
  #   height = 15
  # )
  
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
  # .sct_cluster <- project_sct_region_cluster_allmarkers$sct_cluster[[1]]
  .snn_res <- glue::glue("SCT_snn_res.{.res}")
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
  pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")
  
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


fn_gene_dotplot <- function(.sct_cluster, .marker, .n = 3) {
  
  .marker %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = .n, order_by = avg_log2FC)  ->
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


#
# load data ---------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")


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

# individuals <- c(
#   "UVB" = "Brain UV",
#   "UVM" = "Meninge UV",
#   "UVS" = "Skull UV",
#   "CB" = "Brain Sham ", 
#   "CM" = "Meninge Sham", 
#   "CS" = "Skull Sham", 
#   "DB" = "Brain MCAO", 
#   "DM" = "Meninge MCAO", 
#   "DS" = "Skull MCAO"
# )

sc_sct_raw<- readr::read_rds(
  file = "data/scuvrda/sc_sham_mcao_uv_sct.rds.gz"
)

sc_sct_raw |> 
  dplyr::mutate(
    region = plyr::revalue(x = project, replace = regions),
    case = plyr::revalue(x = project, replace = cases),
    # individual = plyr::revalue(x = project, replace = individuals)
  ) ->
  project_sct_region


# body --------------------------------------------------------------------


# cluster -----------------------------------------------------------------


project_sct_region |> 
  dplyr::mutate(
    sct_cluster = purrr::map(
      .x = sct,
      .f = fn_cluster
    )
  ) ->
  project_sct_region_cluster


project_sct_region_cluster |> 
  dplyr::mutate(
    allmarkers = purrr::map(
      .x = sct_cluster,
      .f = fn_find_all_markers
    )
  ) ->
  project_sct_region_cluster_allmarkers

readr::write_rds(
  x = project_sct_region_cluster_allmarkers,
  file = "data/scuvrda/project_sct_region_cluster_allmarkers.rds.gz"
)



# Heatmap top10 -----------------------------------------------------------

project_sct_region_cluster_allmarkers |> 
  dplyr::mutate(
    heatmap10 = purrr::pmap(
      .l = list(
        .sct_cluster = sct_cluster,
        .region = region,
        .case = case,
        .allmarkers = allmarkers
      ),
      .f = fn_heatmap10,
      .outdir = "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
    )
  ) ->
  project_sct_region_cluster_allmarkers_heatmap10

# project_sct_region_cluster_allmarkers_heatmap10$heatmap10[[1]]

# annotation --------------------------------------------------------------

gs_list <- fn_gs_list()

project_sct_region_cluster_allmarkers_heatmap10 |> 
  dplyr::mutate(
    sct_cluster_sctype = purrr::map(
      .x = sct_cluster,
      .f = fn_sctype,
      .res = 0.3,
      gs_list = gs_list
    )
  ) ->
  project_sct_region_cluster_allmarkers_heatmap10_sctype


project_sct_region_cluster_allmarkers_heatmap10_sctype |> 
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
  project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne


# dot ---------------------------------------------------------------------

project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne |> 
  dplyr::mutate(
    marker_dotplot = purrr::map2(
      .x = sct_cluster_sctype,
      .y = allmarkers,
      .f = fn_gene_dotplot,
      .n = 3
    )
  ) ->
  project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot



readr::write_rds(
  x = project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot,
  file = "data/scuvrda/project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot.rds.gz"
)



# save plots --------------------------------------------------------------

project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot |> 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .heatmap10 = heatmap10,
        .tsne = sct_cluster_sctype_tsne,
        .umap = sct_cluster_sctype_umap,
        .marker_dotplot = marker_dotplot
      ),
      .f = function(.region, .case, .heatmap10, .tsne, .umap, .marker_dotplot, .outdir) {
        .filename <- glue::glue("{.region}_{.case}")
        
        ggsave(
          filename = glue::glue("{.filename}_heatmap_top10.pdf"),
          plot = .heatmap10,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 15
        )
        
        ggsave(
          filename = glue::glue("{.filename}_umap.pdf"),
          plot = .umap,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
        ggsave(
          filename = glue::glue("{.filename}_tsne.pdf"),
          plot = .tsne,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
        ggsave(
          filename = glue::glue("{.filename}_tsne.pdf"),
          plot = .marker_dotplot,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
      },
      .outdir = "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
    )
  )


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "data/scuvrda/09-individual-tissue.rda.gz")
load(file = "data/scuvrda/09-individual-tissue.rda.gz")
