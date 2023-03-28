# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Tue Mar 14 16:41:18 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(HGNChelper)


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


fn_make_gs_list_l3 <- function() {
  
  im <- readr::read_tsv(file = "/home/liuc9/data/refdata/scref/human_pbmc/celltype.l3.tsv") |> 
    dplyr::mutate(tissueType = "Immune system") |> 
    dplyr::mutate(geneSymbolmore2 = NA_character_) |> 
    dplyr::select(
      tissueType,
      cellName = `Expanded Label`,
      # cellName = Label,
      geneSymbolmore1 = Markers,
      geneSymbolmore2,
      shortName = Label,
      # shortName = `Expanded Label`
    )
  cortex <- readr::read_tsv(file = "/mnt/isilon/xing_lab/liuc9/refdata/scref/human_motorcortex/cluster.tsv") |> 
    dplyr::mutate(tissueType = "Brain") |> 
    dplyr::mutate(geneSymbolmore2 = NA_character_) |> 
    dplyr::select(
      tissueType,
      cellName = `Expanded Label`,
      # cellName = Label,
      geneSymbolmore1 = Markers,
      geneSymbolmore2,
      shortName = Label,
      # shortName = `Expanded Label`
    )
  
  dplyr::bind_rows(im, cortex) |> 
    writexl::write_xlsx(
      path = "/home/liuc9/github/scbrain/data/scuvresult/ScTypeDB_full_manual.xlsx"
    )
  
  
}



fn_sctype <- function(.sct_cluster, .annotissuetype, .res = 0.3, gs_list = NULL) {
  # .sct_cluster <- project_sct_region_cluster_allmarkers$sct_cluster[[1]]
  
  db_full = "/home/liuc9/github/scbrain/data/scuvresult/ScTypeDB_full_manual.xlsx"
  
  gs_list = gene_sets_prepare(db_full, .annotissuetype)
  
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
    # coord_fixed(
    #   ratio = 1,
    # ) +
    .split +
    .labs 
}

# load data ---------------------------------------------------------------

project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot <- readr::read_rds(
  file = "data/scuvrda/project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot.rds.gz"
) 


# body --------------------------------------------------------------------
fn_make_gs_list_l3()
# gs_list <- fn_gs_list_l3()

project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot |> 
  dplyr::mutate(
    annotissuetype = ifelse(region == "Brain", "Brain", "Immune system")
  ) |> 
  dplyr::mutate(
    sct_cluster_sctype_l3 = purrr::map2(
      .x = sct_cluster,
      .y = annotissuetype,
      .f = fn_sctype,
      .res = 0.3
    )
  ) ->
  project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot_anno_l3



project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot_anno_l3 |> 
  dplyr::mutate(
    sct_cluster_sctype_l3_umap = purrr::map(
      .x = sct_cluster_sctype_l3,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "umap"
    )
  ) %>%
  dplyr::mutate(
    sct_cluster_sctype_l3_tsne = purrr::map(
      .x = sct_cluster_sctype_l3,
      .f = fn_plot_umap_tsne,
      .celltype = "sctype",
      .reduction = "tsne"
    )
  ) ->
  project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot_anno_l3_sctype_l3


project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot_anno_l3_sctype_l3 |> 
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .heatmap10 = heatmap10,
        .tsne = sct_cluster_sctype_l3_tsne,
        .umap = sct_cluster_sctype_l3_umap,
        .marker_dotplot = marker_dotplot
      ),
      .f = function(.region, .case, .heatmap10, .tsne, .umap, .marker_dotplot, .outdir) {
        .filename <- glue::glue("{.region}_{.case}")
        
        ggsave(
          filename = glue::glue("{.filename}_l3_umap.pdf"),
          plot = .umap,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
        ggsave(
          filename = glue::glue("{.filename}_l3_tsne.pdf"),
          plot = .tsne,
          device = "pdf",
          path = .outdir,
          width = 12,
          height = 6
        )
        
      },
      .outdir = "/home/liuc9/github/scbrain/data/scuvresult/05-indivivudal-tissue"
    )
  )


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
    fullname = "Ttr+ Ependymal cell",
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

{
  b <- project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot$sct_cluster_sctype[[2]]
  d <- project_sct_region_cluster_allmarkers_heatmap10_sctype_umap_tsne_marker_dot$allmarkers[[2]]
  
  
  d %>%
    dplyr::filter(!grepl("^mt-", gene)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 10, order_by = avg_log2FC) %>%
    # dplyr::filter(avg_log2FC > 3) %>%
    dplyr::filter(cluster == 0) |> 
    dplyr::select(cluster, gene, avg_log2FC)
  
  b@meta.data |> 
    tibble::as_tibble() |> 
    dplyr::select(seurat_clusters, sctype) |> 
    dplyr::distinct() |> 
    dplyr::arrange(seurat_clusters)
  
  intersect(rownames(b), "Itm2a")
  
    FeaturePlot(
      object = b,
      features = c("Itm2a", "Cldn5"),
      cols = c("lightgrey", "#CD0000"),
      order = TRUE,
      reduction = "umap",
      max.cutoff = 3,
    )
    
    VlnPlot(b, features = c("Itm2a", "Cldn5"))
  
}

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------