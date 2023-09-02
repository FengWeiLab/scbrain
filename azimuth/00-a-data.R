# dir("/home/liuc9/data/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022", recursive = F, full.names = T)

resdir <- "/home/liuc9/github/scbrain/test/"
tsnedir <- "/home/liuc9/github/scbrain/test/"
qcdir <- "/home/liuc9/github/scbrain/test/"
cmarkdir <- "/home/liuc9/github/scbrain/test/"
degdir <- "/home/liuc9/github/scbrain/test/"
sodir <- "/home/liuc9/github/scbrain/test/"
outdirs <- c(resdir,tsnedir,qcdir,cmarkdir,degdir,sodir)

library(Seurat) # v3.1.1

dirs <- list(
  "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/femur",
  "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/skull"
)

names(dirs) <- c("femur", "skull")



matlist <- lapply(dirs, Read10X)
for(i in 1:length(matlist)) colnames(matlist[[i]]) <- paste0(names(matlist)[[i]], "__", colnames(matlist[[i]]))


## ------------------------- colors ----------------------------
#beige, mint green, ocean blue, steel blue, rose, purple, army green, brown, grey, bubblegum pink
colpan1 = c('#C2B59B','#81DCAD','#0F77BA','#8192DC','#DC818F','#B347A6','#6FAA45','#967352',
            '#89878C','#EEC3DC')
#cyan, teal, goldenrod, gold, cherry, raspberry, orange, sherbert
colpan2 = c('#2EF9F0','#23BBB4','#F9C92E','#BB9723','#F92E73','#BB2356','#FF7F00','#FFBF80')
#periwinkle, baby blue, navy, kelly green, dark green, light green, lavender, dark purple
colpan3 = c('#73B4FF','#BFE0FF','#234A70','#20A815','#426C3F','#B9E0B5','#D5C0E5','#562699')
##contrasting palette
#lime green, army green, pink, magenta, purple, lavender, red, blue, gold
colpan4 = c('#9EFE42','#6A9242','#EEA3B3','#F60677','#6217DB','#C4A8F2','#29ABD8','#D2B70A')

# ===================================== functions ======================================
## --------------------------------- add QC metrics ---------------------------------
add_mito_data <- function(so_object){
  mito_genes <- grep(pattern = "^mt.", x = rownames(x =so_object@assays$RNA@counts))
  percent_mito <- Matrix::colSums(so_object@assays$RNA@counts[mito_genes,])/Matrix::colSums(so_object@assays$RNA@counts)
  so_object <- AddMetaData(object = so_object, metadata = percent_mito, col.name = "percent_mito")
  return(so_object)}
add_ribo_data <- function(so_object){
  ribo_genes <- grep(pattern = "^Rp[l,s].", x = rownames(x =so_object@assays$RNA@counts))
  percent_ribo <- Matrix::colSums(so_object@assays$RNA@counts[ribo_genes,])/Matrix::colSums(so_object@assays$RNA@counts)
  so_object <- AddMetaData(object = so_object, metadata = percent_ribo, col.name = "percent_ribo")
  return(so_object)}

## ---------------------------------- volcanos ----------------------------------
volcano_plot <- function(data, plottitle, up.color, down.color){
  keyvals <- ifelse(
    data$logFC < 0, down.color,
    ifelse(data$logFC > 0, up.color,
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == up.color] <- 'up'
  names(keyvals)[keyvals == down.color] <- 'down'

  EnhancedVolcano(data, lab = data$gene_id, selectLab = NULL, x = "logFC", y = "Pval",
                  labSize = 5.0, subtitle = NULL, title = plottitle,
                  pCutoff = 5e-2, FCcutoff = 0, axisLabSize = 30, gridlines.minor = FALSE, gridlines.major = FALSE,
                  labFace = "italic", legendLabSize = 0, legendIconSize = 0, colCustom = keyvals)
  return(plot)
}

## --------------------------------- pathways --------------------------------
pathway_up <- function(data, mart, GO.cat){
  annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
    filter="external_gene_name",
    values=data$gene_id,
    uniqueRows=TRUE)
  annotLookup <- data.frame(
    data[match(annotLookup$external_gene_name, data$gene_id),],
    annotLookup)
  annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
  up_matrix <- annotLookup_up$logFC
  names(up_matrix) <- annotLookup_up$entrezgene_id
  go_enrich_up <- enrichGO(gene = names(up_matrix),
                           OrgDb = 'org.Mm.eg.db',
                           readable = T,
                           ont = GO.cat,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.10)
  return(go_enrich_up)
}
pathway_down <- function(data, mart, GO.cat){
  annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
    filter="external_gene_name",
    values=data$gene_id,
    uniqueRows=TRUE)
  annotLookup <- data.frame(
    data[match(annotLookup$external_gene_name, data$gene_id),],
    annotLookup)
  annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)
  down_matrix <- annotLookup_down$logFC
  names(down_matrix) <- annotLookup_down$entrezgene_id
  go_enrich_down <- enrichGO(gene = names(down_matrix),
                             OrgDb = 'org.Mm.eg.db',
                             readable = T,
                             ont = GO.cat,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.10)
  return(go_enrich_down)
}
up_path_plot <- function(result, plottitle, color){
  plot <- ggplot(result[1:10], aes(x=reorder(Description, -pvalue), y=Count, fill=p.adjust)) +
    geom_bar(colour = "black", stat = "identity") +
    coord_flip() +
    scale_fill_continuous(low=color, high="white", guide = guide_colorbar(frame.colour = "black",
                                                                          ticks.colour = "black")) +
    labs(x = NULL, y = "Number of Genes",
         title = plottitle, fill = "p.adjust") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    theme(text = element_text(size=10), axis.title = element_text(size = 16),
          axis.text = element_text(size = 18),
          legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(plot)
}
down_path_plot <- function(result, plottitle, color){
  plot <- ggplot(result[1:10], aes(x=reorder(Description, -pvalue), y=Count, fill=p.adjust)) +
    geom_bar(colour = "black", stat = "identity") +
    coord_flip() +
    scale_fill_continuous(low=color, high="white", guide = guide_colorbar(frame.colour = "black",
                                                                          ticks.colour = "black")) +
    labs(x = NULL, y = "Number of Genes",
         title = plottitle, fill = "p.adjust") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    theme(text = element_text(size=10), axis.title = element_text(size = 16),
          axis.text = element_text(size = 18),
          legend.text = element_text(size = 16), legend.title = element_text(size = 16)) +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(plot)
}

# ------------------- Make Seurat Object / QC --------------------
so.list <- lapply(matlist, function(x){CreateSeuratObject(x, min.cells = 6, names.delim = '__')})
so.list <- lapply(so.list, add_mito_data)
so.list <- lapply(so.list, add_ribo_data)

lapply(so.list, function(x){
  VlnPlot(x, c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"),ncol = 2)

  # ggsave(paste0(qcdir, substr(colnames(GetAssayData(x))[1], 1, regexpr("\\_", colnames(GetAssayData(x))[1])-1), "_QC.pdf"), width = 6, height = 10)
  })

lapply(so.list, dim)
so.list <- lapply(so.list, function(x){subset(x, subset = percent_mito < 0.25 &
                                                nFeature_RNA > 200  &
                                                nFeature_RNA < 7000 &
                                                nCount_RNA < 100000 &
                                                nCount_RNA > 1000)})
lapply(so.list, dim)
lapply(so.list, function(x){
  VlnPlot(x, c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"),ncol = 2)

  # ggsave(paste0(qcdir, substr(colnames(GetAssayData(x))[1], 1, regexpr("\\_", colnames(GetAssayData(x))[1])-1), "_filt.pdf"),width = 6, height = 10)
  })

# ------------------------- Normalization ----------------------
set.seed(701)
allso <- merge(so.list[[1]], y = c(so.list[[2]]))
allsce <- as.SingleCellExperiment(allso)
library(scran)
clusters <- quickCluster(allsce, min.mean=0.1, method="igraph")

allsce <- computeSumFactors(allsce, min.mean=0.1, clusters=clusters)

allsce <- scater::logNormCounts(allsce, log = F)

## log for compatibility and then add normalized values to Seurat object
allnorm <- log1p(allsce@assays@data@listData$normcounts)
allso <- NormalizeData(allso)
allso@assays$RNA@data <- allnorm

# ---------------- Scale / Cluster / Visualize ------------------
allso <- ScaleData(allso, features = rownames(allso), vars.to.regress = c("nFeature_RNA","nCount_RNA","percent_mito"))

allso <- FindVariableFeatures(allso, selection.method = "vst")
allso <- RunPCA(allso, npcs = 50, features = VariableFeatures(allso))
ElbowPlot(allso, 50)

allso <- FindNeighbors(allso, dims = 1:30, reduction = "pca")
allso <- FindClusters(allso, resolution = 1.1)

allso <- RunTSNE(allso, dims = 1:30, verbose = T)

TSNEPlot(allso, cols = c(colpan2, colpan1, colpan3, colpan4))
ggsave(paste0(tsnedir, "tsne_clusters.pdf"), width = 7, height = 6)

TSNEPlot(allso, group.by = "orig.ident", cols = colpan4)
ggsave(paste0(tsnedir, "tsne_bysamp.pdf"), width = 7, height = 6)

TSNEPlot(allso, split.by = "orig.ident", cols = c(colpan2, colpan1, colpan3, colpan4))
ggsave(paste0(tsnedir, "tsne_clusters_bysamp.pdf"), width = 9, height = 4)

# ====================== Characterize Clusters ==========================
# QC plots
FeaturePlot(allso, "nCount_RNA")
ggsave(paste0(qcdir, "nCountRNA.pdf"), height = 6, width = 7)

FeaturePlot(allso, "percent_mito")
ggsave(paste0(qcdir, "mito.pdf"), height = 6, width = 7)

FeaturePlot(allso, "nFeature_RNA")
ggsave(paste0(qcdir, "nFeatureRNA.pdf"), height = 6, width = 7)

FeaturePlot(allso, "Mki67")
ggsave(paste0(qcdir, "mki67.pdf"), height = 6, width = 7)

## get cluster markers
future::plan(future::multisession, workers = 10)
clustmarks <- FindAllMarkers(allso, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
future::plan(future::sequential)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

gl<- c()
for(i in as.character(unique(clustmarks$cluster))){
  temp<- clustmarks[which(clustmarks$cluster == i), ]
  temp<- temp[order(-temp$avg_logFC), c(7, 6, 2:5)]
  cat('cluster ', i, ': ', temp$gene[1:5], '\n\n')
  gl<- c(gl, temp$gene[1:5])
  write.csv(temp, paste0(cmarkdir, "cmarks_cluster_", i, ".csv"), row.names = F)
}

# identify cell types
DotPlot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Igll1','Fcer2a','Cd3e','Siglech','Bst2',
                            'Cd74','Ms4a6c','Ms4a7','Plac8','Ly6c2','Lyz2','Adgre1','Mki67','Ptprc',
                            'Fcer1a','Il3ra','Kit','Hemgn','Car1','Hba-a2','Hbb-bs','Jchain','Sdc1','Mpl','Slamf1')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))
# ggsave("b.pdf", width = 10)



## add proposed annotation
new.clusters <- c(
  '0_Neutrophils',
  '1_Neutrophils',
  '2_Immature B cells',
  '3_Neutrophils',
  '4_Macrophages',
  '5_Neutrophils',
  '6_Mature B cells',
  '7_Proliferating Neutrophils',
  '8_Immature B cells',
  '9_HSCs',
  '10_Pre-B cells',
  '11_pDCs',
  '12_T cells',
  '13_Monocytes',
  '14_HSCs',
  '15_Proliferating Neutrophils',
  '16_Mast cells',
  '17_Neutrophils',
  '18_HSCs',
  '19_Pro-B cells',
  '20_Monocytes',
  '21_Basophils',
  '22_Macrophages',
  '23_Neutrophils',
  '24_NK cells',
  '25_Plasma cells',
  '26_Erythroblasts',
  '27_Olfactory sensory neurons')

fn_marker_gene_dotplot <- function(object, assay = NULL, features, cols = c(
  "lightgrey",
  "blue"
), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 8,
idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE,
scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) {

  object <- allso
  assay = NULL
  features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Vpreb3','Mki67','Rag1','Rag2','Fcer2a','Jchain','Xbp1', 'Cd3e','Klrb1c','Siglech','Bst2', 'Cd74','Ms4a6c','Plac8','Ly6c2','Lyz2','Adgre1','Ptprc', 'Fcer1a','Il3ra','Kit','Aqp1','Hemgn','Car1','Hba-a2','Hbb-bs','Omp','Stoml3')
  cols = c("blue", "red")
  col.min = -2.5
  col.max = 2.5
  dot.min = 0
  dot.scale = 8
  idents = NULL
  group.by = NULL
  split.by = NULL
  cluster.idents = FALSE
  scale = TRUE
  scale.by = "radius"
  scale.min = NA
  scale.max = NA

  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in%
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by,
                       size = scale_size,
                       radius = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.",
              call. = FALSE,
              immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(
    object = object, vars = features,
    cells = cells
  )
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(
      rep(x = id.levels, each = length(x = unique.splits)),
      "_", rep(x = unique(x = splits), times = length(x = id.levels))
    )
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1),
                              drop = FALSE
    ]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(
      X = data.use, MARGIN = 2, FUN = PercentAbove,
      threshold = 0
    )
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(
      X = data.plot,
      FUN = unlist
    ))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results",
            call. = FALSE, immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot ==
                              x, "avg.exp"]
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(
          data = data.use, min = col.min,
          max = col.max
        )
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(
      x = avg.exp.scaled,
      breaks = 20
    ))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0(
        "^((",
        paste(sort(x = levels(x = object), decreasing = TRUE),
              collapse = "|"
        ), ")_)"
      ), replacement = "",
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  # object@meta.data |> dplyr::glimpse()
  object@meta.data |>
    dplyr::select(seurat_clusters) |>
    dplyr::distinct() |>
    dplyr::arrange(seurat_clusters) |>
    ggplot(aes(
      x = 1,
      y = seurat_clusters,
      label = seurat_clusters,
    )) +
    geom_tile(
      aes(
        width = 0.8,
        height = 0.8,
      )
    ) +
    geom_text(
      color = "white",
      size = 4
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    coord_fixed() ->
    p1;p1

  data.plot |>
    dplyr::mutate(
      avg.exp.scaled = ifelse(
        avg.exp.scaled > 2,
        2,
        avg.exp.scaled
      )
    ) |>
    dplyr::mutate(
      pct.exp = ifelse(
        pct.exp > 75,
        75,
        pct.exp
      )
    ) |>
    dplyr::mutate(
      id = as.character(id)
    ) |>
    ggplot(mapping = aes_string(
      x = "features.plot",
      y = "id"
    )) +
    geom_point(mapping = aes_string(
      size = "pct.exp",
      color = color.by
    )) +
    scale.func(
      range = c(0, dot.scale),
      limits = c(0, 75),
      breaks = c(0, 25, 50, 75),
      labels = c(0, 25, 50, 75)
    ) +
    scale_color_gradient(
      low = cols[1],
      high = cols[2],
      limits = c(-1, 2),
      name = "Average Expression"
    ) +
    scale_x_discrete(
      position = "bottom",
    ) +
    scale_y_discrete(
      limits = c(
        # Neutrophils
        0,
        2,
        3,
        5,
        15,
        22,
        6,
        # B cell
        18,
        10,
        1,
        7,
        11,
        24,
        # T
        13,
        # NK
        23,
        # DC
        12,
        # Macrophages
        4,
        21,
        26,
        # Monocytes
        8,
        # Basophils
        20,
        # Mast cells
        14,
        16,
        # HSCs
        9,
        # Erythroblasts
        17,
        25,
        # Sensory neurons
        19
      ) |>
        as.character()
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 10,
        color = "black",
        face = "italic"
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.text.y = element_text(size = 10),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    guides(size = guide_legend(title = "Percent Expressed (%)")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                    yes = "Identity", no = "Split Identity"
    ))

  cowplot::plot_grid(
    plotlist = list( p1, p2),
    align = 'h',
    rel_widths = c(0.2, 2)
  ) ->
    plot

  return(p2)
}


# fn_marker_gene_dotplot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Vpreb3','Mki67','Rag1','Rag2','Fcer2a','Jchain','Xbp1', 'Cd3e','Klrb1c','Siglech','Bst2', 'Cd74','Ms4a6c','Plac8','Ly6c2','Lyz2','Adgre1','Ptprc', 'Fcer1a','Il3ra','Kit','Aqp1','Hemgn','Car1','Hba-a2','Hbb-bs','Omp','Stoml3'))

c(
  # Neutrophils
  0,# Neutrophils
  2,# Neutrophils
  3,# Neutrophils
  5,# Neutrophils
  15,# Neutrophils
  22,# Neutrophils
  6, # Proliferating neutrophils
  # B cells
  18,# Pro-B cells
  10,# Pre-B cells
  1,# Immature B cells
  7,# Mature B cells
  11,# Mature B cells
  24,# Plasma cells
  # T
  13,# T cells
  # NK
  23,# NK cells
  # DC
  12,# DCs cells
  # Macrophages
  4,# Macrophages
  21,# Macrophages
  26,# Macrophages
  # Monocytes
  8,# Monocytes
  # Basophils
  20,# Basophils
  # Mast cells
  14,# Mast cells
  16,# Mast cells
  # HSCs
  9,# HSCs
  # Erythroblasts
  17,# Erythroblasts
  25,# Erythroblasts
  # Sensory neurons
  19# Sensory neurons
)
list("0" = "Neutrophils",  "2" = "Neutrophils",  "3" = "Neutrophils",  "5" = "Neutrophils",  "15" = "Neutrophils",  "22" = "Neutrophils",  "6" = "Proliferating neutrophils",  "18" = "Pro-B cells",  "10" = "Pre-B cells",  "1" = "Immature B cells",  "7" = "Mature B cells",  "11" = "Mature B cells",  "24" = "Plasma cells",  "13" = "T cells",  "23" = "NK cells",  "12" = "DCs cells",  "4" = "Macrophages",  "21" = "Macrophages",  "26" = "Macrophages",  "8" = "Monocytes",  "20" = "Basophils",  "14" = "Mast cells",  "16" = "Mast cells",  "9" = "HSCs",  "17" = "Erythroblasts",  "25" = "Erythroblasts",  "19" = "Sensory neurons") ->
  aa

aa |>
  tibble::enframe() |>
  tidyr::unnest() |>
  dplyr::mutate(cells = glue::glue("{name}_{value}")) |>
  dplyr::mutate(name = as.integer(name)) ->
  aaa

aaa |>
  dplyr::mutate(value = factor(value, unique(aaa$value))) |>
  dplyr::mutate(
    cells = factor(cells, unique(cells))
  ) |>
  dplyr::arrange(name) ->
  aaaa

new.clusters <- c(
  '0_Neutrophils',
  '1_Neutrophils',
  '2_Immature B cells',
  '3_Neutrophils',
  '4_Macrophages',
  '5_Neutrophils',
  '6_Mature B cells',
  '7_Proliferating Neutrophils',
  '8_Immature B cells',
  '9_HSCs',
  '10_Pre-B cells',
  '11_pDCs',
  '12_T cells',
  '13_Monocytes',
  '14_HSCs',
  '15_Proliferating Neutrophils',
  '16_Mast cells',
  '17_Neutrophils',
  '18_HSCs',
  '19_Pro-B cells',
  '20_Monocytes',
  '21_Basophils',
  '22_Macrophages',
  '23_Neutrophils',
  '24_NK cells',
  '25_Plasma cells',
  '26_Erythroblasts',
  '27_Olfactory sensory neurons')




new.clusters <- aaaa$cells |> as.character()

allso <- FindClusters(allso, resolution = 1.1)

names(new.clusters) <- allso$seurat_clusters |> levels()
allso <- RenameIdents(allso, new.clusters)

DotPlot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Vpreb3','Mki67','Rag1','Rag2','Fcer2a','Jchain','Xbp1',
                            'Cd3e','Klrb1c','Siglech','Bst2',
                            'Cd74','Ms4a6c','Plac8','Ly6c2','Lyz2','Adgre1','Ptprc',
                            'Fcer1a','Il3ra','Kit','Aqp1','Hemgn','Car1','Hba-a2','Hbb-bs','Omp','Stoml3') ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))

my_levels <- c('0_Neutrophils','1_Neutrophils','3_Neutrophils','5_Neutrophils','17_Neutrophils','23_Neutrophils',
               '7_Proliferating Neutrophils','15_Proliferating Neutrophils',
               '19_Pro-B cells','10_Pre-B cells','2_Immature B cells','8_Immature B cells','6_Mature B cells','25_Plasma cells',
               '12_T cells','24_NK cells','11_pDCs','4_Macrophages','22_Macrophages','13_Monocytes','20_Monocytes',
               '21_Basophils','16_Mast cells','14_HSCs','18_HSCs','9_HSCs',
               '26_Erythroblasts','27_Olfactory sensory neurons')
my_levels <- levels(aaaa$cells)
levels(allso) <- my_levels


DotPlot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Vpreb3','Mki67','Rag1','Rag2','Fcer2a','Jchain','Xbp1',
                            'Cd3e','Klrb1c','Siglech','Bst2',
                            'Cd74','Ms4a6c','Plac8','Ly6c2','Lyz2','Adgre1','Ptprc',
                            'Fcer1a','Il3ra','Kit','Aqp1','Hemgn','Car1','Hba-a2','Hbb-bs','Omp','Stoml3') ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))


# make refs ---------------------------------------------------------------

Idents(allso) |>
  tibble::enframe() |>
  dplyr::mutate(cl = allso$seurat_clusters) |>
  dplyr::select(
    ident = cl,
    cluster = value,
    cell = name
  ) |>
  readr::write_csv("/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/annot_skull.csv")
