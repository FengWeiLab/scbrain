# dir("/home/liuc9/data/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022", recursive = F, full.names = T)

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

  # ggsave(paste0(qcdir, substr(colnames(GetAssayData(x))[1], 1, regexpr("\\_", colnames(GetAssayData(x))[1])-1), "_QC.eps"), width = 6, height = 10)
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

  # ggsave(paste0(qcdir, substr(colnames(GetAssayData(x))[1], 1, regexpr("\\_", colnames(GetAssayData(x))[1])-1), "_filt.eps"),width = 6, height = 10)
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
allso <- ScaleData(allso, features = rownames(allso),
                   vars.to.regress = c("nFeature_RNA","nCount_RNA","percent_mito"))

allso <- FindVariableFeatures(allso, selection.method = "vst")
allso <- RunPCA(allso, npcs = 50, features = VariableFeatures(allso))
ElbowPlot(allso, 50)

allso <- FindNeighbors(allso, dims = 1:30, reduction = "pca")
allso <- FindClusters(allso, resolution = 1.2)

allso <- RunTSNE(allso, dims = 1:30, verbose = T)

TSNEPlot(allso, cols = c(colpan2, colpan1, colpan3, colpan4))
# ggsave(paste0(tsnedir, "tsne_clusters.eps"), width = 7, height = 6)

TSNEPlot(allso, group.by = "orig.ident", cols = colpan4)
# ggsave(paste0(tsnedir, "tsne_bysamp.eps"), width = 7, height = 6)

TSNEPlot(allso, split.by = "orig.ident", cols = c(colpan2, colpan1, colpan3, colpan4))
# ggsave(paste0(tsnedir, "tsne_clusters_bysamp.eps"), width = 9, height = 4)

# ====================== Characterize Clusters ==========================
# QC plots
FeaturePlot(allso, "nCount_RNA")
# ggsave(paste0(qcdir, "nCountRNA.eps"), height = 6, width = 7)

FeaturePlot(allso, "percent_mito")
# ggsave(paste0(qcdir, "mito.eps"), height = 6, width = 7)

FeaturePlot(allso, "nFeature_RNA")
# ggsave(paste0(qcdir, "nFeatureRNA.eps"), height = 6, width = 7)

FeaturePlot(allso, "Mki67")
# ggsave(paste0(qcdir, "mki67.eps"), height = 6, width = 7)

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
  # write.csv(temp, paste0(cmarkdir, "cmarks_cluster_", i, ".csv"), row.names = F)
}

# identify cell types
DotPlot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Igll1','Fcer2a','Cd3e','Siglech','Bst2',
                            'Cd74','Ms4a6c','Ms4a7','Plac8','Ly6c2','Lyz2','Adgre1','Mki67','Ptprc',
                            'Fcer1a','Il3ra','Kit','Hemgn','Car1','Hba-a2','Hbb-bs','Jchain','Sdc1','Mpl','Slamf1')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))

## add proposed annotation
new.clusters <- c('0_Neutrophils','1_Neutrophils','2_Immature B cells','3_Neutrophils','4_Macrophages',
                  '5_Neutrophils','6_Mature B cells','7_Proliferating Neutrophils','8_Immature B cells','9_HSCs',
                  '10_Pre-B cells','11_pDCs','12_T cells','13_Monocytes','14_HSCs',
                  '15_Proliferating Neutrophils','16_Mast cells','17_Neutrophils','18_HSCs','19_Pro-B cells',
                  '20_Monocytes','21_Basophils','22_Macrophages','23_Neutrophils','24_NK cells',
                  '25_Plasma cells','26_Erythroblasts','27_Olfactory sensory neurons')
names(new.clusters) <- levels(allso)
allso <- RenameIdents(allso, new.clusters)
my_levels <- c('0_Neutrophils','1_Neutrophils','3_Neutrophils','5_Neutrophils','17_Neutrophils','23_Neutrophils',
               '7_Proliferating Neutrophils','15_Proliferating Neutrophils',
               '19_Pro-B cells','10_Pre-B cells','2_Immature B cells','8_Immature B cells','6_Mature B cells','25_Plasma cells',
               '12_T cells','24_NK cells','11_pDCs','4_Macrophages','22_Macrophages','13_Monocytes','20_Monocytes',
               '21_Basophils','16_Mast cells','14_HSCs','18_HSCs','9_HSCs',
               '26_Erythroblasts','27_Olfactory sensory neurons')
levels(allso) <- my_levels

## collapse clusters
new.clusters <- c('Neutrophils','Neutrophils','Neutrophils','Neutrophils','Neutrophils','Neutrophils',
                  'Proliferating Neutrophils','Proliferating Neutrophils',
                  'Pro-B cells','Pre-B cells','Immature B cells','Immature B cells','Mature B cells','Plasma cells',
                  'T cells','NK cells','pDCs','Macrophages','Macrophages','Monocytes','Monocytes',
                  'Basophils','Mast cells','HSCs','HSCs','HSCs',
                  'Erythroblasts','Olfactory sensory neurons')
names(new.clusters) <- levels(allso)
allso <- RenameIdents(allso, new.clusters)

DotPlot(allso, features = c('Ly6g','S100a8','Mmp8','Ngp','Cd19','Vpreb1','Vpreb2','Vpreb3','Mki67','Rag1','Rag2','Fcer2a','Jchain','Xbp1',
                            'Cd3e','Klrb1c','Siglech','Bst2',
                            'Cd74','Ms4a6c','Plac8','Ly6c2','Lyz2','Adgre1','Ptprc',
                            'Fcer1a','Il3ra','Kit','Aqp1','Hemgn','Car1','Hba-a2','Hbb-bs','Omp','Stoml3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))
ggsave(paste0(cmarkdir,"dotplot.eps"), width = 12, height = 5)

## cluster proportion barplot
prop <- data.frame(prop.table(table(Idents(allso), allso$orig.ident), margin = 2))
colnames(prop) <- c("cluster","sample", "percent")

ggplot(data=prop, aes(x=sample, y=percent, fill=cluster)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=c(colpan2,colpan1,colpan3))+
  theme_minimal()
ggsave(paste0(resdir,"cluster_props.eps"), width = 5, height = 7)

## absolute cluster counts barplot
prop <- data.frame(table(Idents(allso), allso$orig.ident))
colnames(prop) <- c("cluster","sample", "count")
write.csv(prop, paste0(resdir, "cluster_counts.csv"))

ggplot(data=prop, aes(x=sample, y=count, fill=cluster)) +
  geom_bar(stat = "identity") + scale_fill_manual(values=c(colpan2,colpan1))+
  theme_minimal()
ggsave(paste0(resdir,"cluster_counts.eps"), width = 4, height = 5)

TSNEPlot(allso, cols = c(colpan2, colpan1, colpan3))
ggsave(paste0(tsnedir, "tsne.annot.eps"), width = 8, height = 6)

TSNEPlot(allso, cols = c(colpan2, colpan1, colpan3), split.by = 'orig.ident')
ggsave(paste0(tsnedir, "tsne_annot_bysamp.eps"), width = 11, height = 4)

allso$celltype <- allso@active.ident
