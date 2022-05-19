# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May 19 13:59:07 2022
# @DESCRIPTION: data integration

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(rlang)
library(patchwork)

# Load data ---------------------------------------------------------------

data_dir <- "data/raw/singlecell/report"

tibble::tibble(
  dir_path = list.dirs(
    file.path(data_dir, "1_Cellranger"),
    recursive = FALSE
  )
) %>% 
  dplyr::mutate(project = basename(dir_path)) ->
  project_path

# Function ----------------------------------------------------------------

fn_load_sc_10x <- function(.x) {
  # .x <- project_path$dir_path[[1]]
  # Create all Seurat objects
  # QC and filter each object
  # Calculate percent.mt on raw counts
  # SCTransform(do.scale=FALSE, do.center=FALSE) each Seurat object
  # CellCycleScoring() (while active assay is still SCT Assay)
  # Find cc score differences
  # Change defaultAssay <- 'RNA'
  # SCTransform(vars.to.regress = c("percent.mt", "CC.Difference"), do.scale=TRUE, do.center=TRUE) each Seurat object
  # Merge all Seurat objects (if I want to merge them for pooled analysis)
  # RunPCA()
  # RunUMAP()
  # FindNeighbors()
  # FindCluster()
  # Switch to RNA Assay
  # Perform DE analysis
  
  .regions <- c("CB" = "Brain", "CM" = "Meninge", "CS" = "Skull", "DB" = "Brain", "DM" = "Meninge", "DS" = "Skull")
  .cases <- c("CB" = "Sham", "CM" = "Sham", "CS" = "Sham", "DB" = "MCAO", "DM" = "MCAO", "DS" = "MCAO")
  
  .project <- basename(.x)
  .case <- plyr::revalue(x = .project, replace = .cases)
  .region <- plyr::revalue(x = .project, replace = .regions)
  
  .dir <- file.path(.x, "filtered_feature_bc_matrix")
  .counts <- Seurat::Read10X(data.dir = .dir)
  .sc <- Seurat::CreateSeuratObject(
    counts = .counts,
    project = .project,
    min.cells = 3, 
    min.features = 200
  )
  
  .sc$tissue <- .project
  .sc$case <- .case
  .sc$region <- .region
  
  .sc <- Seurat::PercentageFeatureSet(
    object = .sc, 
    pattern = "^mt-", 
    col.name = "percent.mt"
  )
  
  .vlnplot <- Seurat::VlnPlot(
    .sc, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    ncol = 3
  )
  ggsave(
    filename = "{.project}-vlnplot.pdf" %>% glue::glue(),
    plot = .vlnplot,
    device = "pdf",
    path = "data/result/"
  )
  .plot1 <- Seurat::FeatureScatter(
    object = .sc, 
    feature1 = "nCount_RNA", 
    feature2 = "percent.mt"
  )
  .plot2 <- Seurat::FeatureScatter(
    object = .sc, 
    feature1 = "nCount_RNA", 
    feature2 = "nFeature_RNA"
  )
  .plot <- Seurat::CombinePlots(plots = list(.plot1, .plot2))
  
  .sc_sub <- subset(
    x = .sc, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5
  )
  
  .sc_sub_sct <- Seurat::SCTransform(
    object = .sc_sub,
    do.scale = FALSE,
    do.center = FALSE
  )
  
  .sc_sub_sct <- Seurat::CellCycleScoring(
    object = .sc_sub_sct,
    g2m.features = Seurat::cc.genes$g2m.genes,
    s.features = Seurat::cc.genes$s.genes
  )
  
  .sc_sub_sct$CC.Difference <- .sc_sub_sct$S.Score - .sc_sub_sct$G2M.Score
  
  Seurat::DefaultAssay(.sc_sub_sct) <- "RNA"
  
  .sc_sub_sct_sct <- Seurat::SCTransform(
    object = .sc_sub_sct,
    # method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "CC.Difference"),
    do.scale = TRUE,
    do.center = TRUE
  )
  
  .sc_sub_sct_sct
}


# Load single cell data ---------------------------------------------------

project_sc <- project_path %>% 
  dplyr::mutate(
    sc = purrr::map(
      .x = dir_path,
      .f = fn_load_sc_10x
    )
  )

readr::write_rds(
  x = project_sc, 
  file = "analysis/2022-02-08-single-cell/rda/project_sc.rds.gz"
)
project_sc <- readr::read_rds(file ="analysis/2022-02-08-single-cell/rda/project_sc.rds.gz")
# Integration -------------------------------------------------------------

sc_list <- project_sc$sc
names(sc_list) <- project_sc$project

features <- Seurat::SelectIntegrationFeatures(
  object.list = sc_list, 
  nfeatures = 3000
)
sc_list <- Seurat::PrepSCTIntegration(
  object.list = sc_list, 
  anchor.features = features
)
anchors <- Seurat::FindIntegrationAnchors(
  object.list = sc_list,
  normalization.method = "SCT",
  anchor.features = features
)
sc_sct <- Seurat::IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

readr::write_rds(
  x = sc_sct, 
  file = "analysis/2022-02-08-single-cell/rda/sc_sct.rds.gz"
)

# Cluster -----------------------------------------------------------------


sc_sct %>% 
  Seurat::RunPCA() %>% 
  Seurat::RunUMAP(reduction = "pca", dims = 1:30) %>% 
  Seurat::RunTSNE(reduction = "pca", dims = 1:30) %>% 
  Seurat::FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  Seurat::FindClusters(resolution = 0.4) ->
  sc_sct_cluster

readr::write_rds(
  x = sc_sct_cluster, 
  file = "analysis/2022-02-08-single-cell/rda/sc_sct_cluster.rds.gz"
)


Seurat::DimPlot(
  sc_sct_cluster,
  reduction = "tsne",
  label = TRUE,
  repel = TRUE
)



# Cell annotation ---------------------------------------------------------


# sc-Type -----------------------------------------------------------------


library(HGNChelper)
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

sctype_scores

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sc_sct_cluster@meta.data$customclassif <- ""
for(j in unique(sctype_scores$cluster)) {
  cl_type = sctype_scores[sctype_scores$cluster == j, ]
  sc_sct_cluster@meta.data$customclassif[sc_sct_cluster@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Seurat::DimPlot(
  sc_sct_cluster, 
  reduction = "tsne", 
  label = TRUE, 
  # repel = TRUE, 
  group.by = 'customclassif'
)  +
  theme(
    legend.position = "bottom"
  )



# scMCA -------------------------------------------------------------------

mca <- scMCA::scMCA(
  scdata = exp(sc_sct_cluster[["SCT"]]@scale.data),
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

save.image(file = "analysis/2022-02-08-single-cell/rda/04-sctransform.rda")
