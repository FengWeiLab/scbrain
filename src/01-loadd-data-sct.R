# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu May 19 13:59:07 2022
# @DESCRIPTION: load data and sctransform

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

# Function ----------------------------------------------------------------

fn_load_sc_10x <- function(.x) {
  # .x <- project_path$dir_path[[1]]
  .regions <- c(
    "CB" = "Brain", 
    "CM" = "Meninge", 
    "CS" = "Skull", 
    "DB" = "Brain", 
    "DM" = "Meninge", 
    "DS" = "Skull"
  )
  .cases <- c(
    "CB" = "Sham", 
    "CM" = "Sham", 
    "CS" = "Sham", 
    "DB" = "MCAO", 
    "DM" = "MCAO", 
    "DS" = "MCAO"
  )
  
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
  .sc <- Seurat::PercentageFeatureSet(
    object = .sc, 
    # pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", 
    pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",
    col.name = "percent.ribo"
  )
  apply(
    .sc@assays$RNA@counts,
    2,
    function(x) (100 * max(x)) / sum(x)
  ) ->
    .sc$Percent.Largest.Gene
  

  .vlnplot <- Seurat::VlnPlot(
    .sc, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
    ncol = 4
  )
  ggsave(
    filename = "vlnplot-{.project}.pdf" %>% glue::glue(),
    plot = .vlnplot,
    device = "pdf",
    path = "data/result/01-basic",
    width = 5,
    height = 5
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
  ggsave(
    filename = "qc-{.project}.pdf" %>% glue::glue(),
    plot = .plot,
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
  )
  
  .sc@meta.data %>% 
    dplyr::arrange(percent.mt) %>% 
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    geom_point() +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle("Example of plotting QC metrics") +
    geom_hline(yintercept = 500) +
    geom_hline(yintercept = 3000) +
    theme_bw() ->
    .metrics
  ggsave(
    plot = .metrics,
    filename = "{.project}-metrics.pdf" %>% glue::glue(),
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
  )
  
  .sc@meta.data %>% 
    dplyr::arrange(percent.ribo) %>% 
    ggplot(aes(nCount_RNA, nFeature_RNA, color = percent.ribo)) +
    geom_point() +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle("Example of plotting QC metrics") +
    geom_hline(yintercept = 500) +
    geom_hline(yintercept = 3000) +
    theme_bw() ->
    .metrics_ribo
  ggsave(
    plot = .metrics_ribo,
    filename = "{.project}-metrics-ribo.pdf" %>% glue::glue(),
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
  )
  
  .sc@meta.data %>% 
    ggplot(aes(percent.mt)) +
    geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage Mitochondrion") +
    geom_vline(xintercept = 10) +
    theme_bw() ->
    .percent_mt
  ggsave(
    plot = .percent_mt,
    filename = "{.project}-percent-mt.pdf" %>% glue::glue(),
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
  )
  
  .sc@meta.data %>% 
    ggplot(aes(percent.ribo)) +
    geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage Mitochondrion") +
    geom_vline(xintercept = 20) +
    theme_bw() ->
    .percent_ribo
  ggsave(
    plot = .percent_mt,
    filename = "{.project}-percent-ribo.pdf" %>% glue::glue(),
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
  )
  
  .sc@meta.data %>%
    ggplot(aes(Percent.Largest.Gene)) + 
    geom_histogram(binwidth = 0.7, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage Largest Gene") +
    geom_vline(xintercept = 20) +
    theme_bw() ->
    .percent_largest_gene
  ggsave(
    plot = .percent_largest_gene,
    filename = "{.project}-percent-largest-gene.pdf" %>% glue::glue(),
    device = "pdf",
    path = "data/result/01-basic",
    width = 7,
    height = 5
    )
  
  .sc
}

fn_filter_sct <- function(.sc) {
  .sc <- project_sc$sc[[2]]
  .sc_sub <- subset(
    x = .sc, 
    subset = nFeature_RNA > 500 & 
      nFeature_RNA < 6000 &
      percent.mt < 25 &
      percent.ribo < 50 &
      Percent.Largest.Gene < 40
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
    method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "percent.ribo", "CC.Difference"),
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
  file = "data/rda/project_sc.rds.gz"
)

project_sc %>% 
  dplyr::mutate(
    sct = purrr::map(
      .x = sc,
      .f = fn_filter_sct
    )
  ) %>% 
  dplyr::mutate(
    ratio = purrr::map2_dbl(
      .x = sc,
      .y = sct,
      .f = function(.x, .y) {
        ncol(.y)/ncol(.x)
      }
    )
  ) ->
  project_sc_sct

readr::write_rds(
  x = project_sc_sct, 
  file = "data/rda/project_sc_sct.rds.gz"
)

# save image --------------------------------------------------------------
save.image(file = "data/rda/01-sctransform.rda")
