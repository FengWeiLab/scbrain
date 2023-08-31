library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
library(ggplot2)

countdir <- "/home/liuc9/data/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/skull"
anno_file <- ""
ref.dir <- "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_skull"

# Load count  -------------------------------------------------------------


counts <- Seurat::Read10X(
  data.dir = countdir
)

# transform ---------------------------------------------------------------

sc <- Seurat::CreateSeuratObject(
  counts = counts,
  project = "Skull",
)

sct <- Seurat::SCTransform(object = sc)

sct |>
  Seurat::RunPCA() |>
  Seurat::RunUMAP(dims = 1:30, return.model = TRUE) |>
  Seurat::FindNeighbors(dims = 1:30) |>
  Seurat::FindClusters(resolution = c(0.8, 2)) ->
  sctu
