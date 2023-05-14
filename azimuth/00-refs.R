library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
library(ggplot2)
# args <- commandArgs(trailingOnly = TRUE)



# Dura --------------------------------------------------------------------


# countdir <- "/home/liuc9/data/refdata/brainimmuneatlas/filtered_gene_bc_matrices_dura/mm10"
# anno_dura_file <- "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/annot_dura.csv"
# ref.dir <- "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_dura"
#

# enrSDM --------------------------------------------------------------------


# countdir <- "/home/liuc9/data/refdata/brainimmuneatlas/filtered_gene_bc_matrices_enrSDM/mm10"
# anno_dura_file <- "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/annot_enrSDM.csv"
# ref.dir <- "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_enrSDM"


# wholeBrain -------------------------------------------------------------------

countdir <- "/home/liuc9/data/refdata/brainimmuneatlas/filtered_gene_bc_matrices_wholeBrain/mm10"
anno_dura_file <- "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/annot_wholeBrain.csv"
ref.dir <- "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_wholeBrain"


annotations <- readr::read_csv(
  file = anno_dura_file
) |>
  dplyr::select(cell, cluster) |>
  tibble::deframe()

annotations_names <- glue::glue("{names(annotations)}-1")
names(annotations) <- annotations_names

counts <- Seurat::Read10X(data.dir = countdir)

sc <- Seurat::CreateSeuratObject(
  counts = counts,
  project = "Dura",
)

sct <- Seurat::SCTransform(object = sc)


sct |>
  Seurat::RunPCA() |>
  Seurat::RunUMAP(dims = 1:30, return.model = TRUE) |>
  Seurat::FindNeighbors(dims = 1:30) |>
  Seurat::FindClusters(resolution = c(0.8, 2)) ->
  sctu

annotations
Idents(sctu) |>
  names() ->
  cells


cells_anno <- annotations[cells] |> tidyr::replace_na("remove")
names(cells_anno) <- cells

cells_anno <- factor(cells_anno)


Idents(object = sctu) <- cells_anno

ref <- sctu

if ("remove" %in% levels(x = ref)) {
  ref <- subset(x = ref, idents = "remove", invert = TRUE)
  ref <- RunPCA(object = ref, verbose = FALSE)
}
ref$annotation.l1 <- Idents(object = ref)
ref <- RunUMAP(object = ref, dims = 1:30, return.model = TRUE)
full.ref <- ref
colormap <- list(annotation.l1 = CreateColorMap(object = ref, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)


SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))
saveRDS(object = full.ref, file = file.path(ref.dir, "fullref.Rds"))
