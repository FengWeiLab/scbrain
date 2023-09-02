library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
library(ggplot2)

countdir <- "/home/liuc9/data/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/skull"
anno_file <- "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/annot_skull.csv"
ref.dir <- "/home/liuc9/data/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/azimuth_skull"

dirs <- list(
  "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/femur",
  "/mnt/isilon/xing_lab/liuc9/refdata/brainimmuneatlas/Mazzitelli_Nat_Neurosci_2022/skull"
)

names(dirs) <- c("femur", "skull")

matlist <- lapply(dirs, Read10X)
for(i in 1:length(matlist)) colnames(matlist[[i]]) <- paste0(names(matlist)[[i]], "__", colnames(matlist[[i]]))

so.list <- lapply(matlist, function(x){CreateSeuratObject(x, min.cells = 6, names.delim = '__')})


# Load count  -------------------------------------------------------------
annotations <- readr::read_csv(file = anno_file) |>
  dplyr::select(cell, cluster) |>
  dplyr::mutate(
    cluster = gsub(
      pattern = "[0-9]+_",
      replacement = "",
      x = cluster
    )
  ) |>
  tibble::deframe()

# counts <- Seurat::Read10X(
#   data.dir = countdir
# )

# sc <- merge(so.list[[1]], y = c(so.list[[2]]))


# transform ---------------------------------------------------------------

# sc <- Seurat::CreateSeuratObject(
#   counts = counts,
#   project = "Skull",
# )

sc <- merge(so.list[[1]], y = c(so.list[[2]]))

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

cells_anno <- annotations[cells]
names(cells_anno) <- cells

cells_anno <- factor(cells_anno)


Idents(object = sctu) <- cells_anno

# sctu

sctu <- RenameCells(
  object = sctu,
  new.names = unname(obj = sapply(
    X = Seurat::Cells(x = sctu),
    FUN = function(.s) {
      paste0("DLCJ", .s)
    }
  ))
)


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
