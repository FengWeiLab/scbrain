library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

args_1 <- "/mnt/isilon/xing_lab/liuc9/refdata/azimuth_reference/pbmc_multimodal.h5seurat"
args_2 <- "/mnt/isilon/xing_lab/liuc9/refdata/azimuth_reference/pbmc_multimodal_reference_plotting_cells.txt"
args_3 <- "/mnt/isilon/xing_lab/liuc9/refdata/azimuth_reference/pbmc_multimodal_reference_mapping_cells.txt"


ref.dir <- "/home/liuc9/data/refdata/azimuth_reference"


ob <- LoadH5Seurat(
  file = args_1,
  assays = list(
    "SCT" = c("counts"),
    "ADT" = c("counts", "data")
  )
)

mapping.cells <- readLines(con = args[2])
plotting.cells <- readLines(con = args[3])
mapping.ob <- subset(x = ob, cells = mapping.cells)
plotref <- subset(x = ob[["wnn.umap"]], cells = plotting.cells)
set.seed(28)
ids <- sample(
  x = levels(as.factor(mapping.ob$celltype.l2)),
  size = length(levels(as.factor(mapping.ob$celltype.l2)))
)
colormap.ids <- hue_pal()(n = length(x = ids))
names(x = colormap.ids) <- ids
colormap.ids <- colormap.ids[sort(x = names(x = colormap.ids))]

colormap <- list(
  celltype.l1 = CreateColorMap(ids = unique(x = ob$celltype.l1)),
  celltype.l2 = colormap.ids
)

ob.ref <- AzimuthReference(
  object = mapping.ob,
  refUMAP = "wnn.umap",
  refDR = "spca",
  refAssay = "SCT",
  plotref = plotref,
  assays = "ADT",
  metadata = c("celltype.l1", "celltype.l2", "celltype.l3"),
  dims = 1:50,
  k.param = 31,
  plot.metadata = data.frame(
    celltype.l1 = ob$celltype.l1[plotting.cells],
    celltype.l2 = ob$celltype.l2[plotting.cells],
    celltype.l3 = ob$celltype.l3[plotting.cells]
  ),
  ori.index = match(x = Cells(x = mapping.ob), table = Cells(x = ob)),
  colormap = colormap,
  reference.version = "1.0.0"
)

ValidateAzimuthReference(object = ob.ref)
SaveAnnoyIndex(object = ob.ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ob.ref, file = file.path(ref.dir, "ref.Rds"))



d <- readr::read_rds(
  file = "/mnt/isilon/xing_lab/liuc9/refdata/azimuth_reference/panc/data/panc8_3.0.2.rds"
)
