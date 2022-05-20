

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