# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Feb  6 21:06:20 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
sc_sham_mcao_uv_scn_integrated <- readr::read_rds(
  file = "data/scuvrda/sc_sham_mcao_uv_scn_integrated.rds.gz"
)

Seurat::DefaultAssay(sc_sham_mcao_uv_scn_integrated) <- "integrated"

# body --------------------------------------------------------------------

sc_scale <- sc_sham_mcao_uv_scn_integrated %>% 
  Seurat::ScaleData() %>% 
  Seurat::RunPCA(
    npcs = 10,
    verbose = FALSE
  )


sc_scale@reductions$pca@cell.embeddings %>% 
  as.data.frame() %>% 
  dplyr::bind_cols(sc_scale@meta.data) ->
  v_m





v_m %>% 
  dplyr::mutate(
    batch = plyr::revalue(
      x = tissue,
      replace = c(
        "CB" = "Seq1",
        "DB" = "Seq1",
        "CM" = "Seq2",
        "DM" = "Seq2",
        "CS" = "Seq2",
        "DS" = "Seq2",
        "UVB" = "Seq3",
        "UVM" = "Seq3",
        "UVS" = "Seq3"
      )
    )
  ) ->
  v_m_d

v_m_d %>% 
  dplyr::group_by(batch) %>% 
  dplyr::count()

v_m_d %>% 
  ggplot(aes(
    x = PC_1,
    y = PC_2
  )) +
  geom_point(
    aes(
      color = batch
    )
  ) +
  ggsci::scale_color_jco(
    name = "Batch seq"
  ) +
  theme_classic() ->
  p

# writexl::write_xlsx(
#   x = project_stat, 
#   path = "data/scuvresult/01-basic/reads-stat.xlsx"
# )

ggsave(
  filename = "before-batch-removal.pdf",
  plot = p,
  device = "pdf",
  path = "data/scuvresult/02-batch/",
  width = 6,
  height = 5
)


library(plotly)

plot_ly(
  v_m_d,
  x = ~PC_1,
  y = ~PC_2,
  z = ~PC_3,
  color = ~batch
)

sc_sham_mcao_uv_scn_integrated$batch <- 
  plyr::revalue(
    x = v_m$tissue,
    replace = c(
      "CB" = "Seq1",
      "DB" = "Seq1",
      "CM" = "Seq2",
      "DM" = "Seq2",
      "CS" = "Seq2",
      "DS" = "Seq2",
      "UVB" = "Seq3",
      "UVM" = "Seq3",
      "UVS" = "Seq3"
    )
  )

sc_scale1 <- sc_sham_mcao_uv_scn_integrated %>% 
  Seurat::ScaleData() %>% 
  Seurat::RunPCA(
    npcs = 30,
    verbose = FALSE
  ) %>% 
  harmony::RunHarmony(
    group.by.vars = "batch"
  )




sc_scale1@reductions$harmony@cell.embeddings %>%  
  as.data.frame() %>% 
  dplyr::bind_cols(sc_scale1@meta.data) ->
  v_m1

v_m1 %>% 
  ggplot(aes(
  x = harmony_1,
  y = harmony_2
)) +
  geom_point(
    aes(
      color = batch
    )
  ) +
  ggsci::scale_color_jco(
    name = "Batch seq"
  ) +
  theme_classic() ->
  p

# writexl::write_xlsx(
#   x = project_stat, 
#   path = "data/scuvresult/01-basic/reads-stat.xlsx"
# )

ggsave(
  filename = "after-batch-removal.pdf",
  plot = p,
  device = "pdf",
  path = "data/scuvresult/02-batch/",
  width = 6,
  height = 5
)



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------