#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon Nov 20 15:45:01 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
library(data.table)
#library(rlang)
library(Seurat)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich <-
  readr::read_rds(
    file = "/home/liuc9/github/scbrain/data/azimuth/azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich.rds.gz"
  )
# body --------------------------------------------------------------------
the_case_color <- tibble::tibble(
  case = c("Sham", "MCAO", "UV"),
  label = c("Sham", "tMCAO", "tMCAO+UVB"),
  color = c("#1B1919FF", "#0099B4FF", "#FDAF91FF") |> prismatic::color()
)

# from cell paper.
pro_inflammatory_score_genes <- c("Il6", "Il1a", "Il1b", "Ifng", "Il11", "Il7d", "Il7f", "Il18", "Tnf")

anti_inflammatory_score_genes <- c("Il1rn", "Tgfb1", "Il4", "Il10", "Il12a", "Il13")

#

anti_inflammatory_score_genes <- c("Arg1", "Fgl2", "Adora2b", "Ccl2", "Il4ra", "Msr1", "Il18", "Stat3", "Il1rm", "Spp1", "Mrc1", "Il18bp", "Lgal3", "Tgfb1", "Il4")

brain_repair <- c("Ccr5", "Tgm2", "Kdr", "Vegfa", "Ptprz1", "Fgf1", "Inhbb", "Lif", "Timp1", "Gdnf", "Lifr", "Il1b", "Ptn", "Gm", "Cspg4")

astrocyte_activation <- readr::read_tsv(
  "/scr1/users/liuc9/tmp/GO_term_summary_20231120_165348.txt"
) |>
  dplyr::filter(`Annotated Term` == "astrocyte activation involved in immune response") |>
  dplyr::pull(Symbol) |>
  unique()

angiogenesis <- readr::read_tsv(
  "/scr1/users/liuc9/tmp/GO_term_summary_20231120_165919.txt"
)


angiogenesis_pos <- angiogenesis |>
  dplyr::filter(`Annotated Term` == "positive regulation of vascular wound healing") |>
  dplyr::pull(Symbol) |>
  unique()

angiogenesis_neg <- angiogenesis |>
  dplyr::filter(`Annotated Term` == "negative regulation of vascular wound healing") |>
  dplyr::pull(Symbol) |>
  unique()

inflam <- readr::read_tsv(
  "/scr1/users/liuc9/tmp/GO_term_summary_20231121_153547.txt"
)

inflam_gene <- inflam |>
  dplyr::filter(`Annotated Term` == "negative regulation of inflammatory response to wounding") |>
  dplyr::pull(Symbol) |>
  unique()

immune <- readr::read_tsv(
  "/home/liuc9/tmp/GO_term_summary_20231121_154452.txt"
)

immune_gene <- immune |>
  dplyr::filter(`Annotated Term` == "innate immune response") |>
  dplyr::pull(Symbol) |>
  unique()

azimuth_ref_sunburst_cell_merge_norm_de_change_nn_volcano_enrich$norm[[3]] -> b

a <- UCell::AddModuleScore_UCell(
  b,
  features = list(
    angiogenesis_pos = angiogenesis_pos,
    anti_inflammatory_score_genes = anti_inflammatory_score_genes,
    pro_inflammatory_score_genes = pro_inflammatory_score_genes,
    brain_repair = brain_repair,
    inflam_gene = inflam_gene,
    immune_gene = immune_gene
  )
)


# a <- Seurat::AddModuleScore(
#   b,
#   features = list(
#     anti_inflammatory_score_genes = anti_inflammatory_score_genes,
#     pro_inflammatory_score_genes,
#     brain_repair,
#     astrocyte_activation,
#     angiogenesis,
#     angiogenesis_pos,
#     angiogenesis_neg
#   ),
#   assay = "RNA",
#   name = c("anti", "pro", "repair", "aa", "angio", "angio_pos", "angio_neg"),
#   ctrl = 5
# )

a@meta.data |> colnames()
a@meta.data |>
  dplyr::filter(cell2 == "Astrocyte") |>
  ggplot(
    aes(
      x = case,
      y = anti_inflammatory_score_genes_UCell
    )
  ) +
  geom_boxplot()

a@meta.data |>
  dplyr::filter(cell2 == "Astrocyte") |>
  ggplot(
    aes(
      x = case,
      y = pro_inflammatory_score_genes_UCell
    )
  ) +
  geom_boxplot()

a@meta.data |>
  dplyr::filter(cell2 == "Astrocyte") |>
  ggplot(
    aes(
      x = case,
      y = immune_gene_UCell
    )
  ) +
  geom_boxplot()

a@meta.data -> d

ggpubr::ggboxplot(
  d, x = "case", y = "immune_gene_UCell",
  color = "case",
  # add = "jitter",
  order = c("Sham", "MCAO", "UV")
) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("MCAO", "Sham"),
      c("UV", "MCAO")
    )
  ) +
  ggpubr::stat_compare_means(label.y = 0.4)

ggpubr::ggviolin(
  d |> dplyr::filter(cell2 == "Astrocyte"),
  x = "case", y = "angiogenesis_pos_UCell",
  fill = "case",
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  # color = "case",
  add = "boxplot",
  add.params = list(fill = "white"),
  order = c("Sham", "MCAO", "UV"),
) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("MCAO", "Sham"),
      c("UV", "MCAO")
    ),
    label = "p.signif"
  ) +
  ggpubr::stat_compare_means(label.y = 0.4)


ggpubr::ggviolin(
  d |> dplyr::filter(cell3 == "Astrocyte Aqp4_Slc7a10"),
  x = "case", y = "angiogenesis_pos_UCell",
  fill = "case",
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  # color = "case",
  add = "boxplot",
  add.params = list(fill = "white"),
  order = c("Sham", "MCAO", "UV"),
) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("MCAO", "Sham"),
      c("UV", "MCAO")
    ),
    # label = "p.signif"
  ) +
  ggpubr::stat_compare_means(label.y = 0.4)





ggpubr::ggviolin(
  d |> dplyr::filter(cell3 == "Astrocyte Aqp4_Gfap"),
  x = "case", y = "angiogenesis_pos_UCell",
  fill = "case",
  palette = c("#1B1919FF", "#0099B4FF", "#FDAF91FF"),
  # color = "case",
  add = "boxplot",
  add.params = list(fill = "white"),
  order = c("Sham", "MCAO", "UV"),
) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("MCAO", "Sham"),
      c("UV", "MCAO")
    ),
    # label = "p.signif"
  ) +
  ggpubr::stat_compare_means(label.y = 0.4)

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
