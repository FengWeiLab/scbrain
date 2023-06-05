# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun May  7 17:25:21 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)
library(Seurat)
library(Azimuth)


# args --------------------------------------------------------------------

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(dplyr::desc(cancer_types))
gtex <- yaml::yaml.load_file(input = "https://raw.githubusercontent.com/chunjie-sam-liu/sqtl/main/colors/color_gtexv8_tissue.yaml?token=GHSAT0AAAAAACDDSCQTTVKGIEH55YAEUH7MZD5ECZQ") |>
  tibble::enframe() |>
  tidyr::spread(key = name, value = value) |>
  tidyr::unnest()
# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------





# load data ---------------------------------------------------------------

refcelllevel <- tibble::tibble(
  region = c("Brain", "Meninge", "Skull"),
  ref = c("mousecortexref", "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_dura", "bonemarrowref"),
  celllevel = c("predicted.cluster", "predicted.annotation.l1", "predicted.celltype.l2"),
  supercelllevel = c("predicted.subclass", "predicted.annotation.l1", "predicted.celltype.l1"),
)




project_sc_azimuth <- readr::read_rds(
  # file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/project_sc_azimuth.rds"
  # file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/project_sc_azimuth_newref.rds"
  file = "/mnt/isilon/xing_lab/liuc9/projnet/2022-02-08-single-cell/azimuth/project_sc_azimuth_update.rds"

)

# body --------------------------------------------------------------------

{
  project_sc_azimuth |>
    dplyr::select(-ref) ->
    project_sc_azimuth_ref


  # Ratio -------------------------------------------------------------------

  project_sc_azimuth_ref |>
    dplyr::mutate(
      sunburst = purrr::pmap(
        .l = list(
          .anno = anno,
          .celllevel = celllevel,
          .supercelllevel = supercelllevel
        ),
        .f = function(.anno, .celllevel, .supercelllevel ) {
          # .anno <- project_sc_azimuth_ref$anno[[2]]
          # .celllevel <- project_sc_azimuth_ref$celllevel[[2]]
          # .supercelllevel <- project_sc_azimuth_ref$supercelllevel[[2]]

          .anno@meta.data |>
            dplyr::select(lowres = .supercelllevel, highres = .celllevel) |>
            dplyr::count(lowres, highres) ->
            .d

          .d |>
            dplyr::mutate(highres_r = n/sum(n)) |>
            dplyr::group_by(lowres) |>
            dplyr::mutate(lowres_r = sum(highres_r)) |>
            dplyr::ungroup() |>
            dplyr::mutate(
              lowres = glue::glue("{lowres} {round(lowres_r * 100, 2)}%"),
              highres =  glue::glue("{highres} {round(highres_r * 100, 2)}%")
            ) |>
            dplyr::select(1, 2, 3) |>
            plotme::count_to_sunburst(
              # fill_by_n = TRUE,
              # sort_by_n = TRUE
            )
        }
      )
    ) ->
    project_sc_azimuth_ref_sunburst


  # save image --------------------------------------------------------------



  project_sc_azimuth_ref_sunburst |>
    dplyr::mutate(
      a = purrr::pmap(
        .l = list(
          .region = region,
          .case = case,
          .p = sunburst
        ),
        .f = function(.region, .case, .p, .outdir) {
          # .region <- project_sc_azimuth_ref_sunburst$region[[1]]
          # .case <- project_sc_azimuth_ref_sunburst$case[[1]]
          # .p <- project_sc_azimuth_ref_sunburst$sunburst[[1]]

          reticulate::py_run_string("import sys")

          .filename <- glue::glue("Sunburst_{.region}_{.case}")
          plotly::save_image(
            p = .p,
            file = file.path(
              .outdir,
              glue::glue("{.filename}.pdf")
            ),
            width = 800,
            height = 800,
            device = "pdf"
          )

          htmlwidgets::saveWidget(
            .p,
            file = file.path(
              .outdir,
              glue::glue("{.filename}.html")
            )
          )

        },
        .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel4"
      )
    )

}




# ref sunburst ------------------------------------------------------------


SeuratData::LoadData(ds = "mousecortexref", type = "azimuth") -> mousecortexref
mousecortexref$plot@meta.data |>
  dplyr::select(class, subclass, cluster,) |>
  dplyr::arrange(class, subclass, cluster,) |>
  dplyr::distinct() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    predicted.cluster = cluster
  ) ->
  mousecortexref_cell

mousecortexref_cell |>
  dplyr::count(class, subclass, cluster) |>
  plotme::count_to_sunburst()


SeuratData::LoadData(ds = "pbmcref", type = "azimuth") -> pbmcref
pbmcref$plot@meta.data |>
  dplyr::select(celltype.l1, celltype.l2, ) |>
  dplyr::arrange(celltype.l1, celltype.l2, ) |>
  dplyr::distinct() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    predicted.celltype.l2 = celltype.l2
  ) ->
  pbmcref_cell

pbmcref_cell |>
  dplyr::count(celltype.l1, celltype.l2, ) |>
  plotme::count_to_sunburst()


SeuratData::LoadData(ds = "bonemarrowref", type = "azimuth") -> bonemarrowref
bonemarrowref$plot@meta.data |>
  dplyr::select(celltype.l1, celltype.l2) |>
  dplyr::arrange(celltype.l1, celltype.l2) |>
  dplyr::distinct() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    predicted.celltype.l2 = celltype.l2
  ) ->
  bonemarrowref_cell

bonemarrowref_cell |>
  dplyr::count(celltype.l1, celltype.l2) |>
  plotme::count_to_sunburst()


refcelllevel <- tibble::tibble(
  region = c("Brain", "Meninge", "Skull"),
  refs = c("mousecortexref", "/home/liuc9/data/refdata/brainimmuneatlas/azimuth_dura", "bonemarrowref"),
  celllevel = c("celltype2", "predicted.annotation.l1", "predicted.celltype.l2"),
  supercelllevel = c("celltype1", "predicted.annotation.l1", "predicted.celltype.l1")
) |>
  dplyr::mutate(
    realcell = c(list(mousecortexref_cell), list(pbmcref_cell), list(bonemarrowref_cell))
  )


project_sc_azimuth |>
  dplyr::select(-ref, -refs, -celllevel, -supercelllevel) |>
  dplyr::left_join(
    refcelllevel,
    by = "region"
  ) ->
  project_sc_azimuth_ref_realcell


project_sc_azimuth_ref_realcell |>
  dplyr::mutate(
    sunburst = purrr::pmap(
      .l = list(
        .anno = anno,
        .celllevel = celllevel,
        .realcell = realcell,
        .region = region
      ),
      .f = function(.anno, .celllevel, .realcell, .region ) {
        # .anno <- project_sc_azimuth_ref_realcell$anno[[9]]
        # .celllevel <- project_sc_azimuth_ref_realcell$celllevel[[9]]
        # .realcell <- project_sc_azimuth_ref_realcell$realcell[[9]]
        # .region <- project_sc_azimuth_ref_realcell$region[[9]]



        .d <- if (.region == "Meninge") {
          .anno@meta.data |>
            dplyr::mutate(
              celltype.l1 = predicted.annotation.l1,
              celltype.l2 = predicted.annotation.l1
            )
        } else if (.region == "Brain") {
          .anno@meta.data |>
            dplyr::mutate(
              celltype.l1 = celltype1,
              celltype.l2 = celltype2
            )
        } else {
          .anno@meta.data |>
            dplyr::left_join(
              .realcell,
              by = .celllevel
            )
        }

        # .dd <- if(.region == "Brain") {
        #   nnc <- c("L6 IT_1", "Meis2", "Peri", "Meis2_Top2a")
        #   .d |>
        #     dplyr::filter(!cluster %in% nnc) |>
        #     dplyr::count(class, subclass, cluster) |>
        #     dplyr::mutate(cluster_r = n / sum(n)) |>
        #     dplyr::group_by(subclass) |>
        #     dplyr::mutate(subclass_r = sum(cluster_r)) |>
        #     dplyr::ungroup() |>
        #     dplyr::group_by(class) |>
        #     dplyr::mutate(class_r = sum(cluster_r)) |>
        #     dplyr::ungroup() |>
        #     dplyr::mutate(
        #       class = glue::glue("{class} {round(class_r * 100, 2)}%"),
        #       subclass = glue::glue("{subclass} {round(subclass_r * 100, 2)}%"),
        #       cluster = glue::glue("{cluster} {round(cluster_r * 100, 2)}%")
        #     ) |>
        #     dplyr::select(2,3,4)
        # } else {
        #   .d |>
        #     dplyr::count(celltype.l1, celltype.l2) |>
        #     dplyr::mutate(celltype.l2_r = n / sum(n)) |>
        #     dplyr::group_by(celltype.l1) |>
        #     dplyr::mutate(celltype.l1_r = sum(celltype.l2_r)) |>
        #     dplyr::ungroup() |>
        #     dplyr::mutate(
        #       celltype.l1 = glue::glue("{celltype.l1} {round(celltype.l1_r * 100, 2)}%"),
        #       celltype.l2 = glue::glue("{celltype.l2} {round(celltype.l2_r * 100, 2)}%"),
        #     ) |>
        #     dplyr::select(1,2,3)
        # }
        .d |>
          dplyr::count(celltype.l1, celltype.l2) |>
          dplyr::mutate(celltype.l2_r = n / sum(n)) |>
          dplyr::group_by(celltype.l1) |>
          dplyr::mutate(celltype.l1_r = sum(celltype.l2_r)) |>
          dplyr::ungroup() |>
          dplyr::mutate(
            celltype.l1 = glue::glue("{celltype.l1} {round(celltype.l1_r * 100, 2)}%"),
            celltype.l2 = glue::glue("{celltype.l2} {round(celltype.l2_r * 100, 2)}%"),
          ) |>
          dplyr::select(1,2,3) ->
          .dd

        .dd |>
          plotme::count_to_sunburst() ->
          .p

        .new_anno <- .anno
        .new_anno@meta.data$celltype.l1 <- .d$celltype.l1
        .new_anno@meta.data$celltype.l2 <- .d$celltype.l2

        tibble::tibble(
          sunburst = list(.p),
          anno_new = list(.new_anno)
        )

      }
    )
  ) |>
  tidyr::unnest(cols = sunburst) ->
  project_sc_azimuth_ref_realcell_sunburst



# save image --------------------------------------------------------------



project_sc_azimuth_ref_realcell_sunburst |>
  dplyr::mutate(
    a = purrr::pmap(
      .l = list(
        .region = region,
        .case = case,
        .p = sunburst
      ),
      .f = function(.region, .case, .p, .outdir) {

        dir.create(path = .outdir, showWarnings = F, recursive = TRUE)

        reticulate::py_run_string("import sys")

        .filename <- glue::glue("Sunburst_{.region}_{.case}")
        plotly::save_image(
          p = .p,
          file = file.path(
            .outdir,
            glue::glue("{.filename}.pdf")
          ),
          width = 800,
          height = 800,
          device = "pdf"
        )

        htmlwidgets::saveWidget(
          .p,
          file = file.path(
            .outdir,
            glue::glue("{.filename}.html")
          )
        )

      },
      .outdir = "/home/liuc9/github/scbrain/scuvresult/06-azimuth-celllevel10"
    )
  )

readr::write_rds(
  x = project_sc_azimuth_ref_realcell_sunburst,
  file = "/home/liuc9/github/scbrain/data/azimuth/project_sc_azimuth_ref_realcell_sunburst.rds"
)

# Diff --------------------------------------------------------------------



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
save.image(file = "data/azimuth/03-individual-tissue-supercell.rda")
