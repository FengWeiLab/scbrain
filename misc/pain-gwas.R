# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Jul  7 16:06:19 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(prismatic)
#library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


# body --------------------------------------------------------------------

conn <- DBI::dbConnect(
  drv = duckdb::duckdb(),
  "/home/liuc9/data/refdata/gwas/gwas_catalog_v1.0.2-associations_e107_r2022-10-08-ftp-efo.duckdb"
)

DBI::dbListTables(conn)

gwas <- dplyr::tbl(
  conn,
  "gwas_asso_efo_e107_r2022-10-08"
) |>
  dplyr::collect()


DBI::dbDisconnect(conn, shutdown = T)


gwas |> dplyr::glimpse()


gwas |>
  dplyr::filter(
    grepl(
      pattern = "pain",
      x = `EFO term`
    )
  ) |>
  writexl::write_xlsx(
    path = "/home/liuc9/github/scbrain/scuvresult/pain/gwas-pain.xlsx"
  )


# Disease Net -------------------------------------------------------------

da <- vroom::vroom(
  file = "/mnt/isilon/xing_lab/liuc9/refdata/disgenet/disease_associations.tsv.gz"
)



ga <- vroom::vroom(
  file = "/mnt/isilon/xing_lab/liuc9/refdata/disgenet/gene_associations.tsv.gz"
)

da |>
  dplyr::filter(
    grepl(
      pattern = "pain",
      x = diseaseName,
      ignore.case = T
    )
  )

conn <- DBI::dbConnect(
  drv = RSQLite::SQLite(),
  "/mnt/isilon/xing_lab/liuc9/refdata/disgenet/disgenet_2020.db"
)

DBI::dbListTables(conn)

dplyr::tbl(
  conn,
  "geneDiseaseNetwork"
) |>
  dplyr::glimpse()
dplyr::tbl(
  conn,
  "diseaseClass"
) |>
  dplyr::glimpse()
dplyr::tbl(
  conn,
  "disease2class"
) |>
  dplyr::glimpse()

dplyr::tbl(
  conn,
  "diseaseAttributes"
) |>
  dplyr::glimpse()

dplyr::tbl(
  conn,
  "geneDiseaseNetwork"
) |>
  dplyr::left_join(
    dplyr::tbl(
      conn,
      "diseaseAttributes"
    ),
    by = "diseaseNID"
  ) |>
  dplyr::collect() ->
  disgenet

disgenet |>
  dplyr::filter(
    grepl(
      pattern = "pain",
      x = diseaseName,
      ignore.case = T
    )
  ) |>
  writexl::write_xlsx(
    "/home/liuc9/github/scbrain/scuvresult/pain/disgenet-pain.xlsx"
  )

DBI::dbDisconnect(
  conn,
  shutdown = T
)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
