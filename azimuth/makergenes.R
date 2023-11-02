d <- readxl::read_xlsx(
  path = "/mnt/isilon/xing_lab/liuc9/refdata/cellmarker/Cell_marker_Mouse.xlsx"
)



d 

d |> 
  dplyr::count(tissue_class) |> 
  dplyr::arrange(-n) |> 
  print(n = Inf)

d |> 
  dplyr::select(cellontology_id, cell_name, Symbol) |> 
  dplyr::filter(!is.na(cellontology_id)) |> 
  dplyr::filter(!is.na(Symbol)) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cellontology_id) ->
  m

m |> 
  dplyr::count(cellontology_id, cell_name) |> 
  dplyr::arrange(-n) |> 
  print(n = Inf)

d |> 
  dplyr::filter(cellontology_id == "CL_0000127")

d |> dplyr::glimpse()
d |> 
  dplyr::count(cell_name)

d |> 
  dplyr::filter(
    tissue_class %in% c(
      "Brain",
      "Bone marrow",
      "Blood",
      "Blood vessel",
      "Bone",
      "Lymph node",
      "Nerve",
      "Lymphoid tissue",
      "Epithelium"
    )
  ) ->
  m


cellmarker <- readxl::read_xlsx(
  path = "/mnt/isilon/xing_lab/liuc9/refdata/cellmarker/Cell_marker_Mouse.xlsx"
)
cellmarker_brain <- cellmarker |>
  dplyr::filter(
    tissue_class %in% c(
      "Brain",
      "Bone marrow",
      "Blood",
      "Lymph node",
      "Nerve",
      "Lymphoid tissue",
      "Epithelium"
    )
  ) 

cellmarker_immune <- cellmarker |>
  dplyr::filter(
    tissue_class %in% c(
      "Bone marrow",
      "Blood",
      "Lymph node",
      "Lymphoid tissue",
      "Epithelium"
    )
  )


cellmarker_immune |> 
  dplyr::select(cell_name) |> 
  dplyr::arrange(cell_name) |> 
  dplyr::distinct() |> 
  # print(n = Inf)
  dplyr::filter(grepl("neuron", x = cell_name, ignore.case = T))
  

# skull -------------------------------------------------------------------

dplyr::bind_rows(
cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "T cell", x = cell_name)) |> 
  dplyr::mutate(
    cell3 =  "T cells"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "B cell", x = cell_name)) |> 
  dplyr::mutate(cell3 = "B cells"),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "B cell", x = cell_name)) |> 
  dplyr::mutate(cell3 = "Mature B cells"),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "B cell", x = cell_name)) |> 
  dplyr::mutate(cell3 = "Immature B cells"),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "B cell", x = cell_name)) |> 
  dplyr::mutate(cell3 = "Pre-B cells"),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "B cell", x = cell_name)) |> 
  dplyr::mutate(cell3 = "Pro-B cells"),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Natural killer", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "NK cells"
  ),


cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Monocyte", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Monocytes"
  ),
cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Monocyte", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "CD14 Monocytes"
  ),
cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Monocyte", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "CD16 Monocytes"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Macrophage|macrophage", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Macrophage"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Dendritic|dendritic", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "DC"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Dendritic|dendritic", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "mDC"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Dendritic|dendritic", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "pDC"
  ),



cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Neutrophil|neutrophil", x = cell_name)) |>
  dplyr::mutate(
    cell3 = "Neutrophils"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Hematopoietic", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "HSPC"
  ),

cellmarker_immune |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Erythroid", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Erythroid"
  ),



cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "sensory", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Sensory neurons"
  )) |> 
  dplyr::select(cell3, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::filter(!is.na(Symbol)) |> 
  dplyr::filter(!grepl(
    "H2-",
    Symbol
  )) ->
  acells


# Brain -------------------------------------------------------------------

cellmarker_brain |> 
  dplyr::select(cell_name) |> 
  dplyr::arrange(cell_name) |> 
  dplyr::distinct() |> 
  # print(n = Inf)
  dplyr::filter(grepl("vascular", x = cell_name, ignore.case = T))

astro <- tibble::tibble(
  cell_name = "Astro",
  Symbol = c("Gpc5", "Myoc", "Slc1a2", "Fxyd6", "Slc1a3", "Nrp2", "Slc6a6", "Apoe", "Aqp4", "Gfap", "Lsamp", "Gpc5", "Slc1a2", "Luzp2", "Slc1a3", "Apoe", "Cadm2", "Rora", "Aqp4", "Slc7a10", "Igfbpl1", "Gpc5", "Mki67", "Slc1a2", "Myt1", "Slc1a3", "Apoe", "E2f2", "E2f7", "Top2a") |> unique()
)


dplyr::bind_rows(
 
# cellmarker_brain |> 
#   dplyr::select(cell_name, Symbol) |> 
#   dplyr::distinct() |> 
#   dplyr::arrange(cell_name) |> 
#   dplyr::filter(grepl(pattern = "Astrocyte|astrocyte", x = cell_name)) |> 
#   dplyr::mutate(
#     cell3 = "Astrocyte Aqp4_Slc7a10"
#   ),
# 
# cellmarker_brain |> 
#   dplyr::select(cell_name, Symbol) |> 
#   dplyr::distinct() |> 
#   dplyr::arrange(cell_name) |> 
#   dplyr::filter(grepl(pattern = "Astrocyte|astrocyte", x = cell_name)) |> 
#   dplyr::mutate(
#     cell3 = "Astrocyte Aqp4_Gfap"
#   ),
  
  astro |> 
      dplyr::mutate(
        cell3 = "Astrocyte Aqp4_Slc7a10"
      ),

  astro |> 
      dplyr::mutate(
        cell3 = "Astrocyte Aqp4_Gfap"
      ),
  
cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Microglia|microglia", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Microglia"
  ),


cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Oligodendrocyte", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "OLG"
  ),

cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Oligodendrocyte", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "OPC"
  ),
cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(grepl(pattern = "Endothelial|endothelial", x = cell_name)) |> 
  dplyr::mutate(
    cell3 = "Endothelial"
  )) |> 
  dplyr::select(cell3, Symbol) |> 
  dplyr::filter(!is.na(Symbol)) |> 
  dplyr::distinct() ->
  bcells
  

vlmcgenes <- c("Cubn", "Ptgds", "Bnc2", "Slc7a14", "Sptssb", "Cped1", "Slc7a11", "Dapl1", "Igfbp3", "Bmp6", "Ptgds", "Prg4", "Bnc2", "Tspan8", "Cped1", "Msln", "Slc7a11", "Dapl1", "Bmp6", "Fn1", "Slc47a1", "Ptgds", "Bnc2", "Slc4a10", "Cped1", "Pcdh7", "Slc7a11", "Tmeff2", "Erc2", "Bmp6", "Ptgds", "Sh3gl3", "Bnc2", "Aifm3", "Ppp1r1a", "Cped1", "Gm27151", "Slc7a11", "Atp13a5", "Bmp6", "Pkhd1", "Ptgds", "Rspo2", "Bnc2", "Cftr", "Cped1", "Slc7a11", "Gm2115", "Bmp6", "Atp13a5", "Aox3", "Ptgds", "Bnc2", "Plxna2", "Slc7a11", "Cped1", "Prex2", "Vcam1", "Bmp6", "Coch", "Ptgds", "Gm30624", "Bnc2", "Arhgap22", "Cped1", "Gm47865", "Slc7a11", "Col15a1", "Bmp6", "Adgrl3") |> unique()

cellmarker_brain |> 
  dplyr::select(cell_name, Symbol) |> 
  dplyr::distinct() |> 
  dplyr::arrange(cell_name) |> 
  dplyr::filter(Symbol %in% vlmcgenes) |> 
  dplyr::mutate(
    cell3 = "VLMC"
  ) |> 
  dplyr::select(cell3, Symbol) |> 
  dplyr::filter(!is.na(Symbol)) |> 
  dplyr::distinct() ->
  ccells

dplyr::bind_rows(
  acells, bcells, ccells
) |> 
  dplyr::group_by(cell3) |> 
  tidyr::nest() |> 
  dplyr::ungroup() ->
  candidate_markers

candidate_markers |> 
  # dplyr::filter(cell3 == "mDC")
  tidyr::unnest(cols = data) |> 
  dplyr::filter(Symbol == "H2-Ab1")
  # print(n = Inf)
  dplyr::filter(Symbol == "Aqp1")

readr::write_rds(
  candidate_markers,
  "/home/liuc9/data/refdata/cellmarker/candidate_markers.rds.gz"
)


