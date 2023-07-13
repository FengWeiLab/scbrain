paletteer::paletteer_d(
  "ggsci::schwifty_rickandmorty"
)[c(-9)] ->
  brain_region

paletteer::paletteer_d(
  "ggsci::lanonc_lancet"
)[c(8, 1, 2)] ->
  case
names(case) <- c("C", "A", "I")

case |>
  tibble::enframe()
