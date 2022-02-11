#' Make empty geometry
#'
#' For default arguments in some functions with the right type.
#'
#' @return An \code{sf} data frame whose geometry is an empty GEOMETRYCOLLECTION.
#' @importFrom sf st_sf st_sfc st_geometrycollection
#' @export
make_empty_geometry <- function() {
  st_sf(st_sfc(st_geometrycollection()))
}
