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

.getfun <- function(MARGIN) {
  switch (MARGIN,
          1 = int_elementMetadata,
          2 = int_colData
  )
}
.setfun <- function(MARGIN) {
  switch (MARGIN,
          1 = `int_elementMetadata<-`,
          2 = `int_colData<-`
  )
}
.xdimstr <- function(MARGIN) {
  switch (MARGIN,
          1 = "nrow",
          2 = "ncol"
  )
}
.xdimfun <- function(MARGIN) {
  switch (MARGIN,
    1 = nrow,
    2 = ncol
  )
}
.dg_key <- function(MARGIN) {
  switch (MARGIN,
    1 = "rowGeometries",
    2 = "colGeometries"
  )
}
.unnamed <- "unnamed"
# Modified from SCE to generalize to both rows and columns
.check_dimgeo_names <- function(reference, incoming, MARGIN, withDimnames,
                                fun='dimGeometry', vname='value') {
  if (!is.null(incoming)) {
    rni <- rownames(incoming)
    cnr <- dimnames(reference)[[MARGIN]]
    fun_show <- switch (MARGIN,
      1 = "rownames",
      2 = "colnames"
    )
    if (withDimnames && !is.null(rni)) {
      if (!identical(cnr, rni)) {
        msg <- paste0("non-NULL 'rownames(", vname, ")' should be the same as '",
                      fun_show, "(x)' for '", fun,
                      "<-'. This will be an error in the next release of Bioconductor.")
        warning(paste(strwrap(msg), collapse="\n"))
      }
    }
  }
  incoming
}
.get_internal_all <- SingleCellExperiment:::.get_internal_all
.set_internal_all <- SingleCellExperiment:::.set_internal_all
.get_internal_integer <- SingleCellExperiment:::.get_internal_integer
.get_internal_names <- SingleCellExperiment:::.get_internal_names
.get_internal_missing <- SingleCellExperiment:::.get_internal_missing
.get_internal_character <- SingleCellExperiment:::.get_internal_character
.set_internal_names <- SingleCellExperiment:::.set_internal_names
.set_internal_missing <- SingleCellExperiment:::.set_internal_missing
.set_internal_numeric <- SingleCellExperiment:::.set_internal_numeric
.set_internal_character <- SingleCellExperiment:::.set_internal_character
