#' Simple geometry predicates
#'
#' Unlike functions in \code{sf} like \code{st_intersects}, this function simply
#' returns a logical vector indicating whether each geometry in \code{x}
#' intersects (or returns \code{TRUE} from other predicates) anything in
#' \code{y}, preferably when \code{y} only contains a small number of geometries
#' or is one single MULTI geometry. This is useful when cropping or subsetting
#' an SFE object with a geometry, such as tissue boundary or histological region
#' polygons or a bounding box.
#'
#' @param x An object of class \code{sf}, \code{sfc}, or \code{sfg}.
#' @param y Another object of class \code{sf}, \code{sfc}, or \code{sfg}.
#' @param pred A geometric binary predicate function, such as \code{\link{st_intersects}}.
#' It should return an object of class \code{sgbp}, for sparse predicates.
#' @return A logical vector indicating whether each geometry in \code{x} intersects
#' (or other predicates such as is covered by) anything in \code{y}. Simplified
#' from the \code{sgbp} results which indicate which item in \code{y} each item
#' in \code{x} intersects, which might not always be relevant.
#' @export
#' @importFrom sf st_intersects
st_any_pred <- function(x, y, pred) lengths(pred(x, y)) > 0L

#' @rdname st_any_pred
#' @export
st_any_intersects <- function(x, y) st_any_pred(x, y, st_intersects)

#' Crop an SFE object with a geometry
#'
#' Returns an SFE object whose specified \code{colGeometry} returns \code{TRUE}
#' with a geometric predicate function (usually intersects) with another
#' geometry of interest. This can be used to subset an SFE object with a tissue
#' boundary or histological region polygon, or crop away empty spaces.
#'
#' @param x An SFE object.
#' @param y An object of class \code{sf}, \code{sfc}, or \code{sfg} with which
#'   to crop the SFE object. Optional if \code{xmin}, \code{xmax}, \code{ymin},
#'   and \code{ymax} are specified for a bounding box.
#' @param colGeometryName Column geometry to used to indicate which cells/spots
#'   to keep.
#' @param pred A geometric binary predicate function to indicate which
#'   cells/spots to keep, defaults to \code{\link{st_intersects}}.
#' @param xmin Minimum x coordinate of bounding box. Ignored if \code{y} is
#'   specified.
#' @param xmax Maximum x coordinate of bounding box.
#' @param ymin Minimum y coordinate of bounding box.
#' @param ymax Maximum y coordinate of bounding box.
#' @return An SFE object.
#' @export
crop <- function(x, y = NULL, colGeometryName = 1L, pred = st_intersects,
                 xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL) {
  if (is.null(y)) {
    y <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                           crs = NA))
  }
  cg <- colGeometry(x, type = colGeometryName, sample_id = "all")
  x[, st_any_pred(cg, y, pred)]
}
