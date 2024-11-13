# Split-------------

#' Split SFE object with categorical vector or geometry
#'
#' The \code{split} methods for SFE split an SFE object into multiple SFE
#' objects by geometries (all cells/spots intersecting with each geometry will
#' become a separate SFE object). The \code{splitSamples} function splits the
#' SFE object by \code{sample_id} so each sample will become a separate SFE
#' object. The \code{splitContiguity} function splits the SFE object by
#' contiguity of an \code{annotGeometry}, which by default is "tissueBoundary".
#'
#' @inheritParams crop
#' @param x An SFE object
#' @param f It can be a \code{sf} data frame or \code{sfc} to split by geometry.
#'   Each row of the \code{sf} data frame or each element in the \code{sfc} will
#'   correspond to a new SFE object. The \code{sf} data frame must have a column
#'   \code{sample_id} when splitting multiple samples. Can also be a list of
#'   \code{sfc} whose names correspond to \code{sample_id}s to split.
#' @param sample_id Which samples to split.
#' @param colGeometryName Which \code{colGeometry} to use to determine which
#'   cells or spots should belong to which new SFE object when splitting by
#'   \code{sf} or \code{sfc}. Default to the first one.
#' @param annotGeometryName Name of \code{annotGeometry} to use to split by
#'   contiguity.
#' @param min_area Minimum area in the same unit as the geometry coordinates
#'   (squared) for each piece to be considered a separate piece when splitting
#'   by contiguity. Only pieces that are large enough are considered.
#' @return A list of SFE objects.
#' @concept Geometric operations
#' @name splitByCol
#' @aliases splitByCol
#' @examples
#' # example code
NULL

#' @rdname splitByCol
#' @export
setMethod("splitByCol", c("SpatialFeatureExperiment", "sf"),
          function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
              sample_id <- .check_sample_id(x, sample_id, one = FALSE)
              if (!"sample_id" %in% names(f) && length(sample_id) > 1L)
                  stop("f must have a column sample_id when multiple samples are specified.")
              l <- split(st_geometry(f), f$sample_id)
              splitByCol(x, l, sample_id = sample_id, colGeometryName = colGeometryName, cover = cover)
          })

#' @rdname splitByCol
#' @export
setMethod("splitByCol", c("SpatialFeatureExperiment", "sfc"),
          function(x, f, sample_id = 1L, colGeometryName = 1L, cover = FALSE) {
              sample_id <- .check_sample_id(x, sample_id)
              x <- x[, x$sample_id == sample_id]
              lapply(f, function(g) {
                  crop(x, g, colGeometryName = colGeometryName,
                       keep_whole = "col", cover = cover)
              })
          })

#' @rdname splitByCol
#' @export
setMethod("splitByCol", c("SpatialFeatureExperiment", "list"),
          function(x, f, sample_id = "all", colGeometryName = 1L, cover = FALSE) {
              sample_id <- .check_sample_id(x, sample_id, one = FALSE)
              if (!any(sample_id %in% names(f)))
                  stop("None of the geometries correspond to sample_id")
              f <- f[intersect(sample_id, names(f))]
              out <- lapply(sample_id, function(s) {
                  splitByCol(x, f[[s]], sample_id = s, colGeometryName = colGeometryName,
                        cover = cover)
              })
              names(out) <- sample_id
              unlist(out, recursive = FALSE)
          })

#' @rdname splitByCol
#' @export
splitSamples <- function(x) {
    ss <- sampleIDs(x)
    out <- lapply(ss, function(s) x[, x$sample_id == s])
    names(out) <- ss
    out
}

#' @rdname splitByCol
#' @export
#' @importFrom sf st_collection_extract
splitContiguity <- function(x, colGeometryName = 1L,
                            annotGeometryName = "tissueBoundary",
                            min_area = 0, cover = FALSE) {
    ag <- annotGeometry(x, annotGeometryName)
    gt <- st_geometry_type(ag, by_geometry = FALSE)
    # Will I allow points and linestrings in the future?
    if (!gt %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRY"))
        stop("The geometries must be POLYGON, MULTIPOLYGON, or GEOMETRY.")
    if (gt == "GEOMETRY") {
        tryCatch(ag <- st_collection_extract(ag, type = c("POLYGON", "MULTIPOLYGON")),
                 warning = function(w) stop("None of the geometries are POLYGON or MULTIPOLYGON"))
    }
    ag_union <- st_union(ag) # Merge contiguous pieces
    ag_union <- st_cast(ag_union, "POLYGON", warn = FALSE)
    if (min_area > 0) {
        ag_union$area <- st_area(ag_union)
        ag_union <- ag_union[ag_union$area > min_area,]
    }
    splitByCol(x, ag_union, colGeometryName = colGeometryName, cover = cover)
}
