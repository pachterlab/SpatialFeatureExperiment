#' Annotation geometry methods
#'
#' "Annotation geometry" refers to Simple Feature (\code{sf}) geometries NOT
#' associated with rows (features, genes) or columns (cells or spots) of the
#' gene count matrix in the \code{SpatialFeatureExperiment} object. So there can
#' be any number of rows in the \code{sf} data frame specifying the geometry.
#' Examples of such geometries are tissue boundaries, pathologist annotation of
#' histological regions, and objects not characterized by columns of the gene
#' count matrix (e.g. nuclei segmentation in a Visium dataset where the columns
#' are Visium spots). This page documents getters and setters for the annotation
#' geometries. Internally, annotation geometries are stored in
#' \code{int_metadata}.
#'
#' @inheritParams dimGeometries
#' @param value Value to set. For \code{annotGeometry}, must be a \code{sf} data
#'   frame, or an ordinary data frame that can be converted to a \code{sf} data
#'   frame (see \code{\link{df2sf}}). For \code{annotGeometries}, must be a list
#'   of such \code{sf} or ordinary data frames.
#' @name annotGeometries
NULL

#' @rdname annotGeometries
#' @export
setMethod("annotGeometries", "SpatialFeatureExperiment",
          function(x) int_metadata(x)$annotGeometries)

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometries", "SpatialFeatureExperiment",
                 function(x, ..., value) {
                   value <- .df2sf_list(value, ...)
                   int_metadata(x)$annotGeometries <- value
                   return(x)
                 })

#' @rdname annotGeometries
#' @export
setMethod("annotGeometryNames", "SpatialFeatureExperiment",
          function(x) names(annotGeometries(x)))

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometryNames", c("SpatialFeatureExperiment", "character"),
                 function(x, value) {
                   names(annotGeometries(x)) <- value
                   return(x)
                 })

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", c("SpatialFeatureExperiment", "missing"),
          function(x, type) annotGeometry(x, 1L))

.ag <- function(x, type) annotGeometries(x)[[type]]

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", c("SpatialFeatureExperiment", "numeric"), .ag)

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", c("SpatialFeatureExperiment", "character"), .ag)

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometry", c("SpatialFeatureExperiment", "missing"),
          function(x, type, value) annotGeometry(x, 1L) <- value)

.ag_r <- function(x, type, ..., value) {
  value <- .df2sf_in_list(value, ...)
  int_metadata(x)$annotGeometries[[type]] <- value
  return(x)
}

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometry", c("SpatialFeatureExperiment", "numeric"),
                 .ag_r)

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometry", c("SpatialFeatureExperiment", "character"),
                 .ag_r)
