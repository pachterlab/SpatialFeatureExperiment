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
#' Wrapper for getter and setter of special geometry:
#' \describe{
#' \item{tisseuBoundary}{Boundary of the tissue of interest, including holes.
#' This is usually of geometry type MULTIPOLYGON, though geometries in
#' \code{annotGeometries} can have any type supported by \code{sf}.}
#' }
#'
#' @inheritParams dimGeometries
#' @param value Value to set. For \code{annotGeometry}, must be a \code{sf} data
#'   frame, or an ordinary data frame that can be converted to a \code{sf} data
#'   frame (see \code{\link{df2sf}}). For \code{annotGeometries}, must be a list
#'   of such \code{sf} or ordinary data frames. There must be a column
#'   \code{sample_id} to indicate the sample the geometries are for, and the
#'   \code{sample_id} must also appear in \code{colData}.
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
                   m <- .check_annotgeometries(x)
                   if (length(m)) stop(m)
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
          function(x, type, sample_id = NULL) annotGeometry(x, 1L, sample_id))

.ag <- function(x, type, sample_id = NULL) {
  out <- int_metadata(x)$annotGeometries[[type]]
  if (!is.null(sample_id)) {
    out <- out[out$sample_id %in% sample_id,]
  }
  return(out)
}

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", c("SpatialFeatureExperiment", "numeric"), .ag)

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", c("SpatialFeatureExperiment", "character"), .ag)

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometry", c("SpatialFeatureExperiment", "missing"),
          function(x, type, sample_id, value) annotGeometry(x, 1L, sample_id) <- value)

.ag_r <- function(x, type, sample_id = NULL, ..., value) {
  value <- .df2sf_in_list(value, ...)
  if (!is.null(sample_id)) {
    existing <- int_metadata(x)$annotGeometries[[type]]
    if (!is.null(existing)) {
      value <- .reconcile_cols(existing, value)
      if (sample_id %in% existing$sample_id) {
        # This is a replacement method, so do replace
        existing <- existing[!existing$sample_id %in% sample_id,]
      }
      value <- rbind(existing, value)
    }
  }
  int_metadata(x)$annotGeometries[[type]] <- value
  m <- .check_annotgeometries(x)
  if (length(m)) stop(m)
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

#' @rdname annotGeometries
#' @export
tissueBoundary <- function(x, sample_id = NULL) {
  annotGeometry(x, "tissueBoundary", sample_id)
}

#' @rdname annotGeometries
#' @export
`tissueBoundary<-` <- function(x, sample_id, ..., value) {
  annotGeometry(x, "tissueBoundary", sample_id, ...) <- value
  x
}
