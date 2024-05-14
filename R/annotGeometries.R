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
#' Wrapper for getter and setter of special geometry: \describe{
#' \item{tisseuBoundary}{Boundary of the tissue of interest, including holes.
#' This is usually of geometry type MULTIPOLYGON, though geometries in
#' \code{annotGeometries} can have any type supported by \code{sf}.} }
#'
#' @inheritParams dimGeometries
#' @param value Value to set. For \code{annotGeometry}, must be a \code{sf} data
#'   frame, or an ordinary data frame that can be converted to a \code{sf} data
#'   frame (see \code{\link{df2sf}}). For \code{annotGeometries}, must be a list
#'   of such \code{sf} or ordinary data frames. There must be a column
#'   \code{sample_id} to indicate the sample the geometries are for, and the
#'   \code{sample_id} must also appear in \code{colData}.
#' @return Getters for multiple geometries return a named list. Getters for
#'   names return a character vector of the names. Getters for single geometries
#'   return an \code{sf} data frame. Setters return an SFE object.
#' @name annotGeometries
#' @aliases annotGeometries<- annotGeometry annotGeometry<- annotGeometryNames
#'   annotGeometryNames<-
#' @concept Getters and setters
#' @examples
#' # Example dataset
#' library(SFEData)
#' sfe_small <- McKellarMuscleData(dataset = "small")
#'
#' # Get all annotation geometries, returning a named list
#' annotGeometries(sfe_small)
#'
#' # Set all annotation geometries, in a named list
#' toy <- readRDS(system.file("extdata/sfe_toy.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' ag <- readRDS(system.file("extdata/ag.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' annotGeometries(toy) <- list(hull = ag)
#'
#' # Get names of annotation geometries
#' annotGeometryNames(sfe_small)
#'
#' # Set names of annotation geometries
#' annotGeometryNames(toy) <- "foo"
#'
#' # Get a specific annotation geometry by name
#' # sample_id is optional when there is only one sample present
#' nuclei <- annotGeometry(sfe_small, type = "nuclei", sample_id = "Vis5A")
#'
#' # Get a specific annotation geometry by index
#' tb <- annotGeometry(sfe_small, type = 1L)
#'
#' # Set a specific annotation geometry
#' annotGeometry(sfe_small, type = "nuclei2") <- nuclei
#'
#' # Special convenience function for tissue boundaries
#' # Getter
#' tb <- tissueBoundary(sfe_small, sample_id = "Vis5A")
#' # Setter
#' tissueBoundary(sfe_small, sample_id = "Vis5A") <- tb
NULL

#' @rdname annotGeometries
#' @export
setMethod(
    "annotGeometries", "SpatialFeatureExperiment",
    function(x) int_metadata(x)$annotGeometries
)

#' @rdname annotGeometries
#' @export
setReplaceMethod(
    "annotGeometries", "SpatialFeatureExperiment",
    function(x, translate = TRUE, ..., value) {
        value <- .df2sf_list(value, ...)
        value <- lapply(value, .rm_empty_geometries, MARGIN = 3)
        value <- lapply(value, .translate_value,
            x = x,
            translate = translate
        )
        int_metadata(x)$annotGeometries <- value
        m <- .check_annotgeometries(x)
        if (length(m)) stop(m)
        return(x)
    }
)

#' @rdname annotGeometries
#' @export
setMethod(
    "annotGeometryNames", "SpatialFeatureExperiment",
    function(x) names(annotGeometries(x))
)

#' @rdname annotGeometries
#' @export
setReplaceMethod(
    "annotGeometryNames", c(
        "SpatialFeatureExperiment",
        "character"
    ),
    function(x, value) {
        names(annotGeometries(x)) <- value
        return(x)
    }
)

.ag <- function(x, type = 1L, sample_id = NULL) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    out <- int_metadata(x)$annotGeometries[[type]]
    if (!is.null(sample_id)) {
        out <- out[out$sample_id %in% sample_id, ]
    }
    if (is.null(out))
        stop("annotGeometry ", type, " is absent.")
    return(out)
}

#' @rdname annotGeometries
#' @export
setMethod("annotGeometry", "SpatialFeatureExperiment", .ag)

.ag_r <- function(x, type = 1L, sample_id = NULL, translate = TRUE, ..., value) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (length(sampleIDs(x)) == 1L && !"sample_id" %in% names(value)) {
        value$sample_id <- sample_id
    }
    value <- .df2sf_in_list(value, ...)
    value <- .rm_empty_geometries(value, MARGIN = 3)
    if (!is.null(sample_id) && any(!sampleIDs(x) %in% sample_id)) {
        existing <- int_metadata(x)$annotGeometries[[type]]
        if (!is.null(existing)) {
            existing <- .reconcile_cols(existing, value)
            value <- .reconcile_cols(value, existing)
            value <- value[, names(existing)]
            if (sample_id %in% existing$sample_id) {
                # This is a replacement method, so do replace
                existing <- existing[!existing$sample_id %in% sample_id, ]
            }
            value <- rbind(existing, value)
        }
    }
    value <- .translate_value(x, translate, value)
    int_metadata(x)$annotGeometries[[type]] <- value
    m <- .check_annotgeometries(x)
    if (length(m)) stop(m)
    return(x)
}

#' @rdname annotGeometries
#' @export
setReplaceMethod("annotGeometry", "SpatialFeatureExperiment", .ag_r)

#' @rdname annotGeometries
#' @export
tissueBoundary <- function(x, sample_id = 1L) {
    annotGeometry(x, "tissueBoundary", sample_id)
}

#' @rdname annotGeometries
#' @export
`tissueBoundary<-` <- function(x, sample_id = 1L, translate = TRUE, ...,
                               value) {
    annotGeometry(x, "tissueBoundary", sample_id, translate, ...) <- value
    x
}
