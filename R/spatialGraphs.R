# Spatial graphs Another item in int_metadata. Since I used spdep a lot, I'll use
# spdep's nb and listw. Burning question: Shall I store the nb or listw? The
# listw actually contains nb, and normally the W style is used (I haven't used
# any other style). So I think I'll store listw. Also, what to do with singletons?
# And the graph will be reconstructed when the SFE object is subsetted. Then I need
# to store sufficient info about how the graph was constructed. I think, to do so,
# I need to write wrappers of all the spatial graph functions in spdep to easily
# store the construction info, which would be unavailable if the user calls those
# functions separately. I'll issue a warning if such info is unavailable because
# the wrappers were not used.

#' Spatial graph methods
#'
#' Spatial neighborhood graphs as \code{spdep}'s \code{listw} objects are stored
#' in the \code{int_metadata} of the SFE object. The \code{listw} class is used
#' because \code{spdep} has many useful methods that rely on the neighborhood
#' graph as \code{listw}.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @name spatialGraphs
NULL

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", "SpatialFeatureExperiment",
          function(x) int_metadata(x)$spatialGraphs)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", "SpatialFeatureExperiment",
                 function(x, value) {
                   m <- .check_graphs(value, ncol(x))
                   if (!isTRUE(m)) stop(m)
                   int_metadata(x)$spatialGraphs <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "missing"),
          function(x, type) spatialGraph(x, 1L))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric"),
          function(x, type) spatialGraphs(x)[[type]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character"),
          function(x, type) spatialGraphs(x)[[type]])

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "missing", "listw"),
                 function(x, type, value) spatialGraph(x, 1L) <- value)

.sg_r <- function(x, type, value) {
  if (length(value$neighbours) != ncol(x)) {
    stop("The neighbours field of `value` must be the same as the number of columns",
         " in the gene count matrix.")
  }
  int_metadata(x)$spatialGraphs[[type]] <- value
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric", "listw"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "listw"),
                 .sg_r)

# To do: 1. store info to reconstruct the graph after subsetting
# 2. Eventually deal with sample_id, but assume one sample per object for now
# Might change spatialGraph into a DataFrame with a list column for the listw
