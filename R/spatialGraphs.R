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
#' @param value A \code{listw} object (\code{*Graph}), or a named list of
#'   \code{listw} objects (\code{*Graphs}).
#' @param type An integer specifying the index or string specifying the name of
#'   the *Graph to query or replace. If missing, then the first item in the
#'   *Graph will be returned or replaced.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns. In addition, 3 stands for spatial neighborhood graphs that
#'   correspond to \code{annotGeometries}.
#' @name spatialGraphs
#' @docType methods
NULL

.margin_name <- function(MARGIN) switch (MARGIN, "row", "col", "annot")
.margin_len_fun <- function(MARGIN) switch(MARGIN, nrow, ncol, function(x) return(NA))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", "SpatialFeatureExperiment",
          function(x, MARGIN) int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]])

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", "SpatialFeatureExperiment",
                 function(x, MARGIN, value) {
                   m <- .check_graphs(value, .margin_len_fun(MARGIN)(x),
                                      .margin_name(MARGIN))
                   if (!isTRUE(m)) stop(m)
                   int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]] <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "missing"),
          function(x, type, MARGIN) spatialGraph(x, 1L, MARGIN))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric"),
          function(x, type, MARGIN) spatialGraphs(x, MARGIN)[[type]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character"),
          function(x, type, MARGIN) spatialGraphs(x, MARGIN)[[type]])

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "missing", "numeric", "listw"),
                 function(x, type, MARGIN, value) spatialGraph(x, 1L, MARGIN) <- value)

.sg_r <- function(x, type, MARGIN, value) {
  if (MARGIN < 3 && length(value$neighbours) != .margin_len_fun(MARGIN)(x)) {
    stop("The neighbours field of `value` must be the same as n",
         .margin_name(MARGIN), "s of the gene count matrix.")
  }
  int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[type]] <- value
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric", "numeric", "listw"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "numeric", "listw"),
                 .sg_r)

# To do: 1. store info to reconstruct the graph after subsetting
# 2. Eventually deal with sample_id, but assume one sample per object for now
# Might change spatialGraph into a DataFrame with a list column for the listw