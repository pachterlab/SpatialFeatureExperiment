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
#' @param value A \code{listw} object (\code{*Graph}), or a named list of list
#'   of \code{listw} objects (\code{*Graphs}) where the names of the top level
#'   list are \code{sample_id}s when adding graphs for all samples in the margin
#'   of interest, or a list of \code{listw} objects when adding graphs for one
#'   sample in one margin.
#' @param type An integer specifying the index or string specifying the name of
#'   the *Graph to query or replace. If missing, then the first item in the
#'   *Graph will be returned or replaced.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns. In addition, 3 stands for spatial neighborhood graphs that
#'   correspond to \code{annotGeometries}.
#' @param sample_id Name of the sample the graph is associated with. This is
#'   useful when multiple pieces of tissues are in the same SFE object (say for
#'   a joint dimension reduction and clustering) and the spatial neighborhood is
#'   only meaningful within the same piece of tissue. See the \code{sample_id}
#'   argument in \code{\link{SpatialExperiment}}.
#' @name spatialGraphs
#' @docType methods
NULL

.margin_name <- function(MARGIN) switch (MARGIN, "row", "col", "annot")
.margin_len_fun <- function(MARGIN) switch(MARGIN, nrow, ncol, function(x) return(NA))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "missing"),
          function(x, MARGIN, sample_id)
            int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character"),
          function(x, MARGIN, sample_id)
            int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]])

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "missing"),
                 function(x, MARGIN, sample_id, value) {
                   m <- .check_graphs(value, .margin_len_fun(MARGIN)(x),
                                      .margin_name(MARGIN))
                   if (!isTRUE(m)) stop(m)
                   int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]] <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character"),
                 function(x, MARGIN, sample_id, value) {
                   m <- .check_graphs_sample(value, .margin_len_fun(MARGIN)(x),
                                             .margin_name(MARGIN), sample_id)
                   if (!isTRUE(m)) stop(m)
                   int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]] <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "missing"),
          function(x, type, MARGIN, sample_id) spatialGraph(x, 1L, MARGIN, sample_id))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric"),
          function(x, type, MARGIN, sample_id)
            spatialGraphs(x, MARGIN)[[sample_id]][[type]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character"),
          function(x, type, MARGIN, sample_id)
            spatialGraphs(x, MARGIN)[[sample_id]][[type]])

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "missing",
                                   "numeric", "character", "listw"),
                 function(x, type, MARGIN, sample_id, value)
                   spatialGraph(x, 1L, MARGIN, sample_id) <- value)

.sg_r <- function(x, type, MARGIN, sample_id, value) {
  if (MARGIN < 3 && length(value$neighbours) != .margin_len_fun(MARGIN)(x)) {
    stop("The neighbours field of `value` must be the same as n",
         .margin_name(MARGIN), "s of the gene count matrix.")
  }
  int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]][[type]] <- value
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric",
                                   "numeric", "character", "listw"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character",
                                   "numeric", "character", "listw"),
                 .sg_r)

# To do: 1. store info to reconstruct the graph after subsetting
# 2. Put the barcode row and column info in this package for Visium. Might write
# getters and setters for grid rows and columns of barcodes. Use additional
# columns in int_colData.

#' Find spatial neighborhood graphs for Visium spots
#'
#' Visium spots are arranged in a hexagonal grid. This function uses the known
#' locations of the Visium barcodes to construct a neighborhood graph, so
#' adjacent spots are connected by edges. Since the known rows and columns of
#' the spots are used, the unit the spot centroid coordinates are in does not
#' matter.
#'
#' @inheritParams spdep::nb2listw
#' @param x A \code{SpatialFeatureExperiment} object with Visium data. Column
#'   names of the gene count matrix must be Visium barcodes, which may have a
#'   numeric suffix to distinguish between samples (e.g. "AAACAACGAATAGTTC-1").
#' @param sample_id String specifying the sample, as in \code{sample_id} in
#'   \code{colData}, for which the graph is to be constructed.
#' @return A \code{listw} object for the neighborhood graph. The node index will
#' be the order of barcodes within the sample in the gene count matrix.
#' @importFrom spdep dnearneigh nb2listw
#' @export
findVisiumGraph <- function(x, sample_id, style = "W", zero.policy = NULL) {
  bcs_use <- colnames(x)[colData(x)$sample_id == sample_id]
  bcs_use2 <- sub("-\\d+$", "", bcs_use)
  coords_use <- visium_row_col[match(bcs_use2, visium_row_col), c("col", "row")]
  # So adjacent spots are equidistant
  coords_use$row <- coords_use$row * sqrt(3)
  g <- dnearneigh(as.matrix(coords_use), d1 = 1.9, d2 = 2.1, row.names = bcs_use)
  out <- nb2listw(g, style = style, zero.policy = zero.policy)
  return(out)
}
