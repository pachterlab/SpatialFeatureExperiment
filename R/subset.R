#' Subsetting SpatialFeatureExperiment objects
#'
#' The method for SFE reconstructs the spatial graphs when the SFE object is
#' subsetted as the \code{listw} objects encodes the nodes with indices which
#' are no longer valid after subsetting as some nodes are no longer present.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param i Row indices for subsetting.
#' @param j column indices for subsetting.
#' @param reconstruct_graph Logical, whether to reconstruct graph. If \code{FALSE},
#' then the old graphs will be kept and a warning will be issued that the node
#' indices in the graphs are no longer valid.
#' @importFrom methods callNextMethod
#' @return A subsetted \code{SpatialFeatureExperiment} object.
#' @name SpatialFeatureExperiment-subset
#' @export
setMethod("[", c("SpatialFeatureExperiment", "ANY", "ANY"),
          function(x, i, j, ..., reconstruct_graphs = TRUE, drop = FALSE) {
            x <- callNextMethod()
            if (!is.null(spatialGraphs(x))) {
              graphs_sub <- spatialGraphs(x)
              graphs_sub <- graphs_sub[names(graphs_sub) %in% sampleIDs(x)]
              if (reconstruct_graphs) {
                for (s in seq_along(graphs_sub)) {
                  for (g in seq_along(graphs_sub[[s]])) {
                    method_info <- attr(graphs_sub[[s]][[g]], "method")
                    if (!is.null(method_info)) {
                      graphs_sub[[s]][[g]] <- do.call(method_info$FUN,
                                                      c(x = x, method_info$args))
                    }
                  }
                }
              } else {
                warning("Node indices in the graphs are no longer valid after subsetting.")
              }
              spatialGraphs(x) <- graphs_sub
            }
            return(x)
          })

# To do:
# 1. Add sample_id to *Geometry getters as optional argument
# 2. I might make sample_id optional in the *Graphs functions as I expect most uses of SFE to only have one sample per object
# 3. Method to crop the SFE object with any bbox or (multi)polygon.
# 4. Deal with annotGeometries and Graphs when concatenating SFE objects.
