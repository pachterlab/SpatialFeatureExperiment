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
#' indices in the graphs are no longer valid. At present, this only works with
#' the wrapper functions in this package that take in SFE objects and records
#' the info required to reconstruc the graphs.
#' @importFrom methods callNextMethod
#' @return A subsetted \code{SpatialFeatureExperiment} object.
#' @name SpatialFeatureExperiment-subset
#' @export
setMethod("[", c("SpatialFeatureExperiment", "ANY", "ANY"),
          function(x, i, j, ..., reconstruct_graphs = TRUE, drop = FALSE) {
            x <- callNextMethod()
            # Subset annotGeometries based on sample_id
            if (!is.null(annotGeometries(x))) {
              ag_sub <- annotGeometries(x)
              for (g in seq_along(ag_sub)) {
                ag_ind <- ag_sub[[g]]
                ag_sub[[g]] <- ag_ind[ag_ind$sample_id %in% sampleIDs(x),]
              }
            }
            # Subset *Graphs based on sample_id and reconstruct row and colGraphs
            if (!is.null(spatialGraphs(x))) {
              graphs_sub <- spatialGraphs(x)
              graphs_sub <- graphs_sub[,names(graphs_sub) %in% sampleIDs(x)]
              if (reconstruct_graphs) {
                for (s in seq_along(graphs_sub)) {
                  for (m in 1:2) { # Not reconstructing annotGraphs
                    for (g in seq_along(graphs_sub[[s]][[m]])) {
                      method_info <- attr(graphs_sub[[s]][[m]][[g]], "method")
                      if (is.null(method_info)) {
                        warning("Graph reconstruction info is missing for sample ",
                                names(graphs_sub)[s], " ", .margin_name(m), "Graph ",
                                names(graphs_sub[[s]][[m]])[g], ". ",
                                "Not reconstructing graph. ",
                                "Node indices in the graphs are no longer valid after subsetting.")
                      } else {
                        # To do: check that the package used is present
                        graphs_sub[[s]][[m]][[g]] <- do.call(method_info$FUN,
                                                             c(x = x, method_info$args))
                      }
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
# 3. Method to crop the SFE object with any bbox or (multi)polygon.
# 4. Deal with annotGeometries and Graphs when concatenating SFE objects.
