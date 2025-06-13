#' Identify cells in small bits outside the main piece of tissue
#'
#' This is used for quality control (QC), as the small bits are likely to be low
#' quality technical artifacts and are not informative to spatial analyses.
#' Please confirm the quality of those cells by checking the histology image
#' with ImageJ or QuPath.
#'
#' How this function works: A distance-based spatial neighborhood graph is
#' computed, with \code{distance_cutoff} as the distance cutoff. Then disjoint
#' connected subgraphs are found. Cells in subgraphs with \code{max_cells} or
#' fewer cells are considered debris.
#'
#' @inheritParams findSpatialNeighbors
#' @param x Object with spatial coordinates of cells. Can be a
#'   \code{SpatialExperiment} or \code{SpatialFeatureExperiment} object, a
#'   matrix with 2 columns for x and y coordinates of cells, or a \code{sf} or
#'   \code{sfc} object with cell geometries.
#' @param max_cells Maximum number of cells for a clump of cells to be
#'   considered debris.
#' @param distance_cutoff Minimum distance of cell to the tissue for it to be
#'   considered debris, in the same unit as in \code{x}.
#' @return Depends on the method:
#' \describe{
#' \item{Spatial(Feature)Experiment}{The same object with a logical column 
#' "is_debris" added to \code{colData}.}
#' \item{Matrix and \code{sf(c)}}{A logical vector indicating whether each cell
#' is debris.}
#' }
#' @name findDebrisCells
NULL

#' @rdname findDebrisCells
#' @export
setMethod("findDebrisCells", "matrix", 
          function(x, max_cells = 5, distance_cutoff = 50,
                   BNPARAM = NULL, BPPARAM = SerialParam()) {
              g <- .dnn_bioc(x, distance_cutoff, BNPARAM = BNPARAM,
                             BPPARAM = BPPARAM)
              comps <- spdep::n.comp.nb(g)
              n_cells <- table(comps$comp.id)
              if (any(n_cells <= max_cells)) {
                  comps_use <- names(n_cells)[n_cells <= max_cells]
                  return(comps$comp.id %in% comps_use)
              } else return(rep(FALSE, nrow(x)))
          })

#' @rdname findDebrisCells
#' @export
setMethod("findDebrisCells", "sfc",
          function(x, max_cells = 5, distance_cutoff = 50,
                   BNPARAM = NULL, BPPARAM = SerialParam()) {
              if (st_geometry_type(x, by_geometry = FALSE) != "POINT")
                  x <- st_centroid(x)
              findDebrisCells(st_coordinates(x)[,1:2], max_cells, distance_cutoff,
                              BNPARAM, BPPARAM)
          })

#' @rdname findDebrisCells
#' @export
setMethod("findDebrisCells", "sf", 
          function(x, max_cells = 5, distance_cutoff = 50,
                   BNPARAM = NULL, BPPARAM = SerialParam()) {
              findDebrisCells(st_geometry(x), max_cells, distance_cutoff,
                              BNPARAM, BPPARAM)
          })

#' @rdname findDebrisCells
#' @export
setMethod("findDebrisCells", "SpatialExperiment",
          function(x, max_cells = 5, distance_cutoff = 50,
                   BNPARAM = NULL, BPPARAM = SerialParam()) {
              sfe$is_debris <- findDebrisCells(spatialCoords(x), max_cells, distance_cutoff,
                                               BNPARAM, BPPARAM)
              sfe
          })
