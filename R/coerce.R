#' SpatialFeatureExperiment coercion methods
#'
#' The \code{SpatialFeatureExperiment} class inherits from
#' \code{SpatialExperiment}, which in turn inherits from
#' \code{SingleCellExperiment}. A \code{SpatialExperiment} object with
#' geometries in \code{colGeometries} in the \code{int_colData},
#' \code{rowGeometries} in the \code{int_elementMetadata}, or
#' \code{annotGeometries} in the \code{int_metadata} can be directly converted
#' to \code{SpatialFeatureExperiment} with \code{as(spe,
#' "SpatialFeatureExperiment")}. A \code{SpatialExperiment} object without the
#' geometries can also be converted; the coordinates in the \code{spatialCoords}
#' field will be used to make POINT geometries named "centroids" to add to
#' \code{colGeometries}. The geometries can also be supplied separately when
#' using \code{toSpatialFeatureExperiment}. For now coercion only works for
#' \code{SpatialExperiment}. I'll deal with \code{Seurat} and
#' \code{SingleCellExperiment} later.
#'
#' @inheritParams SpatialFeatureExperiment
#' @param x A \code{SpatialExperiment} object to be coerced to a
#' \code{SpatialFeatureExperiment} object.
#' @importFrom S4Vectors make_zero_col_DFrame
#' @importFrom SpatialExperiment spatialCoords
#' @name SpatialFeatureExperiment-coercion
#' @aliases toSpatialFeatureExperiment
NULL

.sc2cg <- function(coords_use) {
  cg_sfc <- apply(coords_use, 1, st_point, simplify = FALSE)
  st_sf(geometry = cg_sfc, row.names = rownames(coords_use))
}
setAs(from = "SpatialExperiment", to = "SpatialFeatureExperiment",
      function(from) {
        cg <- int_colData(from)[["colGeometries"]]
        if (is.null(cg)) {
          coords_use <- spatialCoords(from)
          cg <- .sc2cg(coords_use)
          int_colData(from)[["colGeometries"]] <- make_zero_col_DFrame(nrow(int_colData(from)))
          int_colData(from)$colGeometries$centroids <- cg
          from
        }
        .spe_to_sfe(from, int_colData(from)[["colGeometries"]],
                    int_elementMetadata(from)[["rowGeometries"]],
                    int_metadata(from)[["annotGeometries"]],
                    spatialCoordsNames(from), "POLYGON",
                    int_metadata(from)[["spatialGraphs"]],
                    int_metadata(from)[["unit"]])
      })

#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod("toSpatialFeatureExperiment", "SpatialExperiment",
          function(x, colGeometries = NULL, rowGeometries = NULL,
                   annotGeometries = NULL, spatialCoordsNames = c("x", "y"),
                   annotGeometryType = "POLYGON",
                   spatialGraphs = NULL, unit = NULL) {
            if (is.null(colGeometries)) {
              colGeometries <- int_colData(x)$colGeometries
            }
            if (is.null(rowGeometries)) {
              rowGeometries <- int_elementMetadata(x)$rowGeometries
            }
            if (is.null(annotGeometries)) {
              annotGeometries <- int_metadata(x)$annotGeometries
            }
            if (is.null(spatialGraphs)) {
              spatialGraphs <- int_metadata(x)$spatialGraphs
            }
            .spe_to_sfe(x, colGeometries, rowGeometries, annotGeometries,
                        spatialCoordsNames, annotGeometryType,
                        spatialGraphs, unit)
          })
