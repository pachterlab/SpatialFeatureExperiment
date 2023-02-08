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
#'   \code{SpatialFeatureExperiment} object.
#' @param BPPARAM Passed to \code{\link{df2sf}}, to parallelize the conversion
#'   of centroid spatial coordinates in the SPE object to \code{sf} point
#'   geometry.
#' @return An SFE object
#' @importFrom S4Vectors make_zero_col_DFrame
#' @importFrom SpatialExperiment spatialCoords toSpatialExperiment
#' @name SpatialFeatureExperiment-coercion
#' @aliases toSpatialFeatureExperiment
#' @examples
#' library(SpatialExperiment)
#' example(read10xVisium)
#' # There can't be duplicate barcodes
#' colnames(spe) <- make.unique(colnames(spe), sep = "-")
#' rownames(spatialCoords(spe)) <- colnames(spe)
#' sfe <- toSpatialFeatureExperiment(spe)
NULL

.sc2cg <- function(coords_use, spotDiameter = NA, BPPARAM = SerialParam()) {
    colnames(coords_use) <- c("x", "y")
    cg_sfc <- df2sf(coords_use, spotDiameter = spotDiameter, BPPARAM = BPPARAM,
                    geometryType = "POINT")
    rownames(cg_sfc) <- rownames(coords_use)
    cg_sfc
}
setAs(
    from = "SpatialExperiment", to = "SpatialFeatureExperiment",
    function(from) {
        cg <- int_colData(from)[["colGeometries"]]
        if (is.null(cg)) {
            coords_use <- spatialCoords(from)
            if (is.null(rownames(coords_use))) {
                rownames(coords_use) <- colnames(from)
            }
            cg <- .sc2cg(coords_use)
            int_colData(from)[["colGeometries"]] <-
                make_zero_col_DFrame(nrow(int_colData(from)))
            int_colData(from)$colGeometries$centroids <- cg
            from
        }
        .spe_to_sfe(from, int_colData(from)[["colGeometries"]],
            int_elementMetadata(from)[["rowGeometries"]],
            int_metadata(from)[["annotGeometries"]],
            spatialCoordsNames(from), "POLYGON",
            int_metadata(from)[["spatialGraphs"]],
            spotDiameter = NA,
            int_metadata(from)[["unit"]],
            BPPARAM = SerialParam()
        )
    }
)

setAs(from = "SingleCellExperiment", to = "SpatialFeatureExperiment",
      function(from) {
          spe <- as(from, "SpatialExperiment")
          as(spe, "SpatialFeatureExperiment")
      })

#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod(
    "toSpatialFeatureExperiment", "SpatialExperiment",
    function(x, colGeometries = NULL, rowGeometries = NULL,
             annotGeometries = NULL, spatialCoordsNames = c("x", "y"),
             annotGeometryType = "POLYGON",
             spatialGraphs = NULL, spotDiameter = NA, unit = NULL,
             BPPARAM = SerialParam()) {
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
        .spe_to_sfe(
            x, colGeometries, rowGeometries, annotGeometries,
            spatialCoordsNames, annotGeometryType,
            spatialGraphs, spotDiameter, unit, BPPARAM
        )
    }
)

#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod("toSpatialFeatureExperiment", "SingleCellExperiment",
          function(x, sample_id="sample01",
                   spatialCoordsNames = c("x", "y"),
                   spatialCoords=NULL,
                   colGeometries = NULL, rowGeometries = NULL,
                   annotGeometries = NULL,
                   annotGeometryType = "POLYGON",
                   spatialGraphs = NULL, spotDiameter = NA,
                   scaleFactors=1,
                   imageSources=NULL,
                   image_id=NULL,
                   loadImage=TRUE,
                   imgData=NULL,
                   unit = NULL,
                   BPPARAM = SerialParam()) {
              spe <- toSpatialExperiment(x, sample_id=sample_id,
                                         spatialCoordsNames=spatialCoordsNames,
                                         spatialCoords=spatialCoords,
                                         scaleFactors=scaleFactors,
                                         imageSources=imageSources,
                                         image_id=image_id,
                                         loadImage=loadImage,
                                         imgData=imgData)
              toSpatialFeatureExperiment(spe, colGeometries = colGeometries,
                                         rowGeometries = rowGeometries,
                                         annotGeometries = annotGeometries,
                                         spatialCoordsNames = spatialCoordsNames,
                                         annotGeometryType = annotGeometryType,
                                         spatialGraphs = spatialGraphs,
                                         spotDiameter = spotDiameter,
                                         unit = unit,
                                         BPPARAM = BPPARAM)
          })
