#' The SpatialFeatureExperiment class
#'
#' This class inherits from the \code{\link{SpatialExperiment}} (SPE) class,
#' which in turn inherits from \code{\link{SingleCellExperiment}} (SCE).
#' \code{SpatialFeatureExperiment} stores geometries of spots or cells in
#' \code{sf} objects which form columns of a \code{DataFrame} which is in turn a
#' column of the \code{int_colData} \code{DataFrame} of the underlying SCE
#' object, just like \code{reducedDim} in SCE. Geometries of the tissue outline,
#' pathologist annotations, and objects (e.g. nuclei segmentation in a Visium
#' dataset) are stored in \code{sf} objects in a named list called
#' \code{geometries} in \code{int_metadata}.
#'
#' @rdname SpatialFeatureExperiment-class
#' @include utils.R
#' @importFrom methods setClass new
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @exportClass SpatialFeatureExperiment
setClass("SpatialFeatureExperiment", contains = "SpatialExperiment")

#' Constructor of SpatialFeatureExperiment object
#'
#' Create a \code{SpatialFeatureExperiment} object.
#'
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @inheritParams SpatialExperiment::SpatialExperiment
#' @param colGeometry Geometry of the entities that correspond to the columns of
#'   the gene count matrix, such as cells and Visium spots. It must be a named
#'   list of be one of the following: \describe{ \item{An \code{sf} data
#'   frame}{The geometry column specifies the geometry of the entities.}
#'   \item{An ordinary data frame specifying centroids}{Column names for the
#'   coordinates are specified in the \code{spatialCoordsNames} argument. For
#'   Visium and ST, in addition to the centroid coordinate data frame, the spot
#'   diameter in the same unit as the coordinates can be specified in the
#'   \code{spotDiamter} argument.} \item{An ordinary data frame specifying
#'   polygons}{Also use \code{spatialCoordsNames}. There should an additional
#'   column "ID" to specify which vertices belong to which polygon. The
#'   coordinates should not be in list columns. Rather, the data frame should
#'   look like it is passed to \code{ggplot2::geom_polygon}. If there are holes,
#'   then there must also be a column "subID" that differentiates between the
#'   outer polygon and the holes.} \item{Path to a GeoJSON file}{For QuPath
#'   annotations.} } In all cases, the data frame should specify the same number
#'   of geometries as the number of columns in the gene count matrix. If the
#'   column "barcode" is present, then it will be matched to column names of the
#'   gene count matrix. Otherwise, the geometries are assumed to be in the same
#'   order as columns in the gene count matrix. If the geometries are specified
#'   in an ordinary data frame, then it will be converted into \code{sf}
#'   internally. Named list of data frames because each entity can have multiple
#'   geometries, such as whole cell and nuclei segmentations. The geometries are
#'   assumed to be POINTs for centroids and POLYGONs for segmentations. If
#'   polygons are specified in an ordinary data frame, then anything with fewer
#'   than 3 vertices will be removed. For anything other than POINTs, attributes
#'   of the geometry will be ignored, but they can be added later with the
#'   \code{AddGeometryAttr} function.
#' @param rowGeometry Geometry associated with genes or features, which
#'   correspond to rows of the gene count matrix.
#' @param annotGeometry Geometry of entities that do not correspond to columns
#'   or rows of the gene count matrix, such as tissue boundary and pathologist
#'   annotations of histological regions, and nuclei segmentation in a Visium
#'   dataset. Also a named list as in \code{primaryGeometry}. The ordinary data
#'   frame may specify POINTs, POLYGONs, or LINESTRINGs, or their MULTI
#'   versions. Each data frame can only specify one type of geometry. For MULTI
#'   versions, there must be a column "group" to identify each MULTI geometry.
#' @param annotGeometryType Character vector specifying geometry type of each
#'   element of the list if \code{annotGeometry} is specified. Each element of
#'   the vector must be one of POINT, LINESTRING, POLYGON, MULTIPOINT,
#'   MULTILINESTRING, and MULTIPOLYGON. Must be either length 1 (same for all
#'   elements of the list) or the same length as the list. Ignored if the
#'   corresponding element is an \code{sf} object.
#' @param spotDiameter Spot diameter for technologies with arrays of spots of
#'   fixed diameter per slide, such as Visium, ST, DBiT-seq, and slide-seq. The
#'   diameter must be in the same unit as the coordinates in the *Geometry
#'   arguments.
#' @param unit Unit the coordinates are in. I'm thinking about using some custom
#'   engineering CRS's which can convert units and invert the y axis for
#'   Cartesian vs. image orientations. Units are also helpful when plotting
#'   scale bars.
#' @param ... Additional arguments passed to the \code{\link{SpatialExperiment}}
#'   and \code{\link{SingleCellExperiment}} constructors.
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment int_colData int_elementMetadata int_metadata
#' @importFrom sf st_point st_sfc st_sf st_polygon st_buffer st_linestring
#'   st_multipoint st_multilinestring st_multipolygon st_is st_coordinates
#'   st_centroid
#' @export
SpatialFeatureExperiment <- function(assays, colGeometry,
                                     rowGeometry = NULL, annotGeometry = NULL,
                                     colData = DataFrame(), rowData = NULL,
                                     spatialCoordsNames = c("x", "y"),
                                     sample_id = "sample01", spotDiameter = NA_real_,
                                     annotGeometryType = "POLYGON",
                                     unit = "full_res_image_pixels",
                                     ...) {
  colGeometry <- .df2sf_list(colGeometry, spatialCoordsNames, spotDiameter, "POLYGON")
  if (!is.null(rowGeometry)) {
    rowGeometry <- .df2sf_list(rowGeometry, spatialCoordsNames, spotDiameter = NA,
                               geometryType = "POLYGON")
  }
  if (!is.null(annotGeometry)) {
    annotGeometry <- .df2sf_list(annotGeometry, spatialCoordsNames,
                                 spotDiameter = NA,
                                 geometryType = annotGeometryType)
  }
  spe_coords <- st_coordinates(st_centroid(colGeometry[[1]]))
  spe <- SpatialExperiment(assays = assays, colData = colData,
                           rowData = rowData, sample_id = sample_id,
                           spatialCoords = spe_coords)
  sfe <- .spe_to_sfe(spe, primaryGeometry, objectGeometry, annotationGeometry,
                     unit)
  return(sfe)
}

.spe_to_sfe <- function(spe, primaryGeometry, objectGeometry,
                        annotationGeometry, unit) {
  old <- S4Vectors:::disableValidity()
  if (!isTRUE(old)) {
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old))
  }
}
# To do: unit test, validity, getters and setters, show, subsetting, cropping with geometry
