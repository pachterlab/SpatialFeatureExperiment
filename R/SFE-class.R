#' The SpatialFeatureExperiment class
#'
#' This class inherits from the \code{\link{SpatialExperiment}} class, which in
#' turn inherits from \code{\link{SingleCellExperiment}}, adding slots for
#' Simple Feature geometries. The additional slots are documented here.
#'
#' @slot primaryGeometry A named list of data frames specifying the geometry of
#'   the items that are the columns of the gene count matrix, such as cells
#'   (e.g. for MERFISH) or Visium spots. Multiple sf data frames are used when
#'   multiple geometries are associated with each item, such as whole cell
#'   segmentation, nuclei segmentation, and transcript localization. All of the
#'   sf data frames must have the same number of rows and a column called "ID"
#'   to identify each primary object such as the Visium spot barcode, the the ID
#'   column must be in the same order in all these data frames. The data frames
#'   should also have the same number and order of rows as there are columns in
#'   the gene count matrix. All geometries in the \code{*Geometry} slots must be
#'   valid geometries (see \code{\link{st_is_valid}}), as this is required for
#'   geometry operations such as intersections and unions.
#' @slot objectGeometry A named list of sf data frames specifying the geometry
#'   of objects other than those in the \code{PrimaryGeometry}, such as nuclei
#'   segmentation from the H&E image associated with a Visium dataset, where the
#'   primary geometry would be for the Visium spots. The number of objects in
#'   this slot can be different from that of \code{PrimaryGeometry}. Each sf
#'   data frame here must also have an "ID" column to identify the objects. When
#'   geometries of intersections between a \code{PrimaryGeometry} and a
#'   \code{ObjectGeometry} or between different \code{ObjectGeometry}s are
#'   computed, the results will be stored in this slot as new elements of the
#'   list.
#' @slot annotationGeometry A named list of sf data frames specifying the
#'   geometry of spatial annotations of the tissue, such as the tissue boundary
#'   and pathologist annotated tissue regions.
#' @rdname SpatialFeatureExperiment-class
#' @include utils.R
#' @importFrom methods setClass new
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @exportClass SpatialFeatureExperiment
setClass("SpatialFeatureExperiment",
         contains = "SpatialExperiment",
         slots = c(primaryGeometry = "list",
                   objectGeometry = "list",
                   annotationGeometry = "list"),
         prototype = list(primaryGeometry = list(spot = make_empty_geometry()),
                          objectGeometry = list(object = make_empty_geometry()),
                          annotationGeometry = list(region = make_empty_geometry())))

#' Constructor of SpatialFeatureExperiment object
#'
#' Create a \code{SpatialFeatureExperiment} object.
#'
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @inheritParams SpatialExperiment::SpatialExperiment
#' @param primaryGeometry Geometry of the entities that correspond to the
#'   columns of the gene count matrix, such as cells and Visium spots. It must
#'   be a named list of be one of the following: \describe{ \item{An \code{sf}
#'   data frame}{The geometry column specifies the geometry of the entities.}
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
#' @param objectGeometry Geometry of objects other than those in
#'   \code{primaryGeometry}, such as nuclei segmentation when the
#'   \code{primaryGeometry} is Visium spots. It can be specified in the same
#'   manner as \code{primaryGeometry}, as a named list. The data frames in the
#'   list don't have to specify the same number of geometries, unlike in
#'   \code{primaryGeometry}. The ordinary data frame must specify either POINTs
#'   or POLYGONs. A warning is issued if the bounding box of
#'   \code{objectGeometry} does not overlap that of \code{primaryGeometry}.
#' @param annotationGeometry Geometry of annotations of the tissue, such as
#'   tissue boundary and pathologist annotations of histological regions. Also a
#'   named list as in \code{primaryGeometry}. The ordinary data frame may
#'   specify POINTs, POLYGONs, or LINESTRINGs, or their MULTI versions. Each
#'   data frame can only specify one type of geometry. For MULTI versions, there
#'   must be a column "group" to identify each MULTI geometry.
#' @param annotationGeometryType Character vector specifying geometry type of
#'   each ordinary data frame in the list if \code{annotationGeometry} is
#'   specified as such. Each element of the vector must be one of POINT,
#'   LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, and MULTIPOLYGON. Must be
#'   the same length as the list of data frames.
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
#' @importFrom sf st_point st_sfc st_sf st_polygon st_buffer st_linestring
#'   st_multipoint st_multilinestring st_multipolygon st_is_empty st_is
#'   st_coordinates st_centroid
#' @export
SpatialFeatureExperiment <- function(assays, primaryGeometry,
                                     objectGeometry = list(object = make_empty_geometry()),
                                     annotationGeometry = list(region = make_empty_geometry()),
                                     colData = DataFrame(),
                                     rowData = NULL, spatialCoordsNames = c("x", "y"),
                                     sample_id = "sample01", spotDiameter = NA_real_,
                                     annotationGeometryType = "MULTIPOLYGON",
                                     unit = "full_res_image_pixels",
                                     ...) {
  # 1. Check the *Geometry arguments, convert them to sf if they are ordinary data frames
  # 2. Pass everything else to the SpatialExperiment constructor
  if (!st_is_empty(annotationGeometry[[1]])) {
    .annot_geom_allowed <- c("POINT", "LINESTRING", "POLYGON",
                             "MULTIPOINT", "MULTILINESTRING", "MULTIPOLYGON")
    annotationGeometryType <- match.arg(annotationGeometryType,
                                        choices = .annot_geom_allowed,
                                        several.ok = TRUE)
    if (length(annotationGeometryType) == 1L) {
      annotationGeometryType <- rep(annotationGeometryType,
                                    length(annotationGeometry))
    } else if (length(annotationGeometryType) != length(annotationGeometry)) {
      stop("annotationGeometryType must be either length 1 or the same length ",
           "as annotationGeometry.")
    }
  }
  primaryGeometry <- lapply(primaryGeometry, .df2sf_list,
                            spatialCoordsNames = spatialCoordsNames,
                            spotDiameter = spotDiameter,
                            geometryType = "POLYGON")
  if (!st_is_empty(objectGeometry[[1]])) {
    objectGeometry <- lapply(objectGeometry, .df2sf_list,
                             spatialCoordsNames = spatialCoordsNames,
                             spotDiameter = spotDiameter,
                             geometryType = "POLYGON")
  }
  if (!st_is_empty(annotationGeometry[[1]])) {
    annotationGeometry <- mapply(.df2sf_list, df = annotationGeometry,
                                 geometryType = annotationGeometryType,
                                 MoreArgs = list(spatialCoordsNames = spatialCoordsNames,
                                                 spotDiameter = NA))
  }
  spe_coords <- st_coordinates(st_centroid(primaryGeometry[[1]]))
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
# To do: unit test, validity, getters and setters, subsetting, cropping with geometry
