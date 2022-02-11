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
#' @importFrom methods setClass
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importClassesFrom sf sf
#' @exportClass
setClass("SpatialFeatureExperiment",
         contains = "SpatialExperiment",
         slots = c(primaryGeometry = "list",
                   objectGeometry = "list",
                   annotationGeometry = "list"),
         prototype = list(primaryGeometry = list(spot = make_empty_geometry()),
                          objectGeometry = list(nuclei = make_empty_geometry()),
                          annotationGeometry = list(tissue = make_empty_geometry())))

#' Constructor of SpatialFeatureExperiment object
#'
#' Create a \code{SpatialFeatureExperiment} object.
#'
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @inheritParams SpatialExperiment::SpatialExperiment
#' @param primaryGeometry Geometry of the entities that correspond to the columns
#' of the gene count matrix, such as cells and Visium spots. It can be one of the
#' following:
#' \describe{
#' \item{An \code{sf} data frame}{The geometry column specifies the geometry of
#' the entities.}
#' \item{An ordinary data frame specifying centroids}{Column names for the
#' coordinates are specified in the \code{spatialCoordsNames} argument. For
#' Visium and ST, in addition to the centroid coordinate data frame, the spot
#' diameter in the same unit as the coordinates can be specified in the
#' \code{spotDiamter} argument.}
#' \item{An ordinary data frame specifying polygons}{Also use
#' \code{spatialCoordsNames}. There should an additional column "ID" to specify
#' which vertices belong to which polygon.}
#' }
#' In all cases, the data frame should specify the same number of geometries as
#' the number of columns in the gene count matrix. If the column "barcode" is
#' present, then it will be matched to column names of the gene count matrix.
#' Otherwise, the geometries are assumed to be in the same order as columns in
#' the gene count matrix. If the geometries are specified in an ordinary data
#' frame, then it will be converted into \code{sf} internally.
#' @param objectGeometry Geometry of objects other than those in
#' \code{primaryGeometry}, such as nuclei segmentation when the
#' \code{primaryGeometry} is Visium spots. It can be specified in the same manner
#' as \code{primaryGeometry}.
#' @param ... Additional arguments passed to the
#'   \code{\link{SpatialExperiment}} and \code{\link{SingleCellExperiment}}
#'   constructors.
#' @importFrom S4Vectors SimpleList
#' @export
SpatialFeatureExperiment <- function(assays, primaryGeometry,
                                     objectGeometry = make_empty_geometry(),
                                     annotationGeometry = make_empty_geometry(),
                                     colData = DataFrame(),
                                     rowData = NULL, spatialCoordsNames = c("x", "y"),
                                     sample_id = "sample01",
                                     ...) {

}
