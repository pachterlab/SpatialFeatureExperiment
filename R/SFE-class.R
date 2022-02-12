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
#'   look like it is passed to \code{ggplot2::geom_polygon}}. \item{Path to a
#'   GeoJSON file}{For QuPath annotations.} } In all cases, the data frame
#'   should specify the same number of geometries as the number of columns in
#'   the gene count matrix. If the column "barcode" is present, then it will be
#'   matched to column names of the gene count matrix. Otherwise, the geometries
#'   are assumed to be in the same order as columns in the gene count matrix. If
#'   the geometries are specified in an ordinary data frame, then it will be
#'   converted into \code{sf} internally. Named list of data frames because each
#'   entity can have multiple geometries, such as whole cell and nuclei
#'   segmentations. The geometries are assumed to be POINTs for centroids and
#'   POLYGONs for segmentations. If polygons are specified in an ordinary data
#'   frame, then anything with fewer than 3 vertices will be removed. For anything
#'   other than POINTs, attributes of the geometry will be ignored, but they can
#'   be added later with the \code{AddGeometryAttr} function.
#' @param objectGeometry Geometry of objects other than those in
#'   \code{primaryGeometry}, such as nuclei segmentation when the
#'   \code{primaryGeometry} is Visium spots. It can be specified in the same
#'   manner as \code{primaryGeometry}, as a named list. The data frames in the
#'   list don't have to specify the same number of geometries, unlike in
#'   \code{primaryGeometry}. The ordinary data frame must specify either POINTs
#'   or POLYGONs. A warning is issued if an \code{objectGeometry} does not
#'   overlap \code{primaryGeometry}.
#' @param annotationGeometry Geometry of annotations of the tissue, such as
#'   tissue boundary and pathologist annotations of histological regions. Also a
#'   named list as in \code{primaryGeometry}. The ordinary data frame may
#'   specify POINTs, POLYGONs, or LINESTRINGs. In the ordinary data frame, there
#'   can be a column "geometry_type" specifying the geometry type. If that
#'   column is absent, then any geometry with one vertice would be made a POINT,
#'   any geometry with 2 vertices would be made a LINESTRING, and any geometry
#'   with at least 3 vertices would be made a POLYGON.
#' @param spotDiameter Spot diameter for technologies with arrays of spots of
#'   fixed diameter per slide, such as Visium, ST, DBiT-seq, and slide-seq. The
#'   diameter must be in the same unit as the coordinates in the *Geometry
#'   arguments.
#' @param unit Unit the coordinates are in. I'm thinking about using some custom
#'   engineering CRS's which can convert units and invert the y axis for
#'   Cartesian vs. image orientations.
#' @param ... Additional arguments passed to the \code{\link{SpatialExperiment}}
#'   and \code{\link{SingleCellExperiment}} constructors.
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom sf st_point st_sfc
#' @export
SpatialFeatureExperiment <- function(assays, primaryGeometry,
                                     objectGeometry = list(object = make_empty_geometry()),
                                     annotationGeometry = list(region = make_empty_geometry()),
                                     colData = DataFrame(),
                                     rowData = NULL, spatialCoordsNames = c("x", "y"),
                                     sample_id = "sample01", spotDiameter = NA_real_,
                                     unit = "full_res_image_pixels",
                                     ...) {
  # SCE will check the validity of assays

}

# From ordinary data frame to sf to construct SFE object
.df2sf <- function(df, spatialCoordsNames, spotDiameter, allPolygons) {
  if (any(!spatialCoordsNames %in% names(df))) {
    cols_absent <- setdiff(spatialCoordsNames, names(df))
    if (length(cols_absent) > 1L) {
      stop("Columns ", paste(cols_absent, collapse = ", "), " are absent.")
    } else {
      stop("Column ", cols_absent, " is absent.")
    }
  }
  if (!"ID" %in% names(df) || !anyDuplicated(df$ID)) {
    # Case 1: centroids, use POINT
    df$geometry <- lapply(seq_len(nrow(df)), function(i) {
      st_point(c(df[[spatialCoordsNames[1]]], df[[spatialCoordsNames[2]]]))
    })
    df$geometry <- st_sfc(df$geometry)
    df <- st_sf(df, sf_column_name = "geometry")
  } else {
    # Case 2: polygons
    if (!"ID" %in% names(df)) {
      stop("Column 'ID' must be present when the geometry is other than POINT.")
    }
    if (any(!names(df) %in% c("ID", spatialCoordsNames))) {
      warning("Geometry attributes are ignored.")
    }
    n_vertices <- table(df$ID)
    ids <- names(n_vertices)
    if (allPolygons) {
      ids_rm <- names(n_vertices[n_vertices < 3])
      df <- df[!df$ID %in% ids_rm,]
      if (!nrow(df)) {
        stop("All geometries have 2 or fewer vertices. Cannot construct polygons.")
      }
    }
    df_split <- split(df, df$ID)
    geometry_use <- lapply(ids, function(i) {

    })
  }

}
