#' The SpatialFeatureExperiment class
#'
#' This class inherits from the \code{\link[SpatialExperiment]{SpatialExperiment}} (SPE) class,
#' which in turn inherits from \code{\link[SingleCellExperiment]{SingleCellExperiment}} (SCE).
#' \code{SpatialFeatureExperiment} stores geometries of spots or cells in
#' \code{sf} objects which form columns of a \code{DataFrame} which is in turn a
#' column of the \code{int_colData} \code{DataFrame} of the underlying SCE
#' object, just like \code{reducedDim} in SCE. Geometries of the tissue outline,
#' pathologist annotations, and objects (e.g. nuclei segmentation in a Visium
#' dataset) are stored in \code{sf} objects in a named list called
#' \code{annotGeometries} in \code{int_metadata}.
#'
#' @rdname SpatialFeatureExperiment-class
#' @include utils.R
#' @importFrom methods setClass new setAs setMethod setGeneric setReplaceMethod callNextMethod is
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @exportClass SpatialFeatureExperiment
#' @concept SpatialFeatureExperiment class
setClass("SpatialFeatureExperiment", contains = "SpatialExperiment")

#' Constructor of SpatialFeatureExperiment object
#'
#' Create a \code{SpatialFeatureExperiment} object.
#'
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @inheritParams SpatialExperiment::SpatialExperiment
#' @param colGeometries Geometry of the entities that correspond to the columns
#'   of the gene count matrix, such as cells and Visium spots. It must be a
#'   named list or vector of one or more of the following:
#'   \describe{
#'   \item{An \code{sf} data frame}{The geometry column specifies the geometry of the entities.}
#'   \item{File names of GeoParquet files}{These will be used to represent geometries
#'   out of memory with \code{DuckDBSpatial}.}
#'   \item{An ordinary data frame specifying centroids or polygons}{Deprecated}
#'   }
#' @param rowGeometries Geometry associated with genes or features, which
#'   correspond to rows of the gene count matrix.
#' @param annotGeometries Geometry of entities that do not correspond to columns
#'   or rows of the gene count matrix, such as tissue boundary and pathologist
#'   annotations of histological regions, and nuclei segmentation in a Visium
#'   dataset. Also a named list as in \code{colGeometries}.
#' @param spatialGraphs A named list of \code{listw} objects (see \code{spdep})
#'   for spatial neighborhood graphs.
#' @param annotGeometryType Deprecated. Ignored if the corresponding element is 
#'   an \code{sf} object or file path of a GeoParquet file.
#' @param spatialCoordsNames A \code{character} vector of column names if
#'   \code{*Geometries} arguments have ordinary data frames, to identify the
#'   columns in the ordinary data frames that specify the spatial coordinates.
#'   If \code{colGeometries} is not specified, then this argument will behave as
#'   in \code{\link[SpatialExperiment]{SpatialExperiment}}, but \code{colGeometries} will be given
#'   precedence if provided.
#' @param spatialCoords A numeric matrix containing columns of spatial
#'   coordinates, as in \code{\link[SpatialExperiment]{SpatialExperiment}}. The coordinates are
#'   centroids of the entities represented by the columns of the gene count
#'   matrix. If \code{colGeometries} is also specified, then it will be given
#'   priority and a warning is issued. Otherwise, the \code{sf} representation
#'   of the centroids will be stored in the \code{colGeometry} called
#'   \code{centroids}.
#' @param spotDiameter Spot diameter for technologies with arrays of spots of
#'   fixed diameter per slide, such as Visium, ST, DBiT-seq, and slide-seq. The
#'   diameter must be in the same unit as the coordinates in the *Geometry
#'   arguments. Ignored for geometries that are not POINT or MULTIPOINT.
#' @param unit Unit the coordinates are in, either microns or pixels in full
#'   resolution image.
#' @param ... Additional arguments passed to the
#'   \code{\link[SpatialExperiment]{SpatialExperiment}}
#'   and \code{\link[SingleCellExperiment]{SingleCellExperiment}} constructors.
#' @return A SFE object. If neither \code{colGeometries} nor \code{spotDiameter}
#'   is specified, then a \code{colGeometry} called "centroids" will be made,
#'   which is essentially the spatial coordinates as sf POINTs. If
#'   \code{spotDiameter} is specified, but not \code{colGeometries}, then the
#'   spatial coordinates will be buffered by half the diameter to get spots with
#'   the desired diameter, and the resulting \code{colGeometry} will be called
#'   "spotPoly", for which there's a convenience getter and setter,
#'   \code{\link{spotPoly}}.
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords<-
#' @importFrom SingleCellExperiment int_colData int_elementMetadata int_metadata
#'   int_metadata<- int_elementMetadata<- int_colData<-
#' @importFrom sf st_point st_sfc st_sf st_polygon st_buffer st_linestring
#'   st_multipoint st_multilinestring st_multipolygon st_coordinates st_centroid
#'   st_geometry_type st_geometry st_is_valid st_geometrycollection
#' @importFrom S4Vectors DataFrame SimpleList
#' @concept SpatialFeatureExperiment class
#' @export
#' @examples
#' library(Matrix)
#' data("visium_row_col")
#' coords1 <- visium_row_col[visium_row_col$col < 6 & visium_row_col$row < 6, ]
#' coords1$row <- coords1$row * sqrt(3)
#' cg <- df2sf(coords1[, c("col", "row")], c("col", "row"), spotDiameter = 0.7)
#'
#' set.seed(29)
#' col_inds <- sample(seq_len(13), 13)
#' row_inds <- sample(seq_len(5), 13, replace = TRUE)
#' values <- sample(seq_len(5), 13, replace = TRUE)
#' mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
#' colnames(mat) <- coords1$barcode
#' rownames(mat) <- sample(LETTERS, 5)
#' rownames(cg) <- colnames(mat)
#'
#' sfe <- SpatialFeatureExperiment(list(counts = mat),
#'     colData = coords1,
#'     spatialCoordsNames = c("col", "row"),
#'     spotDiameter = 0.7
#' )
#' sfe2 <- SpatialFeatureExperiment(list(counts = mat),
#'     colGeometries = list(foo = cg)
#' )
SpatialFeatureExperiment <- function(assays,
                                     colData = DataFrame(), rowData = NULL,
                                     sample_id = "sample01",
                                     spatialCoordsNames = c("x", "y"),
                                     spatialCoords = NULL,
                                     colGeometries = NULL, rowGeometries = NULL,
                                     annotGeometries = NULL,
                                     spotDiameter = NA_real_,
                                     annotGeometryType = deprecated(),
                                     spatialGraphs = NULL,
                                     unit = c("full_res_image_pixel", "micron"),
                                     ...) {
    unit <- match.arg(unit)
    if (is_present(annotGeometryType)) {
        deprecate_warn("1.16.0", "SpatialFeatureExperiment(annotGeometryType)",
                       details = "Using non-spatial data frames in SFE constructor is deprecated.")
    } else annotGeometryType <- "POLYGON"
    if (!length(colData)) {
        colData <- make_zero_col_DFrame(nrow = ncol(assays[[1]]))
    }
    if (any(!spatialCoordsNames %in% names(colData))) {
        scn_use <- NULL
    } else {
        scn_use <- spatialCoordsNames
    }
    if (is.null(colGeometries)) {
        spe <- SpatialExperiment(
            assays = assays, colData = colData,
            rowData = rowData, sample_id = sample_id,
            spatialCoords = spatialCoords,
            spatialCoordsNames = scn_use, ...
        )
    } else {
        if (!is.null(spatialCoords)) {
            warning("Ignoring spatialCoords; coordinates are specified in colGeometries.")
        }
        colGeometries <- .df2sf_list(
            colGeometries, spatialCoordsNames,
            spotDiameter, "POLYGON"
        )
        spe_coords <- st_coordinates(st_centroid(st_geometry(colGeometries[[1]])))
        spe <- SpatialExperiment(
            assays = assays, colData = colData,
            rowData = rowData, sample_id = sample_id,
            spatialCoords = spe_coords,
            spatialCoordsNames = NULL, ...
        )
    }
    rownames(spatialCoords(spe)) <- colnames(assays[[1]])
    sfe <- .spe_to_sfe(
        spe, colGeometries, rowGeometries, annotGeometries,
        spatialCoordsNames, annotGeometryType,
        spatialGraphs, spotDiameter, unit
    )
    return(sfe)
}

.st_geometry_type_show <- function(x) standardGeneric(".st_geometry_type_show")
.st_geometry_type_show.sf <- function(x) as.character(st_geometry_type(x, by_geometry = FALSE))
.st_geometry_type_show.DuckDBColumn <- function(x) {
    check_installed("DuckDBSpatial")
    df <- DuckDBSpatial::st_geometry_type(head(x)) |> as.data.frame()
    unique(df[,1]) |> as.character()
}

.names_types <- function(l) {
    types <- vapply(l, .st_geometry_type_show, FUN.VALUE = character(1))
    paste(paste0(names(l), " (", types, ")"), collapse = ", ")
}

#' Get unit of a SpatialFeatureExperiment
#'
#' Length units can be microns or pixels in full resolution image in SFE
#' objects.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @return A string for the name of the unit. At present it's merely a
#'   string and \code{udunits} is not used.
#' @export
#' @aliases unit
#' @concept SpatialFeatureExperiment class
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' SpatialFeatureExperiment::unit(sfe)
setMethod("unit", "SpatialFeatureExperiment",
          function(x) int_metadata(x)$unit)

#' Print method for SpatialFeatureExperiment
#'
#' Printing summaries of \code{colGeometries}, \code{rowGeometries}, and
#' \code{annotGeometries} in addition to what's shown for
#' \code{SpatialExperiment}. Geometry names and types are printed.
#'
#' @param object A \code{SpatialFeatureExperiment} object.
#' @return None (invisible \code{NULL}).
#' @concept SpatialFeatureExperiment class
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe # The show method is implicitly called
setMethod(
    "show", "SpatialFeatureExperiment",
    function(object) {
        callNextMethod()
        cat("\nunit:", unit(object))
        cg_names <- names(int_colData(object)$colGeometries)
        rg_names <- names(int_elementMetadata(object)$rowGeometries)
        ag_names <- names(int_metadata(object)$annotGeometries)
        skip_geometries <- length(cg_names) < 1 & length(rg_names) < 1 &
            is.null(ag_names)
        if (!skip_geometries) {
            cat("\nGeometries:\n")
            if (length(cg_names)) {
                cat("colGeometries:", .names_types(colGeometries(object, withDimnames = FALSE)), "\n")
            }
            if (length(rg_names)) {
                cat("rowGeometries:", .names_types(rowGeometries(object, withDimnames = FALSE)), "\n")
            }
            if (!is.null(ag_names)) {
                cat(
                    "annotGeometries:", .names_types(annotGeometries(object)),
                    "\n"
                )
            }
        }
        # What to do with the graphs?
        if (!is.null(int_metadata(object)$spatialGraphs)) {
            cat("\nGraphs:")
            df <- int_metadata(object)$spatialGraphs
            for (s in colnames(df)) {
                cat("\n", s, ": ", sep = "")
                out <- vapply(rownames(df), function(r) {
                    l <- df[r, s][[1]]
                    if (length(l)) {
                        paste0(r, ": ", paste(names(l), collapse = ", "))
                    } else {
                        NA_character_
                    }
                }, FUN.VALUE = character(1))
                out <- out[!is.na(out)]
                cat(paste(out, collapse = "; "))
            }
        }
    }
)
