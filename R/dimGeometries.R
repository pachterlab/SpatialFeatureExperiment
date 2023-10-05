#' Dimension geometry methods
#'
#' "Dimension geometry" refers to Simple Feature (\code{sf}) geometries
#' associated with rows (features, genes) or columns (cells or spots) of the
#' gene count matrix in the \code{SpatialFeatureExperiment} object. For each
#' dimension, the number of rows in the \code{sf} data frame specifying the
#' geometries must match the size of the dimension of interest. For example,
#' there must be the same number of rows in the \code{sf} data frame describing
#' cells as there are cells in the gene count matrix. This page documents
#' getters and setters for the dimension geometries. The getters and setters are
#' implemented in a way similar to those of \code{reducedDims} in
#' \code{SingleCellExperiment}.
#'
#' These are convenience wrappers for getters and setters of special geometries:
#' \describe{ \item{colGeometry/ies}{dimGeometry/ies with MARGIN = 2, for
#' geometries associated with columns of the gene count matrix (cells/Visium
#' spots/samples).} \item{rowGeometry/ies}{dimGeometry/ies with MARGIN = 1, for
#' geometries associated with rows of the gene count matrix (genes/features).}
#' \item{spotPoly}{Polygons of spots from technologies such as Visium, ST, and
#' slide-seq, which do not correspond to cells. Centroids of the polygons are
#' stored in \code{spatialCoords} of the underlying \code{SpatialExperiment}
#' object.} \item{ROIPoly}{Polygons of regions of interest (ROIs) from
#' technologies such as laser capture microdissection (LCM) and GeoMX DSP. These
#' should correspond to columns of the gene count matrix.} \item{cellSeg}{Cell
#' segmentation polygons. If the columns of the gene count matrix are single
#' cells, then this is stored in \code{colGeometries}. Otherwise, this is stored
#' in \code{\link{annotGeometries}}.} \item{nucSeg}{Similar to \code{cellSeg},
#' but for nuclei rather than whole cell.} \item{txSpots}{POINT or MULTIPOINT
#' geometries of transcript spots of single molecular resolution technologies,
#' stored in \code{rowGeometries}.} }
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns.
#' @param sample_id Sample ID to get or set geometries.
#' @param withDimnames Logical. If \code{TRUE}, then the dimnames (colnames or
#'   rownames) of the gene count matrix should correspond to row names of the
#'   \code{sf} data frames of interest.
#' @param type An integer specifying the index or string specifying the name of
#'   the *Geometry to query or replace. If missing, then the first item in the
#'   *Geometries will be returned or replaced.
#' @param translate Logical. Only used if \code{\link{removeEmptySpace}} has
#'   been run of the SFE object. If that's the case, this argument indicates
#'   whether the new value to be assigned to the geometry is in the coordinates
#'   prior to removal of empty space so it should be translated to match the new
#'   coordinates after removing empty space. Default to \code{TRUE}.
#' @param value Value to set. For \code{dimGeometry}, must be a \code{sf} data
#'   frame with the same number of rows as size in the dimension of interest, or
#'   an ordinary data frame that can be converted to such a \code{sf} data frame
#'   (see \code{\link{df2sf}}). For \code{dimGeometries}, must be a list of such
#'   \code{sf} or ordinary data frames.
#' @param ... \code{spatialCoordsNames, spotDiameter, geometryType} passed to
#'   \code{\link{df2sf}}. Defaults are the same as in \code{\link{df2sf}}. For
#'   \code{dimGeometries<-} only: \code{geometryType} can be a character vector
#'   of the geometry type of each data frame in the list of the same length as
#'   the list if the data frames specify different types of geometries.
#' @return Getters for multiple geometries return a named list. Getters for names
#'   return a character vector of the names. Getters for single geometries
#'   return an \code{sf} data frame. Setters return an SFE object.
#' @concept Column or row geometries
#' @name dimGeometries
#' @aliases dimGeometries<- dimGeometry dimGeometry<- dimGeometryNames
#'   dimGeometryNames<-
#' @docType methods
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#'
#' # Get all column geometries as a named list
#' # Use MARGIN = 1 or rowGeometry/ies for rowGeometries
#' cgs <- dimGeometries(sfe, MARGIN = 2)
#' # Or equivalently
#' cgs <- colGeometries(sfe)
#'
#' # Set all column geometries with a named list
#' dimGeometries(sfe, MARGIN = 2) <- cgs
#' # Or equivalently
#' colGeometries(sfe) <- cgs
#'
#' # Get names of column geometries
#' cgns <- dimGeometryNames(sfe, MARGIN = 2)
#' cgns <- colGeometryNames(sfe)
#'
#' # Set column geometry names
#' dimGeometryNames(sfe, MARGIN = 2) <- cgns
#' colGeometryNames(sfe) <- cgns
#'
#' # Get a specific column geometry by name
#' spots <- dimGeometry(sfe, "spotPoly", MARGIN = 2)
#' spots <- colGeometry(sfe, "spotPoly")
#' # Or equivalently, the wrapper specifically for Visium spot polygons,
#' # for the name "spotPoly"
#' spots <- spotPoly(sfe)
#' # Other colGeometry wrappers for specific names:
#' # ROIPoly (for LCM and GeoMX DSP), cellSeg and nucSeg (for MERFISH; would
#' # query annotGeometries for Visium)
#' # rowGeometry wrappers for specific names: txSpots (MERFISH transcript spots)
#' # By index
#' spots <- colGeometry(sfe, 1L)
#'
#' # Multiple samples, only get geometries for one sample
#' sfe2 <- McKellarMuscleData("small2")
#' sfe_combined <- cbind(sfe, sfe2)
#' spots1 <- colGeometry(sfe, "spotPoly", sample_id = "Vis5A")
#' spots2 <- spotPoly(sfe_combined, sample_id = "sample02")
#' # Get geometries for multiple samples
#' spots3 <- spotPoly(sfe_combined, sample_id = c("Vis5A", "sample02"))
#' # All samples
#' spots3 <- spotPoly(sfe_combined, sample_id = "all")
#'
#' # Set specific column geometry by name
#' colGeometry(sfe, "foobar") <- spots
#' # Or use wrapper
#' spotPoly(sfe) <- spots
#' # Specify sample_id
#' colGeometry(sfe_combined, "foobar", sample_id = "Vis5A") <- spots1
#' # Only entries for the specified sample are set.
#' foobar <- colGeometry(sfe_combined, "foobar", sample_id = "sample02")
NULL

#' @rdname dimGeometries
#' @export
setMethod(
    "dimGeometries", "SpatialFeatureExperiment",
    function(x, MARGIN = 2, withDimnames = TRUE) {
        .get_intdimdata_all(x, MARGIN, withDimnames,
            getfun = .getfun(MARGIN),
            key = .dg_key(MARGIN)
        )
    }
)

#' @rdname dimGeometries
#' @export
setReplaceMethod(
    "dimGeometries", "SpatialFeatureExperiment",
    function(x, MARGIN, withDimnames = TRUE, translate = TRUE, ...,
             value) {
        .set_intdimdata_all(x, MARGIN,
            withDimnames = withDimnames,
            translate = translate, sf = TRUE,
            getfun = .getfun(MARGIN),
            setfun = .setfun(MARGIN),
            key = .dg_key(MARGIN),
            xdimfun = .xdimfun(MARGIN),
            funstr = "dimGeometries",
            xdimstr = .xdimstr(MARGIN), value, ...
        )
    }
)

#' @rdname dimGeometries
#' @export
setMethod(
    "dimGeometryNames", "SpatialFeatureExperiment",
    function(x, MARGIN) {
        .get_internal_names(x,
            getfun = .getfun(MARGIN),
            key = .dg_key(MARGIN)
        )
    }
)

#' @rdname dimGeometries
#' @export
setReplaceMethod(
    "dimGeometryNames",
    c("SpatialFeatureExperiment", "numeric", "character"),
    function(x, MARGIN, value) {
        .set_internal_names(x, value,
            getfun = .getfun(MARGIN),
            setfun = .setfun(MARGIN),
            key = .dg_key(MARGIN)
        )
    }
)

#' @rdname dimGeometries
#' @export
setMethod(
    "dimGeometry", "SpatialFeatureExperiment",
    function(x, type = 1L, MARGIN, sample_id = NULL, withDimnames = TRUE) {
        .get_internal_id(x, type, MARGIN, sample_id, withDimnames,
            .get_internal,
            getfun = .getfun(MARGIN),
            key = .dg_key(MARGIN), funstr = "dimGeometry",
            substr = "type", namestr = .dg_key2(MARGIN)
        )
    }
)

#' @rdname dimGeometries
#' @export
setReplaceMethod(
    "dimGeometry", "SpatialFeatureExperiment",
    function(x, type = 1L, MARGIN, sample_id = NULL, withDimnames = TRUE,
             translate = TRUE, ..., value) {
        .set_internal_id(x, type, MARGIN, sample_id, withDimnames,
            translate,
            sf = TRUE,
            .get_all_fun = dimGeometries,
            .set_all_fun = `dimGeometries<-`,
            .set_internal_fun = .set_internal,
            getfun = .getfun(MARGIN),
            setfun = .setfun(MARGIN),
            key = .dg_key(MARGIN),
            xdimfun = .xdimfun(MARGIN),
            funstr = "dimGeometry",
            xdimstr = .xdimstr(MARGIN),
            substr = "type", value, ...
        )
    }
)

#' @rdname dimGeometries
#' @export
colGeometry <- function(x, type = 1L, sample_id = NULL, withDimnames = TRUE) {
    dimGeometry(x, type,
        MARGIN = 2, sample_id = sample_id,
        withDimnames = withDimnames
    )
}

#' @rdname dimGeometries
#' @export
`colGeometry<-` <- function(x, type = 1L, sample_id = NULL, withDimnames = TRUE,
                            translate = TRUE, value) {
    `dimGeometry<-`(x, type,
        MARGIN = 2, sample_id = sample_id,
        withDimnames = withDimnames, translate = translate,
        value = value
    )
}

#' @rdname dimGeometries
#' @export
colGeometries <- function(x, withDimnames = TRUE) {
    dimGeometries(x, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`colGeometries<-` <- function(x, withDimnames = TRUE, translate = TRUE, value) {
    `dimGeometries<-`(x,
        MARGIN = 2, withDimnames = withDimnames,
        translate = translate, value = value
    )
}

#' @rdname dimGeometries
#' @export
colGeometryNames <- function(x) {
    dimGeometryNames(x, MARGIN = 2)
}

#' @rdname dimGeometries
#' @export
`colGeometryNames<-` <- function(x, value) {
    dimGeometryNames(x, MARGIN = 2) <- value
    x
}

#' @rdname dimGeometries
#' @export
rowGeometry <- function(x, type = 1L, withDimnames = TRUE) {
    dimGeometry(x, type,
        MARGIN = 1, sample_id = "all",
        withDimnames = withDimnames
    )
}

#' @rdname dimGeometries
#' @export
`rowGeometry<-` <- function(x, type = 1L, withDimnames = TRUE,
                            translate = TRUE, value) {
    `dimGeometry<-`(x, type,
        MARGIN = 1, sample_id = "all",
        withDimnames = withDimnames, translate = translate,
        value = value
    )
    x
}

#' @rdname dimGeometries
#' @export
rowGeometries <- function(x, withDimnames = TRUE) {
    dimGeometries(x, MARGIN = 1, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`rowGeometries<-` <- function(x, withDimnames = TRUE, translate = TRUE, value) {
    dimGeometries(x,
        MARGIN = 1, withDimnames = withDimnames,
        translate = translate
    ) <- value
    x
}

#' @rdname dimGeometries
#' @export
rowGeometryNames <- function(x) {
    dimGeometryNames(x, MARGIN = 1)
}

#' @rdname dimGeometries
#' @export
`rowGeometryNames<-` <- function(x, value) {
    dimGeometryNames(x, MARGIN = 1) <- value
    x
}

#' @rdname dimGeometries
#' @export
spotPoly <- function(x, sample_id = NULL, withDimnames = TRUE) {
    colGeometry(x, "spotPoly", withDimnames, sample_id = sample_id)
}

#' @rdname dimGeometries
#' @export
`spotPoly<-` <- function(x, sample_id = NULL, withDimnames = TRUE,
                         translate = TRUE, value) {
    colGeometry(x, "spotPoly", withDimnames,
        sample_id = sample_id,
        translate = translate
    ) <- value
    x
}

#' @rdname dimGeometries
#' @export
centroids <- function(x, sample_id = NULL, withDimnames = TRUE) {
    colGeometry(x, "centroids", withDimnames, sample_id = sample_id)
}

#' @rdname dimGeometries
#' @export
`centroids<-` <- function(x, sample_id = NULL, withDimnames = TRUE,
                         translate = TRUE, value) {
    colGeometry(x, "centroids", withDimnames,
                sample_id = sample_id,
                translate = translate
    ) <- value
    x
}

#' Add Visium spot polygons to colGeometry
#'
#' For adding the spot polygons to SFE objects converted from SPE.
#'
#' @inheritParams dimGeometries
#' @inheritParams SpatialFeatureExperiment
#' @return A SFE object with a new colGeometry called spotPoly, which has
#' polygons of the spots.
#' @importFrom SpatialExperiment spatialCoordsNames
#' @concept Column or row geometries
#' @export
#' @examples
#' library(SpatialExperiment)
#' example(read10xVisium)
#' # There can't be suplicate barcodes
#' colnames(spe) <- make.unique(colnames(spe), sep = "-")
#' rownames(spatialCoords(spe)) <- colnames(spe)
#' sfe <- toSpatialFeatureExperiment(spe)
#' # A hypothetical spot diameter; check the scalefactors_json.json file for
#' # actual diameter in pixels in full resolution image.
#' sfe <- addVisiumSpotPoly(sfe, spotDiameter = 80)
addVisiumSpotPoly <- function(x, spotDiameter) {
    df <- as.data.frame(spatialCoords(x))
    rownames(df) <- colnames(x)
    spotPoly(x, sample_id = "all", translate = FALSE) <-
        df2sf(df, names(df), spotDiameter = spotDiameter, geometryType = "POINT")
    x
}

#' @rdname dimGeometries
#' @export
ROIPoly <- function(x, sample_id = NULL, withDimnames = TRUE) {
    colGeometry(x, "ROIPoly", withDimnames, sample_id = sample_id)
}

#' @rdname dimGeometries
#' @export
`ROIPoly<-` <- function(x, sample_id = NULL, withDimnames = TRUE,
                        translate = TRUE, value) {
    colGeometry(x, "ROIPoly", withDimnames,
        sample_id = sample_id,
        translate = translate
    ) <- value
    x
}

.get_col_then_annot <- function(x, name, sample_id, withDimnames) {
    if (name %in% colGeometryNames(x)) {
        colGeometry(x, name, sample_id, withDimnames)
    } else {
        annotGeometry(x, name, sample_id)
    }
}

.set_col_then_annot <- function(x, name, sample_id, withDimnames, translate,
                                value) {
    if (nrow(value) == ncol(x)) {
        colGeometry(x, name, sample_id, withDimnames, translate) <- value
    } else {
        annotGeometry(x, name, sample_id, translate) <- value
    }
    return(x)
}

#' @rdname dimGeometries
#' @export
cellSeg <- function(x, sample_id = NULL, withDimnames = TRUE) {
    .get_col_then_annot(x, "cellSeg", sample_id, withDimnames)
}

#' @rdname dimGeometries
#' @export
`cellSeg<-` <- function(x, sample_id = NULL, withDimnames = TRUE,
                        translate = TRUE, value) {
    .set_col_then_annot(x, "cellSeg", sample_id, withDimnames,
        translate = translate, value
    )
}

#' @rdname dimGeometries
#' @export
nucSeg <- function(x, sample_id = NULL, withDimnames = TRUE) {
    .get_col_then_annot(x, "nucSeg", sample_id, withDimnames)
}

#' @rdname dimGeometries
#' @export
`nucSeg<-` <- function(x, sample_id = NULL, withDimnames = TRUE,
                       translate = TRUE, value) {
    .set_col_then_annot(x, "nucSeg", sample_id, withDimnames, translate, value)
}

#' @rdname dimGeometries
#' @export
txSpots <- function(x, withDimnames = TRUE) {
    rowGeometry(x, "txSpots", withDimnames)
}

#' @rdname dimGeometries
#' @export
`txSpots<-` <- function(x, withDimnames = TRUE, translate = TRUE, value) {
    rowGeometry(x, "txSpots", withDimnames, translate) <- value
    x
}
