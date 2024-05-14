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
#' @concept Getters and setters
#' @name dimGeometries
#' @seealso [colGeometries()], [rowGeometries()]
#' @aliases dimGeometries<- dimGeometry dimGeometry<- dimGeometryNames
#'   dimGeometryNames<-
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
    function(x, type = 1L, MARGIN, sample_id = 1L, withDimnames = TRUE) {
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
    function(x, type = 1L, MARGIN, sample_id = 1L, withDimnames = TRUE,
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

#' Column geometry getters and setters
#'
#' \code{colGeometries} are geometries that correspond to columns of the gene
#' count matrix, such as Visium spots or cells. Same as \code{dimGeometry(x,
#' MARGIN = 2L, ...)}, with convenience wrappers for getters and setters of
#' special geometries:
#' \describe{
#' \item{spotPoly}{Polygons of spots from technologies such as Visium, ST, and
#' slide-seq, which do not correspond to cells. Centroids of the polygons are
#' stored in \code{spatialCoords} of the underlying \code{SpatialExperiment}
#' object.}
#' \item{ROIPoly}{Polygons of regions of interest (ROIs) from
#' technologies such as laser capture microdissection (LCM) and GeoMX DSP. These
#' should correspond to columns of the gene count matrix.}
#' \item{cellSeg}{Cell segmentation polygons. If the columns of the gene count
#' matrix are single cells, then this is stored in \code{colGeometries}.
#' Otherwise, this is stored in \code{\link{annotGeometries}}.}
#' \item{nucSeg}{Similar to \code{cellSeg}, but for nuclei rather than whole
#' cell.}}
#'
#' @inheritParams dimGeometries
#' @name colGeometries
#' @concept Getters and setters
#' @seealso [dimGeometries()], [rowGeometries()]
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' cgs <- colGeometries(sfe)
#' spots <- spotPoly(sfe)
NULL

#' @rdname colGeometries
#' @export
colGeometry <- function(x, type = 1L, sample_id = 1L, withDimnames = TRUE) {
    dimGeometry(x, type,
        MARGIN = 2, sample_id = sample_id,
        withDimnames = withDimnames
    )
}

#' @rdname colGeometries
#' @export
`colGeometry<-` <- function(x, type = 1L, sample_id = 1L, withDimnames = TRUE,
                            translate = TRUE, value) {
    `dimGeometry<-`(x, type,
        MARGIN = 2, sample_id = sample_id,
        withDimnames = withDimnames, translate = translate,
        value = value
    )
}

#' @rdname colGeometries
#' @export
colGeometries <- function(x, withDimnames = TRUE) {
    dimGeometries(x, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname colGeometries
#' @export
`colGeometries<-` <- function(x, withDimnames = TRUE, translate = TRUE, value) {
    `dimGeometries<-`(x,
        MARGIN = 2, withDimnames = withDimnames,
        translate = translate, value = value
    )
}

#' @rdname colGeometries
#' @export
colGeometryNames <- function(x) {
    dimGeometryNames(x, MARGIN = 2)
}

#' @rdname colGeometries
#' @export
`colGeometryNames<-` <- function(x, value) {
    dimGeometryNames(x, MARGIN = 2) <- value
    x
}

#' @rdname colGeometries
#' @export
spotPoly <- function(x, sample_id = 1L, withDimnames = TRUE) {
    colGeometry(x, "spotPoly", withDimnames, sample_id = sample_id)
}

#' @rdname colGeometries
#' @export
`spotPoly<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
                         translate = TRUE, value) {
    colGeometry(x, "spotPoly", withDimnames,
                sample_id = sample_id,
                translate = translate
    ) <- value
    x
}

#' @rdname colGeometries
#' @export
centroids <- function(x, sample_id = 1L, withDimnames = TRUE) {
    colGeometry(x, "centroids", withDimnames, sample_id = sample_id)
}

#' @rdname colGeometries
#' @export
`centroids<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
                          translate = TRUE, value) {
    colGeometry(x, "centroids", withDimnames,
                sample_id = sample_id,
                translate = translate
    ) <- value
    x
}

#' @rdname colGeometries
#' @export
ROIPoly <- function(x, sample_id = 1L, withDimnames = TRUE) {
    colGeometry(x, "ROIPoly", withDimnames, sample_id = sample_id)
}

#' @rdname colGeometries
#' @export
`ROIPoly<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
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

#' @rdname colGeometries
#' @export
cellSeg <- function(x, sample_id = 1L, withDimnames = TRUE) {
    .get_col_then_annot(x, "cellSeg", sample_id, withDimnames)
}

#' @rdname colGeometries
#' @export
`cellSeg<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
                        translate = TRUE, value) {
    .set_col_then_annot(x, "cellSeg", sample_id, withDimnames,
                        translate = translate, value
    )
}

#' @rdname colGeometries
#' @export
nucSeg <- function(x, sample_id = 1L, withDimnames = TRUE) {
    .get_col_then_annot(x, "nucSeg", sample_id, withDimnames)
}

#' @rdname colGeometries
#' @export
`nucSeg<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
                       translate = TRUE, value) {
    .set_col_then_annot(x, "nucSeg", sample_id, withDimnames, translate, value)
}

#' Row geometry getters and setters
#'
#' \code{rowGeometries} are geometries that corresponding to rows of the gene
#' count matrix, such as smFISH transcript spots. The \code{txSpots()} function
#' is a convenience wrapper for transcript spots, although this entirely depends
#' on the \code{rowGeometry} being named \code{txSpots}.
#'
#' When there are multiple samples in the SFE object, \code{rowGeometries} for
#' each sample has the \code{sample_id} appended to the name of the geometry.
#' For example, if the name is \code{txSpots} and the sample ID is
#' \code{sample01}, then the actual name of the \code{rowGeometry} is
#' \code{txSpots_sample01}. In the getter, one can still specify
#' \code{rowGeometry(sfe, "txSpots", sample_id = "sample01")}.
#'
#' Appending the \code{sample_id} is unnecessary when there is only one sample,
#' but \code{sample_id} will be appended when to SFE objects are combined with
#' \code{cbind}. It is necessary to distinguish bewteen different samples
#' because they can have overlapping coordinate values.
#'
#' @inheritParams dimGeometries
#' @param partial In setters, if a \code{rowGeometry} of the same name exists,
#'   whether to only replace the rows present in \code{value}.
#' @name rowGeometries
#' @concept Getters and setters
#' @seealso [dimGeometries()], [colGeometries()]
#' @examples
#' library(SFEData)
#' library(RBioFormats)
#' fp <- tempdir()
#' dir_use <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
#' # RBioFormats issue
#' try(sfe <- readXenium(dir_use, add_molecules = TRUE))
#' sfe <- readXenium(dir_use, add_molecules = TRUE)
#' rowGeometries(sfe)
#' rowGeometryNames(sfe)
#' tx <- rowGeometry(sfe, "txSpots")
#' txSpots(sfe)
#' unlink(dir_use, recursive = TRUE)
NULL

.check_rg_type <- function(type, x, sample_id) {
    rg_names <- rowGeometryNames(x)
    if (is.numeric(type)) {
        if (identical(sample_id, "all") || length(sampleIDs(x)) == 1L) {
            type <- rg_names[type]
        } else {
            type <- rg_names[grepl(paste0(sample_id, "$"), rg_names)][type]
        }
    }
    # Deal with: no sample_id in the name, meaning for all samples
    # Have sample_id, but it's not the sample requested
    # What if there's both name and name_sample?
    # the name in type doesn't have to exist in the SFE object,
    # totally fine for setter and is later checked for getter
    other_samples <- setdiff(sampleIDs(x), sample_id)
    if (length(other_samples)) {
        patterns <- paste0(other_samples, "$")
        has_other <- vapply(patterns, grepl, x = type, FUN.VALUE = logical(1))
        if (any(has_other)) {
            stop("Type does not match sample_id")
        }
    }
    if (length(sampleIDs(x)) > 1L && !grepl(paste0(sample_id, "$"), type)) {
        type <- paste(type, sample_id, sep = "_")
    }
    type
}

.check_rg_setter <- function(type, x, sample_id) {
    if (identical(sample_id, "all")) {
        .check_rg_sample_all(type, x)
        return(type)
    }
    sample_id <- .check_sample_id(x, sample_id, TRUE)
    rg_names <- rowGeometryNames(x)
    if (is.numeric(type) && type[1] > length(rg_names)) {
        stop("subscript out of bound for numeric type")
    }
    if (is.numeric(type)) type <- rg_names[type]
    # What do I want from the single geometry setter?
    # When only one sample exists: can use the name as is
    if (length(sampleIDs(x)) == 1L) return(type)
    # When there're multiple samples in sfe: sample must be specified if not sample == "all"
    other_samples <- setdiff(sampleIDs(x), sample_id)
    if (length(other_samples)) {
        patterns <- paste0(other_samples, "$")
        has_other <- vapply(patterns, grepl, x = type, FUN.VALUE = logical(1))
        if (any(has_other)) {
            stop("Type does not match sample_id")
        }
    }
    if (!grepl(paste0(sample_id, "$"), type)) {
        type <- paste(type, sample_id, sep = "_")
    }
    type
}

.check_rg_sample_all <- function(type, x) {
    # Make sure no sample_ids in the name if identical(sample_id, "all")
    rg_names <- rowGeometryNames(x)
    if (is.numeric(type)) type <- rg_names[type]
    patterns <- paste0(sampleIDs(x), "$")
    has_other <- vapply(patterns, grepl, x = type, FUN.VALUE = logical(1))
    if (any(has_other)) {
        stop("Name of rowGeometry for all samples should not include any sample ID.")
    }
}

.check_rg <- function(type, x, sample_id) {
    if (identical(sample_id, "all")) {
        .check_rg_sample_all(type, x)
    } else if (!identical(sample_id, "all")) {
        sample_id <- .check_sample_id(x, sample_id, TRUE)
        # By convention, should be name_sample to distinguish between samples for
        # rowGeometries of the same name
        type <- .check_rg_type(type, x, sample_id)
    }
    type
}

#' @rdname rowGeometries
#' @export
rowGeometry <- function(x, type = 1L, sample_id = 1L, withDimnames = TRUE) {
    type <- .check_rg(type, x, sample_id)
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    dimGeometry(x, type,
                MARGIN = 1, sample_id = sample_id,
                withDimnames = withDimnames
    )
}

#' @rdname rowGeometries
#' @export
`rowGeometry<-` <- function(x, type = 1L, sample_id = 1L, withDimnames = TRUE,
                            partial = FALSE, translate = TRUE, value) {
    type <- .check_rg_setter(type, x, sample_id)
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (partial && type %in% rowGeometryNames(x)) {
        if (!withDimnames)
            stop("withDimnames must be TRUE for partial replace")
        existing <- rowGeometry(x, type, sample_id)
        cols_use <- intersect(names(existing), names(value))
        rows_use <- intersect(rownames(existing), rownames(value))
        existing[rows_use, cols_use] <- value[rows_use, cols_use]
        value <- existing
    }
    `dimGeometry<-`(x, type,
        MARGIN = 1, sample_id = sample_id,
        withDimnames = withDimnames, translate = translate,
        value = value
    )
}

.get_rg_multi <- function(x, sample_id) {
    if (length(sample_id) > 1L) {
        pattern <- paste0("(", paste(sample_id, collapse = ")|("), ")")
    } else pattern <- paste0(sample_id, "$")
    rgns <- rowGeometryNames(x)
    rgns[grepl(pattern, rgns)]
}

.check_rg_multi_sample <- function(rgns, sample_id) {
    if (length(sample_id) > 1L) {
        pattern <- paste0("(", paste(sample_id, collapse = ")|("), ")")
    } else pattern <- paste0(sample_id, "$")
    has_sample <- grepl(pattern, rgns)
    if (!all(has_sample)) {
        if (length(sample_id) == 1L)
            rgns[!has_sample] <- paste(rgns[!has_sample], sample_id, sep = "_")
        else
            stop("sample_id should be indicated in the names of rowGeometries")
    }
    rgns

}

#' @rdname rowGeometries
#' @export
rowGeometries <- function(x, sample_id = "all", withDimnames = TRUE) {
    out <- dimGeometries(x, MARGIN = 1, withDimnames = withDimnames)
    if (identical(sample_id, "all") || length(sampleIDs(x)) == 1L) return(out)
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    rgns <- .get_rg_multi(x, sample_id)
    out[rgns]
}

#' @rdname rowGeometries
#' @export
`rowGeometries<-` <- function(x, sample_id = "all", withDimnames = TRUE,
                              partial = FALSE, translate = TRUE, value) {
    check_names <- !identical(sample_id, "all") && length(sampleIDs(x)) > 1L
    sample_id0 <- sample_id
    sample_id <- .check_sample_id(x, sample_id, one = FALSE, mustWork = FALSE)
    existing <- rowGeometries(x, sample_id = "all")
    # Set to NULL
    if (is.null(value) && check_names) {
        rgns_rm <- .get_rg_multi(x, sample_id)
        value <- existing[setdiff(names(existing), rgns_rm)]
    } else if (!is.null(value)) {
        if (check_names) names(value) <- .check_rg_multi_sample(names(value), sample_id)
        partial_names <- intersect(names(value), names(existing))
        if (partial && length(partial_names)) {
            for (p in partial_names) {
                rowGeometry(x, type = p, sample_id = sample_id0,
                            partial = TRUE, withDimnames = TRUE) <- value[[p]]
            }
            return(x)
        }
        existing <- existing[setdiff(names(existing), names(value))]
        value <- c(existing, value)
    }
    dimGeometries(x,
        MARGIN = 1, withDimnames = withDimnames,
        translate = translate
    ) <- value
    x
}

#' @rdname rowGeometries
#' @export
rowGeometryNames <- function(x) {
    dimGeometryNames(x, MARGIN = 1)
}

#' @rdname rowGeometries
#' @export
`rowGeometryNames<-` <- function(x, value) {
    dimGeometryNames(x, MARGIN = 1) <- value
    x
}

#' @rdname rowGeometries
#' @export
txSpots <- function(x, sample_id = 1L, withDimnames = TRUE) {
    rowGeometry(x, "txSpots", sample_id, withDimnames)
}

#' @rdname rowGeometries
#' @export
`txSpots<-` <- function(x, sample_id = 1L, withDimnames = TRUE,
                        partial = FALSE, translate = TRUE, value) {
    rowGeometry(x, "txSpots", sample_id, withDimnames, partial, translate) <- value
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
#' @concept Geometric operations
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
