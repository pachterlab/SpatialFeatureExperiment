.initDF <- function(m) {
    rownames_use <- colnames(m)
    fd <- make_zero_col_DFrame(nrow = ncol(m))
    rownames(fd) <- rownames_use
    fd
}

#' @importFrom S4Vectors combineCols
.format_fd <- function(x, MARGIN, value = NULL) {
    fd_name <- "featureData"
    dimData <- switch(MARGIN, rowData, colData)
    if (is.null(value)) fd <- metadata(dimData(x))[[fd_name]] else fd <- value
    if (!is.null(fd)) {
        fd <- fd[intersect(rownames(fd), colnames(dimData(x))),, drop = FALSE]
        empty <- .initDF(dimData(x))
        fd <- combineCols(empty, fd)
    }
    return(fd)
}

#' Get global spatial analysis results and metadata of colData, rowData, and geometries
#'
#' Results of spatial analyses on columns in \code{colData}, \code{rowData}, and
#' geometries are stored in their metadata, which can be accessed by the
#' \code{\link{metadata}} function. The \code{colFeaturedata} function allows
#' the users to more directly access these results.
#'
#' @param sfe An SFE object.
#' @param type Which geometry, can be name (character) or index (integer)
#' @param MARGIN Integer, 1 means rowGeometry, 2 means colGeometry, and 3 means
#'   annotGeometry. Defaults to 2, colGeometry.
#' @param dimred Name of a dimension reduction, can be seen in
#'   \code{\link{reducedDimNames}}.
#' @concept Getters and setters
#' @return A \code{DataFrame}.
#' @seealso getParams
#' @export
#' @name colFeatureData
#' @examples
#' library(SpatialFeatureExperiment)
#' library(SingleCellExperiment)
#' library(SFEData)
#' library(Voyager)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' # Moran's I for colData
#' sfe <- colDataMoransI(sfe, "nCounts")
#' colFeatureData(sfe)
colFeatureData <- function(sfe) {
    .format_fd(sfe, 2L)
}

`colFeatureData<-` <- function(sfe, value) {
    if (!is.null(value)) value <- .format_fd(sfe, 2L, value)
    metadata(colData(sfe))$featureData <- value
    sfe
}

#' @rdname colFeatureData
#' @export
rowFeatureData <- function(sfe) {
    .format_fd(sfe, 1L)
}

#' @rdname colFeatureData
#' @export
geometryFeatureData <- function(sfe, type, MARGIN = 2L) {
    geo_fun <- switch (MARGIN, rowGeometry, colGeometry, annotGeometry)
    df <- geo_fun(sfe, type, sample_id = "all")
    attr(df, "featureData")
}

`geometryFeatureData<-` <- function(sfe, type, MARGIN = 2L, value) {
    geo_fun <- switch (MARGIN, rowGeometry, colGeometry, annotGeometry)
    geo_fun_setter <- switch (MARGIN, `rowGeometry<-`, `colGeometry<-`, `annotGeometry<-`)
    df <- geo_fun(sfe, type, sample_id = "all")
    attr(df, "featureData") <- value
    geo_fun_setter(sfe, type, sample_id = "all", value = df)
}

#' @rdname colFeatureData
#' @export
reducedDimFeatureData <- function(sfe, dimred) {
    attr(reducedDim(sfe, dimred), "featureData")
}

`reducedDimFeatureData<-` <- function(sfe, dimred, value) {
    attr(reducedDim(sfe, dimred), "featureData") <- value
    sfe
}

#' Get parameters used in spatial methods
#'
#' The \code{getParams} function allows users to access the parameters used to
#' compute the results that may be stored in \code{\link{colFeatureData}}.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param name Name used to store the results.
#' @param local Logical, whether the results of interest come from a local
#'   spatial method.
#' @param colData Logical, whether the results were computed for a column of
#'   \code{colData(sfe)}.
#' @param colGeometryName To get results for a \code{colGeometry}.
#' @param annotGeometryName To get results for an \code{annotGeometry};
#'   \code{colGeometry} has precedence so this argument is ignored if
#'   \code{colGeometryName} is specified.
#' @param reducedDimName Name of a dimension reduction, can be seen in
#'   \code{\link{reducedDimNames}}. \code{colGeometryName} and
#'   \code{annotGeometryName} have precedence over \code{reducedDimName}.
#' @return A named list showing the parameters
#' @concept Getters and setters
#' @export
#' @examples
#' library(SFEData)
#' library(scater)
#' library(Voyager)
#' sfe <- McKellarMuscleData("small")
#' colGraph(sfe, "visium") <- findVisiumGraph(sfe)
#' sfe <- colDataMoransI(sfe, "nCounts")
#' getParams(sfe, "moran", colData = TRUE)
getParams <- function(sfe, name, local = FALSE, colData = FALSE,
                      colGeometryName = NULL, annotGeometryName = NULL,
                      reducedDimName = NULL) {
    if (local) {
        if (is.null(colGeometryName)) {
            if (is.null(annotGeometryName)) {
                lr <- int_colData(sfe)$localResults
                if (is.null(lr)) return(NULL)
                return(metadata(lr)$params[[name]])
            } else {
                ag <- annotGeometry(sfe, annotGeometryName, "all")
                lr <- ag$localResults
                return(attr(lr, "params")[[name]])
            }
        } else {
            cg <- colGeometry(sfe, colGeometryName, "all")
            lr <- cg$localResults
            return(attr(lr, "params")[[name]])
        }
    } else {
        if (colData) {
            metadata(colData(sfe))$params[[name]]
        } else if (is.null(colGeometryName)) {
            if (is.null(annotGeometryName) && is.null(reducedDimName)) {
                metadata(rowData(sfe))$params[[name]]
            } else if (is.null(reducedDimName)) {
                attr(annotGeometry(sfe, annotGeometryName, sample_id = "all"), "params")[[name]]
            } else {
                attr(reducedDim(sfe, reducedDimName), "params")[[name]]
            }
        } else {
            attr(colGeometry(sfe, colGeometryName, sample_id = "all"), "params")[[name]]
        }
    }
}
