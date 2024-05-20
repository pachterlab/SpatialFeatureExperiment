#' Functions re-exported from other packages
#'
#' These are some commonly used getters and setters of classes that SFE inherits
#' so you don't have to separately attach those packages to use these functions.
#'
#' @name reexports
#' @inheritParams SummarizedExperiment::`colData<-`
#' @inheritParams SummarizedExperiment::rowData
#' @param object A \code{SingleCellExperiment} object, which includes SFE.
#' @param type Name or numeric index to indicate which \code{reducedDim} to get,
#' such as "PCA". By default the first item in \code{reducedDims}.
NULL

#' @rdname reexports
#' @export
SummarizedExperiment::colData
#' @rdname reexports
#' @export
SummarizedExperiment::rowData
#' @rdname reexports
#' @export
SummarizedExperiment::`colData<-`
#' @rdname reexports
#' @export
SpatialExperiment::spatialCoords
#' @rdname reexports
#' @export
SpatialExperiment::`spatialCoords<-`
#' @rdname reexports
#' @export
SpatialExperiment::spatialCoordsNames
#' @rdname reexports
#' @export
SpatialExperiment::getImg
#' @rdname reexports
#' @export
SpatialExperiment::imgData
#' @rdname reexports
#' @export
SpatialExperiment::rmvImg
#' @rdname reexports
#' @export
SingleCellExperiment::counts
#' @importFrom SingleCellExperiment logcounts reducedDim
#' @rdname reexports
#' @export
SingleCellExperiment::logcounts
#' @rdname reexports
#' @export
SingleCellExperiment::reducedDim
