.dg_key <- "dimGeometries"
#' Dimension geometry methods
#'
#' "Dimension geometry" refers to Simple Feature (\code{sf}) geometries associated
#' with rows (features, genes) or columns (cells or spots) of the gene count
#' matrix in the \code{SpatialFeatureExperiment} object. For each dimension, the
#' number of rows in the \code{sf} data frame specifying the geometries must
#' match the size of the dimension of interest. For example, there must be the
#' same number of rows in the \code{sf} data frame describing cells as there are
#' cells in the gene count matrix. This page documents getters and setters for
#' the dimension geometries. The getters and setters are implemented in a way
#' similar to those of \code{reducedDims} in \code{SingleCellExperiment}.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#' columns.
#' @param withDimnames Logical. If \code{TRUE}, then the dimnames (colnames or
#' rownames) of the gene count matrix should correspond to row names of the
#' \code{sf} data frames of interest.
#' @param value Value to set. For \code{dimGeometry}, must be a \code{sf} data
#' frame with the same number of rows as size in the dimension of interest.
#' For \code{dimGeometries}, must be a list of such \code{sf} data frames.
#' @rdname dimGeometries
#' @export
setMethod("dimGeometries", "SpatialFeatureExperiment",
          function(x, MARGIN = 2, withDimnames = TRUE) {
            value <- .get_internal_all(x,
                                       getfun=.getfun(MARGIN),
                                       key=.dg_key)

            if (withDimnames) {
              for (i in seq_along(value)) {
                rownames(value[[i]]) <- dimnames(x)[[MARGIN]]
              }
            }
            value
          })

#' @rdname dimGeometries
#' @export
setReplaceMethod("dimGeometries", "SpatialFeatureExperiment",
                 function(x, MARGIN, withDimnames = TRUE, value) {
                   if (withDimnames) {
                     for (v in seq_along(value)) {
                       .check_dimgeo_names(x, value[[v]], MARGIN, withDimnames=TRUE,
                                           vname=sprintf("value[[%s]]", v), fun='reducedDims')
                     }
                   }

                   .set_internal_all(x, value,
                                     getfun=.getfun(MARGIN),
                                     setfun=.setfun(MARGIN),
                                     key=.dg_key,
                                     convertfun=NULL,
                                     xdimfun=ncol,
                                     vdimfun=nrow,
                                     funstr="dimGeometries",
                                     xdimstr=.xdimstr(MARGIN),
                                     vdimstr="rows")
                 })
