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
                                       key=.dg_key(MARGIN))

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
                                           vname=sprintf("value[[%s]]", v), fun='dimGeometries')
                     }
                   }

                   .set_internal_all(x, value,
                                     getfun=.getfun(MARGIN),
                                     setfun=.setfun(MARGIN),
                                     key=.dg_key(MARGIN),
                                     convertfun=NULL,
                                     xdimfun=.xdimfun(MARGIN),
                                     vdimfun=nrow,
                                     funstr="dimGeometries",
                                     xdimstr=.xdimstr(MARGIN),
                                     vdimstr="rows")
                 })

#' @rdname dimGeometries
#' @export
setMethod("dimGeometryNames", "SpatialFeatureExperiment",
          function(x, MARGIN) {
            .get_internal_names(x,
                                getfun=.getfun(MARGIN),
                                key=.dg_key(MARGIN))
          })

#' @rdname dimGeometries
#' @export
setReplaceMethod("dimGeometryNames",
                 signature(x = "SpatialFeatureExperiment", value = "character"),
                 function(x, MARGIN, value) {
                   .set_internal_names(x, value,
                                       getfun=.getfun(MARGIN),
                                       setfun=.setfun(MARGIN),
                                       key=.dg_key(MARGIN))
                 })

#' @rdname dimGeometries
#' @export
setMethod("dimGeometry", c("SpatialFeatureExperiment", "missing"),
          function(x, type, MARGIN, withDimnames = TRUE) {
            .get_internal_missing(x,
                                  basefun=dimGeometry,
                                  namefun=dimGeometryNames,
                                  funstr="dimGeometry",
                                  withDimnames=withDimnames,
                                  MARGIN = MARGIN)
          })

#' @rdname dimGeometries
#' @export
setMethod("dimGeometry", c("SpatialFeatureExperiment", "numeric"),
          function(x, type, MARGIN, withDimnames = TRUE) {
            out <- .get_internal_integer(x, type,
                                         getfun=.getfun(MARGIN),
                                         key=.dg_key(MARGIN),
                                         funstr="dimGeometry",
                                         substr="type")

            if (withDimnames) {
              rownames(out) <- colnames(x)
            }
            out
          })

#' @rdname dimGeometries
#' @export
setMethod("dimGeometry", c("SpatialFeatureExperiment", "character"),
          function(x, type, MARGIN, withDimnames = TRUE) {
            out <- .get_internal_character(x, type,
                                           getfun=.getfun(MARGIN),
                                           key=.dg_key(MARGIN),
                                           funstr="dimGeometry",
                                           substr="type",
                                           namestr="dimGeometryNames")

            if (withDimnames) {
              rownames(out) <- colnames(x)
            }
            out
          })

#' @rdname dimGeometries
#' @export
setReplaceMethod("dimGeometry",
                 signature(x = "SpatialFeatureExperiment", type = "missing",
                           value = "sf"),
                 function(x, type, MARGIN, withDimnames=TRUE, value) {
                   .set_internal_missing(x, value,
                                         withDimnames=withDimnames,
                                         MARGIN = MARGIN,
                                         basefun=`dimGeometry<-`,
                                         namefun=dimGeometryNames
                   )
                 })

#' @rdname dimGeometries
#' @export
setReplaceMethod("dimGeometry",
                 signature(x = "SpatialFeatureExperiment", type = "numeric",
                           value = "sf"),
                 function(x, type, MARGIN, withDimnames=TRUE, value) {
                   .check_dimgeo_names(x, value, withDimnames)
                   .set_internal_numeric(x, type, value,
                                         getfun=.getfun(MARGIN),
                                         setfun=.setfun(MARGIN),
                                         key=.dg_key(MARGIN),
                                         convertfun=NULL,
                                         xdimfun=.xdimfun(MARGIN),
                                         vdimfun=nrow,
                                         funstr="dimGeometry",
                                         xdimstr=.xdimstr(MARGIN),
                                         vdimstr="rows",
                                         substr="type")
                 })

#' @rdname dimGeometries
#' @export
setReplaceMethod("dimGeometry",
                 signature(x = "SpatialFeatureExperiment", type = "character",
                           value = "sf"),
                 function(x, type, MARGIN, withDimnames=TRUE, value) {
                   .check_dimgeo_names(x, value, withDimnames)
                   .set_internal_character(x, type, value,
                                           getfun=.getfun(MARGIN),
                                           setfun=.setfun(MARGIN),
                                           key=.dg_key(MARGIN),
                                           convertfun=NULL,
                                           xdimfun=.xdimfun(MARGIN),
                                           vdimfun=nrow,
                                           funstr="dimGeometry",
                                           xdimstr=.xdimstr(MARGIN),
                                           vdimstr="rows",
                                           substr="type")
                 })
# To do: replacement methods for ordinary data frames of a certain format.
# Internally convert to sf

#' @rdname dimGeometries
#' @export
colGeometry <- function(x, type, withDimnames = TRUE) {
  dimGeometry(x, type, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`colGeometry<-` <- function(x, type, withDimnames = TRUE, value) {
  `dimGeometry<-`(x, type, MARGIN = 2, withDimnames = withDimnames, value = value)
}

#' @rdname dimGeometries
#' @export
colGeometries <- function(x, withDimnames = TRUE) {
  dimGeometry(x, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`colGeometries<-` <- function(x, withDimnames = TRUE, value) {
  `dimGeometries<-`(x, MARGIN = 2, withDimnames = withDimnames, value = value)
}

#' @rdname dimGeometries
#' @export
colGeometryNames <- function(x) {
  dimGeometryNames(x, MARGIN = 2)
}

#' @rdname dimGeometries
#' @export
`colGeometryNames<-` <- function(x, value) {
  `dimGeometryNames<-`(x, MARGIN = 2, value = value)
}

#' @rdname dimGeometries
#' @export
rowGeometry <- function(x, type, withDimnames = TRUE) {
  dimGeometry(x, type, MARGIN = 1, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`rowGeometry<-` <- function(x, type, withDimnames = TRUE, value) {
  `dimGeometry<-`(x, type, MARGIN = 1, withDimnames = withDimnames, value = value)
}

#' @rdname dimGeometries
#' @export
rowGeometries <- function(x, withDimnames = TRUE) {
  dimGeometry(x, MARGIN = 1, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`rowGeometries<-` <- function(x, withDimnames = TRUE, value) {
  `dimGeometries<-`(x, MARGIN = 1, withDimnames = withDimnames, value = value)
}

#' @rdname dimGeometries
#' @export
rowGeometryNames <- function(x) {
  dimGeometryNames(x, MARGIN = 1)
}

#' @rdname dimGeometries
#' @export
`rowGeometryNames<-` <- function(x, value) {
  `dimGeometryNames<-`(x, MARGIN = 1, value = value)
}
