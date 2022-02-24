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
#' \describe{
#' \item{colGeometry/ies}{dimGeometry/ies with MARGIN = 2, for geometries
#' associated with columns of the gene count matrix (cells/Visium spots/samples).}
#' \item{rowGeometry/ies}{dimGeometry/ies with MARGIN = 1, for geometries
#' associated with rows of the gene count matrix (genes/features).}
#' \item{spotPoly}{Polygons of spots from technologies such as Visium, ST, and
#' slide-seq, which do not correspond to cells. Centroids of the polygons are stored
#' in \code{spatialCoords} of the underlying \code{SpatialExperiment} object.}
#' \item{ROIPoly}{Polygons of regions of interest (ROIs) from technologies such
#' as laser capture microdissection (LCM) and GeoMX DSP. These should correspond
#' to columns of the gene count matrix.}
#' \item{cellSeg}{Cell segmentation polygons. If the columns of the gene count
#' matrix are single cells, then this is stored in \code{colGeometries}.
#' Otherwise, this is stored in \code{\link{annotGeometries}}.}
#' \item{nucSeg}{Similar to \code{cellSeg}, but for nuclei rather than whole
#' cell.}
#' \item{txSpots}{POINT or MULTIPOINT geometries of transcript spots of single
#' molecular resolution technologies, stored in \code{rowGeometries}.}
#' }
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns.
#' @param withDimnames Logical. If \code{TRUE}, then the dimnames (colnames or
#'   rownames) of the gene count matrix should correspond to row names of the
#'   \code{sf} data frames of interest.
#' @param value Value to set. For \code{dimGeometry}, must be a \code{sf} data
#'   frame with the same number of rows as size in the dimension of interest, or
#'   an ordinary data frame that can be converted to such a \code{sf} data frame
#'   (see \code{\link{df2sf}}). For \code{dimGeometries}, must be a list of such
#'   \code{sf} or ordinary data frames.
#' @param ... \code{spatialCoordsNames, spotDiameter, geometryType} passed to
#'   \code{\link{df2sf}}. For \code{dimGeometries<-} only: \code{geometryType}
#'   can be a character vector of the geometry type of each data frame in the
#'   list of the same length as the list if the data frames specify different
#'   types of geometries.
#' @name dimGeometries
#' @aliases dimGeometry dimGeometries dimGeometryNames colGeometry rowGeometry
#'   colGeometries rowGeometries colGeometryNames rowGeometryNames colGeometry<-
#'   rowGeometry<- colGeometries<- rowGeometries<- colGeometryNames<-
#'   rowGeometryNames<- dimGeometry,SpatialFeatureExperiment,missing-method
#'   dimGeometry,SpatialFeatureExperiment,numeric-method
#'   dimGeometry,SpatialFeatureExperiment,character-method
#'   dimGeometries,SpatialFeatureExperiment-method dimGeometry<- dimGeometries<-
#'   dimGeometryNames<- dimGeometry<-,SpatialFeatureExperiment,missing-method
#'   dimGeometry<-,SpatialFeatureExperiment,numeric-method
#'   dimGeometry<-,SpatialFeatureExperiment,character-method
#'   dimGeometries<-,SpatialFeatureExperiment-method
#'   dimGeometryNames<-,SpatialFeatureExperiment,numeric,character-method
NULL

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
                 function(x, MARGIN, withDimnames = TRUE, ..., value) {
                   value <- .df2sf_list(x, ...)
                   if (withDimnames) {
                     for (v in seq_along(value)) {
                       value[[v]] <- .check_dimgeo_names(x, value[[v]], MARGIN,
                                                         withDimnames=TRUE,
                                                         vname=sprintf("value[[%s]]", v),
                                                         fun='dimGeometries')
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
                 c("SpatialFeatureExperiment", "numeric", "character"),
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
setReplaceMethod("dimGeometry", c("SpatialFeatureExperiment", "missing"),
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
setReplaceMethod("dimGeometry", c("SpatialFeatureExperiment", "numeric"),
                 function(x, type, MARGIN, withDimnames=TRUE, ..., value) {
                   value <- .df2sf_in_list(value, ...)
                   value <- .check_dimgeo_names(x, value, withDimnames)
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
setReplaceMethod("dimGeometry", c("SpatialFeatureExperiment", "character"),
                 function(x, type, MARGIN, withDimnames=TRUE, ..., value) {
                   value <- .df2sf_in_list(value, ...)
                   value <- .check_dimgeo_names(x, value, withDimnames)
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

#' @rdname dimGeometries
#' @export
colGeometry <- function(x, type, withDimnames = TRUE) {
  dimGeometry(x, type, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`colGeometry<-` <- function(x, type, withDimnames = TRUE, value) {
  dimGeometry(x, type, MARGIN = 2, withDimnames = withDimnames) <- value
  x
}

#' @rdname dimGeometries
#' @export
colGeometries <- function(x, withDimnames = TRUE) {
  dimGeometry(x, MARGIN = 2, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`colGeometries<-` <- function(x, withDimnames = TRUE, value) {
  dimGeometries(x, MARGIN = 2, withDimnames = withDimnames) <- value
  x
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
}

#' @rdname dimGeometries
#' @export
rowGeometry <- function(x, type, withDimnames = TRUE) {
  dimGeometry(x, type, MARGIN = 1, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`rowGeometry<-` <- function(x, type, withDimnames = TRUE, value) {
  dimGeometry(x, type, MARGIN = 1, withDimnames = withDimnames) <- value
  x
}

#' @rdname dimGeometries
#' @export
rowGeometries <- function(x, withDimnames = TRUE) {
  dimGeometry(x, MARGIN = 1, withDimnames = withDimnames)
}

#' @rdname dimGeometries
#' @export
`rowGeometries<-` <- function(x, withDimnames = TRUE, value) {
  dimGeometries(x, MARGIN = 1, withDimnames = withDimnames) <- value
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
spotPoly <- function(x, withDimnames = TRUE) {
  colGeometry(x, "spotPoly", withDimnames)
}

#' @rdname dimGeometries
#' @export
`spotPoly<-` <- function(x, withDimnames = TRUE, value) {
  colGeometry(x, "spotPoly", withDimnames) <- value
  x
}

#' @rdname dimGeometries
#' @export
ROIPoly <- function(x, withDimnames = TRUE) {
  colGeometry(x, "ROIPoly", withDimnames)
}

#' @rdname dimGeometries
#' @export
`ROIPoly<-` <- function(x, withDimnames = TRUE, value) {
  colGeometry(x, "ROIPoly", withDimnames) <- value
  x
}

.get_col_then_annot <- function(x, name, withDimnames) {
  if (name %in% colGeometryNames(x)) {
    colGeometry(x, name, withDimnames)
  } else {
    annotGeometry(x, name)
  }
}

.set_col_then_annot <- function(x, name, withDimnames, value) {
  if (name %in% colGeometryNames(x)) {
    colGeometry(x, name, withDimnames) <- value
  } else {
    annotGeometry(x, name) <- value
  }
  return(x)
}

#' @rdname dimGeometries
#' @export
cellSeg <- function(x, withDimnames = TRUE) {
  .get_col_then_annot(x, "cellSeg", withDimnames)
}

#' @rdname dimGeometries
#' @export
`cellSeg<-` <- function(x, withDimnames = TRUE, value) {
  .set_col_then_annot(x, "cellSeg", withDimnames, value)
}

#' @rdname dimGeometries
#' @export
nucSeg <- function(x, withDimnames = TRUE) {
  .get_col_then_annot(x, "nucSeg", withDimnames)
}

#' @rdname dimGeometries
#' @export
`nucSeg<-` <- function(x, withDimnames = TRUE, value) {
  .set_col_then_annot(x, "nucSeg", withDimnames, value)
}

#' @rdname dimGeometries
#' @export
txSpots <- function(x, withDimnames = TRUE) {
  rowGeometry(x, "txSpots", withDimnames)
}

#' @rdname dimGeometries
#' @export
`txSpots<-` <- function(x, withDimnames = TRUE, value) {
  rowGeometry(x, "txSpots", withDimnames) <- value
  x
}
