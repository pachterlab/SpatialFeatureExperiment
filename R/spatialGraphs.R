# Spatial graphs Another item in int_metadata. Since I used spdep a lot, I'll use
# spdep's nb and listw. Burning question: Shall I store the nb or listw? The
# listw actually contains nb, and normally the W style is used (I haven't used
# any other style). So I think I'll store listw. Also, what to do with singletons?
# And the graph will be reconstructed when the SFE object is subsetted. Then I need
# to store sufficient info about how the graph was constructed. I think, to do so,
# I need to write wrappers of all the spatial graph functions in spdep to easily
# store the construction info, which would be unavailable if the user calls those
# functions separately. I'll issue a warning if such info is unavailable because
# the wrappers were not used.

#' Spatial graph methods
#'
#' Spatial neighborhood graphs as \code{spdep}'s \code{listw} objects are stored
#' in the \code{int_metadata} of the SFE object. The \code{listw} class is used
#' because \code{spdep} has many useful methods that rely on the neighborhood
#' graph as \code{listw}.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param value A \code{listw} object (\code{*Graph}), or a named list of list
#'   of \code{listw} objects (\code{*Graphs}) where the names of the top level
#'   list are \code{sample_id}s when adding graphs for all samples in the margin
#'   of interest, or a list of \code{listw} objects when adding graphs for one
#'   sample in one margin.
#' @param type An integer specifying the index or string specifying the name of
#'   the *Graph to query or replace. If missing, then the first item in the
#'   *Graph will be returned or replaced.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns. In addition, 3 stands for spatial neighborhood graphs that
#'   correspond to \code{annotGeometries}.
#' @param sample_id Name of the sample the graph is associated with. This is
#'   useful when multiple pieces of tissues are in the same SFE object (say for
#'   a joint dimension reduction and clustering) and the spatial neighborhood is
#'   only meaningful within the same piece of tissue. See the \code{sample_id}
#'   argument in \code{\link{SpatialExperiment}}.
#' @name spatialGraphs
#' @docType methods
NULL

.margin_name <- function(MARGIN) switch (MARGIN, "row", "col", "annot")
.margin_len_fun <- function(MARGIN) switch(MARGIN, nrow, ncol, function(x) return(NA))
.margin_num <- function(name) switch(name, row = 1, col = 2, annot = 3)

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "numeric"),
          function(x, sample_id, MARGIN)
            int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]])

#' @rdname spatialGraphs
#' @export
setMethod("colGraphs", c("SpatialFeatureExperiment", "missing"),
          function(x, sample_id) int_metadata(x)$spatialGraphs[[.margin_name(2)]])

#' @rdname spatialGraphs
#' @export
setMethod("rowGraphs", c("SpatialFeatureExperiment", "missing"),
          function(x, sample_id) int_metadata(x)$spatialGraphs[[.margin_name(1)]])

#' @rdname spatialGraphs
#' @export
setMethod("annotGraphs", c("SpatialFeatureExperiment", "missing"),
          function(x, sample_id) int_metadata(x)$spatialGraphs[[.margin_name(3)]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "character", "numeric"),
          function(x, sample_id, MARGIN)
            int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]])

#' @rdname spatialGraphs
#' @export
setMethod("colGraphs", c("SpatialFeatureExperiment", "character"),
          function(x, sample_id) spatialGraphs(x, sample_id, 2))

#' @rdname spatialGraphs
#' @export
setMethod("rowGraphs", c("SpatialFeatureExperiment", "character"),
          function(x, sample_id) spatialGraphs(x, sample_id, 1))

#' @rdname spatialGraphs
#' @export
setMethod("annotGraphs", c("SpatialFeatureExperiment", "character"),
          function(x, sample_id) spatialGraphs(x, sample_id, 3))

.set_all_graphs_margin <- function(x, MARGIN, value) {
  if (!is.null(value)) {
    m <- .check_graphs(value, .margin_len_fun(MARGIN)(x),
                       .margin_name(MARGIN))
    if (length(m)) stop(m)
  }
  int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]] <- value
  x
}

.set_all_graphs <- function(x, value) {
  if (is.null(value)) {
    int_metadata(x)$spatialGraphs <- NULL
    return(x)
  }
  margins <- names(value)
  m <- c("row", "col", "annot")
  mar_use <- intersect(margins, m)
  if (any(!margins %in% m)) {
    message("Elements ", paste(setdiff(margins, m), collapse = ", "),
            " are ignored.")
    value <- value[mar_use]
    if (!length(value)) {
      warning("Names of value do not match any of row, col, or annot. ",
              "Not changing spatialGeometries.")
      return(x)
    }
    for (mar in mar_use) {
      x <- .set_all_graphs_margin(x, .margin_num(mar), value[[mar]])
    }
    return(x)
  }
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
                 function(x, sample_id, MARGIN, value) .set_all_graphs(x, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "numeric"),
                 function(x, sample_id, MARGIN, value) .set_all_graphs_margin(x, MARGIN, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("colGraphs", c("SpatialFeatureExperiment", "missing"),
                 function(x, sample_id, value) .set_all_graphs_margin(x, 2, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("rowGraphs", c("SpatialFeatureExperiment", "missing"),
                 function(x, sample_id, value) .set_all_graphs_margin(x, 1, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("annotGraphs", c("SpatialFeatureExperiment", "missing"),
                 function(x, sample_id, value) .set_all_graphs_margin(x, 3, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "character", "numeric"),
                 function(x, sample_id, MARGIN, value) {
                   if (!is.null(value)) {
                     m <- .check_graphs_sample(value, .margin_len_fun(MARGIN)(x),
                                               .margin_name(MARGIN), sample_id)
                     if (length(m)) stop(m)
                   }
                   int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]] <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
setReplaceMethod("colGraphs", c("SpatialFeatureExperiment", "character"),
                 function(x, sample_id, value) `spatialGraphs<-`(x, sample_id, 2, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("rowGraphs", c("SpatialFeatureExperiment", "character"),
                 function(x, sample_id, value) `spatialGraphs<-`(x, sample_id, 1, value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("annotGraphs", c("SpatialFeatureExperiment", "character"),
                 function(x, sample_id, value) `spatialGraphs<-`(x, sample_id, 3, value))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphNames", c("SpatialFeatureExperiment", "character", "numeric"),
          function(x, sample_id, MARGIN) names(spatialGraphs(x, MARGIN, sample_id)))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphNames", c("SpatialFeatureExperiment", "character",
                                        "numeric", "character"),
                 function(x, sample_id, MARGIN, value) {
                   names(spatialGraphs(x, MARGIN, sample_id)) <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
colGraphNames <- function(x, sample_id) spatialGraphNames(x, sample_id, 2)

#' @rdname spatialGraphs
#' @export
rowGraphNames <- function(x, sample_id) spatialGraphNames(x, sample_id, 1)

#' @rdname spatialGraphs
#' @export
annotGraphNames <- function(x, sample_id) spatialGraphNames(x, sample_id, 3)

#' @rdname spatialGraphs
#' @export
`colGraphNames<-` <- function(x, sample_id, value) `spatialGraphNames<-`(x, sample_id, 2, value)

#' @rdname spatialGraphs
#' @export
`rowGraphNames<-` <- function(x, sample_id, value) `spatialGraphNames<-`(x, sample_id, 1, value)

#' @rdname spatialGraphs
#' @export
`annotGraphNames<-` <- function(x, sample_id, value) `spatialGraphNames<-`(x, sample_id, 3, value)

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "missing"),
          function(x, sample_id, type, MARGIN) spatialGraph(x, sample_id, 1L, MARGIN))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "numeric"),
          function(x, sample_id, type, MARGIN)
            spatialGraphs(x, MARGIN)[[sample_id]][[type]])

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "character"),
          function(x, sample_id, type, MARGIN)
            spatialGraphs(x, MARGIN)[[sample_id]][[type]])

#' @rdname spatialGraphs
#' @export
colGraph <- function(x, sample_id, type = 1L) spatialGraph(x, sample_id, type, 2)

#' @rdname spatialGraphs
#' @export
rowGraph <- function(x, sample_id, type = 1L) spatialGraph(x, sample_id, type, 1)

#' @rdname spatialGraphs
#' @export
annotGraph <- function(x, sample_id, type = 1L) spatialGraph(x, sample_id, type, 3)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "missing",
                                   "numeric", "listw"),
                 function(x, sample_id, type, MARGIN, value)
                   spatialGraph(x, 1L, MARGIN, sample_id) <- value)

.sg_r <- function(x, sample_id, type, MARGIN, value) {
  if (!is.null(value)) {
    if (!is(value, "listw")) {
      stop("value must be of class listw.")
    } else if (MARGIN < 3 && length(value$neighbours) != .margin_len_fun(MARGIN)(x)) {
      stop("The neighbours field of `value` must be the same as n",
           .margin_name(MARGIN), "s of the gene count matrix.")
    }
  }
  int_metadata(x)$spatialGraphs[[.margin_name(MARGIN)]][[sample_id]][[type]] <- value
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character",
                                   "numeric", "numeric"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "character",
                                   "character", "numeric"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
`colGraph<-` <- function(x, sample_id, type = 1L, value)
  `spatialGraph<-`(x, sample_id, type, 2, value)

#' @rdname spatialGraphs
#' @export
`rowGraph<-` <- function(x, sample_id, type = 1L, value)
  `spatialGraph<-`(x, sample_id, type, 1, value)

#' @rdname spatialGraphs
#' @export
`annotGraph<-` <- function(x, sample_id, type = 1L, value)
  `spatialGraph<-`(x, sample_id, type, 3, value)

# To do: 1. store info to reconstruct the graph after subsetting
