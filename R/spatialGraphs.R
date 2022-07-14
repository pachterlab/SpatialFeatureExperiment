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
#' @param name Name of the graphs to add to each sample_id; used in the
#'   \code{spatialGraphs} replacement method as it must be character while
#'   \code{type} can be either an integer index or a name.
#' @param MARGIN As in \code{\link{apply}}. 1 stands for rows and 2 stands for
#'   columns. In addition, 3 stands for spatial neighborhood graphs that
#'   correspond to \code{annotGeometries}.
#' @param sample_id Name of the sample the graph is associated with. This is
#'   useful when multiple pieces of tissues are in the same SFE object (say for
#'   a joint dimension reduction and clustering) and the spatial neighborhood is
#'   only meaningful within the same piece of tissue. See the \code{sample_id}
#'   argument in \code{\link{SpatialExperiment}}.
#' @name spatialGraphs
#' @return Getters for multiple graphs return a named list. Getters for
#'   names return a character vector of the names. Getters for single graphs
#'   return a \code{listw} object. Setters return an SFE object.
#' @aliases rowGraphs rowGraphs<- spatialGraph spatialGraph<- spatialGraphNames
#'   colGraphs colGraphs<- spatialGraphNames<- spatialGraphs<- annotGraphs
#'   annotGraphs<-
#' @importFrom methods as validObject show
#' @importFrom stats setNames
#' @docType methods
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' g1 <- findVisiumGraph(sfe)
#' g2 <- findSpatialNeighbors(sfe)
#'
#' # Set all graphs of a margin by a named list
#' spatialGraphs(sfe, MARGIN = 2L, sample_id = "Vis5A") <-
#'   list(tri2nb = g2, visium = g1)
#' # Or equivalently
#' colGraphs(sfe, sample_id = "Vis5A") <- list(tri2nb = g2, visium = g1)
#'
#' # Get all graphs of a margin, returning a named list
#' gs <- spatialGraphs(sfe, MARGIN = 2L)
#' # Or equivalently
#' gs <- colGraphs(sfe)
#'
#' # Set graph of the same name and same margin for multiple samples
#' # Each sample has a separate graph
#' sfe2 <- McKellarMuscleData("small2")
#' sfe_combined <- cbind(sfe, sfe2)
#' colGraphs(sfe_combined, name = "visium", sample_id = "all") <-
#'   findVisiumGraph(sfe_combined, sample_id = "all")
#'
#' # Get graph names
#' spatialGraphNames(sfe, MARGIN = 2L, sample_id = "Vis5A")
#' # Or equivalently (sample_id optional as only one sample is present)
#' colGraphNames(sfe)
#'
#' # Set graph names
#' spatialGraphNames(sfe, MARGIN = 2L) <- c("foo", "bar")
#' colGraphNames(sfe) <- c("tri2nb", "visium")
#'
#' # MARGIN = 1 means rowGraphs; MARGIN = 3 means annotation graphs (annotGraphs)
#' # for both getters and setters
#'
#' # Set single graph by
#' # Spatial graph for myofibers
#' g_myofiber <- findSpatialNeighbors(sfe, type = "myofiber_simplified",
#'                                    MARGIN = 3L)
#' spatialGraph(sfe, type = "myofiber", MARGIN = 3L) <- g_myofiber
#' # Or equivalently
#' annotGraph(sfe, "myofiber") <- g_myofiber
#'
#' # Get a specific graph by name
#' g <- spatialGraph(sfe, "myofiber", MARGIN = 3L)
#' g2 <- spatialGraph(sfe, "visium", MARGIN = 2L)
#' # Or equivalently
#' g <- annotGraph(sfe, "myofiber")
#' g2 <- colGraph(sfe, "visium")
NULL

.margin_name <- function(MARGIN) switch (MARGIN, "row", "col", "annot")
.margin_len_fun <- function(MARGIN) switch(MARGIN, nrow, ncol, function(x) return(NA))
.margin_num <- function(name) switch(name, row = 1, col = 2, annot = 3)

.get_graphs <- function(x, type = c("all", "sample_all", "margin_all", "one"),
                        which = NA) {
  gss <- int_metadata(x)$spatialGraphs
  if (is.null(gss)) return(NULL)
  gss <- gss[c("row", "col", "annot"), , drop = FALSE]
  switch (type,
    all = as.list(gss),
    sample_all = as.list(gss[,which]),
    # Get rid of the layer of list indicating margin
    margin_all = lapply(as.list(gss[which,]), function(o) o[[1]]),
    one = {
      out <- gss[which[[1]], which[[2]]]
      if (length(which[[2]]) > 1L) {
        # Multiple samples
        lapply(as.list(out), function(o) o[[1]])
      } else {
        out[[1]]
      }
    }
  )
}

# For lists of listw objects--------
#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "missing", "missing"),
          function(x, MARGIN, sample_id, name) .get_graphs(x, "all"))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "missing", "missing"),
          function(x, MARGIN, sample_id, name) .get_graphs(x, "margin_all", MARGIN))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "character", "missing"),
          function(x, MARGIN, sample_id, name) {
            sample_id <- .check_sample_id(x, sample_id, one = FALSE)
            .get_graphs(x, "sample_all", sample_id)
          })

#' @rdname spatialGraphs
#' @export
setMethod("colGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
          function(x, sample_id, name) spatialGraphs(x, MARGIN = 2))

#' @rdname spatialGraphs
#' @export
setMethod("rowGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
          function(x, sample_id, name) spatialGraphs(x, MARGIN = 1))

#' @rdname spatialGraphs
#' @export
setMethod("annotGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
          function(x, sample_id, name) spatialGraphs(x, MARGIN = 3))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character", "missing"),
          function(x, MARGIN, sample_id, name) {
            sample_id <- .check_sample_id(x, sample_id, one = FALSE)
            .get_graphs(x, "one", list(MARGIN, sample_id))
          })

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character", "character"),
          function(x, MARGIN, sample_id, name) {
            sample_id <- .check_sample_id(x, sample_id, one = FALSE)
            out <- .get_graphs(x, "one", list(MARGIN, sample_id))
            # Get rid of the layer of list indicating name as it's known
            if (length(sample_id) > 1L) {
              out <- lapply(out, function(o) o[[name]])
              out <- out[vapply(out, function(o) !is.null(o), FUN.VALUE = logical(1))]
            } else {
              out <- setNames(list(out[[name]]), sample_id)
            }
            out
          })

#' @rdname spatialGraphs
#' @export
setMethod("colGraphs", c("SpatialFeatureExperiment", "character", "missing"),
          function(x, sample_id, name) spatialGraphs(x, 2, sample_id))

#' @rdname spatialGraphs
#' @export
setMethod("colGraphs", c("SpatialFeatureExperiment", "character", "character"),
          function(x, sample_id, name) spatialGraphs(x, 2, sample_id, name))

#' @rdname spatialGraphs
#' @export
setMethod("rowGraphs", c("SpatialFeatureExperiment", "character", "missing"),
          function(x, sample_id, name) spatialGraphs(x, 1, sample_id))

#' @rdname spatialGraphs
#' @export
setMethod("rowGraphs", c("SpatialFeatureExperiment", "character", "character"),
          function(x, sample_id, name) spatialGraphs(x, 1, sample_id, name))

#' @rdname spatialGraphs
#' @export
setMethod("annotGraphs", c("SpatialFeatureExperiment", "character", "missing"),
          function(x, sample_id, name) spatialGraphs(x, 3, sample_id))

#' @rdname spatialGraphs
#' @export
setMethod("annotGraphs", c("SpatialFeatureExperiment", "character", "character"),
          function(x, sample_id, name) spatialGraphs(x, 3, sample_id, name))

.fill_missing <- function(l, names_use) {
  to_add <- setdiff(names_use, names(l))
  l2 <- vector("list", length = length(to_add))
  names(l2) <- to_add
  l <- c(l, l2)
  l <- l[names_use]
  l
}

.initialize_spatialGraphs <- function(x) {
  empty_col <- list(row = NULL, col = NULL, annot = NULL)
  samples <- sampleIDs(x)
  args_init <- list(sample1 = I(empty_col), row.names = c("row", "col", "annot"))
  names(args_init)[1] <- samples[1]
  df <- do.call(DataFrame, c(args_init, list(check.names = FALSE)))
  if (length(samples) > 1) {
    for (i in 2:length(samples)) {
      df[[samples[i]]] <- empty_col
    }
  }
  df
}
.set_graphs <- function(x, type = c("all", "sample_all", "margin_all", "one"),
                        which = NA, value) {
  if (is.null(int_metadata(x)$spatialGraphs)) {
    df <- .initialize_spatialGraphs(x)
  } else {
    df <- int_metadata(x)$spatialGraphs[c("row", "col", "annot"),,drop = FALSE]
  }
  if (type == "all") {
    if (is.null(value))
      df <- .initialize_spatialGraphs(x)
    else if (!is(value, "DataFrame")) {
      value <- lapply(value, .fill_missing, names_use = c("row", "col", "annot"))
      df <- DataFrame(lapply(value, I), row.names = c("row", "col", "annot"))
    } else df <- value[c("row", "col", "annot"), , drop = FALSE]
  } else if (type == "sample_all") {
    value <- .fill_missing(value, names_use = c("row", "col", "annot"))
    df[[which]] <- value[c("row", "col", "annot")]
  } else if (type == "margin_all") {
    df <- df[-which,,drop = FALSE] # which should be a number. rows have been reordered
    value <- .fill_missing(value, sampleIDs(x))
    value <- lapply(value, function(v) I(list(v)))
    new_row <- as(value, "DataFrame")
    rownames(new_row) <- .margin_name(which)
    new_row <- new_row[,names(df), drop = FALSE]
    df <- rbind(df, new_row)
    df <- df[c("row", "col", "annot"),,drop = FALSE]
  } else {
    df[which[[1]], which[[2]]] <- I(list(value))
  }
  int_metadata(x)$spatialGraphs <- df
  if (!isTRUE(S4Vectors:::disableValidity())) {
    m <- .check_graphs(x)
    if (length(m)) stop(m)
  }
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "missing", "missing"),
                 function(x, MARGIN, sample_id, name, value)
                   .set_graphs(x, "all", value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "missing", "missing"),
                 function(x, MARGIN, sample_id, name, value)
                   .set_graphs(x, "margin_all", which = MARGIN, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "NULL", "missing"),
                 function(x, MARGIN, sample_id, name, value) {
                   sample_id <- .check_sample_id(x, sample_id)
                   .set_graphs(x, "sample_all", which = sample_id, value = value)
                 })

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "missing", "character", "missing"),
                 function(x, MARGIN, sample_id, name, value)
                   .set_graphs(x, "sample_all", which = sample_id, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("colGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
                 function(x, sample_id, name, value)
                   .set_graphs(x, "margin_all", which = 2, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("rowGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
                 function(x, sample_id, name, value)
                   .set_graphs(x, "margin_all", which = 1, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("annotGraphs", c("SpatialFeatureExperiment", "missing", "missing"),
                 function(x, sample_id, name, value)
                   .set_graphs(x, "margin_all", which = 3, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character", "missing"),
                 function(x, MARGIN, sample_id, name, value)
                   .set_graphs(x, "one", which = list(MARGIN, sample_id),
                               value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("colGraphs", c("SpatialFeatureExperiment", "character", "missing"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 2, sample_id, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("rowGraphs", c("SpatialFeatureExperiment", "character", "missing"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 1, sample_id, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("annotGraphs", c("SpatialFeatureExperiment", "character", "missing"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 3, sample_id, value = value))

# I should write another method to add graphs of the same name to each sample_id
# when multiple sample_ids are specified.
#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphs", c("SpatialFeatureExperiment", "numeric", "character",
                                    "character"),
                 function(x, MARGIN, sample_id, name, value) {
                   sample_id <- .check_sample_id(x, sample_id, one = FALSE)
                   for (s in sample_id) {
                     spatialGraph(x, type = name, MARGIN = MARGIN, sample_id = s) <- value[[s]]
                   }
                   x
                 })

#' @rdname spatialGraphs
#' @export
setReplaceMethod("colGraphs", c("SpatialFeatureExperiment", "character", "character"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 2, sample_id,
                                                                       name = name, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("rowGraphs", c("SpatialFeatureExperiment", "character", "character"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 1, sample_id,
                                                                       name = name, value = value))

#' @rdname spatialGraphs
#' @export
setReplaceMethod("annotGraphs", c("SpatialFeatureExperiment", "character", "character"),
                 function(x, sample_id, name, value) `spatialGraphs<-`(x, 3, sample_id,
                                                                       name = name, value = value))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraphNames", c("SpatialFeatureExperiment", "numeric", "ANY"),
          function(x, MARGIN, sample_id) {
            sample_id <- .check_sample_id(x, sample_id)
            names(spatialGraphs(x, MARGIN, sample_id))
          })

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraphNames", c("SpatialFeatureExperiment",
                                        "numeric", "ANY", "character"),
                 function(x, MARGIN, sample_id = NULL, value) {
                   sample_id <- .check_sample_id(x, sample_id)
                   names(spatialGraphs(x, MARGIN, sample_id)) <- value
                   x
                 })

#' @rdname spatialGraphs
#' @export
colGraphNames <- function(x, sample_id = NULL) spatialGraphNames(x, 2, sample_id)

#' @rdname spatialGraphs
#' @export
rowGraphNames <- function(x, sample_id = NULL) spatialGraphNames(x, 1, sample_id)

#' @rdname spatialGraphs
#' @export
annotGraphNames <- function(x, sample_id = NULL) spatialGraphNames(x, 3, sample_id)

#' @rdname spatialGraphs
#' @export
`colGraphNames<-` <- function(x, sample_id = NULL, value) `spatialGraphNames<-`(x, 2, sample_id, value)

#' @rdname spatialGraphs
#' @export
`rowGraphNames<-` <- function(x, sample_id = NULL, value) `spatialGraphNames<-`(x, 1, sample_id, value)

#' @rdname spatialGraphs
#' @export
`annotGraphNames<-` <- function(x, sample_id = NULL, value) `spatialGraphNames<-`(x, 3, sample_id, value)

# For single listws, not in a list-------
#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "missing", "numeric"),
          function(x, type, MARGIN, sample_id = NULL) spatialGraph(x, 1L, MARGIN, sample_id))

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "numeric", "numeric"),
          function(x, type, MARGIN, sample_id = NULL) {
            sample_id <- .check_sample_id(x, sample_id)
            spatialGraphs(x, MARGIN, sample_id)[[type]]
          })

#' @rdname spatialGraphs
#' @export
setMethod("spatialGraph", c("SpatialFeatureExperiment", "character", "numeric"),
          function(x, type, MARGIN, sample_id = NULL) {
            sample_id <- .check_sample_id(x, sample_id)
            spatialGraphs(x, MARGIN, sample_id)[[type]]
          })

#' @rdname spatialGraphs
#' @export
colGraph <- function(x, type = 1L, sample_id = NULL) spatialGraph(x, type, 2, sample_id)

#' @rdname spatialGraphs
#' @export
rowGraph <- function(x, type = 1L, sample_id = NULL) spatialGraph(x, type, 1, sample_id)

#' @rdname spatialGraphs
#' @export
annotGraph <- function(x, type = 1L, sample_id = NULL) spatialGraph(x, type, 3, sample_id)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment", "missing",
                                   "numeric"),
                 function(x, type, MARGIN, sample_id = NULL, value)
                   `spatialGraph<-`(x, 1L, MARGIN, sample_id, value))

.sg_r <- function(x, type, MARGIN, sample_id = NULL, value) {
  sample_id <- .check_sample_id(x, sample_id)
  if (!is.null(value)) {
    if (!is(value, "listw")) {
      stop("value must be of class listw.")
    } else if (MARGIN == 1L && length(value$neighbours) != nrow(x)) {
      stop("The neighbours field of `value` must be the same as nrows of the gene count matrix.")
    } else if (MARGIN == 2L) {
      right_length <- sum(colData(x)$sample_id %in% sample_id)
      if (length(value$neighbours) != right_length) {
        stop("The neighbours field of `value` must be the same as number of samples in this sample_id.")
      }
    }
  }
  if (is.null(int_metadata(x)$spatialGraphs)) {
    df <- .initialize_spatialGraphs(x)
    int_metadata(x)$spatialGraphs <- df
  }
  existing <- int_metadata(x)$spatialGraphs[.margin_name(MARGIN), sample_id][[1]]
  if (type %in% names(existing) || is.null(value)) {
    existing <- existing[names(existing) != type]
  }
  if (!is.null(value)) {
    value <- setNames(list(value), type)
    replacement <- c(existing, value)
  } else {
    replacement <- existing
  }
  int_metadata(x)$spatialGraphs[.margin_name(MARGIN), sample_id] <- I(list(replacement))
  return(x)
}

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment",
                                   "numeric", "numeric"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
setReplaceMethod("spatialGraph", c("SpatialFeatureExperiment",
                                   "character", "numeric"),
                 .sg_r)

#' @rdname spatialGraphs
#' @export
`colGraph<-` <- function(x, type = 1L, sample_id = NULL, value)
  `spatialGraph<-`(x, type, 2, sample_id, value)

#' @rdname spatialGraphs
#' @export
`rowGraph<-` <- function(x, type = 1L, sample_id = NULL, value)
  `spatialGraph<-`(x, type, 1, sample_id, value)

#' @rdname spatialGraphs
#' @export
`annotGraph<-` <- function(x, type = 1L, sample_id = NULL, value)
  `spatialGraph<-`(x, type, 3, sample_id, value)
