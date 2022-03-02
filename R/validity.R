# Validity of SFE object
# 1. All the geometries must be valid geometries
# 2. All elements of *Geometries must be sf objects
.valid_geometry <- function(g, name_show) {
  # Make message if something in an element of *Geometries is not valid
  v <- st_is_valid(g, reason = TRUE)
  which_invalid <- which(v != "Valid Geometry")
  if (length(which_invalid)) {
    paste0("The following items (index shown) in ", name_show, " are invalid:\n",
           paste0(which_invalid, ": ", v[which_invalid], collapse = "\n"))
  } else {
    return(character(0))
  }
}
.check_geometries <- function(gs, name_show) {
  msgs <- lapply(seq_along(gs), function(i) {
    if (is(gs[[i]], "sf")) {
      m <- .valid_geometry(gs[[i]], names(gs)[i])
    } else {
      paste0("Item ", i, " in ", name_show, " is ", class(gs[[i]])[1], " rather than sf.\n")
    }
  })
  msgs <- unlist(msgs)
  if (length(msgs)) {
    paste0("Issues in ", name_show, ":\n", paste(msgs, collapse = "\n"), "\n")
  } else {
    return(character(0))
  }
}

.check_graphs <- function(gs, ncol) {
  # 1. Must be all listw objects
  # 2. Must all have the same length as there are columns in the gene count matrix
  out <- character(0)
  if (is.null(gs)) return(out)
  is_listw <- vapply(gs, is, class2 = "listw", FUN.VALUE = logical(1))
  if (!all(is_listw)) {
    out <- paste("All elements of spatialGraphs must have class listw.",
                 "Elements", paste(which(!is_listw), collapse = ", "),
                 "are not of class listw.\n")
  } else {
    is_good_length <- vapply(gs, function(g) {
      isTRUE(all.equal(ncol, length(gs$neighbours)))
    }, FUN.VALUE = logical(1))
    if (!all(is_good_length)) {
      out <- paste("The `neighbours` field of all elements of spatialGraphs must have the same length",
                   "as there are columns in the gene count matrix.",
                   "Elements", paste(which(!is_good_length), collapse = ", "),
                   "do not have the right length.")
    }
  }
  return(out)
}
setValidity("SpatialFeatureExperiment", function(object) {
  col_msg <- .check_geometries(int_colData(object)$colGeometries, "colGeometries")
  row_msg <- .check_geometries(int_elementMetadata(object)$rowGeometries, "rowGeometries")
  annot_msg <- .check_geometries(int_metadata(object)$annotGeometries, "annotGeometries")
  graph_msg <- .check_graphs(int_metadata(object)$spatialGraphs, ncol(object))
  outs <- c(col_msg, row_msg, annot_msg, graph_msg)
  if (length(outs)) {
    return(outs)
  } else {
    return(TRUE)
  }
})
