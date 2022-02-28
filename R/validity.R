# Validity of SFE object
# 1. All the geometries must be valid geometries
# 2. All elements of *Geometries must be sf objects
.valid_geometry <- function(g, name_show) {
  # Make message if something in an element of *Geometries is not valid
  v <- st_is_valid(g, reason = TRUE)
  which_invalid <- which(v != "ValidGeometry")
  if (length(which_invalid)) {
    paste0("The following items (index shown) in ", name_show, " are invalid:\n",
           paste0(which_invalid, ": ", v[which_invalid], collapse = "\n"))
  } else {
    return("")
  }
}
.check_geometries <- function(gs, name_show) {
  msgs <- vapply(seq_along(gs), function(i) {
    if (is(gs[[i]], "sf")) {
      m <- .valid_geometry(gs[[i]], names(gs)[i])
    } else {
      paste0("Item ", i, " in ", name_show, " is ", class(gs[[i]])[1], " rather than sf.\n")
    }
  }, FUN.VALUE = character(1))
  if (length(msgs)) {
    paste0("Issues in ", name_show, ":\n", paste(msgs, collapse = "\n"), "\n")
  } else {
    return(character(0))
  }
}
setValidity("SpatialFeatureExperiment", function(x) {
  col_msg <- .check_geometries(int_colData(x)$colGeometries, "colGeometries")
  row_msg <- .check_geometries(int_elementMetadata(x)$rowGeometries, "rowGeometries")
  annot_msg <- .check_geometries(int_metadata(x)$annotGeometries, "annotGeometries")
  outs <- c(col_msg, row_msg, annot_msg)
  if (length(outs)) {
    return(outs)
  } else {
    return(TRUE)
  }
})
