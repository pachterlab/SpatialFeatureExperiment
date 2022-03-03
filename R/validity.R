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
.check_annotgeometries <- function(x) {
  gs <- int_metadata(x)$annotGeometries
  ids_coldata <- as.character(unique(colData(x)$sample_id))
  msgs <- lapply(seq_along(gs), function(i) {
    if (!"sample_id" %in% names(gs[[i]])) {
      paste0("Item ", i, " of annotGeometries does not have column sample_id.\n")
    } else {
      samples_seen <- as.character(unique(gs[[i]]$sample_id))
      if (any(!samples_seen %in% ids_coldata)) {
        paste0("Samples ", paste(setdiff(samples_seen, ids_coldata), collapse = ", "),
               " in item ", i, " of annotGeometries are absent from colData.\n")
      } else {
        character(0)
      }
    }
  })
  return(unlist(msgs))
}
.check_graphs_sample <- function(gs, right_length, mar, sample_id) {
  # 1. Must be all listw objects
  # 2. colGraphs: Must all have the same length as there are columns in the gene count matrix
  # rowGraphs: Same length as number of rows
  # annotGraphs: Length is not checked but will thorw error if length doesn't match
  # that of the annotGeometry of interest.
  # 3. spatialGraphs is a list with elements col, row, and annot.
  # 4. Each of col, row, and annot are lists whose elements are sample_id
  # 5. Each sample_id is a list of the listw graphs.
  # 6. I think I should add sample_id to annotGeometries as well.
  out <- character(0)
  if (is.null(gs)) return(out)
  is_listw <- vapply(gs, is, class2 = "listw", FUN.VALUE = logical(1))
  if (!all(is_listw)) {
    out <- paste("All elements of spatialGraphs must have class listw.",
                 "Elements", paste(which(!is_listw), collapse = ", "),
                 "are not of class listw.\n")
  } else if (mar != "annot") {
    is_good_length <- vapply(gs, function(g) {
      isTRUE(all.equal(right_length, length(gs$neighbours)))
    }, FUN.VALUE = logical(1))
    if (!all(is_good_length)) {
      out <- paste0("The `neighbours` field of all elements of ", mar,
                    "Graphs must have the same length as n", mar,
                    "s of the gene count matrix.",
                   "Elements", paste(which(!is_good_length), collapse = ", "),
                   " in sample ", sample_id, " do not have the right length.\n")
    }
  }
  return(out)
}

.check_graphs <- function(gss, right_length, mar) {
  # gss is a list whose elements are sample_id's
  lapply(seq_along(gss), function(i) {
    .check_graphs_sample(gss[[i]], right_length, mar, names(gss)[[i]])
  })
}

.check_graph_sample_id <- function(x) {
  ids_coldata <- as.character(unique(colData(x)$sample_id))
  graph_sample_ids <- lapply(1:3, function(i)
    names(int_metadata(x)$spatialGraphs[[.margin_name(i)]]))
  graph_sample_ids <- unique(unlist(graph_sample_ids))
  if (any(!graph_sample_ids %in% ids_coldata)) {
    paste0("Samples ", paste(setdiff(graph_sample_ids, ids_coldata), collapse = ", "),
           " are in the graphs but not colData.\n")
  } else {
    return(character(0))
  }
}

setValidity("SpatialFeatureExperiment", function(object) {
  col_msg <- .check_geometries(int_colData(object)$colGeometries, "colGeometries")
  row_msg <- .check_geometries(int_elementMetadata(object)$rowGeometries, "rowGeometries")
  annot_msg <- .check_geometries(int_metadata(object)$annotGeometries, "annotGeometries")
  annot_msg2 <- .check_annotgeometries(object)
  colgraph_msg <- .check_graphs(int_metadata(object)$spatialGraphs$col,
                                ncol(object), "col")
  rowgraph_msg <- .check_graphs(int_metadata(object)$spatialGraphs$row,
                                nrow(object), "row")
  annotgraph_msg <- .check_graphs(int_metadata(object)$spatialGraphs$annot,
                                  NA, "annot")
  id_msg <- .check_graph_sample_id(object)
  outs <- c(col_msg, row_msg, annot_msg, annot_msg2, colgraph_msg, rowgraph_msg,
            annotgraph_msg, id_msg)
  outs <- unlist(outs)
  if (length(outs)) {
    return(outs)
  } else {
    return(TRUE)
  }
})
