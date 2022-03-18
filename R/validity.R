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
# I think I'll reimplement the spatial graphs thing, as a DataFrame, with
# sample_ids as columns and margins as rows, for easier subsetting and easier
# implementation of the getters.
.check_graphs_item <- function(gs, right_length, mar, sample_id) {
  # 1. Must be all listw objects
  # 2. colGraphs: Must all have the same length as there are columns in the gene count matrix
  # rowGraphs: Same length as number of rows
  # annotGraphs: Length is not checked but will throw error if length doesn't match
  # that of the annotGeometry of interest.
  # 3. spatialGraphs is a DataFrame with sample_ids as columns
  # 4. And row, col, and annot as rows
  # 5. Each item in the data frame is a list of listw objects that corresponds
  # to the sample and margin
  out <- character(0)
  if (is.null(gs)) return(out)
  is_listw <- vapply(gs, is, class2 = "listw", FUN.VALUE = logical(1))
  if (!all(is_listw)) {
    out <- paste("All non-NULL elements of spatialGraphs must have class listw.",
                 "Elements", paste(which(!is_listw), collapse = ", "),
                 "are not of class listw.\n")
  } else {
    is_good_length <- vapply(seq_along(gs), function(i) {
      if (is.na(right_length)) return(TRUE)
      isTRUE(all.equal(right_length, length(gs[[i]]$neighbours),
                       check.attributes = FALSE))
    }, FUN.VALUE = logical(1))
    if (!all(is_good_length)) {
      out <- paste0("Graphs of ", paste(mar, collapse = ", "),
                   " in sample ", sample_id, " do not have the right length.\n")
    }
  }
  return(out)
}

.check_graphs_sample <- function(gss, right_lengths) {
  # For each column of the spatialGraphs DF, which stands for each sample_id
  mar <- rownames(gss)
  sample_id <- names(gss)
  gsl <- as.list(gss[,1])
  lapply(seq_along(gsl), function(i) {
    .check_graphs_item(gsl[[i]], right_lengths[mar[i]], mar[i], sample_id)
  })
}

.check_graphs <- function(x) {
  sg <- int_metadata(x)$spatialGraphs
  if (is.null(sg)) return(character(0))
  if (is(sg, "DataFrame")) {
    if (!setequal(rownames(sg), c("row", "col", "annot"))) {
      return("Row names of spatialGraphs must be 'row', 'col', and 'annot'.\n")
    }
  } else {
    return("spatialGraphs must be a DataFrame whose columns are sample_ids ",
           "and whose rows are margins (rows, columns, annotation).\n")
  }
  # Check each column of sg, which stands for a sample_id
  right_lengths <- c(row = nrow(x), col = NA_integer_, annot = NA_integer_)
  unlist(lapply(seq_along(sg), function(i) {
    right_lengths[["col"]] <- sum(colData(x)$sample_id == names(sg)[i])
    .check_graphs_sample(sg[, i, drop = FALSE], right_lengths)
  }))
}

.check_graph_sample_id <- function(x) {
  if (is.null(int_metadata(x)$spatialGraphs)) return(character(0))
  ids_coldata <- as.character(unique(colData(x)$sample_id))
  graph_sample_ids <- names(int_metadata(x)$spatialGraphs)
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
  graphs_msg <- .check_graphs(object)
  id_msg <- .check_graph_sample_id(object)
  outs <- c(col_msg, row_msg, annot_msg, annot_msg2, graphs_msg, id_msg)
  outs <- unlist(outs)
  if (length(outs)) {
    return(outs)
  } else {
    return(TRUE)
  }
})
