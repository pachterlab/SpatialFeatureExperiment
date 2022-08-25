#' Get all unique sample IDs
#'
#' The title is self-explanatory.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @return A character vector of all unique entries of the \code{sample_id}
#' column in \code{colData(x)}.
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sampleIDs(sfe)
sampleIDs <- function(sfe) unique(colData(sfe)$sample_id)

#' Change sample IDs
#'
#' Change sample IDs in all fields of the SFE object where sample IDs are
#' present, not just the colData.
#'
#' @inheritParams sampleIDs
#' @param replacement A named character vector whose names are the existing
#'   sample IDs to be changed and whose values are the corresponding
#'   replacements.
#' @return An SFE object.
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
#' sampleIDs(sfe)
changeSampleIDs <- function(sfe, replacement) {
  for (i in seq_along(replacement)) {
    colData(sfe)$sample_id[colData(sfe)$sample_id == names(replacement)[i]] <- replacement[i]
    gs_names <- names(int_metadata(sfe)$spatialGraphs)
    names(int_metadata(sfe)$spatialGraphs)[gs_names == names(replacement)[i]] <- replacement[i]
    if (length(int_metadata(sfe)$annotGeometries)) {
      for (n in names(int_metadata(sfe)$annotGeometries)) {
        ag <- int_metadata(sfe)$annotGeometries[[n]]
        ind <- ag$sample_id == names(replacement)[i]
        int_metadata(sfe)$annotGeometries[[n]]$sample_id[ind] <- replacement[i]
      }
    }
  }
  sfe
}

.check_sample_id <- function(x, sample_id, one = TRUE) {
  if (is.null(sample_id)) {
    sample_id <- sampleIDs(x)
    if (length(sample_id) > 1L) {
      stop("There are more than one sample in this object. sample_id must be specified")
    }
  } else if (identical(sample_id, "all")) {
    sample_id <- sampleIDs(x)
  } else if (!all(sample_id %in% sampleIDs(x))) {
    sample_use <- intersect(sample_id, sampleIDs(x))
    if (!length(sample_use)) {
      stop("None of the samples is present in the SFE object.")
    }
    sample_show <- setdiff(sample_id, sampleIDs(x))
    warning("Sample(s) ", paste(sample_show, sep = ","), " is/are absent from the SFE object.")
    sample_id <- sample_use
  }
  if (one) {
    if (length(sample_id) > 1L)
      stop("Only one sample can be specified at a time.")
  }
  sample_id
}

.rm_empty_geometries <- function(g, MARGIN) {
  empty_inds <- st_is_empty(g)
  if (MARGIN < 3) {
    if (any(empty_inds))
      stop("Empty geometries found in dimGeometry.")
  } else {
    g <- g[!empty_inds,]
  }
  g
}

.translate_value <- function(x, translate, value) {
  if (translate && !is.null(int_metadata(x)$orig_bbox)) {
    value$geometry <- value$geometry - int_metadata(x)$orig_bbox[c("xmin", "ymin")]
  }
  value
}
