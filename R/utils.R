.getfun <- function(MARGIN) switch (MARGIN, int_elementMetadata, int_colData)

.setfun <- function(MARGIN) switch (MARGIN, `int_elementMetadata<-`, `int_colData<-`)

.xdimstr <- function(MARGIN) switch (MARGIN, "nrow", "ncol")

.xdimfun <- function(MARGIN) switch (MARGIN, nrow, ncol)

.dg_key <- function(MARGIN) switch (MARGIN, "rowGeometries", "colGeometries")

.unnamed <- "unnamed"
# Modified from SCE to generalize to both rows and columns
.check_dimgeo_names <- function(reference, incoming, MARGIN, withDimnames,
                                fun='dimGeometry', vname='value') {
  if (!is.null(incoming)) {
    rni <- rownames(incoming)
    cnr <- dimnames(reference)[[MARGIN]]
    fun_show <- switch (MARGIN, "rownames", "colnames")
    if (withDimnames && !is.null(rni)) {
      if (!setequal(cnr, rni)) {
        msg <- paste0("non-NULL 'rownames(", vname, ")' should be the same as '",
                      fun_show, "(x)' for '", fun,
                      "<-'.")
        stop(paste(strwrap(msg), collapse="\n"))
      } else if (!identical(cnr, rni)) {
        # Do the reordering if they have different orders
        rni <- rni[match(cnr, rni)]
        incoming <- incoming[rni,]
      }
    }
  }
  incoming
}
.get_internal_all <- SingleCellExperiment:::.get_internal_all
.set_internal_all <- SingleCellExperiment:::.set_internal_all
.get_internal_integer <- SingleCellExperiment:::.get_internal_integer
.get_internal_names <- SingleCellExperiment:::.get_internal_names
.get_internal_missing <- SingleCellExperiment:::.get_internal_missing
.get_internal_character <- SingleCellExperiment:::.get_internal_character
.set_internal_names <- SingleCellExperiment:::.set_internal_names
.set_internal_missing <- SingleCellExperiment:::.set_internal_missing
.set_internal_numeric <- SingleCellExperiment:::.set_internal_numeric
.set_internal_character <- SingleCellExperiment:::.set_internal_character

#' Get all unique sample IDs
#'
#' The title is self-explanatory.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @return A character vector of all unique entries of the \code{sample_id}
#' column in \code{colData(x)}.
#' @export
#' @importFrom SummarizedExperiment colData colData<-
sampleIDs <- function(x) unique(colData(x)$sample_id)

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

.reconcile_cols <- function(existing, value) {
  # Assume that both existing and value have the geometry column
  if (any(!names(value) %in% names(existing))) {
    names_inter <- intersect(names(value), names(existing))
    value <- value[,names_inter]
  }
  if (any(!names(existing) %in% names(value))) {
    diff_cols <- setdiff(names(existing), names(value))
    additional_cols <- matrix(nrow = nrow(value),
                              ncol = length(diff_cols))
    colnames(additional_cols) <- diff_cols
    rownames(additional_cols) <- rownames(value)
    additional_cols <- as.data.frame(additional_cols)
    value <- cbind(value, additional_cols)
    value <- value[,names(existing)]
  }
  value
}
