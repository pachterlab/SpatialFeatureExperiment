# Internal functions related to int_colData and int_elementMetadata (essentially int_rowData)
# I created this file to generalize internal functions written for dimGeometries
# to localResults and potentially other fields in int_colData.
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

.get_intdimdata_all <- function(x, MARGIN, withDimnames = TRUE, getfun, key) {
  value <- .get_internal_all(x,
                             getfun=getfun,
                             key=key)

  if (withDimnames) {
    for (i in seq_along(value)) {
      rownames(value[[i]]) <- dimnames(x)[[MARGIN]]
    }
  }
  value
}

.set_intdimdata_all <- function(x, MARGIN, withDimnames = TRUE,
                                translate = TRUE, sf = TRUE, getfun, setfun,
                                key, xdimfun, funstr, xdimstr, value, ...) {
  if (sf) value <- .df2sf_list(value, ...)
  if (withDimnames) {
    for (v in seq_along(value)) {
      value[[v]] <- .check_dimgeo_names(x, value[[v]], MARGIN,
                                        withDimnames=withDimnames,
                                        vname=sprintf("value[[%s]]", v),
                                        fun=funstr)
      if (sf) value[[v]] <- .translate_value(x, translate, value[[v]])
    }
  }
  .set_internal_all(x, value,
                    getfun=getfun,
                    setfun=setfun,
                    key=key,
                    convertfun=NULL,
                    xdimfun=xdimfun,
                    vdimfun=nrow,
                    funstr=funstr,
                    xdimstr=xdimstr,
                    vdimstr="rows")
}

.out_intdimdata_id <- function(x, out, MARGIN, sample_id) {
  sample_id <- .check_sample_id(x, sample_id, one = FALSE)
  if (!is.null(sample_id)) {
    if (MARGIN == 1L) {
      # OK, maybe applicable, say to crop the rowGeometries by bbox of a sample
      # I'll consider that later.
      message("sample_id is not applicable to rowGeometries.")
    } else if (MARGIN == 2L) {
      # Somehow I lose the rownames after row subsetting
      inds <- colData(x)$sample_id %in% sample_id
      rns <- rownames(out)[inds]
      out <- out[inds,]
      rownames(out) <- rns
    } else {
      out <- out[out$sample_id %in% sample_id,]
    }
  }
  out
}

.get_internal_id <- function(x, type, MARGIN, sample_id, withDimnames,
                             .get_internal_fun, getfun, key, funstr, substr) {
  out <- .get_internal_integer(x, type,
                               getfun=getfun,
                               key=key,
                               funstr=funstr,
                               substr=substr)

  if (withDimnames) {
    rownames(out) <- dimnames(x)[[MARGIN]]
  }
  .out_intdimdata_id(x, out, MARGIN, sample_id)
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

.intdimdata_partial_replace <- function(existing, value, nrow_full, rownames_full,
                                      all_sample_ids, sample_id, sf = TRUE) {
  if (is.null(existing)) {
    existing <- matrix(nrow = nrow_full, ncol = ncol(value),
                       dimnames = list(rownames_full,
                                       colnames(value)))
    existing <- as.data.frame(existing)
    if (sf) {
      existing$geometry <- st_sfc(lapply(seq_len(nrow(existing)),
                                         function(t) st_geometrycollection()))
      existing <- st_sf(existing, sf_column_name = "geometry",
                        row.names = rownames_full)
      existing <- existing[,colnames(value)]
    }
  } else {
    value <- .reconcile_cols(existing, value)
  }
  existing[all_sample_ids %in% sample_id,] <- value
  existing
}

.set_intdimdata_id <- function(x, value, sample_id, type, MARGIN, sf = TRUE) {
  sample_id <- .check_sample_id(x, sample_id, one = FALSE)
  if (!is.null(sample_id) && any(!sampleIDs(x) %in% sample_id)) {
    if (MARGIN == 1L) {
      message("sample_id is not applicable to rowGeometries.")
    } else {
      # Assuming that the order in value is the same as the
      # order of geometries for this sample in colGeometries
      existing <- .getfun(MARGIN)(x)[[.dg_key(MARGIN)]][[type]]
      value <- .intdimdata_partial_replace(existing, value, ncol(x), colnames(x),
                                         colData(x)$sample_id, sample_id, sf = sf)
    }
  }
  value
}

.initialize_intdimdata <- function(x, .get_all_fun, .set_all_fun,
                                   withDimnames = TRUE, ...) {
  if (!length(.get_all_fun(x, ...))) {
    x <- .set_all_fun(x, withDimnames = withDimnames, value = list(), ...)
  }
  x
}

.set_internal_id <- function(x, type, MARGIN, sample_id, withDimnames = TRUE,
                             translate = TRUE, sf = TRUE, .set_internal_fun,
                             getfun, setfun, key, xdimfun, funstr, xdimstr,
                             substr, value, ...) {
  x <- .initialize_intdimdata(x, .get_all_fun = dimGeometries,
                              .set_all_fun = `dimGeometries<-`,
                              MARGIN = MARGIN, withDimnames = withDimnames)
  if (sf) value <- .df2sf_in_list(value, ...)
  value <- .set_intdimdata_id(x, value, sample_id, type, MARGIN)
  value <- .check_dimgeo_names(x, value, MARGIN = MARGIN,
                               withDimnames = withDimnames, fun = funstr)
  if (sf) value <- .translate_value(x, translate, value)
  .set_internal_fun(x, type, value,
                    getfun=getfun,
                    setfun=setfun,
                    key=key,
                    convertfun=NULL,
                    xdimfun=xdimfun,
                    vdimfun=nrow,
                    funstr=funstr,
                    xdimstr=xdimstr,
                    vdimstr="rows",
                    substr=substr)
}
