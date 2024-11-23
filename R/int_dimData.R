# Internal functions related to int_colData and int_elementMetadata
# (essentially int_rowData)
# I created this file to generalize internal functions written for dimGeometries
# to localResults and potentially other fields in int_colData.
.getfun <- function(MARGIN) {
    switch(MARGIN,
        int_elementMetadata,
        int_colData
    )
}

.setfun <- function(MARGIN) {
    switch(MARGIN,
        `int_elementMetadata<-`,
        `int_colData<-`
    )
}

.xdimstr <- function(MARGIN) {
    switch(MARGIN,
        "nrow",
        "ncol"
    )
}

.xdimfun <- function(MARGIN) {
    switch(MARGIN,
        nrow,
        ncol
    )
}

.dg_key <- function(MARGIN) {
    switch(MARGIN,
        "rowGeometries",
        "colGeometries"
    )
}

.dg_key2 <- function(MARGIN) {
    switch (MARGIN, "rowGeometry", "colGeometry", "annotGeometry")
}

.unnamed <- "unnamed"
# Modified from SCE to generalize to both rows and columns
.check_dimgeo_names <- function(reference, incoming, MARGIN, withDimnames,
                                fun = "dimGeometry", vname = "value") {
    if (!is.null(incoming)) {
        rni <- rownames(incoming)
        cnr <- dimnames(reference)[[MARGIN]]
        fun_show <- switch(MARGIN,
            "rownames",
            "colnames"
        )
        if (withDimnames && !is.null(rni)) {
            if (!all(rni %in% cnr)) {
                msg <- paste0(
                    "non-NULL 'rownames(", vname, ")' should all be in '",
                    fun_show, "(x)' for '", fun,
                    "<-'."
                )
                stop(strwrap(msg))
            }
            if (!identical(cnr, rni)) {
                # Do the reordering if they have different orders
                rni <- rni[match(cnr, rni)]
                incoming <- incoming[rni, ]
                if (anyNA(rownames(incoming))) {
                    rownames(incoming) <- cnr
                }
            }
        }
    }
    incoming
}
.get_internal_all <- SingleCellExperiment:::.get_internal_all
.set_internal_all <- SingleCellExperiment:::.set_internal_all
.get_internal_names <- SingleCellExperiment:::.get_internal_names
.set_internal_names <- SingleCellExperiment:::.set_internal_names

.get_internal <- function(x, index, getfun, key, funstr, substr, namestr) {
    x <- updateObject(x)
    internals <- getfun(x)[[key]]

    tryCatch({
        internals[, index]
    }, error=function(err) {
        if (is.numeric(index)) {
            stop("invalid subscript '", substr, "' in '", funstr, "(<", class(x), ">, type=\"numeric\", ...)':\n  ",
                 conditionMessage(err))
        } else if (is.character(index)) {
            stop("invalid subscript '", substr, "' in '", funstr, "(<", class(x), ">, type=\"character\", ...)':\n  ",
                 "'", index, "' not in '", namestr, "(<", class(x), ">)'")
        }
    })
}

.set_internal <- function(x, type, value, getfun, setfun, key,
                          convertfun, xdimfun, vdimfun, funstr, xdimstr, vdimstr, substr)
{
    x <- updateObject(x)
    if (!is.null(value)) {
        if (!is.null(convertfun)) {
            value <- convertfun(value)
        }
        if (!identical(vdimfun(value), xdimfun(x))) {
            stop("invalid 'value' in '", funstr, "(<", class(x), ">) <- value':\n  ",
                 "'value' should have number of ", vdimstr, " equal to '", xdimstr, "(x)'")
        }
    }

    internals <- getfun(x)
    if (is.numeric(type) && type[1] > ncol(internals[[key]])) {
        stop("'", substr, "' out of bounds in '", funstr, "(<", class(x), ">, type='numeric')")
    }
    internals[[key]][[type]] <- value
    x <- setfun(x, internals)
    x
}

.get_intdimdata_all <- function(x, MARGIN, withDimnames = TRUE, getfun, key) {
    value <- .get_internal_all(x,
        getfun = getfun,
        key = key
    )

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
    for (v in seq_along(value)) {
        value[[v]] <- .check_dimgeo_names(x, value[[v]], MARGIN,
                                          withDimnames = withDimnames,
                                          vname = sprintf("value[[%s]]", v),
                                          fun = funstr
        )
        if (sf) value[[v]] <- .translate_value(x, translate, value[[v]])
    }
    .set_internal_all(x, value,
        getfun = getfun,
        setfun = setfun,
        key = key,
        convertfun = NULL,
        xdimfun = xdimfun,
        vdimfun = nrow,
        funstr = funstr,
        xdimstr = xdimstr,
        vdimstr = "rows"
    )
}

.out_intdimdata_id <- function(x, out, MARGIN, sample_id) {
    if (MARGIN == 1L) {
        return(out)
    }
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    samples <- sampleIDs(x)
    if (!is.null(sample_id)) {
        if (MARGIN == 2L) {
            if (length(samples) > 1L) {
                # Somehow I lose the rownames after row subsetting sf data frames
                inds <- colData(x)$sample_id %in% sample_id
                rns <- rownames(out)[inds]
                out <- out[inds, , drop = FALSE]
                rownames(out) <- rns
            }
        } else { # annotGeometry
            if (length(samples) > 1L)
                out <- out[out$sample_id %in% sample_id, , drop = FALSE]
        }
    }
    out
}

.get_internal_id <- function(x, type, MARGIN, sample_id, withDimnames,
                             .get_internal_fun, getfun, key, funstr, substr,
                             ...) {
    out <- .get_internal_fun(x, type,
        getfun = getfun,
        key = key,
        funstr = funstr,
        substr = substr,
        ...
    )

    if (withDimnames) {
        rownames(out) <- dimnames(x)[[MARGIN]]
    }
    .out_intdimdata_id(x, out, MARGIN, sample_id)
}

.make_na_df <- function(template, nrow, rownames, sf) {
    out <- lapply(template, function(x) {
        if (is.vector(x)) {
            out <- rep(NA, nrow)
        } else if (inherits(x, "sfc")) {
            out <- st_sfc(lapply(seq_len(nrow), function(t) st_geometrycollection()))
        } else {
            # Only deal with columns that are themselves simple matrices
            # or data frames
            # Not when the columns are data frames with columns that are matrices
            # or data frames
            out <- matrix(
                nrow = nrow, ncol = ncol(x),
                dimnames = list(rownames, colnames(x))
            )
            if (is.data.frame(x)) out <- as.data.frame(out)
            if (inherits(x, "DFrame")) out <- DataFrame(out)
            out <- I(out)
        }
        out
    })
    names(out) <- names(template)
    S4 <- inherits(template, "DFrame")
    if (sf && S4) stop("Please use S3 data.frame for sf.")
    df_fun <- if (S4) DataFrame else data.frame
    out <- df_fun(out)
    if (sf) {
        out <- st_sf(out, sf_column_name = "geometry")
    }
    rownames(out) <- rownames
    out
}

.reconcile_cols <- function(existing, value, sf = FALSE) {
    # Assume that both existing and value have the geometry column if sf = TRUE
    # Don't remove columns from `existing` no matter what
    # Add columns only present in `value` to `existing`
    if (any(!colnames(value) %in% colnames(existing))) {
        diff_cols <- setdiff(colnames(value), colnames(existing))
        additional_cols <- .make_na_df(value[, diff_cols, drop = FALSE],
            nrow = nrow(existing),
            rownames = rownames(existing), sf = sf
        )
        existing <- cbind(existing, additional_cols)
    }
    existing
}

.intdimdata_partial_replace <- function(existing, value, nrow_full,
                                        rownames_full, all_sample_ids,
                                        sample_id, sf = TRUE) {
    if (is.null(existing)) {
        existing <- .make_na_df(value,
            nrow = nrow_full, rownames = rownames_full,
            sf = sf
        )
    } else {
        # Should add columns in value but not in existing
        existing <- .reconcile_cols(existing, value)
    }
    existing[all_sample_ids %in% sample_id, names(value)] <- value
    rownames(existing) <- rownames_full
    existing
}

.set_intdimdata_id <- function(x, value, sample_id, type, MARGIN, sf = TRUE,
                               getfun, key) {
    if (MARGIN == 1L) {
        return(value)
    }
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (!is.null(sample_id) && any(!sampleIDs(x) %in% sample_id)) {
        # Assuming that the order in value is the same as the
        # order of geometries for this sample in colGeometries
        existing <- getfun(x)[[key]][[type]]
        value <- .intdimdata_partial_replace(existing, value, ncol(x),
                                             colnames(x), colData(x)$sample_id,
                                             sample_id,
                                             sf = sf
        )
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
                             translate = TRUE, sf = TRUE,
                             .get_all_fun, .set_all_fun,
                             .set_internal_fun,
                             getfun, setfun, key, xdimfun, funstr, xdimstr,
                             substr, value, ...) {
    x <- .initialize_intdimdata(x,
        .get_all_fun = .get_all_fun,
        .set_all_fun = .set_all_fun,
        MARGIN = MARGIN, withDimnames = withDimnames
    )
    if (sf) value <- .df2sf_in_list(value, ...)
    value <- .set_intdimdata_id(x, value, sample_id, type, MARGIN,
        sf = sf,
        getfun, key
    )
    value <- .check_dimgeo_names(x, value,
        MARGIN = MARGIN,
        withDimnames = withDimnames, fun = funstr
    )
    if (sf) value <- .translate_value(x, translate, value, sample_id)
    .set_internal_fun(x, type, value,
        getfun = getfun,
        setfun = setfun,
        key = key,
        convertfun = NULL,
        xdimfun = xdimfun,
        vdimfun = nrow,
        funstr = funstr,
        xdimstr = xdimstr,
        vdimstr = "rows",
        substr = substr
    )
}

.get_internal_feature <- function(x, type, feature, colGeometryName,
                                  annotGeometryName, sample_id, withDimnames,
                                  .get_internal_fun, simplify = TRUE,
                                  swap_rownames = NULL,...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (missing(type)) type <- 1L
    feature <- .symbol2id(x, feature, swap_rownames)
    if (is.null(colGeometryName) && is.null(annotGeometryName)) {
        lr <- .get_internal_id(x,
            type = type, sample_id = sample_id,
            withDimnames = withDimnames,
            .get_internal_fun = .get_internal_fun, ...
        )
    } else if (!is.null(colGeometryName)) {
        cg <- colGeometry(x, type = colGeometryName, sample_id = sample_id)
        lr <- cg$localResults[[type]]
    } else if (!is.null(annotGeometryName)) {
        ag <- annotGeometry(x, type = annotGeometryName, sample_id = sample_id)
        lr <- ag$localResults[[type]]
    }
    if (is.null(feature)) {
        feature_use <- names(lr)
    } else {
        feature_use <- intersect(feature, names(lr))
    }
    if (!length(feature_use)) {
        message_use <- paste0("None of the features are present in ", type)
        if (!is.null(colGeometryName)) {
            message_use <- paste0(
                message_use, " in colGeometry ",
                colGeometryName
            )
        } else if (!is.null(annotGeometryName)) {
            message_use <- paste0(
                message_use, " in annotGeometry ",
                annotGeometryName
            )
        }
        stop(message_use)
    }
    if (length(feature_use) < length(feature)) {
        warning(
            "Features ", paste(setdiff(feature, feature_use), collapse = ", "),
            " are absent in ", type
        )
    }
    out <- lr[, feature_use, drop = FALSE]
    out <- as.list(out)
    if (simplify) out <- out[[1]]
    out
}

.set_geometry_localResults <- function(x, lr_type, feature, sample_id,
                                       colGeometryName, annotGeometryName,
                                       value) {
    if (!is.null(colGeometryName)) {
        get_geom_fun <- colGeometry
        set_geom_fun <- `colGeometry<-`
        geometry_type <- colGeometryName
    } else {
        get_geom_fun <- annotGeometry
        set_geom_fun <- `annotGeometry<-`
        geometry_type <- annotGeometryName
    }

    if (inherits(value, "DFrame")) value <- as.list(value)
    value <- .value2df(value, TRUE, feature = feature)

    g <- get_geom_fun(x, type = geometry_type, sample_id = "all")
    if (is.null(annotGeometryName) && !is.null(colGeometryName)) {
        sample_index <- colData(x)$sample_id %in% sample_id
        all_sample_ids <- colData(x)$sample_id
    } else {
        sample_index <- g$sample_id %in% sample_id
        all_sample_ids <- g$sample_id
    }
    init_lr <- is.null(g$localResults)
    init_type <- is.null(g$localResults[[lr_type]])
    init <- data.frame(..1 = rep(NA, nrow(g)))
    if (init_lr) g$localResults <- init
    if (init_type) g$localResults[[lr_type]] <- init
    lr <- g$localResults[[lr_type]][sample_index, , drop = FALSE]
    lr[, feature] <- value
    if (any(!all_sample_ids %in% sample_id)) {
        g$localResults[[lr_type]] <-
            .intdimdata_partial_replace(
                existing = g$localResults[[lr_type]],
                value = lr,
                nrow_full = nrow(g),
                rownames_full = rownames(g),
                all_sample_ids = all_sample_ids,
                sample_id = sample_id, sf = FALSE
            )
    } else {
        g$localResults[[lr_type]] <- lr
    }
    if (init_type) g$localResults[[lr_type]][["..1"]] <- NULL
    if (init_lr) g$localResults[["..1"]] <- NULL
    x <- set_geom_fun(x, type = geometry_type, sample_id = "all", value = g,
                      translate = FALSE)
    x
}

.set_internal_feature <- function(x, type, feature, colGeometryName,
                                  annotGeometryName, sample_id, withDimnames,
                                  .set_internal_fun, value, ...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (missing(type)) type <- 1L
    if (any(class(value) %in% c("list", "DFrame", "data.frame"))) {
        if (is.null(feature)) {
            feature <- names(value)
        } else {
            feature <- intersect(feature, names(value))
            value <- value[feature]
        }
    }
    not_geometry <- is.null(colGeometryName) && is.null(annotGeometryName)
    value <- .value2df(value, !not_geometry, feature = feature)

    if (not_geometry) {
        lr_type <- localResults(x)[[type]]
        if (is.null(lr_type)) {
            lr_type <- make_zero_col_DFrame(nrow = nrow(value))
        } else {
            lr_type <- lr_type[colData(x)$sample_id %in% sample_id, , drop = FALSE]
        }
        lr_type[, feature] <- value
        .set_internal_id(x,
            withDimnames = withDimnames, type = type,
            sample_id = sample_id, value = lr_type,
            .set_internal_fun = .set_internal_fun, ...
        )
    } else if (!is.null(colGeometryName)) {
        .set_geometry_localResults(x,
            lr_type = type, feature = feature,
            sample_id = sample_id,
            colGeometryName = colGeometryName,
            annotGeometryName = NULL, value = value
        )
    } else {
        .set_geometry_localResults(x,
            lr_type = type, feature = feature,
            sample_id = sample_id,
            colGeometryName = NULL,
            annotGeometryName = annotGeometryName,
            value = value
        )
    }
}
