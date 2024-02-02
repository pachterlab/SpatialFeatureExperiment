.concat_annotgeometries <- function(x, y) {
    # Check the names of the geometries in each object
    # For geometries with the same name, check column names for rbind
    # If the columns are different, then append sample_id to the geometry names
    ag_names_x <- annotGeometryNames(x)
    ag_names_y <- annotGeometryNames(y)
    if (is.null(x) && is.null(y)) {
        return(NULL)
    }
    names_intersect <- intersect(ag_names_x, ag_names_y)
    if (length(names_intersect)) {
        # A bit of repetition, just to get the names right
        which_rbind <- vapply(names_intersect, function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            setequal(names(ag_x), names(ag_y))
        }, FUN.VALUE = logical(1))

        out_rbind <- lapply(names_intersect[which_rbind], function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            ag_y <- ag_y[, names(ag_x)]
            rbind(ag_x, ag_y)
        })
        out_append <- lapply(names_intersect[!which_rbind], function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            setNames(list(ag_x, ag_y), paste(n,
                c(ag_x$sample_id[1], ag_y$sample_id[1]),
                sep = "_"
            ))
        })
        out_append <- unlist(out_append, recursive = FALSE)
        names(out_rbind) <- names_intersect[which_rbind]
        out_x <- annotGeometries(x)[setdiff(ag_names_x, ag_names_y)]
        out_y <- annotGeometries(y)[setdiff(ag_names_y, ag_names_x)]
    } else {
        out_rbind <- out_append <- NULL
        out_x <- annotGeometries(x)
        out_y <- annotGeometries(y)
    }
    out <- c(out_rbind, out_append, out_x, out_y)
    return(out)
}

.rg_sample_rename <- function(nms, sample_id) {
    if (is.null(nms)) return(NULL)
    # Append sample_id to items that don't already have sample_id
    # For one sample
    inds <- !str_detect(nms, paste0(sample_id, "$"))
    nms[inds] <- paste(nms[inds], sample_id, sep = "_")
    nms
}

.rg_xy_rename <- function(nms, sample_id, overlap, xy) {
    if (is.null(nms)) return(NULL)
    inds <- str_detect(nms, paste0(sample_id, "$"))
    nms[inds] <- str_remove(nms[inds], paste0("_", sample_id, "$"))
    which_overlap <- which(nms %in% overlap)
    nms[which_overlap] <- paste(nms[which_overlap], xy, sep = "_")
    nms <- paste(nms, sample_id, sep = "_")
    nms
}

.has_any_sample <- function(nms, sample_ids) {
    patterns <- paste0(sample_ids, "$")
    out <- vapply(patterns, str_detect, string = nms, FUN.VALUE = logical(length(nms)))
    # out is a matrix with patterns as columns
    apply(out, 1, any)
}

.add_xy <- function(rgns, sample_ids, xy) {
    pattern <- paste0("_", paste0(sample_ids, collapse = "|"), "$")
    base <- str_remove(rgns, pattern)
    samples <- str_extract(rgns, pattern)
    samples[is.na(samples)] <- ""
    paste0(base, "_", xy, samples)
}

.disambiguate <- function(rgns, samples, other, xy) {
    if (is.null(rgns)) return(NULL)
    has_samples <- .has_any_sample(rgns, samples)
    to_add_x <- has_samples & rgns %in% other
    rgns[to_add_x] <- .add_xy(rgn1[to_add_x], samples, xy)
    rgns
}

.rg_xy_multi <- function(rgn1, rgn2, samples1, samples2) {
    # a. items that overlap, with the same sample_ids, add x or y
    samples <- intersect(samples1, samples2)
    rgn1 <- .disambiguate(rgn1, samples, rgn2, "x")
    rgn2 <- .disambiguate(rgn2, samples, rgn1, "y")
    # b. items with different sample_ids, no need to change
    list(x = rgn1, y = rgn2)
}

.combine_rowgeos <- function(x, y) {
    rg_names_x <- rowGeometryNames(x)
    rg_names_y <- rowGeometryNames(y)
    if (is.null(rg_names_x) && is.null(rg_names_y))
        return(NULL)

    samples_x <- sampleIDs(x)
    samples_y <- sampleIDs(y)
    # Scenarios:
    # 1. x and y each has one sample, with different sample_ids. In this case
    # I'll simply append the sample_id if not already present
    if (length(samples_x) == 1L && length(samples_y) == 1L) {
        if (samples_x != samples_y) {
            rg_names_x <- .rg_sample_rename(rg_names_x, samples_x)
            rg_names_y <- .rg_sample_rename(rg_names_y, samples_y)
        } else {
            # 2. x and y have different sample_ids and rowGeometry with the same name
            # What to do? Maybe I'll just append "y" on the one from y. Sometimes they're
            # meant to be the same sample but sometimes they're not and the user forgot
            # to change it. What to do?
            # Don't want to break my rule to put sample_id by the end of the name
            # So it will become name_x/y_sample
            overlap <- intersect(rg_names_x, rg_names_y)
            rg_names_x <- .rg_xy_rename(rg_names_x, samples_x, overlap, "x")
            rg_names_y <- .rg_xy_rename(rg_names_y, samples_y, overlap, "y")
        }
    } else {
        # 3. At least one of x and y have multiple samples, some of which are for all
        # samples in the object and some are sample-specific. What to do? For those
        # that are not sample-specific or have overlapping names in x and y, I'll
        # append x or y to the name.
        new_rgns <- .rg_xy_multi(rg_names_x, rg_names_y, samples_x, samples_y)
        rg_names_x <- new_rgns[[1]]
        rg_names_y <- new_rgns[[2]]

        # Then what if we split the SFE object to form two samples? Might not need to split rowGeometries?
        # Let me not deal with this case (1 and 2 here) right now, until I implement the tree.
        # 1. items that overlap, without sample_id, i.e. for all samples in the object
        # This system of _sample needs to generalize when I implement the tree of sample relations
        # 2. items that don't overlap, without sample_id => now what, it's not sample specific but doesn't apply for the entire SFE object?
        # So that's connected to the tree.
    }
    rgsx <- rowGeometries(x)
    rgsy <- rowGeometries(y)

    names(rgsx) <- rg_names_x
    names(rgsy) <- rg_names_y
    return(c(rgsx, rgsy))
    # Output: new list of the combined rowGeometries
}

#' Concatenate SpatialFeatureExperiment objects
#'
#' On top of the \code{cbind} method of \code{SpatialExperiment}, this method is
#' needed to properly merge the \code{spatialGraphs} field in the different SFE
#' objects.
#'
#' @param ... SFE objects to cbind.
#' @param deparse.level See \code{?\link[base]{rbind}}.
#' @return A combined SFE object.
#' @importFrom BiocGenerics cbind
#' @export
#' @concept Non-spatial operations
#' @examples
#' library(SFEData)
#' sfe_small <- McKellarMuscleData(dataset = "small")
#' sfe_small2 <- McKellarMuscleData(dataset = "small2")
#' sfe2 <- cbind(sfe_small, sfe_small2)
setMethod(
    "cbind", "SpatialFeatureExperiment",
    function(..., deparse.level = 1) {
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }
        args <- list(...)
        new_rg <- Reduce(.combine_rowgeos, args)
        args <- lapply(args, function(x) {
            rowGeometries(x) <- NULL
            x
        })
        out <- do.call(
            callNextMethod,
            c(args, list(deparse.level = deparse.level))
        )
        colnames(out) <- make.unique(colnames(out), sep = "-")
        rowGeometries(out) <- new_rg
        # Combine the annotGeometries
        has_ag <- vapply(args, function(a) length(annotGeometries(a)) > 0L,
            FUN.VALUE = logical(1)
        )
        if (any(has_ag)) {
            new_ag <- Reduce(.concat_annotgeometries, args)
            int_metadata(out)[["annotGeometries"]] <- NULL
            int_metadata(out)[["annotGeometries"]] <- new_ag
        }
        # Combine the spatialGraphs
        # TODO: benchmark spatial methods. Consider switching to sparse matrices
        has_sg <- vapply(args, function(a) !is.null(spatialGraphs(a)),
            FUN.VALUE = logical(1)
        )
        if (any(has_sg)) {
            sgs_use <- lapply(args, function(a) {
                if (is.null(spatialGraphs(a))) {
                    .initialize_spatialGraphs(a)
                } else {
                    int_metadata(a)[["spatialGraphs"]]
                }
            })
            new_sgs <- do.call(cbind, sgs_use)
            int_metadata(out)[["spatialGraphs"]] <- NULL
            int_metadata(out)$spatialGraphs <- new_sgs
        }

        # Clean up duplicate version field in int_metadata
        # Bug in SCE
        im_names <- unique(names(int_metadata(out)))
        names(int_metadata(out)) <- make.unique(names(int_metadata(out)))
        int_metadata(out) <- int_metadata(out)[im_names]
        out
    }
)

# TODO: merge function, kind of like cbind but allow for different number of rows
# Error when the row names don't overlap, warning when the overlap is too small
# Add 0's to the matrices and empty geometries to rowGeometries. Kind of like full_join
