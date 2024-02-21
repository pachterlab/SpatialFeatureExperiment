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
            ag_x <- annotGeometry(x, n, sample_id = "all")
            ag_y <- annotGeometry(y, n, sample_id = "all")
            setequal(names(ag_x), names(ag_y))
        }, FUN.VALUE = logical(1))

        out_rbind <- lapply(names_intersect[which_rbind], function(n) {
            ag_x <- annotGeometry(x, n, sample_id = "all")
            ag_y <- annotGeometry(y, n, sample_id = "all")
            ag_y <- ag_y[, names(ag_x)]
            rbind(ag_x, ag_y)
        })
        out_append <- lapply(names_intersect[!which_rbind], function(n) {
            ag_x <- annotGeometry(x, n, sample_id = "all")
            ag_y <- annotGeometry(y, n, sample_id = "all")
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

#' Concatenate SpatialFeatureExperiment objects
#'
#' On top of the \code{cbind} method of \code{SpatialExperiment}, this method is
#' needed to properly merge the \code{spatialGraphs} field in the different SFE
#' objects. \code{rowGeometries} and \code{annotGeometries} also need to be
#' combined properly.
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
        args <- list(...)
        if (length(args) < 2L) return(args[[1]]) # just like cbind for matrix
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }
        # Like in SPE, make sure sample_ids are unique
        sids_list <- lapply(args, sampleIDs)
        sids <- unlist(sids_list)
        sid_df <- data.frame(id = sids,
                             object = rep(seq_along(sids_list), times = lengths(sids_list)))
        if (length(sids) > length(unique(sids))) {
            message(
                "'sample_id's are duplicated across",
                " 'SpatialFeatureExperiment' objects to cbind;",
                " appending sample indices.")
            # Might not be fast but since I only expect a small number of SFE
            # objects it should be fine
            for (i in seq_along(args)) {
                ind <- sid_df$object < i
                for (s in sids_list[[i]]) {
                    n_existing <- sum(s == sid_df$id[ind])
                    if (n_existing == 0L) next
                    id_new <- paste(s, n_existing, sep = "_")
                    names(id_new) <- s
                    args[[i]] <- changeSampleIDs(args[[i]], id_new)
                }
            }
        }
        # Rename rowGeometries that are not sample specific
        for (i in seq_along(args)) {
            if (length(sids_list[[i]]) == 1L) {
                rgns <- rowGeometryNames(args[[i]])
                if (is.null(rgns)) next
                s <- sampleIDs(args[[i]])
                to_add <- !grepl(paste0(s, "$"), rgns)
                rgns[to_add] <- paste(rgns[to_add], s, sep = "_")
                rowGeometryNames(args[[i]]) <- rgns
            }
        }
        # Multiple samples, with rowGeometries that are not sample specific
        rgns_list <- lapply(args, rowGeometryNames)
        rgns <- unlist(rgns_list)
        if (length(rgns) > length(unique(rgns))) {
            # Presumably some of the SFE objects have rowGeometries
            rgn_df <- data.frame(id = rgns,
                                 object = rep(seq_along(rgns_list), times = lengths(rgns_list)))
            for (i in seq_along(args)) {
                if (is.null(rgns_list[[i]])) next
                ind <- rgn_df$object < i
                rgns_use <- rgns_list[[i]]
                rgns_old <- rgns_use
                for (j in seq_along(rgns_use)) {
                    s <- rgns_use[j]
                    n_existing <- sum(s == rgn_df$id[ind])
                    if (n_existing == 0L) next
                    rgns_use[j] <- paste(s, n_existing, sep = "_")
                    # This is tentative. I'll implement the tree thing to make it clearer
                }
                if (!identical(rgns_use, rgns_old))
                    rowGeometryNames(args[[i]]) <- rgns_use
            }
        }
        # Still need to combine rowGeometries separately
        if (length(rgns)) {
            rgs <- lapply(args, rowGeometries)
            rgs <- do.call(c, rgs)
            args <- lapply(args, function(x) {
                rowGeometries(x) <- NULL
                x
            })
        }
        out <- do.call(
            callNextMethod,
            c(args, list(deparse.level = deparse.level))
        )
        colnames(out) <- make.unique(colnames(out), sep = "-")
        if (length(rgns)) {
            rowGeometries(out) <- rgs
        }
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
