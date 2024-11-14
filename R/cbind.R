#' @importFrom S4Vectors combineCols
.combine_fd <- function(fd_x, fd_y) {
    if (is.null(fd_x) && is.null(fd_y)) {
        return(NULL)
    }
    # Edge case: what if fd_x and fd_y have different numbers of rows
    empty <- make_zero_col_DFrame(nrow(fd_x) %||% nrow(fd_y))
    rns <- rownames(empty) <- rownames(fd_x) %||% rownames(fd_y)
    fd_x <- fd_x %||% empty
    fd_y <- fd_y %||% empty
    combineCols(fd_x, fd_y)
}

.combine_attr_fd <- function(x, y) {
    # x and y should have the same number of columns and rows
    o <- rbind(x, y)
    if (!is.null(attr(x, "featureData")) ||
        !is.null(attr(y, "featureData"))) {
        fd_x <- attr(x, "featureData")
        fd_y <- attr(y, "featureData")
        fd <- .combine_fd(fd_x, fd_y)
        attr(o, "featureData") <- fd
    }
    o
}
#' @importFrom SingleCellExperiment reducedDims reducedDims<- reducedDimNames
.combine_reduceddim_fd <- function(x, y) {
    # x and y are named lists of featureData
    # By this time x and y should have the same reducedDim or colGeometry names
    if (is.null(x) && is.null(y)) {
        return(NULL)
    }
    out <- lapply(names(x), function(n) {
        .combine_fd(x[[n]], y[[n]])
    })
    names(out) <- names(x)
    out
}

.concat_annotgeometries <- function(x, y) {
    # x, y are lists of annotGeometries from each SFE object
    # Check the names of the geometries in each object
    # For geometries with the same name, check column names for rbind
    # If the columns are different, then append sample_id to the geometry names
    if (is.null(x) && is.null(y)) {
        return(NULL)
    }
    ag_names_x <- names(x)
    ag_names_y <- names(y)
    names_intersect <- intersect(ag_names_x, ag_names_y)
    if (length(names_intersect)) {
        # A bit of repetition, just to get the names right
        which_rbind <- vapply(names_intersect, function(n) {
            setequal(names(x[[n]]), names(y[[n]]))
        }, FUN.VALUE = logical(1))

        out_rbind <- lapply(names_intersect[which_rbind], function(n) {
            .combine_attr_fd(x[[n]], y[[n]])
        })
        out_append <- lapply(names_intersect[!which_rbind], function(n) {
            ag_x <- x[[n]]
            ag_y <- y[[n]]
            setNames(list(ag_x, ag_y), paste(n,
                c(ag_x$sample_id[1], ag_y$sample_id[1]),
                sep = "_"
            ))
        })
        out_append <- unlist(out_append, recursive = FALSE)
        names(out_rbind) <- names_intersect[which_rbind]
        out_x <- x[setdiff(ag_names_x, ag_names_y)]
        out_y <- y[setdiff(ag_names_y, ag_names_x)]
    } else {
        out_rbind <- out_append <- NULL
        out_x <- x
        out_y <- y
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
#' @importFrom S4Vectors metadata<- metadata
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
            rgs <- lapply(rgs, as.list) # S4 SimpleList won't unlist
            names(rgs) <- NULL
            rgs <- unlist(rgs, recursive = FALSE)
            args <- lapply(args, function(x) {
                rowGeometries(x) <- NULL
                x
            })
        }
        # Remove rowData's metadata; I might regret it but I don't care so much
        # about the params
        args <- lapply(args, function(x) {
            metadata(rowData(x)) <- list()
            x
        })
        out <- do.call(
            callNextMethod,
            c(args, list(deparse.level = deparse.level))
        )
        colnames(out) <- make.unique(colnames(out), sep = "-")
        if (length(rgns)) {
            rowGeometries(out) <- rgs
        }
        has_ag <- vapply(args, function(a) length(annotGeometries(a)) > 0L,
            FUN.VALUE = logical(1)
        )
        if (any(has_ag)) {
            ags <- lapply(args, annotGeometries)
            new_ag <- Reduce(.concat_annotgeometries, ags)
            int_metadata(out)[["annotGeometries"]] <- NULL
            int_metadata(out)[["annotGeometries"]] <- new_ag
        }
        # Combine reducedDims featureData
        if (!is.null(reducedDimNames(out))) {
            fds <- lapply(args, function(x) {
                rdns <- reducedDimNames(x)
                f <- lapply(rdns, reducedDimFeatureData, sfe = x)
                names(f) <- rdns
                f
            })
            new_rd_fds <- Reduce(.combine_reduceddim_fd, fds)
            rdns <- reducedDimNames(out)
            for (n in rdns) {
                reducedDimFeatureData(out, n) <- new_rd_fds[[n]]
            }
        }
        # Combine colGeometry featureData
        if (!is.null(colGeometryNames(out))) {
            fds <- lapply(args, function(x) {
                cgns <- colGeometryNames(x)
                f <- lapply(cgns, geometryFeatureData, sfe = x, MARGIN = 2L)
                names(f) <- cgns
                f
            })
            new_cg_fds <- Reduce(.combine_reduceddim_fd, fds)
            cgns <- colGeometryNames(out)
            for (n in cgns) {
                geometryFeatureData(out, n, MARGIN = 2) <- new_cg_fds[[n]]
            }
        }
        # Combine colData featureData
        cd_fds <- lapply(args, colFeatureData)
        cd_fds_new <- Reduce(.combine_fd, cd_fds)
        colFeatureData(out) <- cd_fds_new

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
            names(sgs_use) <- NULL
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
